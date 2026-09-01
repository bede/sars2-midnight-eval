#!/usr/bin/env python3
"""Validate metadata and committed workflow artifacts used by the analysis."""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import pandas as pd
import zstandard as zstd
from Bio import SeqIO


REPO = Path(__file__).resolve().parents[1]
RESULTS = REPO / "wf-artic-results"
GENOME_LENGTH = 29_903
CONTROL_SAMPLES = {"WHO", "TWIST", "NEG", "NEG1", "NEG2", "POS", "PC", "PC2", "NC", "IQC"}
EXCLUDED_SAMPLES = {"b1da6318c02e644d", "773df88e2e8c8412"}


class ValidationErrors:
    def __init__(self) -> None:
        self.messages: list[str] = []

    def check(self, condition: bool, message: str) -> None:
        if not condition:
            self.messages.append(message)

    def finish(self) -> None:
        if self.messages:
            details = "\n".join(f"- {message}" for message in self.messages)
            raise SystemExit(f"data validation failed:\n{details}")


def open_result(path: Path):
    if path.is_file():
        return path.open("rt", encoding="utf-8")
    compressed = Path(f"{path}.zst")
    if compressed.is_file():
        return zstd.open(compressed, "rt", encoding="utf-8")
    raise FileNotFoundError(path)


def artifact_status(run_dir: Path, barcode: str, fasta_barcodes: set[str], lineage_barcodes: set[str]) -> set[str]:
    missing: set[str] = set()
    if barcode not in fasta_barcodes:
        missing.add("consensus")
    if barcode not in lineage_barcodes:
        missing.add("lineage")
    if not (run_dir / f"barcode{barcode}.depth_J.tsv.zst").is_file():
        missing.add("depth")
    return missing


def main() -> None:
    errors = ValidationErrors()
    runs = pd.read_csv(REPO / "in/runs.tsv", sep="\t")
    barcodes = pd.read_csv(REPO / "in/barcodes.tsv", sep="\t", dtype={"barcode": str})
    ct = pd.read_csv(REPO / "in/ct.tsv", sep="\t")
    exceptions = pd.read_csv(
        REPO / "in/artifact-exceptions.tsv", sep="\t", dtype={"barcode": str}
    )

    for name, frame, keys in (
        ("runs", runs, ["run"]),
        ("barcodes", barcodes, ["run", "barcode"]),
        ("Ct", ct, ["plate", "sample"]),
        ("artifact exceptions", exceptions, ["run", "barcode"]),
    ):
        duplicates = frame.duplicated(keys, keep=False)
        errors.check(not duplicates.any(), f"{name} has duplicate key(s) {keys}")

    errors.check(set(barcodes["run"]) == set(runs["run"]), "run sets differ between runs.tsv and barcodes.tsv")
    errors.check(
        set(zip(exceptions["run"], exceptions["barcode"])).issubset(
            set(zip(barcodes["run"], barcodes["barcode"]))
        ),
        "artifact exception references an unknown run/barcode",
    )

    plate_lineages = runs.groupby("plate")["plate_lineage"].nunique()
    errors.check(plate_lineages.eq(1).all(), "plate_lineage is inconsistent within a plate")

    exception_map = {
        (row.run, row.barcode): set(row.missing_artifacts.split(","))
        for row in exceptions.itertuples(index=False)
    }
    observed_missing: dict[tuple[str, str], set[str]] = {}
    lineage_frames: list[pd.DataFrame] = []
    depth_count = 0
    expected_cov20: dict[tuple[str, str], float] = {}

    for run_row in runs.itertuples(index=False):
        run_dir = RESULTS / run_row.run
        with open_result(run_dir / "all_consensus.fasta") as handle:
            records = list(SeqIO.parse(handle, "fasta"))
        fasta_barcodes = {record.id.removeprefix("barcode") for record in records}
        errors.check(len(fasta_barcodes) == len(records), f"{run_row.run}: duplicate FASTA identifiers")
        for record in records:
            sequence = str(record.seq).upper()
            invalid = set(sequence) - set("ACGTNRYKMSWBDHV-")
            errors.check(not invalid, f"{run_row.run}/{record.id}: invalid consensus symbols {sorted(invalid)}")

        with open_result(run_dir / "lineage_report.csv") as handle:
            lineage = pd.read_csv(handle)
        lineage["barcode"] = (
            lineage["taxon"].astype(str).str.replace("_MN908947.3", "", regex=False).str.removeprefix("barcode")
        )
        lineage["lineage"] = lineage["lineage"].replace("None", pd.NA)
        lineage["run"] = run_row.run
        lineage["plate"] = run_row.plate
        lineage["sample"] = lineage["barcode"].map(
            barcodes.loc[barcodes["run"].eq(run_row.run)].set_index("barcode")["sample"]
        )
        lineage_frames.append(lineage)
        lineage_barcodes = set(lineage["barcode"])
        errors.check(not lineage["barcode"].duplicated().any(), f"{run_row.run}: duplicate lineage barcode")

        expected = barcodes.loc[barcodes["run"].eq(run_row.run), ["barcode"]]
        for barcode in expected["barcode"]:
            missing = artifact_status(run_dir, barcode, fasta_barcodes, lineage_barcodes)
            if missing:
                observed_missing[(run_row.run, barcode)] = missing

        for depth_path in run_dir.glob("barcode*.depth_J.tsv.zst"):
            depth_count += 1
            with zstd.open(depth_path, "rt", encoding="utf-8") as handle:
                depth = pd.read_csv(handle, sep="\t", header=None, usecols=[1, 2], names=["pos", "depth"])
            valid_positions = depth["pos"].equals(pd.Series(range(1, GENOME_LENGTH + 1)))
            errors.check(valid_positions, f"{depth_path.relative_to(REPO)}: positions are not exactly 1..{GENOME_LENGTH}")
            valid_depth = depth["depth"].notna().all() and depth["depth"].ge(0).all() and depth["depth"].mod(1).eq(0).all()
            errors.check(valid_depth, f"{depth_path.relative_to(REPO)}: invalid depth value")
            barcode = depth_path.name.removeprefix("barcode")[:2]
            expected_cov20[(run_row.run, barcode)] = depth["depth"].ge(20).sum() / GENOME_LENGTH * 100

    errors.check(observed_missing == exception_map, f"observed missing artifacts differ from manifest: {observed_missing}")

    positive = ct.loc[
        ct["orf1_ct"].gt(0)
        & ~ct["sample"].isin(CONTROL_SAMPLES | EXCLUDED_SAMPLES),
        ["plate", "sample"],
    ]
    missing_positive = exceptions.merge(barcodes, on=["run", "barcode"], validate="one_to_one").merge(
        runs[["run", "plate"]], on="run", validate="many_to_one"
    ).merge(positive, on=["plate", "sample"], how="inner", validate="many_to_one")
    errors.check(len(missing_positive) == 3, f"expected 3 documented missing positive observations, found {len(missing_positive)}")

    lineages = pd.concat(lineage_frames, ignore_index=True)
    plate4_samples = set(positive.loc[positive["plate"].eq(4), "sample"])
    plate4 = lineages.loc[lineages["plate"].eq(4) & lineages["sample"].isin(plate4_samples)].copy()
    plate4["scorpio"] = plate4["scorpio_call"].str.partition(" ")[0]
    sample_calls = plate4.groupby("sample")["scorpio"].agg(lambda values: next(iter(values.dropna().mode()), pd.NA))
    call_counts = Counter(sample_calls.fillna("uncalled"))

    derived = pd.read_csv(REPO / "cov20.csv", dtype={"barcode": str})
    errors.check(not derived.duplicated(["run", "barcode"]).any(), "cov20.csv has duplicate run/barcode rows")
    errors.check(derived["cov20"].between(0, 100).all(), "cov20.csv contains an out-of-range cov20 value")
    errors.check(derived["acgt"].between(0, 1).all(), "cov20.csv contains an out-of-range acgt value")
    expected_derived_keys = {
        (row.run, row.barcode)
        for row in barcodes.dropna(subset=["sample"]).itertuples(index=False)
        if (row.run, row.barcode) in expected_cov20
    }
    derived_keys = set(zip(derived["run"], derived["barcode"]))
    errors.check(derived_keys == expected_derived_keys, "cov20.csv run/barcode keys differ from available mapped depth results")
    cov20_lookup = derived.set_index(["run", "barcode"])["cov20"]
    errors.check(
        all(abs(cov20_lookup.loc[key] - expected_cov20[key]) < 1e-10 for key in expected_derived_keys),
        "cov20.csv contains a value inconsistent with an archived depth file",
    )

    run_plate = runs.set_index("run")["plate"].to_dict()
    sample_for_key = barcodes.set_index(["run", "barcode"])["sample"].to_dict()
    available_by_run = {
        run: {
            sample_for_key[key]
            for key in expected_cov20
            if key[0] == run and pd.notna(sample_for_key.get(key))
        }
        for run in runs["run"]
    }
    complete_by_plate = {
        plate: set.intersection(*(available_by_run[run] for run in runs.loc[runs["plate"].eq(plate), "run"]))
        for plate in runs["plate"].unique()
    }
    expected_eligibility = [
        row.sample in complete_by_plate[run_plate[row.run]]
        and (run_plate[row.run], row.sample) in set(map(tuple, positive.to_numpy()))
        for row in derived.itertuples(index=False)
    ]
    errors.check(
        derived["analysis_eligible"].tolist() == expected_eligibility,
        "cov20.csv analysis_eligible flags are inconsistent with complete-case availability",
    )

    errors.finish()
    print(f"OK: {len(runs)} runs, {len(barcodes)} barcode mappings, {depth_count} depth files")
    print(f"OK: {len(exception_map)} documented missing observations; {len(missing_positive)} are positive clinical samples")
    print(f"Sample set 4 Scorpio classifications: {dict(call_counts)}")


if __name__ == "__main__":
    main()
