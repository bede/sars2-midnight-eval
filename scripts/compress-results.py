#!/usr/bin/env python3
"""Compress the result files used by the notebook."""

from pathlib import Path

import pandas as pd
import zstandard as zstd


repo = Path(__file__).resolve().parents[1]
results = repo / "wf-artic-results"

runs = pd.read_csv(repo / "in/runs.tsv", sep="\t")["run"].tolist()
barcodes = pd.read_csv(
    repo / "in/barcodes.tsv", sep="\t", dtype={"barcode": str}
)
exceptions = pd.read_csv(
    repo / "in/artifact-exceptions.tsv", sep="\t", dtype={"barcode": str}
)
allowed_missing_depth = {
    (row.run, row.barcode)
    for row in exceptions.itertuples(index=False)
    if "depth" in row.missing_artifacts.split(",")
}

files = [
    results / run / filename
    for run in runs
    for filename in ("all_consensus.fasta", "lineage_report.csv")
]
depth_files = [
    results / row.run / f"barcode{row.barcode}.depth_J.tsv"
    for row in barcodes.itertuples(index=False)
]
available_depth_files = [path for path in depth_files if path.is_file()]
unexpected_missing_depth = [
    path
    for path in depth_files
    if not path.is_file()
    and (path.parent.name, path.name[7:9]) not in allowed_missing_depth
]
if unexpected_missing_depth:
    formatted = "\n".join(f"  {path.relative_to(repo)}" for path in unexpected_missing_depth)
    raise FileNotFoundError(f"unexpected missing depth files:\n{formatted}")
files.extend(available_depth_files)

compressor = zstd.ZstdCompressor(
    level=10,
    write_checksum=True,
    write_content_size=True,
)
written = 0
for source in files:
    destination = Path(f"{source}.zst")
    compressed = compressor.compress(source.read_bytes())
    if not destination.is_file() or destination.read_bytes() != compressed:
        destination.write_bytes(compressed)
        written += 1

print(
    f"{len(files)} artifacts: {written} written; "
    f"{len(depth_files) - len(available_depth_files)} documented depth files absent"
)
