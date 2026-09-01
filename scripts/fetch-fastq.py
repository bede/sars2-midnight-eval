#!/usr/bin/env python3
"""Fetch PRJEB124656"""

import argparse
import concurrent.futures
import csv
import hashlib
import os
import re
import sys
import urllib.parse
import urllib.request
from pathlib import Path

PORTAL = "https://www.ebi.ac.uk/ena/portal/api/filereport"
# Submitted files are named s2se__<run>__barcode<NN>.fastq.gz
PATTERN = re.compile(r"s2se__(.+)__barcode(\d{2})\.fastq\.gz")

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--project", default="PRJEB124656")
parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[1] / "fastq")
parser.add_argument("--jobs", type=int, default=8)
parser.add_argument("--dry-run", action="store_true", help="List destinations without downloading.")


def md5(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch(item: tuple[str, str, Path]) -> None:
    url, expected, target = item
    if target.is_file() and md5(target) == expected:
        return
    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_name(target.name + ".part")
    try:
        with urllib.request.urlopen(url, timeout=300) as response, partial.open("wb") as handle:
            while chunk := response.read(1 << 20):
                handle.write(chunk)
        if md5(partial) != expected:
            raise ValueError(f"MD5 mismatch, expected {expected}")
        os.replace(partial, target)
    except Exception:
        partial.unlink(missing_ok=True)
        raise


def main() -> None:
    args = parser.parse_args()
    output = args.output.resolve()
    query = urllib.parse.urlencode(
        {
            "accession": args.project,
            "result": "read_run",
            "fields": "submitted_ftp,submitted_md5",
            "format": "tsv",
            "download": "true",
        }
    )
    with urllib.request.urlopen(f"{PORTAL}?{query}", timeout=120) as response:
        rows = list(csv.DictReader(response.read().decode().splitlines(), delimiter="\t"))

    work = []
    for row in rows:
        url, checksum = row["submitted_ftp"], row["submitted_md5"]
        match = PATTERN.fullmatch(url.rsplit("/", 1)[-1])
        if not match or not checksum:
            raise SystemExit(f"unexpected ENA record: {row}")
        run, barcode = match.groups()
        work.append((f"https://{url}", checksum, output / run / f"barcode{barcode}" / f"barcode{barcode}.fastq.gz"))
    if not work:
        raise SystemExit(f"no runs found for {args.project}")

    if args.dry_run:
        for _, _, target in sorted(work, key=lambda item: item[2]):
            print(target.relative_to(output))
        return

    failures = []
    with concurrent.futures.ThreadPoolExecutor(args.jobs) as executor:
        futures = {executor.submit(fetch, item): item[2] for item in work}
        for done, future in enumerate(concurrent.futures.as_completed(futures), 1):
            try:
                future.result()
            except Exception as error:
                failures.append(f"{futures[future]}: {error}")
            print(f"\r{done}/{len(work)}", end="", file=sys.stderr)
    print(file=sys.stderr)
    if failures:
        raise SystemExit("failed:\n- " + "\n- ".join(sorted(failures)))
    print(f"restored {len(work)} FASTQs beneath {output}")


if __name__ == "__main__":
    main()
