This notebook enables reproduction of results from the following manuscript:

*Rapid turnaround multiplex sequencing of SARS-CoV-2: comparing tiling amplicon protocol performance* https://www.medrxiv.org/content/10.1101/2021.12.28.21268461.full

Dependencies are pinned in `uv.lock`. Run from the repository root; the notebook reads `in/` and `wf-artic-results/`. The notebook also contains instructions to optionally regenerate results from raw sequence data deposited in [ENA BioProject XXXXX](#).

Sample set 4 is an SGTF-selected December 2021 cohort that is mostly, but not entirely,
Omicron. The `plate_lineage` value in `in/runs.tsv` is retained for compatibility with the
original analysis and denotes the intended cohort, not independently established lineage
truth for every sample. The archived lineage reports currently classify 49 of its 71 positive
samples as Omicron, 15 as Delta, and leave 7 without a Scorpio call; the manuscript reports
48 Omicron assignments from its original Pangolin/Scorpio run.

Run notebook interactively:

```bash
uv run jupyter lab results.ipynb
```

Validate metadata keys, archived artifacts, sequence/depth structure, documented
missingness, and the derived `cov20.csv` table:

```bash
uv run python scripts/validate-data.py
```

Known absent workflow observations are documented in `in/artifact-exceptions.tsv`.
Unexpected missing artifacts cause validation and compression to fail. The three positive
clinical samples absent from the plate-2 Midnight run are excluded only from complete-case
protocol comparisons; they remain biologically positive in the metadata.

Regenerate figures into `results/`:

```bash
uv run jupyter nbconvert --to notebook --execute results.ipynb --output /tmp/out.ipynb
```

Swap `--output /tmp/out.ipynb` for `--inplace` to also refresh the outputs stored in the notebook.

`results/` holds main and supplementary figures.

To create or update the committed artifacts from uncompressed workflow output:

```bash
uv run python scripts/compress-results.py
```
