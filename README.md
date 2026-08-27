This notebook enables reproduction of results from the following manuscript:

*Rapid turnaround multiplex sequencing of SARS-CoV-2: comparing tiling amplicon protocol performance* https://www.medrxiv.org/content/10.1101/2021.12.28.21268461.full

Dependencies are pinned in `uv.lock`. Run from the repository root; the notebook reads `in/` and `wf-artic/results-barcode/`.

Interactive:

```bash
uv run jupyter lab results.ipynb
```

Command line — regenerates every figure into `results/`, exiting non-zero if any cell fails (~20s):

```bash
uv run jupyter nbconvert --to notebook --execute results.ipynb --output /tmp/out.ipynb
```

Swap `--output /tmp/out.ipynb` for `--inplace` to also refresh the outputs stored in the notebook.

`results/` holds `f1.pdf`–`f5.svg` (Figures 1–5), `ct_dist.pdf` and `supp.cov20_vs_acgt.pdf` (Supplementary Figures 1–2), plus `cov20.csv` and `problems.csv`. `f1.1.pdf` and `eccmid.f1.png` are unused variants. `cov20.csv` should reproduce the committed root-level `cov20.csv` exactly.
