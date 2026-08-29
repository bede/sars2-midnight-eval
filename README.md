This notebook enables reproduction of results from the following manuscript:

*Rapid turnaround multiplex sequencing of SARS-CoV-2: comparing tiling amplicon protocol performance* https://www.medrxiv.org/content/10.1101/2021.12.28.21268461.full

Dependencies are pinned in `uv.lock`. Run from the repository root; the notebook reads `in/` and `wf-artic-results/`. The notebook also contains instructions to optionally regenerate results from raw sequence data deposited in [ENA BioProject XXXXX](#).

Run notebook interactively:

```bash
uv run jupyter lab results.ipynb
```

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

