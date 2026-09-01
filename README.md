This notebook enables reproduction of results from the manuscript:

*Rapid turnaround multiplex sequencing of SARS-CoV-2: comparing tiling amplicon protocol performance* https://www.medrxiv.org/content/10.1101/2021.12.28.21268461.full

Plots can be regenerated from intermediate files inside `in/` and `wf-artic-results/`. The notebook also contains instructions to optionally regenerate results from raw sequence data deposited in [ENA BioProject PRJEB124656](https://www.ebi.ac.uk/ena/browser/view/PRJEB124656).

**To run notebook interactively:**

```bash
uv run jupyter lab results.ipynb
```

Validate metadata keys, archived artifacts, sequence/depth structure, documented
missingness, and the derived `cov20.csv` table:

```bash
uv run python scripts/validate-data.py
```

**To regenerate figures into `results/`:**

```bash
uv run jupyter nbconvert --to notebook --execute results.ipynb --output /tmp/out.ipynb
```

Swap `--output /tmp/out.ipynb` for `--inplace` to also refresh the outputs stored in the notebook.

Figures are saved inside `results/`.

**To update committed artifacts from workflow outputs:**

```bash
uv run python scripts/compress-results.py
```
