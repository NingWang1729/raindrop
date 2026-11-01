# Raindrop

A browser viewer for single-cell and spatial transcriptomics data, in the spirit of
[cirro](https://github.com/lilab-bcb/cirrocumulus): every cell is a drop in a 3D point cloud.
It follows the integration method [PIANO](https://github.com/NingWang1729/piano), and is named
after Chopin's "Raindrop" Prelude.

Unlike cirro, Raindrop needs no server program. A Python script converts an `.h5ad` file into a folder of
static files once; the viewer is plain HTML/JavaScript that downloads only the pieces it needs (positions,
the obs columns in use, and genes as they are picked) and does all filtering in the browser. Any static web
host can serve it, and any page can embed it.

## Converting a dataset

```bash
python convert.py data.h5ad data/my_dataset --name "My dataset" \
    --color Subclass --filters Subclass section --axis-labels AP DV LR --up spatial 0 -1 0
```

- Exports every `obs` column (numbers as ranges, everything else as categories), every 2D/3D `obsm`
  entry, and every gene in `X` (dense `uint16` if all values are whole counts up to 65535, else `float32`).
- `--color` and `--filters` choose what the viewer shows when it opens, and `--embedding` the embedding
  (default: `spatial` if present). `--axis-labels` names the three axes. `--up OBSM_KEY X Y Z` sets the direction
  drawn upward for one embedding (default `0 1 0`), since coordinate frames differ (raw CCF has DV increasing
  downward). `--exclude` leaves obs columns out.
- Colors come from scanpy-style `uns['<column>_colors']` or per-cell `obs['<column>_color']` columns.
- See `python convert.py --help` and the top of `convert.py` for the folder layout.

The repository's `data/` folder is ignored by git; datasets can live there or anywhere else the server can reach.

## Development

Install the test dependencies and run the unit tests:

```bash
uv venv .venv && uv pip install --python .venv/bin/python -e '.[test]'
.venv/bin/python -m pytest tests/                 # converter tests
```
