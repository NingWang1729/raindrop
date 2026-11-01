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
    --color Subclass --filters Subclass section region \
    --axis-labels AP DV LR --ends X_ccf anterior posterior dorsal ventral left right --view left
```

- Exports every `obs` column (numbers as ranges, everything else as categories), every 2D/3D `obsm`
  entry, and every gene in `X` (dense `uint16` if all values are whole counts up to 65535, else `float32`).
- `--color` and `--filters` choose what the viewer shows when it opens, and `--embedding` the embedding
  (default: `spatial` if present). `--exclude` leaves obs columns out. `--axis-labels` names the three axes of
  spatial embeddings; an `X_` embedding without `--ends` or `--roles` (UMAP, PCA, ...) names its axes after itself
  (UMAP 1, PC 1, t-SNE 1, or e.g. piano_umap 1).
- `--ends` names, per embedding, the low and high end of each axis; they label the side buttons, and an end named
  `dorsal` is drawn upward. Mouse CCF runs anterior to posterior, dorsal to ventral and left to right, hence
  `--ends X_ccf anterior posterior dorsal ventral left right`. `--view left` opens on the left side
  (anterior to the left). Embeddings with no end named dorsal and no `--roles` can set the upward direction with
  `--up OBSM_KEY X Y Z` instead (default `0 1 0`), since coordinate frames differ (raw CCF has DV increasing downward).
- Axis roles say which data axis is which anatomical axis, whatever the axes' order. When every end of an
  embedding has an anatomical name (anterior/posterior, dorsal/ventral, left/right, in any case), its roles follow
  from them; `--roles OBSM_KEY X Y Z` gives them for an embedding whose end names are not anatomical. A role names
  its axis's low end first: `AP` or `PA`, `DV` or `VD`, `LR` or `RL` (e.g. `--roles spatial LR PA VD` for x left to
  right, y posterior to anterior, z ventral to dorsal). Each anatomical axis appears once, roles must agree with
  anatomical end names, and an end name needs its partner (anterior with posterior, ...). With roles, dorsal is drawn
  upward.
- A mirror-image frame cannot show the left side with dorsal up and anterior to the left; prefer an unmirrored
  embedding for anatomy. A frame is a mirror image when its axes, mapped to the anatomy, are left-handed compared
  with CCF's (anterior to posterior, dorsal to ventral, left to right): an odd number of flipped axes plus swapped
  axis pairs, e.g. CCF with only DV flipped, or CCF with AP and DV swapped.
- `--locked` names two axis labels whose turns are locked when the viewer opens, so drags turn the data only about
  the third (e.g. `--locked AP LR`); the default locks none.
- `--hierarchy PARENT CHILD` groups a column's values by another obs column (e.g. `--hierarchy Class Subclass`);
  `--hierarchy-csv FILE` does the same from a two-column CSV (parent, child), whose parent need not be an obs
  column but which must cover every child value in the data. Each child value must have exactly one parent;
  values without one (missing values) go in an NA group. One hierarchy per column, two levels for now.
- Colors come from scanpy-style `uns['<column>_colors']` or per-cell `obs['<column>_color']` columns.
- `--gene-source NAME H5AD` adds the genes in `X` of another file as gene source NAME (repeatable), e.g. imputed
  genes. Its cells are matched to the converted file's by `obs_names`; a cell it does not name gets no value (NaN,
  drawn grey, the color of missing values), and a cell it names that the converted file lacks, or names twice,
  raises a ValueError. Values are stored as given (the viewer colors genes by log(1 + value)).
- `--log-scaled SOURCE SCALE` marks a gene source (`measured` for `X`, or a `--gene-source` name) whose values are
  already on a log scale, e.g. `--log-scaled measured log2` for log-normalized data: the viewer then colors them
  linearly instead of taking the log again, and legends name the scale (e.g. `Gene1 (log2)`). Repeatable.
- `--one-byte SOURCE` stores a gene source's fractional values in one byte per cell instead of four (float32): each
  gene's values are rounded to 254 steps from 0 to its top value (e.g. steps of 0.04 for log2 values up to 10),
  which colors cannot tell apart, and the viewer decodes them on download. For large datasets on static hosts
  (500 genes of 100,000 cells: 50 MB instead of 200 MB). Values must not be negative; whole counts of a source that
  covers every cell are already stored exactly in two bytes, so the option has no effect on them (a source that misses
  cells is stored as float32, so its counts are rounded to the steps too). Repeatable.
- See `python convert.py --help` and the top of `convert.py` for the folder layout.

The repository's `data/` folder is ignored by git; datasets can live there or anywhere else the server can reach.

## Development

Install the test dependencies and run the unit tests:

```bash
uv venv .venv && uv pip install --python .venv/bin/python -e '.[test]'
.venv/bin/python -m pytest tests/                 # converter tests
```
