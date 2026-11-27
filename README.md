# Raindrop

A browser viewer for single-cell and spatial transcriptomics data, in the spirit of
[cirro](https://github.com/lilab-bcb/cirrocumulus): every cell is a drop in a 3D point cloud.
It follows the integration method [PIANO](https://github.com/NingWang1729/piano), and is named
after Chopin's "Raindrop" Prelude.

Unlike cirro, Raindrop needs no server program. A Python script converts an `.h5ad` file into a folder of
static files once; the viewer is plain HTML/JavaScript that downloads only the pieces it needs (positions,
the obs columns in use, and genes as they are picked) and does all filtering in the browser. Any static web
host can serve it, and any page can embed it.

## Quick start

Install the converter's dependencies (Python 3.10 or later), from a clone of this repository:

```bash
git clone https://github.com/NingWang1729/raindrop.git
cd raindrop
python -m pip install -e .    # or just the dependencies: python -m pip install anndata numpy pandas scipy
```

Convert an `.h5ad` file into a dataset folder (the folder must not exist yet). The file needs at least one 2D or 3D
`obsm` entry to plot (e.g. `spatial` or `X_umap`); genes come from `X`.

```bash
python convert.py my.h5ad data/my_dataset
```

Serve the repository root with any static file server, for example Python's:

```bash
python -m http.server 8000
```

and open <http://localhost:8000/viewer/?data=../data/my_dataset/> in a browser.

Opening `viewer/index.html` straight from disk (a `file://` address) does not work: browsers block the viewer's
JavaScript modules and its downloads of the dataset's files from `file://` pages, so the files must come from a web
server.

**On a remote machine**, serve only to that machine (`python -m http.server 8000 --bind 127.0.0.1`) and forward the
port from your own computer:

```bash
ssh -L 8000:127.0.0.1:8000 <user>@<host>
```

then open the same address, <http://localhost:8000/viewer/?data=../data/my_dataset/>, on your own computer.

## Using the viewer

- **Turning and zooming.** Drag to turn the data like a turntable that can also tip over (left/right spins it about
  its up axis, up/down tilts it); Shift + drag spins or tilts only, whichever way the drag mostly goes. Scroll to
  zoom; right-drag to pan. 2D embeddings (e.g. UMAP) pan on a drag instead. On touch screens one finger turns (pans in
  2D), two fingers zoom and pan.
- **Embedding, sides and angles.** In the side panel's View section, the embedding menu picks what is plotted (e.g.
  `spatial` or `X_umap`), a button per axis end shows that side of the data (e.g. Left, Dorsal), the Turn about menu
  limits drags to one axis of the data, and the rotation table has a slider and an angle box per axis. Reset
  rotation sets every angle to 0 (Reset view for 2D embeddings). Show axes toggles the axes box and labels.
- **Color by** any obs column from the color menu, or by a gene: type its name in the gene box (with suggestions)
  and press Enter. Datasets with several gene sources have a menu to pick the source. Categories get colors from
  the data (see [Colors](#colors)) or generated ones; numbers and genes get a color bar.
- **Filters.** Add a filter panel per obs column from the "+ Add filter..." menu: tick values of a categorical
  column (a search box narrows long lists; drag over rows to tick several), or type a minimum and maximum for a
  numeric one. Filters combine: a cell is shown when it passes all of them. Each panel has **reverse** (show the cells
  it would hide instead), **clear** (show all values) and × (remove the panel), and sorts its values by name or by
  cell count. The count at the top says how many cells are shown.
- **Legend.** The panel to the right of the plots lists the colored column's values with their cell counts. Ticking
  its rows sets the filter on that column, with the same reverse and clear buttons. Grouped columns (see
  [Groups](#groups)) show as a tree in both the filter panel and the legend.
- **Shown and filtered-out cells.** The point size and opacity table sets both for the cells that pass the filters
  (Shown) and for the rest (Filtered out), each with a slider and a number box; "Filtered-out cells keep their colors"
  draws filtered-out cells in their colors rather than grey.
- **Views.** The Views menu shows up to 12 views side by side, each with its own embedding and coloring; click a
  view to select it, and the side panel and legend then apply to it. Linked views share filters and display
  settings, and views on the same embedding also share the camera and rotation; each view's chain button unlinks it
  (or links it again).
- **Hover label.** Pointing at a shown cell (tapping it on a touch screen) shows its value for the column or gene the
  view is colored by.

## Converter options

`python convert.py --help` lists every option. The converter exports every `obs` column (numbers as ranges,
everything else as categories; `<column>_color` columns only give colors, see [Colors](#colors)), every 2D or 3D
`obsm` entry, and every gene in `X` (stored as two-byte whole counts if all values are whole numbers up to 65535,
else as four-byte floats). An example with several options:

```bash
python convert.py my.h5ad data/my_dataset --name "My dataset" \
    --color cell_type --filters cell_type region --select region A B \
    --axis-labels AP DV LR --ends spatial anterior posterior dorsal ventral left right --view left
```

### What the viewer shows at start

- `--color COLUMN` colors by an obs column when the viewer opens (default: the first categorical column).
- `--filters COLUMN ...` opens filter panels for these obs columns.
- `--select COLUMN VALUE ...` ticks values of a `--filters` column when the viewer opens (repeatable, once per
  column). Afterwards the panels work as usual: `clear` shows all values (reload for the starting ticks).
- `--start-style KEY=VALUE ...` sets the starting Shown / Filtered out style: `shown-size` and `filtered-size` as
  multiples of the default point size (1/16 to 8; the slider takes the nearest 0.1 step in log2), `shown-opacity`
  and `filtered-opacity` in percent (0 to 100; the slider takes the nearest whole percent), and `keep-colors` on or
  off; keys left out keep the viewer's own starting values. For example `--filters region --select region A B C
  --start-style shown-size=1 shown-opacity=100 filtered-size=0.5 filtered-opacity=10` shows regions A, B and C large
  and opaque inside the faint rest of the data.
- `--embedding OBSM_KEY` is the embedding shown when the viewer opens (default: `spatial` if present, else the first
  2D or 3D `obsm` entry).
- `--exclude COLUMN ...` leaves obs columns out.
- `--name NAME` is the dataset name shown in the viewer (default: the file name without `.h5ad`).

### Axes and anatomy

- `--axis-labels` names the three axes of spatial embeddings (e.g. `AP DV LR`; default `x y z`). An `X_` embedding
  without `--ends` or `--roles` (UMAP, PCA, ...) names its axes after itself (UMAP 1, PC 1, t-SNE 1, or e.g.
  `X_my_umap` gives my_umap 1).
- `--ends OBSM_KEY LOW_X HIGH_X LOW_Y HIGH_Y LOW_Z HIGH_Z` names, per embedding, the low and high end of each axis;
  they label the side buttons, and an end named `dorsal` is drawn upward. The common coordinate framework (CCF)
  runs anterior to posterior, dorsal to ventral and left to right, so an embedding in CCF coordinates takes
  `--ends spatial anterior posterior dorsal ventral left right`, and `--view left` then opens on the left side
  (anterior to the left, dorsal up).
- `--view END` is the end name of the side shown when the viewer opens (default: the low x end).
- `--up OBSM_KEY X Y Z` sets the direction drawn upward (default `0 1 0`) for an embedding with no end named dorsal
  and no `--roles`, since coordinate frames differ (raw CCF has dorsal to ventral increasing downward:
  `--up spatial 0 -1 0`).
- `--roles OBSM_KEY X Y Z` says which data axis is which anatomical axis, whatever the axes' order. When every end
  of an embedding has an anatomical name (anterior/posterior, dorsal/ventral, left/right, in any case), its roles
  follow from them; `--roles` gives them for an embedding whose end names are not anatomical. A role names its
  axis's low end first: `AP` or `PA`, `DV` or `VD`, `LR` or `RL` (e.g. `--roles spatial LR PA VD` for x left to
  right, y posterior to anterior, z ventral to dorsal). Each anatomical axis appears once, roles must agree with
  anatomical end names, and an end name needs its partner (anterior with posterior, ...). With roles, dorsal is
  drawn upward.
- A mirror-image frame cannot show the left side with dorsal up and anterior to the left; prefer an unmirrored
  embedding for anatomy. A frame is a mirror image when its axes, mapped to the anatomy, are left-handed compared
  with CCF's (anterior to posterior, dorsal to ventral, left to right): an odd number of flipped axes plus swapped
  axis pairs, e.g. CCF with only dorsal-ventral flipped, or CCF with AP and DV swapped.
- `--locked LABEL LABEL` names two axis labels whose turns are locked when the viewer opens, so drags turn the data
  only about the third (e.g. `--locked AP LR`); the default locks none.

### Groups

- `--hierarchy PARENT CHILD` groups a column's values by another obs column (e.g. `--hierarchy cell_class
  cell_type`); filter panels and the legend then show the values as a tree.
- `--hierarchy-csv FILE` does the same from a two-column CSV (parent, child), whose parent need not be an obs column
  but which must cover every child value in the data.
- Each child value must have exactly one parent; values without one (missing values) go in an NA group. One
  hierarchy per column, two levels for now. Both options are repeatable.

### Genes

- `--gene-source NAME H5AD` adds the genes in `X` of another file as gene source NAME (repeatable), e.g. imputed
  genes. Its cells are matched to the converted file's by `obs_names`; a cell it does not name gets no value (drawn
  grey, the color of missing values), and a cell it names that the converted file lacks, or names twice, raises an
  error. Values are stored as given (the viewer colors genes by log(1 + value)).
- `--log-scaled SOURCE SCALE` marks a gene source (`measured` for `X`, or a `--gene-source` name) whose values are
  already on a log scale, e.g. `--log-scaled measured log2` for log-normalized data: the viewer then colors them
  linearly instead of taking the log again, and legends name the scale (e.g. `Gene1 (log2)`). Repeatable.
- `--one-byte SOURCE` stores a gene source's fractional values in one byte per cell instead of four: each gene's
  values are rounded to 254 steps from 0 to its top value (e.g. steps of 0.04 for log2 values up to 10), which
  colors cannot tell apart, and the viewer decodes them on download. This suits large datasets on static hosts (500
  genes of 100,000 cells: 50 MB instead of 200 MB). Values must not be negative; whole counts of a source that covers
  every cell are already stored exactly in two bytes, so the option has no effect on them (a source that misses
  cells is stored as floats, so its counts are rounded to the steps too). Repeatable.

### Colors

Category colors come from scanpy-style `uns['<column>_colors']` (one color per category) or from a per-cell
`obs['<column>_color']` column that gives every cell of a category the same color. Categories without colors get
generated ones.

The dataset folder's layout is described at the top of `convert.py` and in `python convert.py --help`.

## Hosting and embedding

Raindrop consists only of static files, so any static web host can serve it (a lab web server, GitHub Pages, cloud
storage serving a website, ...). Put `viewer/` and the dataset folders on the host, and open
`viewer/?data=<dataset folder>`. The dataset folder must be served alongside the viewer, or be reachable from it:
`?data=` is a URL resolved relative to the viewer page (`viewer/`), so `../data/my_dataset/` is the folder
`data/my_dataset` beside `viewer/`. An absolute URL works too; a dataset on another host must allow cross-origin
requests (CORS), since the viewer downloads its files with JavaScript. The repository's `data/` folder is ignored
by git, so converted datasets there are not committed. A dataset's `manifest.json` records the source of the
`.h5ad` file it was converted from.

### URL parameters

| Parameter | What it does |
| --- | --- |
| `data=<dataset folder URL>` | The dataset to open (required). |
| `embed=1` or `embed=true` | Widget mode for another page's iframe: a bar over the view with a More / Fewer controls toggle and a Full screen button. A plain scroll over the view scrolls the host page (with a short hint); Ctrl + scroll (Cmd + scroll on a Mac) or a trackpad pinch zooms. `embed=0`, `embed=false` or no `embed` opens the full viewer. Values are read in any case (`embed=TRUE`); any other value (an empty one included) stops the page with a message naming the allowed values. |
| `controls=compact` or `controls=full` | In widget mode, the controls the widget opens with, in any case (`controls=Full`). `compact` (the default) hides the side panel and puts the cell count, the color menu, the gene box and the reset button in the bar, which suits narrow frames (the legend still filters by ticking values); `full` shows the side panel with every control. In widget mode any other value stops the page with a message naming the allowed values. Without widget mode it has no effect, whatever its value: the viewer opens as usual and warns in the browser console and, once the dataset has opened, in the status line. |

### Embedding in a page

```html
<iframe src="https://<host>/viewer/?data=../data/my_dataset/&embed=1"
        width="100%" height="600" allowfullscreen></iframe>
```

Give the iframe `allowfullscreen` so the browser lets the widget go full screen; where the browser does not allow
full screen in the frame, the Full screen button is left out.

## Development

Install the test dependencies and run the unit tests:

```bash
uv venv .venv && uv pip install --python .venv/bin/python -e '.[test]'
.venv/bin/python -m pytest tests/                 # converter tests
node tests/test_viewer.mjs                        # viewer filtering, color and position logic (Node.js)
```

three.js is vendored in `viewer/vendor/` (MIT license, see `THREE_LICENSE`), so the viewer has no build step and no
CDN dependency.
