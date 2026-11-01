"""Convert an AnnData .h5ad file into a Raindrop dataset folder.

A dataset folder is static files that the viewer downloads on demand:
    manifest.json             names, types, categories, colors and file paths
    embeddings/<i>.bin        float32 positions, cells x dims (row major)
    fields/<i>.bin            one obs column: uint8/uint16 category codes or float32 values
    genes/measured/<i>.bin    one gene of X: dense uint16 counts or float32 values, one value per cell
Files are numbered rather than named after columns or genes, so any name is safe in a URL.
All arrays are little endian, as typed arrays in the browser expect.
"""
import argparse
import json
import os
import warnings

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

FORMAT_VERSION = 1
# Category name given to missing values (NaN or None) of a categorical obs column
MISSING_CATEGORY = 'NA'
# Name of the gene source written from X
MEASURED = 'measured'
DEFAULT_AXIS_LABELS = ('x', 'y', 'z')
# Direction drawn upward on screen for an embedding without --up: +y
DEFAULT_UP = [0, 1, 0]
# Suffix of obs columns that hold per-cell colors for another column (e.g. Class_color)
COLOR_SUFFIX = '_color'
# Largest values the category code and count dtypes hold
UINT8_MAX = np.iinfo(np.uint8).max
UINT16_MAX = np.iinfo(np.uint16).max


def write_array(path, array, dtype):
    """Write a numpy array as raw little-endian bytes, creating its folder."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    np.ascontiguousarray(array, dtype=np.dtype(dtype).newbyteorder('<')).tofile(path)


def color_columns(obs):
    """Per-cell color columns: '<column>_color' where '<column>' is also in obs."""
    return {c for c in obs.columns if c.endswith(COLOR_SUFFIX) and c[:-len(COLOR_SUFFIX)] in obs}


def exported_columns(obs, exclude):
    """obs columns to export, in obs order: all except excluded ones and color columns (which color others)."""
    skipped = color_columns(obs) | set(exclude)
    return [c for c in obs.columns if c not in skipped]


def embedding_keys(adata):
    """obsm keys that hold 2D or 3D positions, in obsm order."""
    return [k for k in adata.obsm if np.asarray(adata.obsm[k]).ndim == 2 and adata.obsm[k].shape[1] in (2, 3)]


def default_embedding(embeddings):
    """Embedding shown when the viewer opens if none is chosen."""
    # h5ad files list obsm keys alphabetically, so prefer the spatial coordinates when present
    return 'spatial' if 'spatial' in embeddings else embeddings[0]


def default_ends(axis_labels):
    """Generic names for the low and high end of each axis: '-x', '+x', ..."""
    return [[f'-{label}', f'+{label}'] for label in axis_labels]


def check_exclude(obs, exclude):
    """Raise if --exclude names a column that is not in obs."""
    unknown = [c for c in exclude if c not in obs]
    if unknown:
        raise ValueError(f'--exclude columns not in obs: {unknown}')


def check_columns(columns, color, filters):
    """Raise if --color or --filters names a column that is not exported."""
    for column in [color, *filters]:
        if column is not None and column not in columns:
            raise ValueError(f'{column} is not an exported obs column; choose from {columns}')


def check_embeddings(embeddings, embedding, up):
    """Raise if there is nothing to show, or --embedding or --up names an obsm entry that is not 2D or 3D."""
    if not embeddings:
        raise ValueError('adata.obsm has no 2D or 3D embedding to show')
    if embedding is not None and embedding not in embeddings:
        raise ValueError(f'{embedding} is not a 2D or 3D obsm entry; choose from {embeddings}')
    unknown = [k for k in up if k not in embeddings]
    if unknown:
        raise ValueError(f'--up names obsm entries that are not 2D/3D embeddings: {unknown}; choose from {embeddings}')


def check_out_dir(out_dir):
    """Raise if the output folder exists: a dataset folder is never overwritten."""
    if os.path.exists(out_dir):
        raise FileExistsError(f'{out_dir} already exists; remove it or choose another output folder')


def category_colors(obs, uns, column, categories):
    """One color per category, None where the data has none.

    Uses scanpy-style uns['<column>_colors'] first (when it has one color per category), then a per-cell
    obs['<column>_color'] column, which must give the same color to every cell of a category (an empty color cell
    becomes the text 'None' or 'nan').
    """
    if f'{column}_colors' in uns and len(uns[f'{column}_colors']) == len(categories):
        return [str(c) for c in uns[f'{column}_colors']]
    if f'{column}{COLOR_SUFFIX}' not in obs:
        return [None] * len(categories)
    per_category = obs.groupby(column, observed=True)[f'{column}{COLOR_SUFFIX}'].agg(lambda c: set(c.astype(str)))
    if (per_category.map(len) > 1).any():
        warnings.warn(f'{column}{COLOR_SUFFIX} has several colors for one category; the viewer will generate colors')
        return [None] * len(categories)
    # Categories without cells get no color; the viewer generates one
    return [next(iter(per_category[c])) if c in per_category.index else None for c in categories]


def category_codes(values):
    """Integer code per cell and the category names; missing values get an 'NA' category at the end."""
    values = values.astype('category')
    categories = [str(c) for c in values.cat.categories]
    # pandas stores codes as int8/int16 sized to the categories; int64 leaves room for the extra NA code
    codes = values.cat.codes.to_numpy().astype(np.int64)
    # pandas marks missing values with code -1
    if (codes < 0).any():
        codes[codes < 0] = len(categories)
        categories.append(MISSING_CATEGORY)
    return codes, categories


def code_dtype(n_categories, column):
    """Unsigned integer type of the category codes: uint8 for at most UINT8_MAX categories, uint16 for at most
    UINT16_MAX, else an error. (Each limit is one below what the type could hold.)"""
    if n_categories > UINT16_MAX:
        raise ValueError(f'obs column {column} has {n_categories} categories; at most {UINT16_MAX} are supported')
    return 'uint8' if n_categories <= UINT8_MAX else 'uint16'


def is_numeric(values):
    """Numbers are exported as numbers; everything else (strings, bools, categoricals) as categories."""
    return pd.api.types.is_numeric_dtype(values) and not pd.api.types.is_bool_dtype(values)


def finite_range(values):
    """[min, max] of the finite values, or [None, None] if there are none."""
    finite = values[np.isfinite(values)]
    return [float(finite.min()), float(finite.max())] if len(finite) else [None, None]


def gene_dtype(matrix):
    """uint16 if every stored value is a whole number from 0 to UINT16_MAX (counts, stored exactly), else float32."""
    data = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
    whole = np.all(data == np.round(data))
    in_range = data.size == 0 or (data.min() >= 0 and data.max() <= UINT16_MAX)
    return 'uint16' if whole and in_range else 'float32'


def write_embeddings(adata, keys, up, axis_labels, out_dir):
    """Write positions (float32) and return their manifest entries.

    Each entry names its axes, the low and high end of every axis (generic '-x', '+x', ... names) and the
    direction drawn upward on screen (from --up, else +y).
    """
    entries = []
    for i, key in enumerate(keys):
        positions = np.asarray(adata.obsm[key], dtype=np.float32)
        write_array(f'{out_dir}/embeddings/{i}.bin', positions, 'float32')
        entries.append({'name': key, 'file': f'embeddings/{i}.bin', 'dims': positions.shape[1],
                        'labels': list(axis_labels), 'ends': default_ends(axis_labels),
                        'up': list(up.get(key, DEFAULT_UP))})
    return entries


def write_fields(adata, columns, out_dir):
    """Write obs columns and return their manifest entries.

    Numbers are float32 with their finite range; everything else becomes category codes with names and colors.
    """
    entries = []
    for i, column in enumerate(columns):
        values = adata.obs[column]
        entry = {'name': column, 'file': f'fields/{i}.bin'}
        if is_numeric(values):
            numbers = values.to_numpy(dtype=np.float32)
            value_range = finite_range(numbers)
            entry.update(kind='numeric', dtype='float32', min=value_range[0], max=value_range[1])
            write_array(f"{out_dir}/{entry['file']}", numbers, 'float32')
        else:
            codes, categories = category_codes(values)
            colors = category_colors(adata.obs, adata.uns, column, categories)
            entry.update(kind='category', dtype=code_dtype(len(categories), column), categories=categories,
                         colors=colors)
            write_array(f"{out_dir}/{entry['file']}", codes, entry['dtype'])
        entries.append(entry)
    return entries


def write_genes(matrix, genes, name, folder, out_dir):
    """Write one dense file per gene of a cells x genes matrix and return the gene source's manifest entry.

    The entry's path has an {index} placeholder that the viewer fills with the gene's position in 'genes'.
    """
    # A CSC copy makes each gene's column a contiguous slice
    X = sparse.csc_matrix(matrix) if sparse.issparse(matrix) else np.asarray(matrix)
    dtype = gene_dtype(X)
    for i in range(len(genes)):
        column = X[:, i].toarray().ravel() if sparse.issparse(X) else X[:, i]
        write_array(f'{out_dir}/genes/{folder}/{i}.bin', column, dtype)
    return {'name': name, 'dtype': dtype, 'path': f'genes/{folder}/{{index}}.bin', 'genes': [str(g) for g in genes]}


def convert(h5ad_path, out_dir, *, name=None, color=None, filters=(), embedding=None,
            axis_labels=DEFAULT_AXIS_LABELS, up=None, exclude=()):
    """Write the dataset folder and return its manifest.

    `name` is shown in the viewer (default: the file name without extension); `color` and `filters` are the obs
    columns colored by and given filter panels when the viewer opens; `embedding` is the obsm key shown first
    (default: 'spatial' if present); `axis_labels` names the three axes; `up` maps obsm keys to the direction
    drawn upward on screen (default +y); `exclude` lists obs columns to leave out.
    Every request is checked before anything is written.
    """
    up = up or {}
    name = name or os.path.splitext(os.path.basename(h5ad_path))[0]

    # Check the request against the data
    adata = ad.read_h5ad(h5ad_path)
    check_exclude(adata.obs, exclude)
    columns = exported_columns(adata.obs, exclude)
    embeddings = embedding_keys(adata)
    check_columns(columns, color, filters)
    check_embeddings(embeddings, embedding, up)
    check_out_dir(out_dir)

    # Write the arrays, then the manifest that lists them
    manifest = {
        'format_version': FORMAT_VERSION,
        'name': name,
        'source': os.path.abspath(h5ad_path),
        'n_cells': int(adata.n_obs),
        'axes': {'labels': list(axis_labels)},
        'defaults': {'embedding': embedding or default_embedding(embeddings), 'color': color,
                     'filters': list(filters)},
        'embeddings': write_embeddings(adata, embeddings, up, axis_labels, out_dir),
        'fields': write_fields(adata, columns, out_dir),
        'gene_sources': [write_genes(adata.X, adata.var_names, MEASURED, MEASURED, out_dir)],
    }
    with open(f'{out_dir}/manifest.json', 'w') as f:
        json.dump(manifest, f, indent=1)
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('h5ad', help='input AnnData file')
    parser.add_argument('out_dir', help='dataset folder to create (must not exist)')
    parser.add_argument('--name', help='dataset name shown in the viewer (default: the file name)')
    parser.add_argument('--color', help='obs column to color by when the viewer opens')
    parser.add_argument('--filters', nargs='*', default=[],
                        help='obs columns that get a filter panel when the viewer opens')
    parser.add_argument('--embedding', help="obsm key shown when the viewer opens (default: 'spatial' if present)")
    parser.add_argument('--axis-labels', nargs=3, default=list(DEFAULT_AXIS_LABELS),
                        help='names of the three axes, e.g. AP DV LR')
    # Coordinate frames differ between embeddings (e.g. raw CCF has DV increasing downward), so each gets its own
    parser.add_argument('--up', nargs=4, action='append', default=[], metavar=('OBSM_KEY', 'X', 'Y', 'Z'),
                        help='direction drawn upward for one embedding, e.g. X_ccf 0 -1 0 '
                             '(default 0 1 0); repeatable')
    parser.add_argument('--exclude', nargs='*', default=[], help='obs columns to leave out')
    args = parser.parse_args()

    up = {key: [float(x), float(y), float(z)] for key, x, y, z in args.up}
    manifest = convert(args.h5ad, args.out_dir, name=args.name, color=args.color, filters=args.filters,
                       embedding=args.embedding, axis_labels=args.axis_labels, up=up, exclude=args.exclude)
    genes = ', '.join(f"{len(s['genes'])} {s['name']}" for s in manifest['gene_sources'])
    print(f"Wrote {args.out_dir}: {manifest['n_cells']:,} cells, {len(manifest['fields'])} obs columns, "
          f"genes: {genes}", flush=True)


if __name__ == '__main__':
    main()
