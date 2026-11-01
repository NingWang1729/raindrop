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

# An axis end with this name (in any case) is drawn upward on screen, as an alternative to --up
UP_END = 'dorsal'
# Anatomical end names (in any case) and their letters. An axis's role names its anatomical axis by these letters,
# low end first: AP runs from anterior to posterior, PA from posterior to anterior
END_LETTERS = {'anterior': 'A', 'posterior': 'P', 'dorsal': 'D', 'ventral': 'V', 'left': 'L', 'right': 'R'}
# The six roles and the anatomical axis each names (both directions of an axis count as the same axis)
AXIS_OF_ROLE = {'AP': 'AP', 'PA': 'AP', 'DV': 'DV', 'VD': 'DV', 'LR': 'LR', 'RL': 'LR'}
# Axis names of common non-spatial embeddings, by obsm key without 'X_', lowercased
METHOD_AXIS_NAMES = {'umap': 'UMAP', 'tsne': 't-SNE', 'pca': 'PC'}


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


def embedding_axis_labels(key, axis_labels, ends, roles):
    """Names of an embedding's three axes.

    An 'X_' key (UMAP, PCA, ...) without anatomical end names (--ends) or roles (--roles) is not spatial: its axes
    are named after the method ('UMAP 1', 'PC 1', 't-SNE 1'; other keys by their name without 'X_', e.g.
    'piano_umap 1'). Every other embedding takes the dataset's axis labels (--axis-labels).
    """
    if not key.startswith('X_') or key in ends or key in roles:
        return list(axis_labels)
    method = key[2:]
    name = METHOD_AXIS_NAMES.get(method.lower(), method)
    return [f'{name} {k}' for k in (1, 2, 3)]


def has_up_end(ends):
    """Whether one of an embedding's axis ends is named 'dorsal' (in any case)."""
    return any(name.lower() == UP_END for pair in ends for name in pair)


def up_end_direction(ends):
    """Unit vector toward the end named 'dorsal', e.g. [0, -1, 0] when dorsal is the low end of y."""
    for axis, pair in enumerate(ends):
        # The low end of an axis points toward -1, the high end toward +1 (in the embedding's own data coordinates)
        for side, name in zip((-1, 1), pair):
            if name.lower() == UP_END:
                direction = [0, 0, 0]
                direction[axis] = side
                return direction
    raise ValueError(f'no axis end is named {UP_END}')


def has_anatomical_ends(pair):
    """Whether either end of an axis has an anatomical name."""
    return any(name.lower() in END_LETTERS for name in pair)


def role_from_ends(key, pair):
    """Role of an axis from its end names, e.g. ['posterior', 'anterior'] -> 'PA'."""
    # A non-anatomical name becomes '?', so a half-anatomical pair (e.g. anterior, top) is not a role
    role = ''.join(END_LETTERS.get(name.lower(), '?') for name in pair)
    if role not in AXIS_OF_ROLE:
        raise ValueError(f'--ends {key}: {pair} are not the two ends of one anatomical axis; name both ends '
                         f'anatomically (e.g. anterior posterior) or neither')
    return role


def has_roles(key, roles, names):
    """Whether an embedding's axes have anatomical roles: from --roles, or from anatomical names of all its ends."""
    return key in roles or all(has_anatomical_ends(pair) for pair in names)


def axis_roles(key, roles, names):
    """Anatomical role of each axis (e.g. ['LR', 'PA', 'VD']): from --roles, else from the end names."""
    return roles[key] if key in roles else [role_from_ends(key, pair) for pair in names]


def role_end_direction(roles, letter):
    """Unit vector toward one anatomical end, e.g. 'D' (dorsal) -> [0, 0, 1] for roles LR PA VD."""
    # The checks give each of AP, DV and LR to exactly one of the three axes, so exactly one role holds the letter
    axis = next(i for i, role in enumerate(roles) if letter in role)
    # A role names its low end first, so the end's letter in first place points toward -1 (in the embedding's own
    # data coordinates; the viewer turns the camera so this vector points up on screen)
    direction = [0, 0, 0]
    direction[axis] = -1 if roles[axis][0] == letter else 1
    return direction


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


def check_axes_once(option, key, roles):
    """Raise if roles name an anatomical axis more than once (e.g. AP and PA)."""
    axes = [AXIS_OF_ROLE[role] for role in roles]
    repeated = sorted({axis for axis in axes if axes.count(axis) > 1})
    if repeated:
        raise ValueError(f'{option} {key} names the {repeated[0]} axis more than once; give each of AP, DV and LR '
                         f'to one axis')


def check_end_pairs(key, names):
    """Raise if an embedding's anatomical end names do not pair up into distinct anatomical axes."""
    check_axes_once('--ends', key, [role_from_ends(key, pair) for pair in names if has_anatomical_ends(pair)])


def check_ends(embeddings, ends, up):
    """Raise if --ends names an obsm entry that is not 2D or 3D, sets the up direction twice (an end named dorsal
    and --up), or names ends that are not pairs of one anatomical axis."""
    unknown = [k for k in ends if k not in embeddings]
    if unknown:
        raise ValueError(f'--ends names obsm entries that are not 2D/3D embeddings: {unknown}; choose from {embeddings}')
    both = [k for k in ends if k in up and has_up_end(ends[k])]
    if both:
        raise ValueError(f'{both}: an end named {UP_END} already sets the up direction; drop --up for these embeddings')
    for key, names in ends.items():
        check_end_pairs(key, names)


def check_role_list(key, roles, names):
    """Raise if one embedding's --roles are unknown, repeat an anatomical axis or disagree with its end names."""
    unknown = [role for role in roles if role not in AXIS_OF_ROLE]
    if unknown:
        raise ValueError(f'--roles {key}: {unknown} are not axis roles; give each axis one of {list(AXIS_OF_ROLE)} '
                         f'(its low end first, e.g. PA for posterior to anterior)')
    check_axes_once('--roles', key, roles)
    for axis, (role, pair) in enumerate(zip(roles, names)):
        if has_anatomical_ends(pair) and role_from_ends(key, pair) != role:
            raise ValueError(f'--roles {key}: axis {axis} is {role} but its --ends are {pair}; make them agree or '
                             f'drop one of them')


def check_roles(embeddings, roles, up, ends, axis_labels):
    """Raise if --roles names an obsm entry that is not 2D or 3D, sets the up direction twice (roles and --up), or
    gives roles that check_role_list rejects. Embeddings without --ends are checked against generic end names."""
    unknown = [k for k in roles if k not in embeddings]
    if unknown:
        raise ValueError(f'--roles names obsm entries that are not 2D/3D embeddings: {unknown}; choose from {embeddings}')
    both = [k for k in roles if k in up]
    if both:
        raise ValueError(f'{both}: --roles already sets the up direction (toward dorsal); drop --up for these embeddings')
    for key, given in roles.items():
        check_role_list(key, given, ends.get(key, default_ends(axis_labels)))


def check_view(view, embedding, embeddings, ends, axis_labels):
    """Raise if --view is not an end name of the embedding shown when the viewer opens.

    An embedding without --ends is checked against the generic names of the dataset's axis labels ('-x', ...),
    also for an 'X_' embedding, whose written ends use its method's axis names instead.
    """
    if view is None:
        return
    shown = embedding or default_embedding(embeddings)
    names = [name for pair in ends.get(shown, default_ends(axis_labels)) for name in pair]
    if view not in names:
        raise ValueError(f'--view {view} is not an end of embedding {shown}; choose from {names}')


def check_locked(locked, axis_labels):
    """Raise if --locked names an unknown axis label, or locks other than two axes or none."""
    for label in locked:
        if label not in axis_labels:
            raise ValueError(f'--locked {label} is not an axis label; choose from {list(axis_labels)}')
    # The viewer turns the data about any axis (none locked) or about the one axis left free (two locked)
    if len(set(locked)) not in (0, 2):
        raise ValueError(f'--locked {" ".join(locked)}: lock two axes (drags then turn about the third) or none')


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


def write_embeddings(adata, keys, up, ends, roles, axis_labels, out_dir):
    """Write positions (float32) and return their manifest entries.

    Each entry names its axes (embedding_axis_labels), the low and high end of every axis (for the side buttons;
    from --ends, else generic '-x', '+x', ... after the axis names), the anatomical role of each axis (None if
    unknown; the viewer orients side views by them) and the direction drawn upward on screen: toward dorsal if the
    axes have roles or an end is named 'dorsal', else from --up, else +y.
    """
    entries = []
    for i, key in enumerate(keys):
        positions = np.asarray(adata.obsm[key], dtype=np.float32)
        write_array(f'{out_dir}/embeddings/{i}.bin', positions, 'float32')
        labels = embedding_axis_labels(key, axis_labels, ends, roles)
        names = ends.get(key, default_ends(labels))
        if has_roles(key, roles, names):
            known = axis_roles(key, roles, names)
            direction = role_end_direction(known, 'D')
        else:
            known = None
            direction = up_end_direction(names) if has_up_end(names) else list(up.get(key, DEFAULT_UP))
        entries.append({'name': key, 'file': f'embeddings/{i}.bin', 'dims': positions.shape[1], 'labels': labels,
                        'ends': names, 'roles': known, 'up': direction})
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


def hierarchy_pairs_from_obs(obs, columns, parent, child):
    """Distinct (parent, child) value pairs of the cells, missing values dropped, for --hierarchy PARENT CHILD."""
    for column in (parent, child):
        if column not in columns or is_numeric(obs[column]):
            raise ValueError(f'--hierarchy {parent} {child}: {column} is not an exported obs column with categories')
    return obs[[parent, child]].dropna().astype(str).drop_duplicates()


def hierarchy_pairs_from_csv(path, obs, columns):
    """Distinct (parent, child) pairs from a two-column CSV, for --hierarchy-csv. The CSV's second column names the
    child obs column, which must have categories, and the CSV must cover all of its values in the data."""
    pairs = pd.read_csv(path, dtype=str)
    if pairs.shape[1] != 2:
        raise ValueError(f'{path}: a hierarchy CSV needs exactly two columns (parent, child), got {list(pairs.columns)}')
    child = pairs.columns[1]
    if child not in columns or is_numeric(obs[child]):
        raise ValueError(f'{path}: its child column {child} is not an exported obs column with categories')
    missing = sorted(set(obs[child].dropna().astype(str)) - set(pairs[child].dropna()))
    if missing:
        raise ValueError(f'{child} values missing from {path}: {missing[:10]}')
    return pairs.dropna().drop_duplicates()


def check_nesting(pairs, source):
    """Raise if a child value has more than one parent."""
    parent, child = pairs.columns
    n_parents = pairs.groupby(child)[parent].nunique()
    several = sorted(n_parents[n_parents > 1].index)
    if several:
        raise ValueError(f'{source}: {child} values with more than one {parent}: {several[:10]}')


def hierarchy_entry(pairs, adata, columns):
    """Manifest entry of a two-level hierarchy.

    'groups' lists the parent values in order of first use by the child's categories; 'colors' gives each group's
    color when the parent is an exported obs column with colors (else None); 'parent_of' gives each child category's
    group index. Child categories without a parent (missing values, values without cells) go in an 'NA' group.
    """
    parent, child = pairs.columns
    _, categories = category_codes(adata.obs[child])
    parent_by_child = dict(zip(pairs[child], pairs[parent]))
    names = [parent_by_child.get(c, MISSING_CATEGORY) for c in categories]
    # dict.fromkeys keeps the first occurrence of each name, in order
    groups = list(dict.fromkeys(names))
    if parent in columns:
        _, parent_categories = category_codes(adata.obs[parent])
        colors = dict(zip(parent_categories, category_colors(adata.obs, adata.uns, parent, parent_categories)))
    else:
        colors = {}
    return {'parent': parent, 'child': child, 'groups': groups, 'colors': [colors.get(g) for g in groups],
            'parent_of': [groups.index(n) for n in names]}


def build_hierarchies(adata, columns, hierarchies, hierarchy_csvs):
    """Manifest entries for --hierarchy PARENT CHILD pairs and --hierarchy-csv files, checked before anything is
    written: one hierarchy per child column, and every child value with exactly one parent."""
    sources = [(f'--hierarchy {p} {c}', hierarchy_pairs_from_obs(adata.obs, columns, p, c)) for p, c in hierarchies]
    sources += [(path, hierarchy_pairs_from_csv(path, adata.obs, columns)) for path in hierarchy_csvs]
    children = [pairs.columns[1] for _, pairs in sources]
    repeated = sorted({c for c in children if children.count(c) > 1})
    if repeated:
        raise ValueError(f'{repeated[0]} has more than one hierarchy; give one --hierarchy or --hierarchy-csv per column')
    for source, pairs in sources:
        check_nesting(pairs, source)
    return [hierarchy_entry(pairs, adata, columns) for _, pairs in sources]


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
            axis_labels=DEFAULT_AXIS_LABELS, up=None, exclude=(), ends=None, roles=None, view=None, locked=(),
            hierarchies=(), hierarchy_csvs=()):
    """Write the dataset folder and return its manifest.

    `name` is shown in the viewer (default: the file name without extension); `color` and `filters` are the obs
    columns colored by and given filter panels when the viewer opens; `embedding` is the obsm key shown first
    (default: 'spatial' if present); `axis_labels` names the three axes; `up` maps obsm keys to the direction
    drawn upward on screen (default +y); `exclude` lists obs columns to leave out.
    `ends` maps obsm keys to [low, high] end names per axis; `roles` maps obsm keys to the anatomical role of each
    axis (e.g. ['LR', 'PA', 'VD'], low end first) for embeddings whose end names are not anatomical; `view` is the
    end name of the side shown first (default: the low x end); `locked` lists two axis labels (drags then turn the
    data only about the third when the viewer opens) or none; `hierarchies` lists (parent, child) obs column pairs
    and `hierarchy_csvs` two-column CSV files that group a child column's values.
    Every request is checked before anything is written.
    """
    up, ends, roles = up or {}, ends or {}, roles or {}
    name = name or os.path.splitext(os.path.basename(h5ad_path))[0]

    # Check the request against the data
    adata = ad.read_h5ad(h5ad_path)
    check_exclude(adata.obs, exclude)
    columns = exported_columns(adata.obs, exclude)
    embeddings = embedding_keys(adata)
    check_columns(columns, color, filters)
    check_embeddings(embeddings, embedding, up)
    check_ends(embeddings, ends, up)
    check_roles(embeddings, roles, up, ends, axis_labels)
    check_view(view, embedding, embeddings, ends, axis_labels)
    check_locked(locked, axis_labels)
    check_out_dir(out_dir)
    hierarchy_entries = build_hierarchies(adata, columns, hierarchies, hierarchy_csvs)

    # Write the arrays, then the manifest that lists them
    manifest = {
        'format_version': FORMAT_VERSION,
        'name': name,
        'source': os.path.abspath(h5ad_path),
        'n_cells': int(adata.n_obs),
        'axes': {'labels': list(axis_labels)},
        'defaults': {'embedding': embedding or default_embedding(embeddings), 'color': color,
                     'filters': list(filters), 'view': view, 'locked': list(locked)},
        'embeddings': write_embeddings(adata, embeddings, up, ends, roles, axis_labels, out_dir),
        'fields': write_fields(adata, columns, out_dir),
        'gene_sources': [write_genes(adata.X, adata.var_names, MEASURED, MEASURED, out_dir)],
        'hierarchies': hierarchy_entries,
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
    # Names of both ends of each axis label the side buttons (e.g. Left shows the left side facing the viewer)
    parser.add_argument('--ends', nargs=7, action='append', default=[],
                        metavar=('OBSM_KEY', 'LOW_X', 'HIGH_X', 'LOW_Y', 'HIGH_Y', 'LOW_Z', 'HIGH_Z'),
                        help="names of each axis' low and high end for one embedding, e.g. spatial anterior posterior "
                             "dorsal ventral left right; an end named dorsal sets the up direction; repeatable")
    # Roles say which data axis is which anatomical axis when the end names do not (anatomical end names imply them)
    parser.add_argument('--roles', nargs=4, action='append', default=[], metavar=('OBSM_KEY', 'X', 'Y', 'Z'),
                        help='anatomical axis of x, y and z for one embedding, each named low end first from AP/PA, '
                             'DV/VD and LR/RL, e.g. spatial LR PA VD; sets the up direction (dorsal); repeatable')
    parser.add_argument('--view', help='end name of the side shown when the viewer opens, e.g. left '
                                       '(default: the low x end)')
    parser.add_argument('--locked', nargs='*', default=[],
                        help='two axis labels locked when the viewer opens, so drags turn the data only about the '
                             'third, e.g. AP LR (default: none, drags turn it like a turntable)')
    # A hierarchy groups a column's values (e.g. subclasses by class) in its filter panel and legend
    parser.add_argument('--hierarchy', nargs=2, action='append', default=[], metavar=('PARENT', 'CHILD'),
                        help='obs columns where each CHILD value belongs to one PARENT value, e.g. Class Subclass; '
                             'repeatable')
    parser.add_argument('--hierarchy-csv', action='append', default=[], metavar='CSV',
                        help='two-column CSV (parent, child) mapping each value of the child obs column to a group; '
                             'the parent need not be an obs column; repeatable')
    parser.add_argument('--exclude', nargs='*', default=[], help='obs columns to leave out')
    args = parser.parse_args()

    up = {key: [float(x), float(y), float(z)] for key, x, y, z in args.up}
    # --ends gives six names per embedding: low and high end of x, then y, then z
    ends = {key: [names[0:2], names[2:4], names[4:6]] for key, *names in args.ends}
    roles = {key: [role.upper() for role in given] for key, *given in args.roles}
    manifest = convert(args.h5ad, args.out_dir, name=args.name, color=args.color, filters=args.filters,
                       embedding=args.embedding, axis_labels=args.axis_labels, up=up, exclude=args.exclude,
                       ends=ends, roles=roles, view=args.view, locked=args.locked,
                       hierarchies=[tuple(pair) for pair in args.hierarchy], hierarchy_csvs=args.hierarchy_csv)
    genes = ', '.join(f"{len(s['genes'])} {s['name']}" for s in manifest['gene_sources'])
    print(f"Wrote {args.out_dir}: {manifest['n_cells']:,} cells, {len(manifest['fields'])} obs columns, "
          f"genes: {genes}", flush=True)


if __name__ == '__main__':
    main()
