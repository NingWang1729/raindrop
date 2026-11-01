"""Tests for convert.py: build a tiny AnnData, convert it, and read the files back."""
import json
import os
import subprocess
import sys

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from convert import convert, embedding_axis_labels

CONVERT_PY = os.path.join(os.path.dirname(__file__), '..', 'convert.py')


def make_adata():
    """Five cells, two genes (and a third, too large for uint16, in X), a categorical column with colors and a
    missing value, a string column, a number with NaN and a bool; 3D, 2D and 10D embeddings."""
    obs = pd.DataFrame({
        'Subclass': pd.Categorical(['a', 'b', 'a', None, 'c']),
        'Subclass_color': ['#ff0000', '#00ff00', '#ff0000', None, '#0000ff'],
        'section': ['s1', 's1', 's2', 's2', 's2'],
        'depth': [0.5, 1.5, np.nan, 2.5, 3.5],
        'passed': [True, False, True, True, False],
    }, index=[f'cell{i}' for i in range(5)])
    X = sparse.csr_matrix(np.array([[0, 1, 2], [3, 0, 0], [0, 0, 70000], [1, 1, 1], [0, 5, 0]], dtype=np.float32))
    adata = ad.AnnData(X=X[:, :2], obs=obs, var=pd.DataFrame(index=['g1', 'g2']))
    adata.obsm['spatial'] = np.arange(15, dtype=np.float64).reshape(5, 3)
    adata.obsm['X_umap'] = np.arange(10, dtype=np.float64).reshape(5, 2)
    adata.obsm['X_pca'] = np.zeros((5, 10))
    return adata, X


def read(out, entry, dtype):
    """Read the array file of a manifest entry."""
    return np.fromfile(os.path.join(out, entry['file']), dtype=np.dtype(dtype).newbyteorder('<'))


def run(tmp_path, adata, **kwargs):
    """Convert `adata` with test defaults (color Subclass, a section filter, spatial drawn with -y up)."""
    path = tmp_path / 'in.h5ad'
    adata.write_h5ad(path)
    args = dict(name='test', color='Subclass', filters=['section'], up={'spatial': [0, -1, 0]})
    args.update(kwargs)
    out = str(tmp_path / 'out')
    return out, convert(str(path), out, **args)


def test_manifest_and_arrays(tmp_path):
    adata, _ = make_adata()
    out, manifest = run(tmp_path, adata)
    assert json.load(open(f'{out}/manifest.json')) == manifest
    assert manifest['n_cells'] == 5
    assert manifest['defaults'] == {'embedding': 'spatial', 'color': 'Subclass', 'filters': ['section'], 'view': None,
                                    'locked': []}
    assert manifest['hierarchies'] == []

    # Only 2D/3D embeddings are kept, with values intact
    assert [(e['name'], e['dims']) for e in manifest['embeddings']] == [('X_umap', 2), ('spatial', 3)]
    np.testing.assert_array_equal(read(out, manifest['embeddings'][1], 'float32').reshape(5, 3), adata.obsm['spatial'])

    # Each embedding has its own up direction (unlisted ones default to +y) and generic axis end names; an X_ key
    # (UMAP, PCA, ...) names its axes after itself, others take the dataset's axis labels
    assert [e['up'] for e in manifest['embeddings']] == [[0, 1, 0], [0, -1, 0]]
    umap, spatial = manifest['embeddings']
    assert umap['labels'] == ['UMAP 1', 'UMAP 2', 'UMAP 3']
    assert umap['ends'] == [['-UMAP 1', '+UMAP 1'], ['-UMAP 2', '+UMAP 2'], ['-UMAP 3', '+UMAP 3']]
    assert spatial['labels'] == ['x', 'y', 'z'] and spatial['ends'] == [['-x', '+x'], ['-y', '+y'], ['-z', '+z']]

    # The color column is not a field of its own; it colors Subclass
    fields = {f['name']: f for f in manifest['fields']}
    assert list(fields) == ['Subclass', 'section', 'depth', 'passed']
    subclass = fields['Subclass']
    assert subclass['categories'] == ['a', 'b', 'c', 'NA']
    assert subclass['colors'] == ['#ff0000', '#00ff00', '#0000ff', None]
    np.testing.assert_array_equal(read(out, subclass, subclass['dtype']), [0, 1, 0, 3, 2])

    # Strings and bools become categories; numbers keep NaN and report their finite range
    np.testing.assert_array_equal(read(out, fields['section'], 'uint8'), [0, 0, 1, 1, 1])
    assert fields['section']['colors'] == [None, None]
    assert fields['passed']['categories'] == ['False', 'True']
    assert fields['depth']['kind'] == 'numeric' and [fields['depth']['min'], fields['depth']['max']] == [0.5, 3.5]
    np.testing.assert_array_equal(read(out, fields['depth'], 'float32'), adata.obs['depth'].to_numpy(np.float32))

    # Integer counts are stored as dense uint16, one file per gene
    source = manifest['gene_sources'][0]
    assert source['name'] == 'measured' and source['genes'] == ['g1', 'g2'] and source['dtype'] == 'uint16'
    for i in range(2):
        gene = np.fromfile(f"{out}/{source['path'].format(index=i)}", dtype='<u2')
        np.testing.assert_array_equal(gene, adata.X[:, i].toarray().ravel())


def test_float_genes_and_large_counts(tmp_path):
    adata, X = make_adata()

    # Values above 65535 do not fit uint16, so the genes are stored as float32
    big = ad.AnnData(X=X, obs=adata.obs, var=pd.DataFrame(index=['g1', 'g2', 'g3']), obsm=adata.obsm)
    out, manifest = run(tmp_path, big)
    assert manifest['gene_sources'][0]['dtype'] == 'float32'
    np.testing.assert_array_equal(np.fromfile(f'{out}/genes/measured/2.bin', dtype='<f4'), [2, 0, 70000, 1, 0])


def test_chosen_embedding_and_axis_labels(tmp_path):
    adata, _ = make_adata()
    out, manifest = run(tmp_path, adata, embedding='X_umap', axis_labels=['AP', 'DV', 'LR'])
    assert manifest['defaults']['embedding'] == 'X_umap'
    assert manifest['axes']['labels'] == ['AP', 'DV', 'LR']
    spatial = manifest['embeddings'][1]
    assert spatial['labels'] == ['AP', 'DV', 'LR'] and spatial['ends'][0] == ['-AP', '+AP']

    # Without a 'spatial' embedding, the first 2D/3D obsm key opens
    del adata.obsm['spatial']
    (tmp_path / 'no_spatial').mkdir()
    _, manifest = run(tmp_path / 'no_spatial', adata, up={})
    assert manifest['defaults']['embedding'] == 'X_umap'


def test_dense_genes_and_wide_categories(tmp_path):
    adata, _ = make_adata()
    # Dense X with whole counts is stored as uint16 like sparse X
    adata.X = adata.X.toarray()
    # 300 categories need uint16 codes
    names = [f'c{i}' for i in range(300)]
    adata.obs['cell'] = pd.Categorical(names[:5], categories=names)
    out, manifest = run(tmp_path, adata)
    source = manifest['gene_sources'][0]
    assert source['dtype'] == 'uint16'
    np.testing.assert_array_equal(np.fromfile(f'{out}/genes/measured/1.bin', dtype='<u2'), [1, 0, 0, 1, 5])
    cell = next(f for f in manifest['fields'] if f['name'] == 'cell')
    assert cell['dtype'] == 'uint16'
    np.testing.assert_array_equal(read(out, cell, 'uint16'), [0, 1, 2, 3, 4])


def test_scanpy_colors(tmp_path):
    adata, _ = make_adata()
    # uns['<column>_colors'] wins over a color column when it has one color per category
    adata.uns['Subclass_colors'] = ['#111111', '#222222', '#333333']
    adata.uns['section_colors'] = ['#444444', '#555555']
    _, manifest = run(tmp_path, adata)
    fields = {f['name']: f for f in manifest['fields']}
    assert fields['section']['colors'] == ['#444444', '#555555']
    # Subclass has a fourth category (NA), so its three uns colors do not fit: the color column is used instead
    assert fields['Subclass']['colors'] == ['#ff0000', '#00ff00', '#0000ff', None]


def test_anatomical_ends(tmp_path):
    adata, _ = make_adata()
    frame_a = [['anterior', 'posterior'], ['dorsal', 'ventral'], ['left', 'right']]

    # The end named dorsal sets up (-y here); --view picks the first side shown
    _, manifest = run(tmp_path, adata, axis_labels=['AP', 'DV', 'LR'], up={}, ends={'spatial': frame_a}, view='left',
                      locked=['AP', 'LR'])
    spatial = manifest['embeddings'][1]
    assert spatial['ends'] == frame_a and spatial['up'] == [0, -1, 0] and spatial['roles'] == ['AP', 'DV', 'LR']
    assert manifest['defaults']['view'] == 'left' and manifest['defaults']['locked'] == ['AP', 'LR']

    # Dorsal as the high end of x points up along +x; with only one anatomical axis there are no roles
    sideways = [['ventral', 'dorsal'], ['a', 'b'], ['c', 'd']]
    os.makedirs(tmp_path / 'other')
    _, manifest = run(tmp_path / 'other', adata, up={}, ends={'spatial': sideways})
    assert manifest['embeddings'][1]['up'] == [1, 0, 0] and manifest['embeddings'][1]['roles'] is None


def test_axis_roles(tmp_path):
    adata, _ = make_adata()
    for folder in ['none', 'ends', 'given', 'both']:
        os.makedirs(tmp_path / folder)
    # Without anatomical ends or --roles an embedding has no roles
    _, manifest = run(tmp_path / 'none', adata)
    assert [e['roles'] for e in manifest['embeddings']] == [None, None]

    # Roles come from anatomical end names in any axis order and case (e.g. LR, PA, VD), each
    # naming its axis's low end first; up points to dorsal
    frame_b = [['Left', 'Right'], ['posterior', 'anterior'], ['ventral', 'dorsal']]
    _, manifest = run(tmp_path / 'ends', adata, up={}, ends={'spatial': frame_b})
    assert manifest['embeddings'][1]['roles'] == ['LR', 'PA', 'VD'] and manifest['embeddings'][1]['up'] == [0, 0, 1]

    # --roles gives them for an embedding whose end names are not anatomical; up then points to dorsal (-y)
    _, manifest = run(tmp_path / 'given', adata, up={}, roles={'spatial': ['AP', 'DV', 'LR']})
    spatial = manifest['embeddings'][1]
    assert spatial['roles'] == ['AP', 'DV', 'LR'] and spatial['up'] == [0, -1, 0]
    assert spatial['ends'] == [['-x', '+x'], ['-y', '+y'], ['-z', '+z']]

    # --roles that agree with anatomical end names are accepted
    _, manifest = run(tmp_path / 'both', adata, up={}, ends={'spatial': frame_b}, roles={'spatial': ['LR', 'PA', 'VD']})
    assert manifest['embeddings'][1]['roles'] == ['LR', 'PA', 'VD']


def test_axis_role_errors(tmp_path):
    adata, _ = make_adata()
    frame_b = [['left', 'right'], ['posterior', 'anterior'], ['ventral', 'dorsal']]
    with pytest.raises(ValueError, match=r"--roles spatial: \['AX'\] are not axis roles"):
        run(tmp_path, adata, up={}, roles={'spatial': ['AX', 'DV', 'LR']})
    with pytest.raises(ValueError, match=r"--roles spatial names the AP axis more than once; give each of AP, DV and LR"):
        run(tmp_path, adata, up={}, roles={'spatial': ['AP', 'PA', 'LR']})
    with pytest.raises(ValueError, match=r"--roles spatial: axis 2 is DV but its --ends are \['ventral', 'dorsal'\]"):
        run(tmp_path, adata, up={}, ends={'spatial': frame_b}, roles={'spatial': ['LR', 'PA', 'DV']})
    with pytest.raises(ValueError, match='--roles names obsm entries'):
        run(tmp_path, adata, up={}, roles={'X_pca': ['AP', 'DV', 'LR']})
    with pytest.raises(ValueError, match='--roles already sets the up direction'):
        run(tmp_path, adata, roles={'spatial': ['AP', 'DV', 'LR']})

    # End names must name both ends of one anatomical axis, each axis once
    with pytest.raises(ValueError, match=r"--ends spatial: \['anterior', 'dorsal'\] are not the two ends of one anatomical axis"):
        run(tmp_path, adata, up={}, ends={'spatial': [['anterior', 'dorsal'], ['a', 'b'], ['c', 'd']]})
    with pytest.raises(ValueError, match=r"--ends spatial: \['left', 'b'\] are not the two ends of one anatomical axis"):
        run(tmp_path, adata, up={}, ends={'spatial': [['left', 'b'], ['c', 'd'], ['e', 'f']]})
    with pytest.raises(ValueError, match='--ends spatial names the LR axis more than once'):
        run(tmp_path, adata, up={}, ends={'spatial': [['left', 'right'], ['right', 'left'], ['ventral', 'dorsal']]})


def test_hierarchy_from_obs(tmp_path):
    adata, _ = make_adata()
    # Subclasses a, b, a, NA, c: a and b belong to class X, c to Y; the missing subclass goes under NA
    adata.obs['Class'] = pd.Categorical(['X', 'X', 'X', None, 'Y'])
    adata.obs['Class_color'] = ['#111111', '#111111', '#111111', None, '#222222']
    _, manifest = run(tmp_path, adata, hierarchies=[('Class', 'Subclass')])
    assert manifest['hierarchies'] == [{'parent': 'Class', 'child': 'Subclass', 'groups': ['X', 'Y', 'NA'],
                                        'colors': ['#111111', '#222222', None], 'parent_of': [0, 0, 1, 2]}]


def test_hierarchy_from_csv(tmp_path):
    adata, _ = make_adata()
    # The parent level need not be an obs column; extra rows are ignored
    pd.DataFrame({'Group': ['G1', 'G1', 'G2', 'G3'], 'Subclass': ['a', 'b', 'c', 'z']}).to_csv(tmp_path / 'map.csv', index=False)
    _, manifest = run(tmp_path, adata, hierarchy_csvs=[str(tmp_path / 'map.csv')])
    assert manifest['hierarchies'] == [{'parent': 'Group', 'child': 'Subclass', 'groups': ['G1', 'G2', 'NA'],
                                        'colors': [None, None, None], 'parent_of': [0, 0, 1, 2]}]


def test_hierarchy_errors(tmp_path):
    adata, _ = make_adata()
    adata.obs['Class'] = pd.Categorical(['X', 'X', 'Y', None, 'Y'])
    with pytest.raises(ValueError, match=r"Subclass values with more than one Class: \['a'\]"):
        run(tmp_path, adata, hierarchies=[('Class', 'Subclass')])
    with pytest.raises(ValueError, match="--hierarchy Class Foo: Foo is not an exported obs column"):
        run(tmp_path, adata, hierarchies=[('Class', 'Foo')])
    # A numeric column has no categories to group by
    with pytest.raises(ValueError, match="--hierarchy depth Subclass: depth is not an exported obs column with categories"):
        run(tmp_path, adata, hierarchies=[('depth', 'Subclass')])

    # A CSV must cover every value in the data, give each one parent, and have exactly two columns
    csv = tmp_path / 'map.csv'
    pd.DataFrame({'Group': ['G1', 'G1'], 'Subclass': ['a', 'b']}).to_csv(csv, index=False)
    with pytest.raises(ValueError, match=r"Subclass values missing from .*map.csv: \['c'\]"):
        run(tmp_path, adata, hierarchy_csvs=[str(csv)])
    pd.DataFrame({'Group': ['G1', 'G2', 'G1', 'G1'], 'Subclass': ['a', 'a', 'b', 'c']}).to_csv(csv, index=False)
    with pytest.raises(ValueError, match=r"Subclass values with more than one Group: \['a'\]"):
        run(tmp_path, adata, hierarchy_csvs=[str(csv)])
    pd.DataFrame({'A': ['x'], 'B': ['y'], 'Subclass': ['a']}).to_csv(csv, index=False)
    with pytest.raises(ValueError, match='needs exactly two columns'):
        run(tmp_path, adata, hierarchy_csvs=[str(csv)])
    pd.DataFrame({'Group': ['G1'], 'depth': ['1.5']}).to_csv(csv, index=False)
    with pytest.raises(ValueError, match='its child column depth is not an exported obs column with categories'):
        run(tmp_path, adata, hierarchy_csvs=[str(csv)])

    # One hierarchy per column
    pd.DataFrame({'Group': ['G1', 'G1', 'G2'], 'Subclass': ['a', 'b', 'c']}).to_csv(csv, index=False)
    with pytest.raises(ValueError, match='Subclass has more than one hierarchy'):
        run(tmp_path, adata, hierarchies=[('section', 'Subclass')], hierarchy_csvs=[str(csv)])


def test_embedding_axis_labels():
    labels = ['AP', 'DV', 'LR']
    # Known methods get their usual names; other X_ keys their name without X_; others the dataset's labels
    assert embedding_axis_labels('X_umap', labels, {}, {}) == ['UMAP 1', 'UMAP 2', 'UMAP 3']
    assert embedding_axis_labels('X_tsne', labels, {}, {}) == ['t-SNE 1', 't-SNE 2', 't-SNE 3']
    assert embedding_axis_labels('X_pca', labels, {}, {}) == ['PC 1', 'PC 2', 'PC 3']
    assert embedding_axis_labels('X_piano_umap', labels, {}, {}) == ['piano_umap 1', 'piano_umap 2', 'piano_umap 3']
    assert embedding_axis_labels('spatial', labels, {}, {}) == labels
    # An X_ key with anatomical end names or roles is spatial: it keeps the dataset's labels
    ends = {'X_ccf': [['anterior', 'posterior'], ['dorsal', 'ventral'], ['left', 'right']]}
    assert embedding_axis_labels('X_ccf', labels, ends, {}) == labels
    assert embedding_axis_labels('X_other', labels, {}, {'X_other': ['LR', 'PA', 'VD']}) == labels


def test_conflicting_colors_warn(tmp_path):
    adata, _ = make_adata()

    # Two colors for category 'a': warn and leave the colors to the viewer
    adata.obs['Subclass_color'] = ['#ff0000', '#00ff00', '#0000ff', None, '#0000ff']
    with pytest.warns(UserWarning, match='several colors for one category'):
        _, manifest = run(tmp_path, adata)
    assert manifest['fields'][0]['colors'] == [None] * 4


def test_too_many_categories(tmp_path):
    adata, _ = make_adata()
    # Unused categories count too, so five cells can carry 65540 categories
    names = [f'c{i}' for i in range(65540)]
    adata.obs['cell'] = pd.Categorical(names[:5], categories=names)
    with pytest.raises(ValueError, match='obs column cell has 65540 categories'):
        run(tmp_path, adata)


def test_errors(tmp_path):
    adata, _ = make_adata()
    with pytest.raises(ValueError, match='not an exported obs column'):
        run(tmp_path, adata, color='Subclass_color')
    with pytest.raises(ValueError, match='not a 2D or 3D obsm entry'):
        run(tmp_path, adata, embedding='X_pca')
    with pytest.raises(ValueError, match='--up names obsm entries'):
        run(tmp_path, adata, up={'X_pca': [0, 1, 0]})
    frame_a = [['anterior', 'posterior'], ['dorsal', 'ventral'], ['left', 'right']]
    with pytest.raises(ValueError, match='--ends names obsm entries'):
        run(tmp_path, adata, ends={'X_pca': frame_a})
    with pytest.raises(ValueError, match='already sets the up direction; drop --up'):
        run(tmp_path, adata, ends={'spatial': frame_a})
    with pytest.raises(ValueError, match='--view top is not an end of embedding spatial'):
        run(tmp_path, adata, up={}, ends={'spatial': frame_a}, view='top')
    with pytest.raises(ValueError, match=r"--locked AP is not an axis label; choose from \['x', 'y', 'z'\]"):
        run(tmp_path, adata, locked=['AP'])
    # Drags turn about any axis or about the one left free: one or three locked axes have no meaning
    for locked in (['x'], ['x', 'y', 'z']):
        with pytest.raises(ValueError, match=r'lock two axes \(drags then turn about the third\) or none'):
            run(tmp_path, adata, locked=locked)
    with pytest.raises(ValueError, match='--exclude columns not in obs'):
        run(tmp_path, adata, exclude=['missing'])
    with pytest.raises(ValueError, match='no 2D or 3D embedding'):
        run(tmp_path, ad.AnnData(X=adata.X, obs=adata.obs, var=adata.var))

    # Nothing is written when a check fails, and an existing output folder is never overwritten
    assert not os.path.exists(tmp_path / 'out')
    os.makedirs(tmp_path / 'out')
    with pytest.raises(FileExistsError):
        run(tmp_path, adata)


def test_command_line(tmp_path):
    """End to end: the script as a user runs it, with the default name taken from the file name."""
    adata, _ = make_adata()
    adata.write_h5ad(tmp_path / 'example.h5ad')
    result = subprocess.run([sys.executable, CONVERT_PY, str(tmp_path / 'example.h5ad'), str(tmp_path / 'out'),
                             '--color', 'Subclass', '--up', 'spatial', '0', '-1', '0', '--exclude', 'passed'],
                            capture_output=True, text=True, check=True)
    assert 'Wrote' in result.stdout and '5 cells' in result.stdout
    manifest = json.load(open(tmp_path / 'out' / 'manifest.json'))
    assert manifest['name'] == 'example'
    assert [f['name'] for f in manifest['fields']] == ['Subclass', 'section', 'depth']
    assert manifest['embeddings'][1]['up'] == [0.0, -1.0, 0.0]


def test_command_line_anatomy(tmp_path):
    """End to end: --ends, --roles (any case), --view, --locked and --hierarchy as a user types them."""
    adata, _ = make_adata()
    adata.obsm['X_ccf'] = adata.obsm['spatial'].copy()
    # Subclasses a and b belong to class X, c to Y
    adata.obs['Class'] = pd.Categorical(['X', 'X', 'X', None, 'Y'])
    adata.write_h5ad(tmp_path / 'example.h5ad')
    subprocess.run([sys.executable, CONVERT_PY, str(tmp_path / 'example.h5ad'), str(tmp_path / 'out'),
                    '--axis-labels', 'AP', 'DV', 'LR',
                    '--ends', 'spatial', 'anterior', 'posterior', 'dorsal', 'ventral', 'left', 'right',
                    '--roles', 'X_ccf', 'ap', 'dv', 'lr', '--view', 'left', '--locked', 'AP', 'LR',
                    '--hierarchy', 'Class', 'Subclass'], capture_output=True, text=True, check=True)
    manifest = json.load(open(tmp_path / 'out' / 'manifest.json'))
    embeddings = {e['name']: e for e in manifest['embeddings']}
    assert embeddings['spatial']['ends'][2] == ['left', 'right'] and embeddings['spatial']['up'] == [0, -1, 0]
    assert embeddings['X_ccf']['roles'] == ['AP', 'DV', 'LR'] and embeddings['X_ccf']['labels'] == ['AP', 'DV', 'LR']
    assert manifest['defaults']['view'] == 'left' and manifest['defaults']['locked'] == ['AP', 'LR']
    assert manifest['hierarchies'][0]['groups'] == ['X', 'Y', 'NA']
