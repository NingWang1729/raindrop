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

from convert import convert

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
    assert manifest['defaults'] == {'embedding': 'spatial', 'color': 'Subclass', 'filters': ['section']}

    # Only 2D/3D embeddings are kept, with values intact
    assert [(e['name'], e['dims']) for e in manifest['embeddings']] == [('X_umap', 2), ('spatial', 3)]
    np.testing.assert_array_equal(read(out, manifest['embeddings'][1], 'float32').reshape(5, 3), adata.obsm['spatial'])

    # Each embedding has its own up direction (unlisted ones default to +y), the dataset's axis labels and generic
    # axis end names
    assert [e['up'] for e in manifest['embeddings']] == [[0, 1, 0], [0, -1, 0]]
    for entry in manifest['embeddings']:
        assert entry['labels'] == ['x', 'y', 'z'] and entry['ends'] == [['-x', '+x'], ['-y', '+y'], ['-z', '+z']]

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
    for entry in manifest['embeddings']:
        assert entry['labels'] == ['AP', 'DV', 'LR'] and entry['ends'][0] == ['-AP', '+AP']

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
