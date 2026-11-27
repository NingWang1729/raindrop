"""
Tests of how raindrop is packaged and released.

raindrop: a browser viewer for single-cell and spatial data
Copyright (C) 2026 Ning Wang

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
import json
import os
import shutil
import subprocess
import sys
import zipfile

from test_convert import CONVERT_PY, make_adata

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
# What a wheel could be built from (the root convert.py too, to show it is left out); anything else in a clone
# (datasets, environments, build output) is not copied
SOURCE_FILES = ['pyproject.toml', 'README.md', 'LICENSE', 'convert.py']
SOURCE_FOLDERS = ['raindrop', 'viewer']


def test_converter_is_in_the_package():
    import raindrop.convert

    for name in ('main', 'convert', 'add_arguments'):
        assert callable(getattr(raindrop.convert, name)), name


def test_root_convert_py_runs_the_clones_converter_from_any_folder(tmp_path):
    """python convert.py, run from another folder, uses the clone's raindrop package ahead of another one on the path
    (here a decoy whose converter fails), and shows the packaged converter's help under the name convert.py."""
    import raindrop.convert

    decoy = tmp_path / 'decoy' / 'raindrop'
    decoy.mkdir(parents=True)
    (decoy / '__init__.py').write_text('')
    (decoy / 'convert.py').write_text("raise ImportError('the decoy raindrop was imported')")
    # The decoy goes ahead of any search path the environment needs (e.g. to find anndata)
    search_path = [str(tmp_path / 'decoy'), os.environ.get('PYTHONPATH')]
    env = dict(os.environ, PYTHONPATH=os.pathsep.join(filter(None, search_path)))

    script = subprocess.run([sys.executable, CONVERT_PY, '--help'], cwd=tmp_path, env=env, capture_output=True,
                            text=True)
    assert script.returncode == 0, script.stderr
    assert script.stdout.startswith('usage: convert.py')
    assert raindrop.convert.__doc__ in script.stdout

    # A conversion with relative paths, from that folder
    adata, _ = make_adata()
    adata.write_h5ad(tmp_path / 'example.h5ad')
    converted = subprocess.run([sys.executable, CONVERT_PY, 'example.h5ad', 'out'], cwd=tmp_path, env=env,
                               capture_output=True, text=True)
    assert converted.returncode == 0, converted.stderr
    with open(tmp_path / 'out' / 'manifest.json') as manifest:
        assert json.load(manifest)['n_cells'] == 5


def build_wheel(tmp_path):
    """Build a wheel from a copy of the source (so the clone gets no build output); return the wheel's path."""
    source = tmp_path / 'source'
    source.mkdir()
    for name in SOURCE_FILES:
        shutil.copy(os.path.join(ROOT, name), source / name)
    for name in SOURCE_FOLDERS:
        shutil.copytree(os.path.join(ROOT, name), source / name, ignore=shutil.ignore_patterns('__pycache__'))
    # setuptools' own build hook, without pip or an isolated environment, so the test needs no network
    built = subprocess.run([sys.executable, '-c', 'import sys; from setuptools import build_meta; '
                            'print(build_meta.build_wheel(sys.argv[1]))', str(tmp_path / 'wheel')],
                           cwd=source, capture_output=True, text=True)
    assert built.returncode == 0, built.stderr
    return tmp_path / 'wheel' / built.stdout.strip().splitlines()[-1]


def test_wheel_has_the_package_and_no_top_level_convert(tmp_path):
    with zipfile.ZipFile(build_wheel(tmp_path)) as wheel:
        names = wheel.namelist()
    assert 'convert.py' not in names
    for name in ('raindrop/convert.py', 'raindrop/cli.py', 'raindrop/viewer/index.html',
                 'raindrop/viewer/vendor/THREE_LICENSE'):
        assert name in names, name
    assert any(name.endswith('.dist-info/licenses/LICENSE') for name in names)
