"""
Tests of how raindrop is packaged: the converter lives in the package, the clone's convert.py only runs it, the
wheel installs no top-level convert module, and the package's PyPI name and page are right.

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
import configparser
import email.parser
import json
import os
import re
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

from test_convert import CONVERT_PY, make_adata

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
# What a wheel could be built from (the root convert.py too, to show it is left out); anything else in a clone
# (datasets, environments, build output) is not copied
SOURCE_FILES = ['pyproject.toml', 'README.md', 'LICENSE', 'convert.py']
SOURCE_FOLDERS = ['raindrop', 'viewer']
# Link targets in Markdown: inline links and images ([text](target "title"), with or without <>), reference links
# ([name]: target) and HTML links (href="target")
LINK_TARGETS = [re.compile(r'\]\(\s*<?([^)\s>]+)'), re.compile(r'^\s*\[[^\]]+\]:\s*<?([^\s>]+)', re.MULTILINE),
                re.compile(r'href="([^"]+)"')]


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
    subprocess.run([sys.executable, CONVERT_PY, 'example.h5ad', 'out'], cwd=tmp_path, env=env, capture_output=True,
                   check=True)
    assert json.load(open(tmp_path / 'out' / 'manifest.json'))['n_cells'] == 5


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
    names = zipfile.ZipFile(build_wheel(tmp_path)).namelist()
    assert 'convert.py' not in names
    for name in ('raindrop/convert.py', 'raindrop/cli.py', 'raindrop/viewer/index.html',
                 'raindrop/viewer/vendor/THREE_LICENSE'):
        assert name in names, name
    assert any(name.endswith('.dist-info/licenses/LICENSE') for name in names)


def test_wheel_metadata_names_the_pypi_package(tmp_path):
    """The wheel installs as raindrop-viewer (raindrop is another project's name on PyPI), with the raindrop command
    and the README as its PyPI page."""
    wheel = zipfile.ZipFile(build_wheel(tmp_path))
    dist_info = next(name for name in wheel.namelist() if name.endswith('.dist-info/METADATA')).rsplit('/', 1)[0]
    metadata = email.parser.Parser().parsestr(wheel.read(f'{dist_info}/METADATA').decode())
    # Names are compared as pip does (PEP 503: case, '-', '_' and '.' alike)
    assert re.sub(r'[-_.]+', '-', metadata['Name']).lower() == 'raindrop-viewer'
    assert metadata['Description-Content-Type'] == 'text/markdown'
    assert '# Raindrop' in metadata.get_payload()
    assert not [c for c in metadata.get_all('Classifier', []) if c.startswith('License ::')]
    entry_points = configparser.ConfigParser()
    entry_points.read_string(wheel.read(f'{dist_info}/entry_points.txt').decode())
    assert entry_points['console_scripts']['raindrop'] == 'raindrop.cli:main'


def test_readme_links_work_on_pypi():
    """PyPI shows the README without the repository's files and without heading anchors, so every link in it is a
    full address."""
    readme = Path(ROOT, 'README.md').read_text(encoding='utf-8')
    targets = [target for pattern in LINK_TARGETS for target in pattern.findall(readme)]
    assert targets
    assert [target for target in targets if not target.startswith('https://')] == []
