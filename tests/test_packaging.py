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
import configparser
import email.parser
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import zipfile
from pathlib import Path

import raindrop
from test_convert import CONVERT_PY, make_adata

# The repository on GitHub
REPO_URL = 'https://github.com/NingWang1729/raindrop'
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
# What a wheel could be built from (the root convert.py too, to show it is left out); anything else in a clone
# (datasets, environments, build output) is not copied
SOURCE_FILES = ['pyproject.toml', 'README.md', 'LICENSE', 'convert.py']
SOURCE_FOLDERS = ['raindrop', 'viewer']
# What an sdist could be built from besides that: MANIFEST.in and the files it may add
SDIST_FILES = ['MANIFEST.in', 'CHANGELOG.md', 'CITATION.cff']
SDIST_FOLDERS = ['tests']
# Link targets in Markdown: inline links and images ([text](target "title"), with or without <>), reference links
# ([name]: target) and HTML links (href="target")
LINK_TARGETS = [re.compile(r'\]\(\s*<?([^)\s>]+)'), re.compile(r'^\s*\[[^\]]+\]:\s*<?([^\s>]+)', re.MULTILINE),
                re.compile(r'href="([^"]+)"')]
# A release's heading in CHANGELOG.md (## [0.1.1] - 2026-11-27; not [Unreleased])
CHANGELOG_RELEASE = re.compile(r'^## \[(\d[^\]]*)\]', re.MULTILINE)
# The version line of CITATION.cff
CITATION_VERSION = re.compile(r'^version: (\S+)$', re.MULTILINE)
# The top-level doi line of CITATION.cff, and the values in its identifiers list (indented)
CITATION_DOI = re.compile(r'^doi: (\S+)$', re.MULTILINE)
CITATION_IDENTIFIER = re.compile(r'^ +value: (\S+)$', re.MULTILINE)
# A Zenodo DOI (10.5281 is Zenodo's DOI prefix, then zenodo.<record number>)
ZENODO_DOI = re.compile(r'10\.5281/zenodo\.\d+')


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


def copy_source(source, files, folders):
    """Copy those files and folders of the clone that exist into the folder source (missing ones are left for the
    test's checks to notice)."""
    source.mkdir()
    for name in files:
        if os.path.exists(os.path.join(ROOT, name)):
            shutil.copy(os.path.join(ROOT, name), source / name)
    for name in folders:
        shutil.copytree(os.path.join(ROOT, name), source / name, ignore=shutil.ignore_patterns('__pycache__'))


def build(tmp_path, kind, files, folders):
    """Build a wheel or an sdist (kind) from a copy of the source (so the clone gets no build output); return its
    path."""
    source = tmp_path / f'source_{kind}'
    copy_source(source, files, folders)
    # setuptools' own build hook, without pip or an isolated environment, so the test needs no network
    built = subprocess.run([sys.executable, '-c', 'import sys; from setuptools import build_meta; '
                            f'print(build_meta.build_{kind}(sys.argv[1]))', str(tmp_path / kind)],
                           cwd=source, capture_output=True, text=True)
    assert built.returncode == 0, built.stderr
    return tmp_path / kind / built.stdout.strip().splitlines()[-1]


def build_wheel(tmp_path):
    return build(tmp_path, 'wheel', SOURCE_FILES, SOURCE_FOLDERS)


def build_sdist(tmp_path):
    return build(tmp_path, 'sdist', SOURCE_FILES + SDIST_FILES, SOURCE_FOLDERS + SDIST_FOLDERS)


def test_wheel_has_the_package_and_no_top_level_convert(tmp_path):
    with zipfile.ZipFile(build_wheel(tmp_path)) as wheel:
        names = wheel.namelist()
    assert 'convert.py' not in names
    for name in ('raindrop/convert.py', 'raindrop/cli.py', 'raindrop/viewer/index.html',
                 'raindrop/viewer/vendor/THREE_LICENSE'):
        assert name in names, name
    assert any(name.endswith('.dist-info/licenses/LICENSE') for name in names)


def dist_info(wheel):
    """The .dist-info folder of a wheel (a ZipFile)."""
    return next(name for name in wheel.namelist() if name.endswith('.dist-info/METADATA')).rsplit('/', 1)[0]


def read_metadata(wheel):
    """The METADATA of a wheel (a ZipFile) as an email message."""
    return email.parser.Parser().parsestr(wheel.read(f'{dist_info(wheel)}/METADATA').decode())


def test_wheel_installs_as_raindrop_viewer(tmp_path):
    """pip installs the wheel as raindrop-viewer (raindrop is another project's name on PyPI), with the raindrop
    command."""
    with zipfile.ZipFile(build_wheel(tmp_path)) as wheel:
        metadata = read_metadata(wheel)
        entry_points = configparser.ConfigParser()
        entry_points.read_string(wheel.read(f'{dist_info(wheel)}/entry_points.txt').decode())
    # Names are compared as pip does (PEP 503: case, '-', '_' and '.' alike)
    assert re.sub(r'[-_.]+', '-', metadata['Name']).lower() == 'raindrop-viewer'
    assert metadata['Version'] == raindrop.__version__
    assert entry_points['console_scripts']['raindrop'] == 'raindrop.cli:main'


def test_wheel_metadata_fills_the_pypi_page(tmp_path):
    """The README is the PyPI page's description, next to the author, keywords, classifiers and project links."""
    with zipfile.ZipFile(build_wheel(tmp_path)) as wheel:
        metadata = read_metadata(wheel)
    assert metadata['Description-Content-Type'] == 'text/markdown'
    assert '# Raindrop' in metadata.get_payload()
    assert metadata['Author'] == 'Ning Wang'
    assert 'single-cell' in metadata['Keywords']
    classifiers = metadata.get_all('Classifier', [])
    assert 'Intended Audience :: Science/Research' in classifiers
    assert not [c for c in classifiers if c.startswith('License ::')]
    assert f'Issues, {REPO_URL}/issues' in metadata.get_all('Project-URL', [])


def test_readme_links_work_on_pypi():
    """PyPI shows the README without the repository's files and without heading anchors, so every link in it is a
    full address."""
    readme = Path(ROOT, 'README.md').read_text(encoding='utf-8')
    targets = [target for pattern in LINK_TARGETS for target in pattern.findall(readme)]
    assert targets
    assert [target for target in targets if not target.startswith('https://')] == []


def test_changelog_names_the_version():
    """The newest release in CHANGELOG.md is the package's version."""
    changelog = Path(ROOT, 'CHANGELOG.md').read_text(encoding='utf-8')
    assert CHANGELOG_RELEASE.findall(changelog)[0] == raindrop.__version__


def test_pypi_page_links_the_changelog(tmp_path):
    with zipfile.ZipFile(build_wheel(tmp_path)) as wheel:
        metadata = read_metadata(wheel)
    assert f'Changelog, {REPO_URL}/blob/master/CHANGELOG.md' in metadata.get_all('Project-URL', [])


def test_citation_names_the_version():
    """The version in CITATION.cff is the package's version."""
    citation = Path(ROOT, 'CITATION.cff').read_text(encoding='utf-8')
    assert CITATION_VERSION.findall(citation) == [raindrop.__version__]


def test_citation_gives_the_zenodo_doi():
    """CITATION.cff cites one Zenodo DOI, which its identifiers list too, and every identifier is a Zenodo DOI."""
    citation = Path(ROOT, 'CITATION.cff').read_text(encoding='utf-8')
    dois = CITATION_DOI.findall(citation)
    assert len(dois) == 1
    assert ZENODO_DOI.fullmatch(dois[0])
    identifiers = CITATION_IDENTIFIER.findall(citation)
    assert dois[0] in identifiers
    assert all(ZENODO_DOI.fullmatch(identifier) for identifier in identifiers)


def test_readme_badge_links_the_citation_doi():
    """The README's DOI badge shows and links the doi of CITATION.cff."""
    citation = Path(ROOT, 'CITATION.cff').read_text(encoding='utf-8')
    readme = Path(ROOT, 'README.md').read_text(encoding='utf-8')
    dois = CITATION_DOI.findall(citation)
    assert dois
    assert f'[![DOI](https://zenodo.org/badge/DOI/{dois[0]}.svg)](https://doi.org/{dois[0]})' in readme


def test_sdist_runs_its_own_tests(tmp_path):
    """The sdist carries every file its tests need (the packaging tests also read the changelog and the citation);
    its command-line and packaging tests pass from the unpacked sdist (the sdist tests themselves are left out, as
    they would build another sdist)."""
    with tarfile.open(build_sdist(tmp_path)) as sdist:
        top = sdist.getnames()[0].split('/')[0]
        for name in ('tests/conftest.py', 'convert.py', 'tests/test_viewer.mjs', 'CHANGELOG.md', 'CITATION.cff'):
            assert f'{top}/{name}' in sdist.getnames(), name
        # The 'data' filter (Python 3.12, backported to later 3.10 and 3.11 releases) refuses unsafe members and
        # keeps 3.14 from warning; older releases extract as before
        safe = {'filter': 'data'} if hasattr(tarfile, 'data_filter') else {}
        sdist.extractall(tmp_path / 'unpacked', **safe)
    tests = subprocess.run([sys.executable, '-m', 'pytest', '-q', '-p', 'no:cacheprovider', '-k', 'not sdist',
                            '--basetemp', str(tmp_path / 'inner'), 'tests/test_cli.py', 'tests/test_packaging.py'],
                           cwd=tmp_path / 'unpacked' / top, capture_output=True, text=True)
    assert tests.returncode == 0, tests.stdout[-3000:]


def test_license_metadata_declares_three_js(tmp_path):
    """The wheel and the sdist declare both licenses (Raindrop's LGPL, three.js's MIT) and the wheel carries both
    texts."""
    with zipfile.ZipFile(build_wheel(tmp_path)) as wheel:
        metadata = read_metadata(wheel)
        licenses = [f'{dist_info(wheel)}/licenses/{name}' for name in ('LICENSE', 'viewer/vendor/THREE_LICENSE')]
        assert set(licenses) <= set(wheel.namelist())
    assert metadata['License-Expression'] == 'LGPL-3.0-or-later AND MIT'
    with tarfile.open(build_sdist(tmp_path)) as sdist:
        pkg_info = next(name for name in sdist.getnames() if name.count('/') == 1 and name.endswith('/PKG-INFO'))
        sdist_metadata = email.parser.Parser().parsestr(sdist.extractfile(pkg_info).read().decode())
    assert sdist_metadata['License-Expression'] == 'LGPL-3.0-or-later AND MIT'


def test_version():
    # Raised with each release, together with CHANGELOG.md and CITATION.cff (checked above)
    assert raindrop.__version__ == '0.1.1'
