"""
Tests for the raindrop command: `raindrop convert` behaves like `python convert.py`.

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
import os
import subprocess
import sys

import pytest

from raindrop import cli
from test_convert import CONVERT_PY, make_adata

# Options that touch several parts of the converter: defaults, filters, start ticks and style, axes
OPTIONS = ['--color', 'Subclass', '--filters', 'section', '--select', 'section', 's2',
           '--start-style', 'shown-size=2', 'filtered-opacity=10', '--up', 'spatial', '0', '-1', '0',
           '--exclude', 'passed']


def folder_files(folder):
    """Every file in a folder, by path relative to it, with its bytes."""
    files = {}
    for root, _, names in os.walk(folder):
        for name in names:
            path = os.path.join(root, name)
            with open(path, 'rb') as f:
                files[os.path.relpath(path, folder)] = f.read()
    return files


def test_convert_matches_convert_py(tmp_path, capsys):
    """`raindrop convert IN OUT [options]` writes the same files and the same summary as `python convert.py`."""
    adata, _ = make_adata()
    adata.write_h5ad(tmp_path / 'example.h5ad')
    script = subprocess.run([sys.executable, CONVERT_PY, str(tmp_path / 'example.h5ad'), str(tmp_path / 'script'),
                             *OPTIONS], capture_output=True, text=True, check=True)
    cli.main(['convert', str(tmp_path / 'example.h5ad'), str(tmp_path / 'command'), *OPTIONS])

    assert folder_files(tmp_path / 'command') == folder_files(tmp_path / 'script')
    summary = capsys.readouterr().out
    assert summary.replace(str(tmp_path / 'command'), 'OUT') == script.stdout.replace(str(tmp_path / 'script'), 'OUT')


def test_convert_errors_match_convert_py(tmp_path, capsys):
    """A bad option value stops `raindrop convert` with the converter's message and exit status, writing nothing."""
    adata, _ = make_adata()
    adata.write_h5ad(tmp_path / 'example.h5ad')
    bad = ['--start-style', 'shown-opacity=150']
    script = subprocess.run([sys.executable, CONVERT_PY, str(tmp_path / 'example.h5ad'), str(tmp_path / 'script'),
                             *bad], capture_output=True, text=True)
    with pytest.raises(SystemExit) as stop:
        cli.main(['convert', str(tmp_path / 'example.h5ad'), str(tmp_path / 'command'), *bad])

    assert stop.value.code == script.returncode == 2
    message = 'shown-opacity=150: give 0 to 100'
    assert message in capsys.readouterr().err and message in script.stderr
    assert not os.path.exists(tmp_path / 'command')


def test_convert_help_lists_converter_options(capsys):
    """`raindrop convert --help` shows the converter's description and options under the command's name."""
    with pytest.raises(SystemExit) as stop:
        cli.main(['convert', '--help'])
    assert stop.value.code == 0
    text = capsys.readouterr().out
    assert text.startswith('usage: raindrop convert')
    assert 'Convert an AnnData .h5ad file into a Raindrop dataset folder.' in text
    for option in ('--color', '--filters', '--select', '--start-style', '--gene-source', '--one-byte'):
        assert option in text


def test_command_needed(capsys):
    """`raindrop` without a command stops with usage naming the commands."""
    with pytest.raises(SystemExit) as stop:
        cli.main([])
    assert stop.value.code == 2
    assert 'convert' in capsys.readouterr().err
