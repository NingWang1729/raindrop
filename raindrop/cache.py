# raindrop: a browser viewer for single-cell and spatial data
# Copyright (C) 2026 Ning Wang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""Dataset folders for `raindrop view`: where they live, when one can be reused, and converting into one safely.

A folder made by `raindrop view` holds, beside the dataset, a marker file with a key: a hash of the input file's
identity (resolved path, size, modification time), of the other files the options read, of the converter options and
of the raindrop version. A folder whose key matches is reused as it is; one with another key is converted again and
replaced; a folder without the marker (the user's own, or one written by convert.py) is never touched.
"""
import hashlib
import json
import os
import shutil
import tempfile
from pathlib import Path

import platformdirs

import convert
from raindrop import __version__

# Environment variable that sets the cache folder, e.g. to a large scratch disk
CACHE_DIR_VARIABLE = 'RAINDROP_CACHE_DIR'
# Marker file in each dataset folder raindrop view makes; it holds only the key, since folders may be published
MARKER = 'view_key.json'
# Hex digits of the key in a cache folder's name (48 bits: collisions among one user's datasets are negligible, and the
# marker holds the full key anyway)
NAME_KEY_DIGITS = 12


def choose_cache_root(cache_dir, environ):
    """The cache folder and where it was set: --cache-dir, else RAINDROP_CACHE_DIR, else the platform's cache folder
    (on Linux $XDG_CACHE_HOME/raindrop or ~/.cache/raindrop)."""
    if cache_dir:
        return Path(cache_dir).expanduser(), '--cache-dir'
    if environ.get(CACHE_DIR_VARIABLE):
        return Path(environ[CACHE_DIR_VARIABLE]).expanduser(), CACHE_DIR_VARIABLE
    return Path(platformdirs.user_cache_dir('raindrop', appauthor=False)), 'the default cache folder'


def check_writable_folder(folder, source):
    """Create a folder if missing; raise ValueError naming `source` (the option that set it) if it is not a folder
    or cannot be written."""
    # Writing a temporary folder in it is the reliable test on every platform (os.access is not, on Windows)
    try:
        folder.mkdir(parents=True, exist_ok=True)
        tempfile.TemporaryDirectory(dir=folder).cleanup()
    except OSError as error:
        raise ValueError(f'{folder} (from {source}) is not a writable folder ({error}); choose another with '
                         f'--cache-dir, --out or {CACHE_DIR_VARIABLE}') from None


def file_identity(path):
    """What makes a file the same file for the cache: resolved path, size and modification time (nanoseconds)."""
    status = os.stat(path)
    return {'path': str(Path(path).resolve()), 'size': status.st_size, 'mtime_ns': status.st_mtime_ns}


def dataset_key(h5ad, kwargs):
    """Hex SHA-256 key of a conversion: the input file, the other files the options read (--gene-source files and
    --hierarchy-csv files), convert()'s keyword arguments, and the raindrop and dataset format versions."""
    read_files = [path for _, path in kwargs['gene_sources']] + list(kwargs['hierarchy_csvs'])
    record = {'raindrop': __version__, 'format': convert.FORMAT_VERSION, 'input': file_identity(h5ad),
              'files': [file_identity(path) for path in read_files], 'options': kwargs}
    # sort_keys makes the text, and so the key, independent of dict order
    return hashlib.sha256(json.dumps(record, sort_keys=True).encode()).hexdigest()


def cache_folder(root, h5ad, key):
    """Dataset folder of a conversion in the cache: <input file name without extension>-<start of the key>."""
    return root / f'{Path(h5ad).stem}-{key[:NAME_KEY_DIGITS]}'


def has_marker(folder):
    """Whether raindrop view made this folder."""
    return (folder / MARKER).is_file()


def is_reusable(folder, key):
    """Whether a folder holds the conversion with this key."""
    if not has_marker(folder):
        return False
    # A damaged marker (e.g. from an interrupted copy) counts as another key, so the folder is converted again
    try:
        return json.loads((folder / MARKER).read_text())['key'] == key
    except (ValueError, KeyError, TypeError):
        return False


def check_target(folder):
    """Raise FileExistsError if the folder exists but raindrop view did not make it, so it must not be replaced."""
    if folder.exists() and not has_marker(folder):
        raise FileExistsError(f'{folder} exists and was not made by raindrop view, so it is left as it is; remove it '
                              f'or choose another folder with --out')


def build_dataset(folder, h5ad, kwargs, key):
    """Convert into `folder`, replacing a folder raindrop view made from another key (check_target first).

    The conversion goes into a temporary folder beside it, which is then renamed into place, so `folder` is never
    seen half written: a crash or Ctrl-C leaves only the dot-named temporary folder, which is never reused. If another
    raindrop view of the same file and options finishes first, its folder is kept.
    """
    temporary = Path(tempfile.mkdtemp(prefix=f'.{folder.name}.', suffix='.tmp', dir=folder.parent))
    try:
        converted = temporary / 'dataset'
        convert.convert(str(h5ad), str(converted), **kwargs)
        (converted / MARKER).write_text(json.dumps({'key': key}))
        if is_reusable(folder, key):
            return
        # Move a folder from another key out of the way; it is deleted with the temporary folder
        if has_marker(folder):
            os.rename(folder, temporary / 'replaced')
        try:
            os.rename(converted, folder)
        except OSError:
            # Renaming onto an existing non-empty folder fails: another run created it since the check above
            if not is_reusable(folder, key):
                raise RuntimeError(f'{folder} was changed by another program while converting; run raindrop view '
                                   f'again') from None
    finally:
        shutil.rmtree(temporary)
