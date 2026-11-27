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

"""The raindrop command: `raindrop convert IN OUT [options]` runs the converter as `python convert.py` does;
`raindrop view my.h5ad [options]` converts into a reused folder and shows it in the browser."""
import argparse
import functools
import json
import os
import sys
import warnings
from pathlib import Path

import convert
from raindrop import cache, serve

VIEW_DESCRIPTION = f"""Show an .h5ad file in the browser: convert it (or reuse an earlier conversion of the same file
with the same options), serve the viewer and the dataset at 127.0.0.1, and open the browser there. Ctrl-C stops.

Converted folders go in the cache folder: --cache-dir if given, else ${cache.CACHE_DIR_VARIABLE} if set, else the
platform's cache folder (on Linux $XDG_CACHE_HOME/raindrop or ~/.cache/raindrop), and are never deleted
automatically. --out uses one given dataset folder instead.

Where no browser can open (an ssh session, no display, or --no-browser), the address is printed with the ssh command
that forwards the port from another computer."""
DEFAULT_PORT = 8000


def build_parser():
    """The raindrop command's parser, with one subcommand per task."""
    parser = argparse.ArgumentParser(prog='raindrop', description='Raindrop: a browser viewer for single-cell and '
                                                                  'spatial data.')
    commands = parser.add_subparsers(dest='command', required=True, title='commands')
    # The converter's own parser setup and help text, so both ways of running it take the same arguments
    convert_parser = commands.add_parser('convert', help='convert an .h5ad file into a dataset folder',
                                         description=convert.__doc__,
                                         formatter_class=argparse.RawDescriptionHelpFormatter)
    convert.add_arguments(convert_parser)
    # Errors in option values then show this subcommand's usage
    convert_parser.set_defaults(run=functools.partial(convert.run, parser=convert_parser))
    add_view_command(commands)
    return parser


def port_number(text):
    """argparse type of --port: a TCP port number, 0 to 65535."""
    port = int(text)
    if not 0 <= port <= 65535:
        raise argparse.ArgumentTypeError(f'{text} is not a port number (0 to 65535)')
    return port


def add_view_command(commands):
    """Add the `view` subcommand: the input file, the view options and every converter option."""
    parser = commands.add_parser('view', help='show an .h5ad file in the browser', description=VIEW_DESCRIPTION,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('h5ad', help='input AnnData file')
    parser.add_argument('--out', metavar='DIR',
                        help='dataset folder to use instead of the cache: reused if raindrop view converted it from '
                             'the same file and options, converted again if from others; a folder raindrop view did '
                             'not make is never overwritten; missing parent folders are created')
    parser.add_argument('--cache-dir', metavar='DIR',
                        help=f'cache folder for this run, created if missing (default: ${cache.CACHE_DIR_VARIABLE}, '
                             'else the platform cache folder, e.g. ~/.cache/raindrop)')
    parser.add_argument('--port', type=port_number, default=DEFAULT_PORT,
                        help=f'port at 127.0.0.1 to serve on (default {DEFAULT_PORT}; 0 picks a free one); a busy '
                             'port is an error')
    parser.add_argument('--no-browser', action='store_true',
                        help='do not open a browser; print the address and an ssh tunnel command')
    convert.add_options(parser.add_argument_group('converter options (as for raindrop convert)'))
    parser.set_defaults(run=functools.partial(run_view, parser=parser))


def choose_dataset_folder(args, key):
    """The dataset folder to use: --out, else the cache's folder for this key; raises ValueError for a cache or
    --out parent folder that cannot be written."""
    if args.out:
        if args.cache_dir:
            warnings.warn('--cache-dir has no effect with --out', stacklevel=2)
        folder = Path(args.out)
        cache.check_writable_folder(folder.absolute().parent, '--out')
        return folder
    root, source = cache.choose_cache_root(args.cache_dir, os.environ)
    cache.check_writable_folder(root, source)
    return cache.cache_folder(root, args.h5ad, key)


def prepare_dataset(folder, h5ad, kwargs, key):
    """Reuse the folder if it holds this conversion, else convert into it; say which."""
    if cache.is_reusable(folder, key):
        print(f'Reusing {folder} (converted earlier from the same file and options)', flush=True)
        return
    if folder.exists():
        print(f'Replacing {folder} (converted earlier from another version of the file or other options)',
              flush=True)
    else:
        print(f'Converting {h5ad} into {folder}', flush=True)
    cache.build_dataset(folder, h5ad, kwargs, key)
    manifest = json.loads((folder / 'manifest.json').read_text())
    print(f'Converted: {convert.describe(manifest)}', flush=True)


def show(server, no_browser):
    """Print the address, and open the browser or print the ssh command that reaches the address."""
    print(f'Serving {server.name} at {server.url} (Ctrl-C to stop)', flush=True)
    opened = not no_browser and serve.should_try_browser(os.environ, sys.platform) and serve.open_browser(server.url)
    if not opened:
        print('Open the address above in a browser. From another computer, first forward the port, then open the '
              'same address there:', flush=True)
        print(f'    {serve.tunnel_command(server.port)}', flush=True)


def run_view(args, parser):
    """Run `raindrop view`: check the request, take the port, reuse or convert the dataset, then serve until Ctrl-C."""
    # Everything that can be checked quickly is checked before a conversion that may take minutes
    try:
        kwargs = convert.convert_kwargs(args)
        key = cache.dataset_key(args.h5ad, kwargs)
        folder = choose_dataset_folder(args, key)
        cache.check_target(folder)
        server = serve.ViewerServer(args.port, folder.name, folder)
    # OSError: an input file that is missing or cannot be read, a user folder at --out (FileExistsError), or a port
    # that cannot be bound for another reason
    except (ValueError, OSError) as error:
        parser.error(str(error))
    with server:
        try:
            prepare_dataset(folder, args.h5ad, kwargs, key)
            show(server, args.no_browser)
            server.serve_forever()
        except KeyboardInterrupt:
            print('Stopped.', flush=True)


def main(argv=None):
    """Run the raindrop command with the given arguments (default: the command line's)."""
    parser = build_parser()
    args = parser.parse_args(argv)
    args.run(args)
