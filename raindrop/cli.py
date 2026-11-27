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

"""The raindrop command: `raindrop convert IN OUT [options]` runs the converter as `python convert.py` does."""
import argparse
import functools

import convert


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
    return parser


def main(argv=None):
    """Run the raindrop command with the given arguments (default: the command line's)."""
    parser = build_parser()
    args = parser.parse_args(argv)
    args.run(args)
