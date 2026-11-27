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

"""Run the converter (raindrop/convert.py) from a clone of the repository: python convert.py IN OUT [options].

This script is not installed; an installed raindrop runs the same converter as `raindrop convert`. Python puts the
script's own folder first on the module search path, so the clone's raindrop package is used from any working
directory, ahead of an installed one (except under python -P or PYTHONSAFEPATH, which leave that folder out).
"""
from raindrop.convert import main

if __name__ == '__main__':
    main()
