# Changelog

All notable changes to Raindrop are listed here. The format follows [Keep a Changelog](https://keepachangelog.com/),
and versions follow [Semantic Versioning](https://semver.org/).

## [Unreleased]

## [0.1.1] - 2026-11-27

### Changed

- Raindrop is on PyPI as `raindrop-viewer` (`pip install raindrop-viewer`), since the name `raindrop` belongs to
  another project there. The command (`raindrop`) and the Python package (`import raindrop`) keep their names. If you
  installed 0.1.0 (from GitHub or a clone), run `python -m pip uninstall raindrop` first: both versions install the
  `raindrop` package and command.
- The converter is now part of the package, as `raindrop.convert`, and an install no longer adds a top-level
  `convert` module: code that did `import convert` now does `from raindrop import convert`. In a clone,
  `python convert.py ...` works as before, with the same options and help.
- The package's license metadata names three.js's MIT license (for the files in `viewer/vendor/`) next to
  Raindrop's LGPL.
- `raindrop view` converts each dataset again once after the upgrade, since the version is part of its cache key; the
  folders 0.1.0 made stay in the cache folder (e.g. `~/.cache/raindrop`) until you delete them.

## [0.1.0] - 2026-11-27

First release.

### Added

- The converter: an `.h5ad` file becomes a folder of static files that any web host can serve
  (`raindrop convert IN OUT`, or `python convert.py` in a clone). Options set the opening color, filters, selection
  and point style; anatomical axis names and orientation; two-level groupings of categories; extra gene sources,
  log-scaled values and one-byte gene storage for large datasets; category colors from scanpy-style or per-cell color
  columns.
- The viewer, plain HTML and JavaScript with no server program: 3D and 2D embeddings (e.g. spatial coordinates and
  UMAP) as point clouds; color by any obs column or gene; filter panels with search, reverse and clear; a legend that
  filters; separate point size and opacity for shown and filtered-out cells; up to 12 linked views side by side;
  turning controls, side views and axis labels; a hover label with each cell's value. It downloads only what it
  needs and filters in the browser.
- `raindrop view my.h5ad`: converts the file into a reused cache folder, serves it on this machine and opens the
  browser; on a remote machine it prints the address and the `ssh -L` command to reach it.
- Embedding in web pages: `?embed=1` gives a compact widget for an iframe, with a More / Fewer controls toggle, Full
  screen and a layout for phones and tablets; `controls=compact` or `controls=full` picks how it opens.
- Robust loading: a download that finishes after a newer choice is dropped, and a failed one shows an error and can
  be retried.

[Unreleased]: https://github.com/NingWang1729/raindrop/compare/v0.1.1...HEAD
[0.1.1]: https://github.com/NingWang1729/raindrop/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/NingWang1729/raindrop/releases/tag/v0.1.0
