"""
Tests for `raindrop view`: the dataset cache, the local server and the browser or tunnel hint.

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
import http.client
import json
import os
import re
import socket
import sys
import threading
import webbrowser

import pytest

import convert
from raindrop import cli, serve
from test_convert import make_adata

# Environment variables that decide the cache folder and whether a browser can open; each test sets its own
ENVIRONMENT = ('RAINDROP_CACHE_DIR', 'XDG_CACHE_HOME', 'SSH_CONNECTION', 'SSH_CLIENT', 'SSH_TTY', 'DISPLAY',
               'WAYLAND_DISPLAY')


@pytest.fixture(autouse=True)
def clean_environment(monkeypatch, tmp_path):
    """No cache or display settings from the machine running the tests; the default cache lives in tmp_path."""
    for name in ENVIRONMENT:
        monkeypatch.delenv(name, raising=False)
    monkeypatch.setenv('XDG_CACHE_HOME', str(tmp_path / 'xdg'))
    monkeypatch.setattr(serve.getpass, 'getuser', lambda: 'alice')
    monkeypatch.setattr(serve.socket, 'getfqdn', lambda name='': 'host.example')


@pytest.fixture
def h5ad(tmp_path):
    """A five-cell .h5ad file."""
    adata, _ = make_adata()
    path = tmp_path / 'my.h5ad'
    adata.write_h5ad(path)
    return path


@pytest.fixture
def no_serving(monkeypatch):
    """Stop the server as Ctrl-C would, as soon as it starts serving."""
    def interrupt(self, *args, **kwargs):
        raise KeyboardInterrupt
    monkeypatch.setattr(serve.ViewerServer, 'serve_forever', interrupt)


def view(*args):
    """Run `raindrop view` with test settings (a free port, no browser) and the given arguments."""
    cli.main(['view', *map(str, args), '--port', '0', '--no-browser'])


def used_folder(output):
    """The dataset folder `raindrop view` says it used."""
    return re.search(r'^(?:Reusing|Replacing|Converting \S+ into) (\S+)', output, re.MULTILINE).group(1)


def manifest_time(folder):
    return os.stat(os.path.join(folder, 'manifest.json')).st_mtime_ns


# ---- The cache


def test_second_view_reuses_the_folder(h5ad, tmp_path, capsys, no_serving):
    """The same file with the same options converts once, into the default cache; later runs reuse that folder."""
    view(h5ad, '--color', 'Subclass')
    first = capsys.readouterr().out
    folder = used_folder(first)
    assert first.startswith(f'Converting {h5ad} into ')
    assert os.path.dirname(folder) == str(tmp_path / 'xdg' / 'raindrop')
    assert os.path.basename(folder).startswith('my-')
    # The converter's options pass through
    assert json.load(open(os.path.join(folder, 'manifest.json')))['defaults']['color'] == 'Subclass'
    written = manifest_time(folder)

    view(h5ad, '--color', 'Subclass')
    second = capsys.readouterr().out
    assert second.startswith(f'Reusing {folder}') and 'Converting' not in second
    assert manifest_time(folder) == written


def test_changes_make_a_new_conversion(h5ad, tmp_path, capsys, no_serving):
    """A changed input file, other options, or a changed --hierarchy-csv or --gene-source file give another cache
    folder."""
    csv = tmp_path / 'groups.csv'
    csv.write_text('Class,Subclass\nX,a\nX,b\nY,c\n')
    view(h5ad, '--hierarchy-csv', csv)
    first = used_folder(capsys.readouterr().out)

    view(h5ad, '--hierarchy-csv', csv, '--color', 'section')
    other_options = used_folder(capsys.readouterr().out)

    csv.write_text('Class,Subclass\nX,a\nY,b\nY,c\n')
    view(h5ad, '--hierarchy-csv', csv)
    other_csv = used_folder(capsys.readouterr().out)

    adata, _ = make_adata()
    adata.obs['section'] = 's9'
    adata.write_h5ad(h5ad)
    view(h5ad, '--hierarchy-csv', csv)
    other_file = used_folder(capsys.readouterr().out)

    # A --gene-source file counts too
    extra = tmp_path / 'extra.h5ad'
    adata.write_h5ad(extra)
    view(h5ad, '--hierarchy-csv', csv, '--gene-source', 'extra', extra)
    with_source = used_folder(capsys.readouterr().out)
    adata.X = adata.X * 2
    adata.write_h5ad(extra)
    view(h5ad, '--hierarchy-csv', csv, '--gene-source', 'extra', extra)
    other_source = used_folder(capsys.readouterr().out)

    folders = (first, other_options, other_csv, other_file, with_source, other_source)
    assert len(set(folders)) == len(folders)
    assert all(os.path.exists(os.path.join(f, 'manifest.json')) for f in folders)


def test_out_folder(h5ad, tmp_path, capsys, no_serving):
    """--out converts into that folder; a second run reuses it; other options replace it, saying so."""
    out = tmp_path / 'shown'
    view(h5ad, '--out', out)
    assert used_folder(capsys.readouterr().out) == str(out)
    written = manifest_time(out)

    view(h5ad, '--out', out)
    assert capsys.readouterr().out.startswith(f'Reusing {out}')
    assert manifest_time(out) == written

    view(h5ad, '--out', out, '--color', 'section')
    assert capsys.readouterr().out.startswith(f'Replacing {out}')
    assert json.load(open(out / 'manifest.json'))['defaults']['color'] == 'section'
    # Only the dataset folder is left: no temporary folders beside it
    assert sorted(os.listdir(tmp_path)) == ['my.h5ad', 'shown']


@pytest.mark.parametrize('kind', ['user folder', 'convert.py dataset'])
def test_out_never_overwrites_other_folders(h5ad, tmp_path, capsys, no_serving, kind):
    """A folder that raindrop view did not make is refused, untouched, before anything is converted or served."""
    out = tmp_path / 'mine'
    if kind == 'user folder':
        out.mkdir()
        (out / 'notes.txt').write_text('keep me')
    else:
        convert.convert(str(h5ad), str(out))
    before = sorted(os.listdir(out))
    with pytest.raises(SystemExit) as stop:
        view(h5ad, '--out', out)
    assert stop.value.code == 2
    error = capsys.readouterr().err
    assert f'{out} exists and was not made by raindrop view' in error and '--out' in error
    assert sorted(os.listdir(out)) == before


def test_failed_conversion_leaves_nothing(h5ad, tmp_path, monkeypatch, capsys, no_serving):
    """A conversion that stops midway leaves neither the dataset folder nor temporary folders; the next run
    converts again."""
    real_convert = convert.convert

    def fail_midway(h5ad_path, out_dir, **kwargs):
        os.makedirs(out_dir)
        open(os.path.join(out_dir, 'manifest.json'), 'w').close()
        raise RuntimeError('stopped midway')
    monkeypatch.setattr(convert, 'convert', fail_midway)
    out = tmp_path / 'shown'
    with pytest.raises(RuntimeError, match='stopped midway'):
        view(h5ad, '--out', out)
    assert sorted(os.listdir(tmp_path)) == ['my.h5ad']

    monkeypatch.setattr(convert, 'convert', real_convert)
    capsys.readouterr()
    view(h5ad, '--out', out)
    assert capsys.readouterr().out.startswith(f'Converting {h5ad} into {out}')
    assert json.load(open(out / 'manifest.json'))['n_cells'] == 5


def test_concurrent_conversion_is_reused(h5ad, tmp_path, monkeypatch, capsys, no_serving):
    """When another raindrop view of the same file finishes first, its folder is kept and served."""
    out = tmp_path / 'shown'
    real_convert = convert.convert

    def other_finishes_first(h5ad_path, out_dir, **kwargs):
        # The other run converts the same file with the same options into the same folder, then this one finishes
        monkeypatch.setattr(convert, 'convert', real_convert)
        view(h5ad, '--out', out)
        real_convert(h5ad_path, out_dir, **kwargs)
    monkeypatch.setattr(convert, 'convert', other_finishes_first)
    view(h5ad, '--out', out)
    output = capsys.readouterr().out
    assert 'Serving' in output
    assert json.load(open(out / 'manifest.json'))['n_cells'] == 5
    assert sorted(os.listdir(tmp_path)) == ['my.h5ad', 'shown']


def test_damaged_marker_is_converted_again(h5ad, tmp_path, capsys, no_serving):
    """A raindrop view folder whose marker cannot be read is converted again, not reused or refused."""
    out = tmp_path / 'shown'
    view(h5ad, '--out', out)
    (out / 'view_key.json').write_text('{"key": ')
    capsys.readouterr()
    view(h5ad, '--out', out)
    assert capsys.readouterr().out.startswith(f'Replacing {out}')
    assert json.load(open(out / 'manifest.json'))['n_cells'] == 5


def test_missing_input_stops_early(tmp_path, capsys):
    with pytest.raises(SystemExit) as stop:
        view(tmp_path / 'missing.h5ad')
    assert stop.value.code == 2
    assert 'missing.h5ad' in capsys.readouterr().err
    assert not (tmp_path / 'xdg').exists()


# ---- Where the cache lives


def test_cache_dir_precedence(h5ad, tmp_path, monkeypatch, capsys, no_serving):
    """--out, then --cache-dir, then RAINDROP_CACHE_DIR, then XDG_CACHE_HOME/raindrop; missing roots are created."""
    view(h5ad)
    assert os.path.dirname(used_folder(capsys.readouterr().out)) == str(tmp_path / 'xdg' / 'raindrop')

    monkeypatch.setenv('RAINDROP_CACHE_DIR', str(tmp_path / 'env' / 'cache'))
    view(h5ad)
    assert os.path.dirname(used_folder(capsys.readouterr().out)) == str(tmp_path / 'env' / 'cache')

    view(h5ad, '--cache-dir', tmp_path / 'option')
    assert os.path.dirname(used_folder(capsys.readouterr().out)) == str(tmp_path / 'option')

    with pytest.warns(UserWarning, match='--cache-dir has no effect with --out'):
        view(h5ad, '--cache-dir', tmp_path / 'option', '--out', tmp_path / 'shown')
    assert used_folder(capsys.readouterr().out) == str(tmp_path / 'shown')


def test_default_cache_without_xdg(h5ad, tmp_path, monkeypatch, capsys, no_serving):
    """Without XDG_CACHE_HOME the default cache is ~/.cache/raindrop (Linux)."""
    if not sys.platform.startswith('linux'):
        pytest.skip('platform cache folder differs')
    monkeypatch.delenv('XDG_CACHE_HOME')
    monkeypatch.setenv('HOME', str(tmp_path / 'home'))
    view(h5ad)
    assert os.path.dirname(used_folder(capsys.readouterr().out)) == str(tmp_path / 'home' / '.cache' / 'raindrop')


@pytest.mark.parametrize('source', ['--cache-dir', 'RAINDROP_CACHE_DIR'])
@pytest.mark.parametrize('problem', ['a file', 'read-only'])
def test_unusable_cache_dir(h5ad, tmp_path, monkeypatch, capsys, no_serving, source, problem):
    """A cache folder that is a file or cannot be written stops with a message naming where it came from."""
    root = tmp_path / 'root'
    if problem == 'a file':
        root.write_text('')
    else:
        if not hasattr(os, 'geteuid') or os.geteuid() == 0:
            pytest.skip('root can write anywhere')
        root.mkdir()
        root.chmod(0o500)
    args = [h5ad, '--cache-dir', root] if source == '--cache-dir' else [h5ad]
    if source == 'RAINDROP_CACHE_DIR':
        monkeypatch.setenv('RAINDROP_CACHE_DIR', str(root))
    try:
        with pytest.raises(SystemExit) as stop:
            view(*args)
    finally:
        if problem == 'read-only':
            root.chmod(0o700)
    assert stop.value.code == 2
    error = capsys.readouterr().err
    assert f'{root} (from {source}) is not a writable folder' in error
    assert not (tmp_path / 'xdg').exists()


# ---- Port, browser and tunnel hint


def test_busy_port_names_the_option(h5ad, capsys):
    """A port in use stops before converting, naming --port."""
    with socket.socket() as taken:
        taken.bind(('127.0.0.1', 0))
        taken.listen()
        port = taken.getsockname()[1]
        with pytest.raises(SystemExit) as stop:
            cli.main(['view', str(h5ad), '--port', str(port), '--no-browser'])
    assert stop.value.code == 2
    captured = capsys.readouterr()
    assert f'port {port} is in use; choose another with --port' in captured.err
    assert 'Converting' not in captured.out


def test_privileged_port_names_the_option(h5ad, capsys):
    """A port this user may not take stops before converting, naming --port."""
    if not hasattr(os, 'geteuid') or os.geteuid() == 0:
        pytest.skip('ports below 1024 are restricted only for non-root users on Unix')
    with pytest.raises(SystemExit) as stop:
        cli.main(['view', str(h5ad), '--port', '1', '--no-browser'])
    assert stop.value.code == 2
    assert 'port 1 is not allowed for this user; choose another with --port' in capsys.readouterr().err


def test_port_zero_names_the_free_port(h5ad, capsys, no_serving):
    view(h5ad)
    output = capsys.readouterr().out
    port = int(re.search(r'http://127\.0\.0\.1:(\d+)/', output).group(1))
    assert port > 0
    assert f'ssh -L {port}:127.0.0.1:{port} alice@host.example' in output


def browser_calls(monkeypatch, result):
    """Replace webbrowser.open: record the URLs asked for, and return `result` (or raise it)."""
    calls = []

    def fake_open(url):
        calls.append(url)
        if isinstance(result, Exception):
            raise result
        return result
    monkeypatch.setattr(webbrowser, 'open', fake_open)
    return calls


@pytest.mark.parametrize('setting, opens, hint', [
    ({'DISPLAY': ':0'}, True, False),                              # a desktop: the browser opens, no hint
    ({'WAYLAND_DISPLAY': 'wayland-0'}, True, False),
    ({}, False, True),                                             # Linux without a display
    ({'DISPLAY': ':0', 'SSH_CONNECTION': '1.2.3.4 5 6.7.8.9 22'}, False, True),  # an SSH session
    ({'DISPLAY': ':0', 'SSH_TTY': '/dev/pts/1'}, False, True),
])
def test_browser_or_tunnel_hint(h5ad, monkeypatch, capsys, no_serving, setting, opens, hint):
    if not sys.platform.startswith('linux'):
        pytest.skip('display variables matter on Linux only')
    for name, value in setting.items():
        monkeypatch.setenv(name, value)
    calls = browser_calls(monkeypatch, True)
    cli.main(['view', str(h5ad), '--port', '0'])
    output = capsys.readouterr().out
    url = re.search(r'(http://127\.0\.0\.1:\d+/\S+)', output).group(1)
    assert calls == ([url] if opens else [])
    assert ('ssh -L ' in output) == hint


@pytest.mark.parametrize('result', [False, webbrowser.Error('no runnable browser')])
def test_browser_failure_gives_the_hint(h5ad, monkeypatch, capsys, no_serving, result):
    monkeypatch.setenv('DISPLAY', ':0')
    calls = browser_calls(monkeypatch, result)
    cli.main(['view', str(h5ad), '--port', '0'])
    assert len(calls) == 1
    assert 'alice@host.example' in capsys.readouterr().out


def test_no_browser_option(h5ad, monkeypatch, capsys, no_serving):
    monkeypatch.setenv('DISPLAY', ':0')
    calls = browser_calls(monkeypatch, True)
    view(h5ad)
    output = capsys.readouterr().out
    assert calls == []
    assert 'Serving' in output and 'alice@host.example' in output and 'Stopped.' in output


def test_hint_placeholders(h5ad, monkeypatch, capsys, no_serving):
    """Without a user or host name the hint keeps placeholders."""
    def no_user():
        raise OSError('no user name')
    monkeypatch.setattr(serve.getpass, 'getuser', no_user)
    monkeypatch.setattr(serve.socket, 'getfqdn', lambda name='': 'localhost')
    view(h5ad)
    assert re.search(r'ssh -L (\d+):127\.0\.0\.1:\1 <user>@<host>', capsys.readouterr().out)


# ---- The server


@pytest.fixture
def server(tmp_path):
    """A server on a free port for a small dataset folder, with a secret file beside it."""
    folder = tmp_path / 'served'
    folder.mkdir()
    (folder / 'manifest.json').write_text('{"n_cells": 5}')
    (folder / 'sub').mkdir()
    (tmp_path / 'secret.txt').write_text('secret')
    running = serve.ViewerServer(0, 'served', folder)
    thread = threading.Thread(target=running.serve_forever, daemon=True)
    thread.start()
    yield running
    running.shutdown()
    running.server_close()


def get(server, path, host='127.0.0.1'):
    """Status, Content-Type and body of a GET for a raw (unnormalized) path."""
    connection = http.client.HTTPConnection('127.0.0.1', server.server_address[1], timeout=10)
    connection.putrequest('GET', path, skip_host=True)
    connection.putheader('Host', f'{host}:{server.server_address[1]}')
    connection.endheaders()
    response = connection.getresponse()
    result = response.status, response.getheader('Content-Type'), response.read()
    connection.close()
    return result


def test_server_serves_viewer_and_dataset(server):
    assert server.server_address[0] == '127.0.0.1'
    base = re.match(r'http://127\.0\.0\.1:\d+(/[^/]+/)viewer/\?data=\.\./data/served/$', server.url).group(1)
    index = (serve.viewer_dir() / 'index.html').read_bytes()
    assert get(server, f'{base}viewer/')[::2] == (200, index)
    assert get(server, f'{base}viewer/index.html')[::2] == (200, index)
    # ES modules load only with a JavaScript type
    status, kind, _ = get(server, f'{base}viewer/js/app.js')
    assert (status, kind) == (200, 'text/javascript')
    assert get(server, f'{base}data/served/manifest.json')[::2] == (200, b'{"n_cells": 5}')
    # Hosts a browser uses for this machine, with any port (an ssh tunnel may forward another port)
    for host in ('localhost', '[::1]'):
        assert get(server, f'{base}viewer/index.html', host)[0] == 200


def test_server_refuses_everything_else(server, tmp_path):
    base = re.match(r'http://127\.0\.0\.1:\d+(/[^/]+/)', server.url).group(1)
    refused = ['/', '/viewer/index.html', '/data/served/manifest.json', f'{base}convert.py', f'{base}',
               f'{base}viewer/../../secret.txt', f'{base}viewer/%2e%2e/%2e%2e/secret.txt',
               f'{base}data/served/../../../secret.txt', f'{base}data/served/%2e%2e/%2e%2e/secret.txt',
               f'{base}data/other/manifest.json', f'{base}data/served/', f'{base}data/served/sub/',
               f'{base}viewer/js/', f'{base}data/{tmp_path}/secret.txt']
    for path in refused:
        status, _, body = get(server, path)
        assert status == 404, path
        assert b'secret' not in body
    # A page on another site whose name points at this machine (DNS rebinding) cannot read anything
    status, _, body = get(server, f'{base}data/served/manifest.json', 'attacker.example')
    assert status == 403 and b'n_cells' not in body


def test_viewer_files_found():
    assert (serve.viewer_dir() / 'index.html').is_file()


def test_every_viewer_folder_is_packaged():
    """pyproject.toml lists each folder of viewer/ as a package, so the wheel carries all of it."""
    root = os.path.join(os.path.dirname(__file__), '..')
    pyproject = open(os.path.join(root, 'pyproject.toml')).read()
    for folder in os.listdir(os.path.join(root, 'viewer')):
        if os.path.isdir(os.path.join(root, 'viewer', folder)):
            assert f'"raindrop.viewer.{folder}"' in pyproject, folder


def test_view_help_lists_options(capsys):
    with pytest.raises(SystemExit):
        cli.main(['view', '--help'])
    text = capsys.readouterr().out
    for option in ('--out', '--cache-dir', 'RAINDROP_CACHE_DIR', '--port', '--no-browser', '--color', '--select'):
        assert option in text


def test_bad_option_value_stops_early(h5ad, capsys):
    with pytest.raises(SystemExit) as stop:
        view(h5ad, '--start-style', 'shown-opacity=150')
    assert stop.value.code == 2
    assert 'shown-opacity=150: give 0 to 100' in capsys.readouterr().err
