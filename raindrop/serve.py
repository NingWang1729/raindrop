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

"""The local web server of `raindrop view`, and opening a browser on it (or saying how to reach it over ssh).

The server listens on 127.0.0.1 only and serves two folders under a random token, so other users of a shared machine
cannot guess the address: /<token>/viewer/ is the viewer and /<token>/data/<name>/ the dataset folder. Every other
path, directory listings included, is 404.
"""
import errno
import getpass
import os
import secrets
import socket
import webbrowser
from http import HTTPStatus
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import quote, unquote, urlsplit

HOST = '127.0.0.1'
# Host header names of this machine; any other (e.g. a web page's own name pointed at 127.0.0.1, "DNS rebinding") is
# refused. The port is not checked: an ssh tunnel may forward another local port
LOCAL_HOSTS = {'127.0.0.1', 'localhost', '::1'}
# Bytes of randomness in the URL token (16 bytes: 22 URL-safe characters)
TOKEN_BYTES = 16
PACKAGE = Path(__file__).parent
# The viewer is inside the package when installed from a wheel or sdist (pyproject.toml copies viewer/ in), and beside
# it in a clone (an editable install, or tests run from the source tree)
VIEWER_PLACES = (PACKAGE / 'viewer', PACKAGE.parent / 'viewer')
# Environment variables set in an ssh session, where a browser would open on the remote machine, not the user's
SSH_VARIABLES = ('SSH_CONNECTION', 'SSH_CLIENT', 'SSH_TTY')
# Environment variables of a graphical session on Linux and other Unix systems; without one, Python's webbrowser would
# start a text browser (lynx, w3m) in the terminal
DISPLAY_VARIABLES = ('DISPLAY', 'WAYLAND_DISPLAY')
# Platforms whose browsers open without a display variable
DESKTOP_PLATFORMS = ('darwin', 'win32')


def viewer_dir():
    """The viewer's folder (the one with index.html)."""
    for place in VIEWER_PLACES:
        if (place / 'index.html').is_file():
            return place
    raise RuntimeError(f'the viewer is missing from both {VIEWER_PLACES[0]} and {VIEWER_PLACES[1]}; reinstall raindrop')


def url_parts(url_path):
    """The percent-decoded parts of a request's path, without query and fragment ('/a/b/' gives ['a', 'b', ''])."""
    return [unquote(part) for part in urlsplit(url_path).path.split('/')[1:]]


def is_local_host(host):
    """Whether a Host header names this machine (any port); a malformed one does not."""
    try:
        return host is not None and urlsplit(f'//{host}').hostname in LOCAL_HOSTS
    except ValueError:
        return False


class ViewerHandler(SimpleHTTPRequestHandler):
    """Serves only the server's two folders, to requests addressed to this machine."""

    # ES modules load only with a JavaScript type; on Windows, mimetypes may read another type for .js from the registry
    extensions_map = {**SimpleHTTPRequestHandler.extensions_map, '.js': 'text/javascript'}

    def send_head(self):
        if not is_local_host(self.headers.get('Host')):
            self.send_error(HTTPStatus.FORBIDDEN, 'Host not allowed')
            return None
        try:
            self.server.local_path(self.path)
        except LookupError:
            self.send_error(HTTPStatus.NOT_FOUND)
            return None
        return super().send_head()

    def translate_path(self, path):
        return str(self.server.local_path(path))

    def list_directory(self, path):
        self.send_error(HTTPStatus.NOT_FOUND)
        return None

    def log_request(self, code='-', size='-'):
        # One line per downloaded file would bury the address; errors are still logged (log_error)
        pass


class ViewerServer(ThreadingHTTPServer):
    """Web server on 127.0.0.1 for the viewer and one dataset folder (`name` is its name in the URL)."""

    # SO_REUSEADDR lets a restarted server take its port at once on Unix; on Windows it would let two servers share a
    # port, hiding a busy one
    allow_reuse_address = os.name != 'nt'
    daemon_threads = True

    def __init__(self, port, name, folder):
        self.token = secrets.token_urlsafe(TOKEN_BYTES)
        self.name = name
        self.mounts = {(self.token, 'viewer'): viewer_dir(), (self.token, 'data', name): Path(folder)}
        try:
            super().__init__((HOST, port), ViewerHandler)
        except OSError as error:
            # A port in use, or one this user may not take (below 1024 on Unix), is the user's choice to change
            if error.errno == errno.EADDRINUSE:
                problem = 'is in use'
            elif error.errno == errno.EACCES:
                problem = 'is not allowed for this user'
            else:
                raise
            raise ValueError(f'port {port} {problem}; choose another with --port (--port 0 picks a free one)') from None

    @property
    def port(self):
        return self.server_address[1]

    @property
    def url(self):
        """The viewer's address with this dataset: ?data= is relative to /<token>/viewer/."""
        return f'http://{HOST}:{self.port}/{self.token}/viewer/?data=../data/{quote(self.name)}/'

    def local_path(self, url_path):
        """The file or folder a request path names; LookupError if it is outside the two served folders."""
        parts = url_parts(url_path)
        # '..' and separators inside a part could step out of a folder; the check after resolving catches symlinks
        if any(part in ('.', '..') or '/' in part or os.sep in part for part in parts):
            raise LookupError(url_path)
        for prefix, root in self.mounts.items():
            if tuple(parts[:len(prefix)]) == prefix:
                path = root.joinpath(*parts[len(prefix):])
                if not path.resolve().is_relative_to(root.resolve()):
                    raise LookupError(url_path)
                return path
        raise LookupError(url_path)


def should_try_browser(environ, platform):
    """Whether a browser can open on the user's screen: not in an ssh session, and with a display on Unix."""
    if any(environ.get(name) for name in SSH_VARIABLES):
        return False
    return platform in DESKTOP_PLATFORMS or any(environ.get(name) for name in DISPLAY_VARIABLES)


def open_browser(url):
    """Open the address in the default browser; whether it worked."""
    try:
        return webbrowser.open(url)
    except webbrowser.Error:
        return False


def tunnel_command(port):
    """The ssh command that forwards the port from another computer, with this machine's user and host names where
    known (else <user> and <host>)."""
    try:
        user = getpass.getuser()
    except (OSError, KeyError):
        user = '<user>'
    host = socket.getfqdn()
    if not host or host.startswith('localhost'):
        host = '<host>'
    return f'ssh -L {port}:{HOST}:{port} {user}@{host}'
