import os
import time
import logging
from urllib.error import HTTPError
from urllib.parse import urlparse
import hashlib
from astropy.utils.data import download_file
from .hdf import *
from .record import *
from .gracedb import *

logger = logging.getLogger('pycbc.io')

# Backup locations for the GWOSC data files used by the test suite / examples,
# in gwastro/pycbc_data.  Binary frame/ROM files are stored as release assets
# (not Git LFS, so downloads are not rate-limited), with the old Git-LFS media
# location as a second fallback; the small pre-computed GWOSC JSON responses
# are plain files in the repo, keyed by an md5 of the request URL.
base_backup_url = "https://raw.githubusercontent.com/gwastro/pycbc_data/master/{}"
base_lfs_backup_url = (
    "https://github.com/gwastro/pycbc_data/releases/download/ci-data/{}"
)
base_lfs_media_url = (
    "https://media.githubusercontent.com/media/gwastro/pycbc_data/master/{}"
)

# HTTP status codes for which retrying the same URL is pointless; move straight
# on to the next candidate location instead of waiting and trying again.
_NO_RETRY_HTTP_CODES = (400, 401, 403, 404, 410)

# Hosts whose frame/ROM files are mirrored in gwastro/pycbc_data.  GWOSC has
# used several names over the years (gwosc.org, gw-openscience.org,
# losc.ligo.org) and the "cleaned" O2 strain lives on the DCC; all of them are
# regularly unreachable or slow from GitHub Actions.
_MIRRORED_HOSTS = (
    "gwosc.org/",
    "gw-openscience.org/",
    "losc.ligo.org/",
    "dcc.ligo.org/",
)
# ...of those, the ones that also serve the small GWOSC JSON API responses
# cached (by md5 of the request URL) in gwastro/pycbc_data.
_JSON_API_HOSTS = ("gwosc.org/", "gw-openscience.org/")


def _candidate_urls(url):
    """Ordered list of URLs to try for ``url``.

    Accessing GWOSC / the DCC directly from GitHub Actions frequently fails,
    and the files the test suite / examples / tutorials ask for there are
    exactly the ones mirrored in gwastro/pycbc_data.  So *only* inside GitHub
    Actions, and *only* for those hosts, the mirrored copies are tried first
    (both the release asset and the Git-LFS media copy for frame/ROM files)
    with the original URL last.  Everywhere else ``url`` is used as given -- a
    direct request outside Actions is very likely for a file that is not in
    the mirror.
    """
    if (os.getenv("GITHUB_ACTIONS") == "true"
            and any(host in url for host in _MIRRORED_HOSTS)):
        basename = os.path.basename(urlparse(url).path)
        if basename.endswith(('gwf', 'hdf5')):
            candidates = [
                base_lfs_backup_url.format(basename),
                base_lfs_media_url.format(basename),
                url,
            ]
        elif any(host in url for host in _JSON_API_HOSTS):
            hh = hashlib.md5(url.strip().lower().encode('utf-8')).hexdigest()
            candidates = [base_backup_url.format(hh + '.json'), url]
        else:
            candidates = [url]
    else:
        candidates = [url]

    # de-duplicate while preserving order
    seen = set()
    return [c for c in candidates if not (c in seen or seen.add(c))]


def get_file(url, retry=5, backoff_cap=5, **args):
    """Retrieve a file, retrying with back-off and falling back to backups.

    Wraps :func:`astropy.utils.data.download_file`, adding bounded
    exponential back-off between attempts and, for GWOSC / DCC URLs fetched
    from GitHub Actions, automatic fall-back to the mirrored copies in
    gwastro/pycbc_data.  A definitive HTTP error (e.g. a 404) moves straight
    on to the next location instead of retrying, so genuine "not found"
    failures are reported quickly.  See astropy for the full set of keyword
    options.

    Parameters
    ----------
    url : str
        The URL of the file to retrieve.
    retry : int
        Maximum number of attempts per candidate location.
    backoff_cap : int
        Upper bound, in seconds, on the wait between attempts.
    """
    if retry < 1:
        raise ValueError(f"retry must be a positive integer, got {retry}")
    candidates = _candidate_urls(url)
    last_exc = None
    for cand in candidates:
        if cand != url:
            logger.warning("Redirecting %s to backup location %s", url, cand)
        for i in range(1, retry + 1):
            try:
                return download_file(cand, **args)
            except HTTPError as exc:
                last_exc = exc
                if exc.code in _NO_RETRY_HTTP_CODES:
                    logger.warning(
                        "%s returned HTTP %d; not retrying this location",
                        cand, exc.code
                    )
                    break
                logger.warning(
                    "Failed on attempt %d to download %s (HTTP %d)",
                    i, cand, exc.code
                )
            except Exception as exc:
                # Any other failure (URLError, ContentTooShortError, socket
                # timeout, ...) is treated as transient and retried.
                last_exc = exc
                logger.warning(
                    "Failed on attempt %d to download %s (%s)",
                    i, cand, type(exc).__name__
                )
            if i < retry:
                time.sleep(min(backoff_cap, 2 ** i))
    logger.error("Giving up on %s", url)
    raise last_exc
