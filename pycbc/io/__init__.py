import time
import logging
from urllib.error import HTTPError
from astropy.utils.data import download_file
from .hdf import *
from .record import *
from .gracedb import *

logger = logging.getLogger('pycbc.io')

# HTTP status codes for which retrying is pointless -- the server has told us
# the file is not there, or not ours to have.
_NO_RETRY_HTTP_CODES = (400, 401, 403, 404, 410)


def get_file(url, retry=5, backoff_cap=5, **args):
    """Retrieve a file, retrying a transient failure with back-off.

    Wraps :func:`astropy.utils.data.download_file`, adding bounded
    exponential back-off between attempts.  A definitive HTTP error (e.g. a
    404) is raised at once rather than retried, so a genuine "not found" is
    reported quickly.  See astropy for the full set of keyword options.

    Retrying here covers the caller who reaches the origin server directly.
    In CI the request goes through a caching proxy instead (see the
    github-ci-cache action), so a file fetched by an earlier run is served
    from disk and never depends on the origin being reachable at all.

    Parameters
    ----------
    url : str
        The URL of the file to retrieve.
    retry : int
        Maximum number of attempts.
    backoff_cap : int
        Upper bound, in seconds, on the wait between attempts.
    """
    if retry < 1:
        raise ValueError(f"retry must be a positive integer, got {retry}")
    last_exc = None
    for i in range(1, retry + 1):
        try:
            return download_file(url, **args)
        except HTTPError as exc:
            if exc.code in _NO_RETRY_HTTP_CODES:
                raise
            last_exc = exc
            logger.warning("Failed on attempt %d to download %s (HTTP %d)",
                           i, url, exc.code)
        except Exception as exc:
            # Any other failure (URLError, ContentTooShortError, socket
            # timeout, ...) is treated as transient and retried.
            last_exc = exc
            logger.warning("Failed on attempt %d to download %s (%s)",
                           i, url, type(exc).__name__)
        if i < retry:
            time.sleep(min(backoff_cap, 2 ** i))
    logger.error("Giving up on %s", url)
    raise last_exc
