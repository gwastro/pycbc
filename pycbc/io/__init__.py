import logging
from astropy.utils.data import download_file
from .hdf import *
from .record import *
from .gracedb import *

logger = logging.getLogger('pycbc.io')


def get_file(url, retry=5, **args):
    """ Retrieve file with retry upon failure

    Uses the astropy download_file but adds a retry feature for flaky
    connections. See astropy for full options
    """
    # DEBUGGING, will be removed
    if os.getenv("GITHUB_ACTIONS") == "true":
        if "gwosc.org/" in url:
            print("CANNOT ACCESS GWOSC FROM GITHUB CI")
            print(url)
            raise ValueError()
    i = 0
    while True:
        i += 1
        try:
            return download_file(url, **args)
        except Exception as e:
            logger.warning("Failed on attempt %d to download %s", i, url)
            if i >= retry:
                logger.error("Giving up on %s", url)
                raise e
