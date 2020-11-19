from astropy.utils.data import download_file
from .hdf import *
from .record import *


def get_file(url, retry=5, **args):
    """ Retrieve file with retry upon failure

    Uses the astropy download_file but adds a retry feature for flaky
    connections. See astropy for full options
    """
    i = 0
    while 1:
        i += 1
        try:
            return download_file(url, **args)
        except Exception as e:
            print("Failed on attempt {} to download {}".format(i, url))
            if i >= retry:
                print("Giving up on {}".format(url))
                raise e
