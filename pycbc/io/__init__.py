import os
import logging
from astropy.utils.data import download_file
import hashlib
from urllib.parse import urlparse
from .hdf import *
from .record import *
from .gracedb import *

logger = logging.getLogger('pycbc.io')

# Backup URL in case GWOSC fails
base_backup_url = "https://raw.githubusercontent.com/gwastro/pycbc_data/master/{}"
base_lfs_backup_url = "https://media.githubusercontent.com/media/gwastro/pycbc_data/master/{}"


def get_file(url, retry=5, **args):
    """ Retrieve file with retry upon failure

    Uses the astropy download_file but adds a retry feature for flaky
    connections. See astropy for full options
    """
    i = 0
    if os.getenv("GITHUB_ACTIONS") == "true":
        # Accessing GWOSC from GitHub Actions is a pain and often fails.
        # If this is in GitHub Actions we divert the URLs to a backup path
        if "gwosc.org/" in url:
            basename = os.path.basename(urlparse(url).path)
            if basename.endswith('hdf5') or basename.endswith('gwf'):
                # Just download file directly from backup
                new_url = base_lfs_backup_url.format(basename)
                logger.warning(
                    "Redirecting GWOSC URL %s to backup url %s",
                    url,
                    new_url
                )
                url = new_url
            else:
                cleaned_url = url.strip().lower()
                hash_object = hashlib.md5(cleaned_url.encode('utf-8'))
                hh = hash_object.hexdigest()
                new_url = base_backup_url.format(hh + '.json')
                logger.warning(
                    "Redirecting GWOSC URL %s to backup url %s",
                    url,
                    new_url
                )
                url = new_url
    while True:
        i += 1
        try:
            return download_file(url, **args)
        except Exception as e:
            logger.warning("Failed on attempt %d to download %s", i, url)
            if i >= retry:
                logger.error("Giving up on %s", url)
                raise e
