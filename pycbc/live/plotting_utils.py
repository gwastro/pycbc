import logging

from astropy.time import Time
from datetime import time as dtt
import datetime

logger = logging.getLogger('pycbc.live.plotting_utils')

def strip_time(gps_time):
    """
    Strip off the time information from a GPS time, returning the GPS time at the previous midnight and the corresponding date string.

    Parameters
    ----------
    gps_time : float or int
        GPS time to be stripped to midnight (seconds since GPS epoch).

    Returns
    -------
    midnight_gps : float
        GPS time at the previous midnight (00:00:00 UTC).
    date_str : str
        Date string in the format 'YYYY-MM-DD' corresponding to the previous midnight.
    """
    midnight_date = Time(gps_time, format='gps', scale='utc').to_datetime().date()
    midnight_dt = datetime.datetime.combine(midnight_date, dtt.min)
    return Time(midnight_dt, format='datetime').gps, midnight_date.strftime("%Y-%m-%d")


__all__ = [
    'strip_time',
]