import logging

from astropy.time import Time
from datetime import time as dtt
import datetime

logger = logging.getLogger('pycbc.live.plotting_utils')

def strip_time(time):
    """Strip off the time information, give gps time at the previous midnight"""
    midnight_date = Time(time, format='gps', scale='utc').to_datetime().date()
    midnight_dt = datetime.datetime.combine(midnight_date, dtt.min)
    return Time(midnight_dt, format='datetime').gps, midnight_date.strftime("%Y-%m-%d")


__all__ = [
    'strip_time',
]