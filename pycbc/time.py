"""
Module to contain time conversions used in pycbc
"""

from astropy.time import Time
from datetime import timezone


def _ensure_utc_datetime(date):
    """
    Ensure a datetime is timezone-aware in UTC.

    - If `date` is naive, assume it represents UTC and attach UTC tzinfo.
    - If `date` is timezone-aware, convert it to UTC.
    """
    if date.tzinfo is None:
        return date.replace(tzinfo=timezone.utc)
    return date.astimezone(timezone.utc)


def gps_to_utc_datetime(gps):
    """
    Convert a GPS time to a UTC datetime.

    Parameters
    ----------
    gps : float
        The GPS time to convert.

    Returns
    -------
    datetime
    """
    dt = Time(gps, format='gps', scale='utc').to_datetime()
    return _ensure_utc_datetime(dt)


def datetime_to_str(date, format="%Y-%m-%d %H:%M:%S"):
    """
    Convert a datetime to a string.

    Parameters
    ----------
    date : datetime
        The datetime to convert.
    format : str
        The format to use for the string. Default is yyyy-mm-dd HH:MM:SS.
        Supply as a format according to datetime documentation
        https://docs.python.org/3/library/datetime.html#strftime-strptime-behavior
    
    Returns
    -------
    str
    """
    return date.strftime(format)


def gps_to_utc_str(gps, format="%Y-%m-%d %H:%M:%S"):
    """
    Convert a GPS time to a UTC string.

    Parameters
    ----------
    gps : float
        The GPS time to convert.
    format : str

    Returns
    -------
    str
    """
    return datetime_to_str(gps_to_utc_datetime(gps), format=format)


def strip_time_from_date(date):
    """
    Strip the time from a datetime object.

    Parameters
    ----------
    date : datetime
        The datetime object to strip the time from.

    Returns
    -------
    datetime
    """
    stripped = date.replace(hour=0, minute=0, second=0, microsecond=0)
    return _ensure_utc_datetime(stripped)


def strip_time_from_gps(gps, format="%Y-%m-%d"):
    """
    Round a GPS time to the midnight of the same UTC day, and return the date of that day as a string.

    Parameters
    ----------
    gps : float
        The GPS time to strip the time from.

    Returns
    -------
    float, str
    """
    gps_datetime = gps_to_utc_datetime(gps)
    midnight_datetime = strip_time_from_date(gps_datetime)
    return utc_datetime_to_gps(midnight_datetime), datetime_to_str(midnight_datetime, format=format)


def utc_datetime_to_gps(date):
    """
    Convert a UTC datetime to a GPS time.

    Parameters
    ----------
    date : datetime
        The UTC datetime to convert.

    Returns
    -------
    float
    """
    date_utc = _ensure_utc_datetime(date)
    return float(Time(date_utc, format='datetime', scale='utc').gps)


def gps_now():
    """Return the current GPS time as a float using Astropy.

    Returns
    -------
    float
    """

    return float(Time.now().gps)


def gmst_accurate(gps_time):
    """
    Calculate the Greenwich Mean Sidereal Time in radians
    for a given GPS time.

    Parameters
    ----------
    gps_time : float
        The GPS time to calculate the GMST for.

    Returns
    -------
    float
    """
    gmst = Time(gps_time, format='gps', scale='utc',
                location=(0, 0)).sidereal_time('mean').rad
    return gmst

__all__ = [
    'gps_to_utc_datetime',
    'gps_to_utc_str',
    'strip_time_from_date',
    'strip_time_from_gps',
    'utc_datetime_to_gps',
    'gps_now',
    'gmst_accurate',
]