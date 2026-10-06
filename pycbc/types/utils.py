import logging
import numpy as _numpy
from numpy import  float64

from pycbc.libutils import import_optional

logger = logging.getLogger('pycbc.type.utils')

_lal = import_optional('lal')

def determine_epoch(epoch, initial_array):
    """
    Determine what the value should be given the epoch input
    and initial array input to creating an array.
    Errors giving TypeError if the type cannot be determined.

    We gave a nonsensical default value ("") to FrequencySeries
    and TimeSeries epoch so we can test if it has been set.

    If this function receives this default value, then we test
    `initial_array`; if `initial_array` has an 'epoch' attribute,
    we use that, otherwise return zero

    But if the user passed in any value to FrequencySeries or Timeseries
    - even 'None' - then that will take precedence over anything set in
    the initial_array. None values are returned directly, all others
    we try to convert to float64 first.

    Parameters
    ----------
    epoch: 
        float64/number-type, LIGOTimeGPS, None
    initial_array:
        Array - only really matters if this has an _epoch set already

    Returns
    -------
    epoch: float64 or None - see logic above
    """


    if isinstance(epoch, float64) or epoch is None:
        return epoch
    
    if epoch == "":
        # The default has been given, try these:
        try:
            # inherit epoch from initial array
            return initial_array._epoch
        except AttributeError:
            # default epoch given, and we can't grab the epoch
            # from the initial array - fall back to zero
            return float64(0)

    # If we reach here, then the epoch has been given
    # but is not already a float64 or None, so we try to do conversions

    # LIGOTimeGPS is a special case, as numpy.isscalar fails, but
    # it can be converted using float64().
    # We require lal to be imported to do this check
    is_ltg = _lal is not None and isinstance(epoch, _lal.LIGOTimeGPS)

    # It looks like this is an array/list/tuple, so float conversion could
    # succeed, but we shouldn't be trying it
    if not is_ltg and not _numpy.isscalar(epoch):
        # Its not a 
        raise TypeError("epoch must be a number, not array-like")
    
    try:
        # Okay we have gone through the special cases now, just try it and see
        return float64(epoch)
    except TypeError as e:
        # Give something helpful before failing.
        logger.warning(
            "epoch cannot be determined: "
            f"type: {type(epoch)}, value: {epoch}"
        )
        raise e