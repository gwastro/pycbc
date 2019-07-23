#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
""" This module contains functions for calculating expected rates of noise
    and signal coincidences.
"""

import itertools
import logging
import numpy
import pycbc.detector


def multiifo_noise_lograte(log_rates, slop):
    """
    Calculate the expected rate of noise coincidences for multiple
    combinations of detectors

    Parameters
    ----------
    log_rates: dict
        Key: ifo string, Value: sequence of log single-detector trigger rates,
        units assumed to be Hz
    slop: float
        time added to maximum time-of-flight between detectors to account
        for timing error

    Returns
    -------
    expected_log_rates: dict
        Key: ifo combination string
        Value: expected log coincidence rate in the combination, units log Hz
    """
    expected_log_rates = {}

    # Order of ifos must be stable in output dict keys, so sort them
    ifos = sorted(list(log_rates.keys()))
    ifostring = ' '.join(ifos)

    # Calculate coincidence for all-ifo combination
    expected_log_rates[ifostring] = \
        combination_noise_lograte(log_rates, slop)

    # If more than one possible coincidence type exists,
    # calculate coincidence for subsets through recursion
    if len(ifos) > 2:
        # Calculate rate for each 'miss-one-out' detector combination
        subsets = itertools.combinations(ifos, len(ifos) - 1)
        for subset in subsets:
            rates_subset = {}
            for ifo in subset:
                rates_subset[ifo] = log_rates[ifo]
            sub_coinc_rates = multiifo_noise_lograte(rates_subset, slop)
            # add these sub-coincidences to the overall dictionary
            for sub_coinc in sub_coinc_rates:
                expected_log_rates[sub_coinc] = sub_coinc_rates[sub_coinc]

    return expected_log_rates


def combination_noise_rate(rates, slop):
    """
    Calculate the expected rate of noise coincidences for a combination of
    detectors
    WARNING: for high stat values, numerical underflow can occur

    Parameters
    ----------
    rates: dict
        Key: ifo string, Value: sequence of single-detector trigger rates,
        units assumed to be Hz
    slop: float
        time added to maximum time-of-flight between detectors to account
        for timing error

    Returns
    -------
    numpy array
        Expected coincidence rate in the combination, units Hz
    """
    logging.warning('combination_noise_rate() is liable to numerical '
                    'underflows, use combination_noise_lograte '
                    'instead')
    log_rates = {k: numpy.log(r) for (k, r) in rates.items()}
    # exp may underflow
    return numpy.exp(combination_noise_lograte(log_rates, slop))


def combination_noise_lograte(log_rates, slop):
    """
    Calculate the expected rate of noise coincidences for a combination of
    detectors given log of single detector noise rates

    Parameters
    ----------
    log_rates: dict
        Key: ifo string, Value: sequence of log single-detector trigger rates,
        units assumed to be Hz
    slop: float
        time added to maximum time-of-flight between detectors to account
        for timing error

    Returns
    -------
    numpy array
        Expected log coincidence rate in the combination, units Hz
    """
    # multiply product of trigger rates by the overlap time
    allowed_area = multiifo_noise_coincident_area(list(log_rates), slop)
    # list(dict.values()) is python-3-proof
    rateprod = numpy.sum(list(log_rates.values()), axis=0)
    return numpy.log(allowed_area) + rateprod


def multiifo_noise_coincident_area(ifos, slop):
    """
    Calculate the total extent of time offset between 2 detectors,
    or area of the 2d space of time offsets for 3 detectors, for
    which a coincidence can be generated
    Cannot yet handle more than 3 detectors.

    Parameters
    ----------
    ifos: list of strings
        list of interferometers
    slop: float
        extra time to add to maximum time-of-flight for timing error

    Returns
    -------
    allowed_area: float
        area in units of seconds^(n_ifos-1) that coincident values can fall in
    """
    # set up detector objects
    dets = {}
    for ifo in ifos:
        dets[ifo] = pycbc.detector.Detector(ifo)
    n_ifos = len(ifos)

    if n_ifos == 2:
        allowed_area = 2. * \
            (dets[ifos[0]].light_travel_time_to_detector(dets[ifos[1]]) + slop)
    elif n_ifos == 3:
        tofs = numpy.zeros(n_ifos)
        ifo2_num = []

        # calculate travel time between detectors (plus extra for timing error)
        # TO DO: allow for different timing errors between different detectors
        for i, ifo in enumerate(ifos):
            ifo2_num.append(int(numpy.mod(i + 1, n_ifos)))
            det0 = dets[ifo]
            det1 = dets[ifos[ifo2_num[i]]]
            tofs[i] = det0.light_travel_time_to_detector(det1) + slop

        # combine these to calculate allowed area
        allowed_area = 0
        for i, _ in enumerate(ifos):
            allowed_area += 2 * tofs[i] * tofs[ifo2_num[i]] - tofs[i]**2
    else:
        raise NotImplementedError("Not able to deal with more than 3 ifos")

    return allowed_area


def multiifo_signal_coincident_area(ifos):
    """
    Calculate the area in which signal time differences are physically allowed

    Parameters
    ----------
    ifos: list of strings
        list of interferometers

    Returns
    -------
    allowed_area: float
        area in units of seconds^(n_ifos-1) that coincident signals will occupy
    """
    n_ifos = len(ifos)

    if n_ifos == 2:
        det0 = pycbc.detector.Detector(ifos[0])
        det1 = pycbc.detector.Detector(ifos[1])
        allowed_area = 2 * det0.light_travel_time_to_detector(det1)
    elif n_ifos == 3:
        dets = {}
        tofs = numpy.zeros(n_ifos)
        ifo2_num = []
        # set up detector objects
        for ifo in ifos:
            dets[ifo] = pycbc.detector.Detector(ifo)

        # calculate travel time between detectors
        for i, ifo in enumerate(ifos):
            ifo2_num.append(int(numpy.mod(i + 1, n_ifos)))
            det0 = dets[ifo]
            det1 = dets[ifos[ifo2_num[i]]]
            tofs[i] = det0.light_travel_time_to_detector(det1)

        # calculate allowed area
        phi_12 = numpy.arccos((tofs[0]**2 + tofs[1]**2 - tofs[2]**2)
                              / (2 * tofs[0] * tofs[1]))
        allowed_area = numpy.pi * tofs[0] * tofs[1] * numpy.sin(phi_12)
    else:
        raise NotImplementedError("Not able to deal with more than 3 ifos")

    return allowed_area
