#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
""" This modules contains functions for calculating and manipulating
coincident triggers.
"""

import numpy
import itertools
import pycbc.detector


def multiifo_noise_coinc_rate(rates, ifos, slop):
    """
    Calculate the expected rate of coincidences for multiple detectors
    """
    ifos = numpy.array(ifos)
    rates = numpy.array(rates)
    expected_coinc_rates = {}
    n_ifos = len(ifos)
    assert len(rates) == n_ifos
    # Calculate coincidence for all-ifo combination
    # multiply the two rates and by the overlap time
    allowed_area = multiifo_noise_coincident_area(ifos, slop)
    rateprod = [numpy.prod(rs) for rs in zip(*rates)]
    ifostring = ''.join(ifos)
    expected_coinc_rates[ifostring] = allowed_area * numpy.array(rateprod)
    # if more than one possible coicidence type exists,
    # calculate coincidences for subsets
    if n_ifos > 2:
        # Calculate rate for each 'miss-one-out' detector combination
        subsets = itertools.combinations(ifos, n_ifos-1)
        for subset in subsets:
            i_set = [numpy.nonzero(ifo == ifos)[0][0] for ifo in subset]
            # calculate coincidence rates in subsets through iteration
            sub_coinc_rates = multiifo_noise_coinc_rate(rates[i_set],
                                                        ifos[i_set], slop)
            # add these sub-coincidences to the overall dictionary
            for sub_coinc in sub_coinc_rates:
                expected_coinc_rates[sub_coinc] = sub_coinc_rates[sub_coinc]

    return expected_coinc_rates


def multiifo_noise_coincident_area(ifos, slop):
    """
    calculate the multiplicative factor of the individual detector trigger
    rates in order to calculate combined rates
    """
    # TO DO: add in capability for more than 3 detectors
    n_ifos = len(ifos)
    if n_ifos == 2:
        det0 = pycbc.detector.Detector(ifos[0])
        det1 = pycbc.detector.Detector(ifos[1])
        allowed_area = 2*(det0.light_travel_time_to_detector(det1) + slop)
    elif n_ifos == 3:
        dets = {}
        tofs = numpy.zeros(n_ifos)
        ifo2_num = []
        # set up detector objects
        for ifo in ifos:
            dets[ifo] = pycbc.detector.Detector(ifo)

        # calculate travel time between detectors (plus extra for timing error)
        # TO DO: allow for different timing errors between different detectors
        for i, ifo in enumerate(ifos):
            ifo2_num.append(int(numpy.mod(i+1, n_ifos)))
            det0 = dets[ifo]
            det1 = dets[ifos[ifo2_num[i]]]
            tofs[i] = det0.light_travel_time_to_detector(det1) + slop

        # combine these to calculate allowed area
        allowed_area = 0
        for i, _ in enumerate(ifos):
            allowed_area += 2*tofs[i]*tofs[ifo2_num[i]] - tofs[i]**2

    return allowed_area

