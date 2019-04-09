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

import h5py, numpy, logging, pycbc.pnutils, copy, lal, itertools
import pycbc.detector as detector

def multiifo_noise_coinc_rate(rates, ifos, slop):
    """
    Calculate the expected rate of coincidences for multiple detectors
    """
    ifos = numpy.array(ifos)
    rates = numpy.array(rates)
    expected_coinc_rates = {}
    n_ifos = len(ifos)
    assert len(rates)==n_ifos
    # Calculate coincidence for all-ifo combination
    # multiply the two rates and by the overlap time
    allowed_area = multiifo_noise_coincident_area(ifos, slop)
    rateprod = [numpy.prod(rs) for rs in zip(*rates)]
    ifostring = ''.join(ifos)
    expected_coinc_rates[ifostring] = allowed_area * numpy.array(rateprod)
    # if more than one possible coicidence type exists, calculate coincidences for subsets
    if n_ifos > 2:
        # Calculate rate for each 'miss-one-out' detector combination
        subsets = itertools.combinations(ifos, n_ifos-1)
        for subset in subsets:
            i_set = [numpy.nonzero(ifo==ifos)[0][0] for ifo in subset]
            ifostring = ''.join(ifos[i_set])
            # calculate coincidence rates in subsets through iteration
            sub_coincidences = \
                multiifo_noise_coinc_rate(rates[i_set],
                                          ifos[i_set], slop)
            #add these sub-coincidences to the overall dictionary
            for ct in sub_coincidences:
                expected_coinc_rates[ct] = sub_coincidences[ct]

    return expected_coinc_rates

def multiifo_noise_coincident_area(ifos, slop):
    """
    calculate the multiplicative factor of the individual detector trigger
    rates in order to calculate combined rates
    """
    #TODO: add in capability for more than 3 detectors
    if len(ifos) == 2:
        det0, det1 = detector.Detector(ifos[0]), detector.Detector(ifos[1])
        allowed_area = 2*(det0.light_travel_time_to_detector(det1) + slop)
    elif len(ifos) ==3:
        dets = {}
        tofs = numpy.zeros(len(ifos))
        # set up detector objects
        for ifo in ifos:
            dets[ifo] = detector.Detector(ifo)

        # calculate travel time between detectors (plus extra for timing error)
        # TODO: allow for different timing errors between different sets of detectors
        for i in range(0,len(ifos)):
            det0 = dets[ifos[i]]
            det1 = dets[ifos[numpy.mod(i+1,len(ifos))]]
            tofs[i] = det0.light_travel_time_to_detector(det1) + slop

        # combine these to calculate allowed area
        allowed_area = 0
        for i in range(0,len(ifos)):
            allowed_area += 2*tofs[i]*tofs[numpy.mod(i+1,len(ifos))] - tofs[i]**2
    return allowed_area
