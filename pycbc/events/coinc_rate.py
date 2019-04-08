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
    assert len(rates)==len(ifos)
    if len(ifos) == 2:
        # Use 2-ifo coincidence calculation
        # multiply the two rates and by the overlap time
        det0, det1 = detector.Detector(ifos[0]), detector.Detector(ifos[1])
        time_window = det0.light_travel_time_to_detector(det1) + slop
        expected_coinc_rates[''.join(ifos)] = 2 * time_window * \
                                              numpy.multiply(rates[0], rates[1])
    elif len(ifos) == 3:
        # Calculate rate for each two-detector combination
        subsets = itertools.combinations(ifos,2)
        for subset in subsets:
            i_set = [numpy.nonzero(ifo==ifos)[0][0] for ifo in subset]
            ifostring = ''.join(ifos[i_set])
            expected_coinc_rates[ifostring] = \
                multiifo_noise_coinc_rate(rates[i_set],
                                          ifos[i_set], slop)[ifostring]
        # Use 3-ifo coincidence calculation                                                       
        # calculate three-detector coincident rate
        dets = {}
        tofs = numpy.zeros(len(ifos))
        for ifo in ifos:
            dets[ifo] = detector.Detector(ifo)

        for i in range(0,len(ifos)):
            det0 = dets[ifos[i]]
            det1 = dets[ifos[numpy.mod(i+1,len(ifos))]]
            tofs[i] = abs(det0.light_travel_time_to_detector(det1)) + slop

        allowed_area = 0
        for i in range(0,len(ifos)):
            allowed_area += 2*tofs[i]*tofs[numpy.mod(i+1,len(ifos))] - tofs[i]**2

        ifostring = ''.join(ifos)
        rateprod = [r1*r2*r3 for r1,r2,r3 in zip(*rates)]
        expected_coinc_rates[ifostring] = allowed_area * numpy.array(rateprod)

    return expected_coinc_rates



class MultiRingBuffer(object):
    """Dynamic size n-dimensional ring buffer that can expire elements."""

    def __init__(self, num_rings, max_time, dtype):
        """
        Parameters
        ----------
        num_rings: int
            The number of ring buffers to create. They all will have the same
            intrinsic size and will expire at the same time.
        max_time: int
            The maximum "time" an element can exist in each ring.
        dtype: numpy.dtype
            The type of each element in the ring buffer.
        """
        self.max_time = max_time
        self.buffer = []
        self.buffer_expire = []
        for _ in range(num_rings):
            self.buffer.append(numpy.zeros(0, dtype=dtype))
            self.buffer_expire.append(numpy.zeros(0, dtype=int))
        self.time = 0

    @property
    def filled_time(self):
        return min(self.time, self.max_time)

    def num_elements(self):
        return sum([len(a) for a in self.buffer])

    @property
    def nbytes(self):
        return sum([a.nbytes for a in self.buffer])

    def discard_last(self, indices):
        """Discard the triggers added in the latest update"""
        for i in indices:
            self.buffer_expire[i] = self.buffer_expire[i][:-1]
            self.buffer[i] = self.buffer[i][:-1]

    def advance_time(self):
        """Advance the internal time increment by 1, expiring any triggers that
        are now too old.
        """
        self.time += 1

        expired = self.time - self.max_time
        for j, exp in enumerate(self.buffer_expire):
            if (len(exp) > 0) and (exp[0] < expired):
                self.buffer_expire[j] = exp[1:].copy()
                self.buffer[j] = self.buffer[j][1:].copy()

    def add(self, indices, values):
        """Add triggers in 'values' to the buffers indicated by the indices
        """
        for i, v in zip(indices, values):
            self.buffer[i] = numpy.append(self.buffer[i], v)
            self.buffer_expire[i] = numpy.append(self.buffer_expire[i], self.time)
        self.advance_time()

    def expire_vector(self, buffer_index):
        """Return the expiration vector of a given ring buffer """
        return self.buffer_expire[buffer_index]

    def data(self, buffer_index):
        """Return the data vector for a given ring buffer"""
        return self.buffer[buffer_index]


