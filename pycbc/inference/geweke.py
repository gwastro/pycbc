""" Functions for computing the Geweke convergence statistic.
"""

import numpy

def geweke(x, seg_length, seg_stride, end_idx, ref_start,
           ref_end=None, seg_start=0):
    """ Calculates Geweke conervergence statistic for a chain of data.
    This function will advance along the chain and calculate the
    statistic for each step.

    x : numpy.array
        A one-dimensional array of data.
    seg_length : int
        Number of samples to use for each Geweke calculation.
    seg_stride : int
        Number of samples to advance before next Geweke calculation.
    end_idx : int
        Index of last start.
    ref_start : int
        Index of beginning of end reference segment.
    ref_end : int
        Index of end of end reference segment. Default is None which
        will go to the end of the data array.
    seg_start : int
        What index to start computing the statistic. Default is 0 which
        will go to the beginning of the data array.
    """

    # lists to hold statistic and end index
    stats = []
    ends = []

    # get the beginning of all segments
    starts = numpy.arange(seg_start, end_idx, seg_stride)

    # get second segment of data at the end to compare
    x_end = x[ref_start:ref_end]

    # loop over all segments
    for start in starts:

        # find the end of the first segment
        x_start_end = int(start + seg_length)

        # get first segment
        x_start = x[start:x_start_end]

        # compute statistic
        stats.append((x_start.mean() - x_end.mean())
                     / numpy.sqrt(x_start.var() + x_end.var()))

        # store end of first segment
        ends.append(x_start_end)

    return numpy.array(starts), numpy.array(ends), numpy.array(stats)
