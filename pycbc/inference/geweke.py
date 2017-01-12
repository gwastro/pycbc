""" Functions for computing the Geweke convergence statistic.
"""

import numpy

def geweke(x, num_samples, last_start, last, intervals, stride):
    """ Calculates Geweke.

    x : numpy.array
        A one-dimensional array of data.
    num_samples : int
        Number of samples to use for each Geweke calculation.
    last_start : int
        Index of last start.
    last : int
        Index of beginning of end segment.
    intervals : int
        Number of intervals to split data into.
    stride : int
        How far of a step to take for each point.
    """
    stats = []

    end = len(x) - 1

    step = int(last / (intervals - 1))

    starts = numpy.arange(0, last_start, stride)
    ends = []

    for start in starts:

        x_start_end = int(start + num_samples)
        x_end_start = int(end - last)

        x_start = x[start:x_start_end]
        x_end = x[x_end_start:]

        print start, x_start_end, x_end_start

        stats.append((x_start.mean() - x_end.mean())
                     / numpy.sqrt(x_start.var() + x_end.var()))
        ends.append(x_start_end)

    return numpy.array(starts), numpy.array(ends), numpy.array(stats)
