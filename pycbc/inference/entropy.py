""" The module contains functions for calculating the
Kullback-Leibler divergence.
"""

import numpy

from scipy import stats


def kl(samples1, samples2, pdf1=False, pdf2=False,
       bins=30, hist_min=None, hist_max=None):
    """ Computes the Kullback-Leibler divergence for a single parameter
    from two distributions.

    Parameters
    ----------
    samples1 : numpy.array
        Samples or probability density function (must also set `pdf1=True`).
    samples2 : numpy.array
        Samples or probability density function (must also set `pdf2=True`).
    pdf1 : bool
        Set to `True` if `samples1` is a probability density funtion already.
    pdf2 : bool
        Set to `True` if `samples2` is a probability density funtion already.
    bins : int
        Number of bins to use when calculating probability density function
        from a set of samples of the distribution.
    hist_min : numpy.float64
        Minimum of the distributions' values to use.
    hist_max : numpy.float64
        Maximum of the distributions' values to use.

    Returns
    -------
    numpy.float64
        The Kullback-Leibler divergence value.
    """
    hist_range = (hist_min, hist_max)
    if not pdf1:
        samples1, _ = numpy.histogram(samples1, bins=bins,
                                      range=hist_range, normed=True)
    if not pdf2:
        samples2, _ = numpy.histogram(samples2, bins=bins,
                                      range=hist_range, normed=True)
    return stats.entropy(samples1, qk=samples2)
