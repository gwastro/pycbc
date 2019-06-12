""" The module contains functions for calculating the
Kullback-Leibler divergence.
"""

import numpy
from scipy import stats


def check_hist_params(hist_min, hist_max, hist_bins):
    """ Checks that the bound values given for the histogram are consistent,
    returning the range if they are or raising an error if they are not.
    Also checks that if hist_bins is a str, it corresponds to a method
    available in numpy.histogram

    Parameters
    ----------
    hist_min : numpy.float64
        Minimum value for the histogram.
    hist_max : numpy.float64
        Maximum value for the histogram.
    hist_bins: int or str
        If int, number of equal-width bins to use in numpy.histogram. If str,
        it should be one of the methods to calculate the optimal bin width
        available in numpy.histogram: ['auto', 'fd', 'doane', 'scott', 'stone',
        'rice', 'sturges', 'sqrt']. Default is 'fd' (Freedman Diaconis
        Estimator). This option will be ignored if `kde=True`.

    Returns
    -------
    range : tuple
        The bounds (hist_min, hist_max).
    """

    hist_methods = ['auto', 'fd', 'doane', 'scott', 'stone', 'rice',
                     'sturges', 'sqrt']
    if not hist_bins:
        hist_bins = 'fd'
    elif isinstance(hist_bins, str) and hist_bins not in hist_methods:
        raise ValueError('Method for calculating bins width must be one of'
                         ' {}'.format(hist_methods))

    # One of the bounds is missing
    if hist_min and not hist_max:
        raise ValueError('hist_min provided but hist_max missing.')
    elif hist_max and not hist_min:
        raise ValueError('hist_max provided but hist_min missing.')

    # Both bounds given
    if hist_min and hist_max:
        if hist_min >= hist_max:
            raise ValueError('hist_min must be lower than hist_max.')
        else:
            hist_range = (hist_min, hist_max)
    # No bounds given
    elif not hist_min and not hist_max:
        hist_range = None

    return hist_range, hist_bins

def entropy(pdf1, base=numpy.e):
    """ Computes the information entropy for a single parameter
    from one probability density function.

    Parameters
    ----------
    pdf1 : numpy.array
        Probability density function.
    base : {numpy.e, numpy.float64}, optional
        The logarithmic base to use (choose base 2 for information measured
        in bits, default is nats).

    Returns
    -------
    numpy.float64
        The information entropy value.
    """

    return stats.entropy(pdf1, base=base)


def kl(samples1, samples2, pdf1=False, pdf2=False, kde=False,
       bins=None, hist_min=None, hist_max=None, base=numpy.e):
    """ Computes the Kullback-Leibler divergence for a single parameter
    from two distributions.

    Parameters
    ----------
    samples1 : numpy.array
        Samples or probability density function (for the latter must also set
        `pdf1=True`).
    samples2 : numpy.array
        Samples or probability density function (for the latter must also set
        `pdf2=True`).
    pdf1 : bool
        Set to `True` if `samples1` is a probability density funtion already.
    pdf2 : bool
        Set to `True` if `samples2` is a probability density funtion already.
    kde : bool
        Set to `True` if at least one of `pdf1` or `pdf2` is `False` to
        estimate the probability density function using kernel density
        estimation (KDE).
    bins : int or str, optional
        If int, number of equal-width bins to use when calculating probability
        density function from a set of samples of the distribution. If str, it
        should be one of the methods to calculate the optimal bin width
        available in numpy.histogram: ['auto', 'fd', 'doane', 'scott', 'stone',
        'rice', 'sturges', 'sqrt']. Default is 'fd' (Freedman Diaconis
        Estimator). This option will be ignored if `kde=True`.
    hist_min : numpy.float64
        Minimum of the distributions' values to use. This will be ignored if
        `kde=True`.
    hist_max : numpy.float64
        Maximum of the distributions' values to use. This will be ignored if
        `kde=True`.
    base : numpy.float64
        The logarithmic base to use (choose base 2 for information measured
        in bits, default is nats).

    Returns
    -------
    numpy.float64
        The Kullback-Leibler divergence value.
    """
    if pdf1 and pdf2 and kde:
        raise ValueError('KDE can only be used when at least one of pdf1 or '
                         'pdf2 is False.')

    pdfs = {}
    for n, (samples, pdf) in enumerate(((samples1, pdf1), (samples2, pdf2))):
        if pdf:
            pdfs[n] = samples
        elif kde:
            samples_kde = stats.gaussian_kde(samples)
            pdfs[n] = samples_kde.evaluate(samples)
        else:
            hist_range, hist_bins = check_hist_params(hist_min, hist_max, bins)
            pdfs[n], _ = numpy.histogram(samples, bins=hist_bins,
                                         range=hist_range, normed=True)

    return stats.entropy(pdfs[0], qk=pdfs[1], base=base)


def js(samples1, samples2, kde=False, bins=None, hist_min=None, hist_max=None,
       base=numpy.e):
    """ Computes the Jensen-Shannon divergence for a single parameter
    from two distributions.

    Parameters
    ----------
    samples1 : numpy.array
        Samples.
    samples2 : numpy.array
        Samples.
    kde : bool
        Set to `True` to estimate the probability density function using
        kernel density estimation (KDE).
    bins : int or str, optional
        If int, number of equal-width bins to use when calculating probability
        density function from a set of samples of the distribution. If str, it
        should be one of the methods to calculate the optimal bin width
        available in numpy.histogram: ['auto', 'fd', 'doane', 'scott', 'stone',
        'rice', 'sturges', 'sqrt']. Default is 'fd' (Freedman Diaconis
        Estimator). This option will be ignored if `kde=True`.
    hist_min : numpy.float64
        Minimum of the distributions' values to use. This will be ignored if
        `kde=True`.
    hist_max : numpy.float64
        Maximum of the distributions' values to use. This will be ignored if
        `kde=True`.
    base : numpy.float64
        The logarithmic base to use (choose base 2 for information measured
        in bits, default is nats).

    Returns
    -------
    numpy.float64
        The Jensen-Shannon divergence value.
    """
    join_samples = numpy.concatenate((samples1, samples2))
    if kde:
        samplesm_kde = stats.gaussian_kde(join_samples)
        samplesm = samplesm_kde.evaluate(join_samples)
    else:
        hist_range, hist_bins = check_hist_params(hist_min, hist_max, bins)
        samplesm, _ = numpy.histogram(join_samples, bins=hist_bins,
                                      range=hist_range, normed=True)
    samplesm = (1./2) * samplesm

    js_div = 0
    for samples in (samples1, samples2):
        js_div += (1./2) * kl(samples, samplesm, pdf1=False, pdf2=True,
                              kde=kde, bins=bins, hist_min=hist_min,
                              hist_max=hist_max, base=base)

    return js_div
