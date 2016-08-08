"""
Tools for maximum likelihood fits to single trigger statistic values

For some set of values above a threshold, e.g. trigger SNRs, the functions
in this module perform maximum likelihood fits with 1-sigma uncertainties
to various simple functional forms of PDF, all normalized to 1.
You can also obtain the fitted function and its (inverse) CDF and perform
a Kolmogorov-Smirnov test.

Usage:
# call the fit function directly if the threshold is known
alpha, sigma_alpha = fit_exponential(snrs, 5.5)

# apply a threshold explicitly
alpha, sigma_alpha = fit_above_thresh('exponential', snrs, thresh=6.25)

# let the code work out the threshold from the smallest value via the default thresh=None
alpha, sigma_alpha = fit_above_thresh('exponential', snrs)

# or only fit the largest N values, i.e. tail fitting
thresh = tail_threshold(snrs, N=500)
alpha, sigma_alpha = fit_above_thresh('exponential', snrs, thresh)

# obtain the fitted function directly
xvals = numpy.xrange(5.5, 10.5, 20)
exponential_fit = expfit(xvals, alpha, thresh)

# or access function by name
exponential_fit_1 = fit_fn('exponential', xvals, alpha, thresh)

# get the KS test statistic and p-value - see scipy.stats.kstest
ks_stat, ks_pval = KS_test('exponential', snrs, alpha, thresh)

"""

# Copyright T. Dent 2015 (thomas.dent@aei.mpg.de)
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.

from __future__ import division
import numpy
from scipy.stats import kstest

fitalpha_dict = {
    'exponential' : lambda vals, thresh : 1. / (numpy.mean(vals) - thresh),
    'rayleigh'    : lambda vals, thresh : 2. / (numpy.mean(vals**2.) - thresh**2.),
    'power'       : lambda vals, thresh : numpy.mean(numpy.log(vals/thresh))**-1. + 1.
}

# measurement standard deviation = (-d^2 log L/d alpha^2)^(-1/2)
fitstd_dict = {
    'exponential' : lambda vals, alpha : alpha / len(vals)**0.5,
    'rayleigh'    : lambda vals, alpha : alpha / len(vals)**0.5,
    'power'       : lambda vals, alpha : (alpha - 1.) / len(vals)**0.5
}

def fit_above_thresh(distr, vals, thresh=None):
    """
    Maximum likelihood fit for the coefficient alpha

    Fitting a distribution of discrete values above a given threshold.
    Exponential  p(x) = alpha exp(-alpha (x-x_t))
    Rayleigh     p(x) = alpha x exp(-alpha (x**2-x_t**2)/2)
    Power        p(x) = ((alpha-1)/x_t) (x/x_t)**-alpha
    Values below threshold will be discarded.
    If no threshold is specified the minimum sample value will be used.

    Parameters
    ----------
    distr : {'exponential', 'rayleigh', 'power'}
        Name of distribution
    vals : sequence of floats
        Values to fit
    thresh : float
        Threshold to apply before fitting; if None, use min(vals)

    Returns
    -------
    alpha : float
        Fitted value
    sigma_alpha : float
        Standard error in fitted value
    """
    vals = numpy.array(vals)
    if thresh is None:
        thresh = min(vals)
    else:
        vals = vals[vals >= thresh]
    alpha = fitalpha_dict[distr](vals, thresh)
    return alpha, fitstd_dict[distr](vals, alpha)


fitfn_dict = {
    'exponential' : lambda x, alpha, t : alpha * numpy.exp(-alpha * (x - t)),
    'rayleigh' : lambda x, alpha, t : alpha * x * \
                                      numpy.exp(-alpha * (x**2. - t**2.) / 2.),
    'power' : lambda x, alpha, t : (alpha - 1.) * x**(-alpha) * t**(alpha - 1.)
}

def fit_fn(distr, xvals, alpha, thresh):
    """
    The fitted function normalized to 1 above threshold

    To normalize to a given total count multiply by the count.

    Parameters
    ----------
    xvals : sequence of floats
        Values where the function is to be evaluated
    alpha : float
        The fitted parameter
    thresh : float
        Threshold value applied to fitted values

    Returns
    -------
    fit : array of floats
        Fitted function at the requested xvals
    """
    xvals = numpy.array(xvals)
    fit = fitfn_dict[distr](xvals, alpha, thresh)
    # set fitted values below threshold to 0
    numpy.putmask(fit, xvals < thresh, 0.)
    return fit


cum_fndict = {
    'exponential' : lambda x, alpha, t : numpy.exp(-alpha * (x - t)),
    'rayleigh' : lambda x, alpha, t : numpy.exp(-alpha * (x**2. - t**2.) / 2.),
    'power' : lambda x, alpha, t : x**(1. - alpha) * t**(alpha - 1.)
}

def cum_fit(distr, xvals, alpha, thresh):
    """
    Integral of the fitted function above a given value (reverse CDF)

    The fitted function is normalized to 1 above threshold

    Parameters
    ----------
    xvals : sequence of floats
        Values where the function is to be evaluated
    alpha : float
        The fitted parameter
    thresh : float
        Threshold value applied to fitted values

    Returns
    -------
    cum_fit : array of floats
        Reverse CDF of fitted function at the requested xvals
    """
    xvals = numpy.array(xvals)
    cum_fit = cum_fndict[distr](xvals, alpha, thresh)
    # set fitted values below threshold to 0
    numpy.putmask(cum_fit, xvals < thresh, 0.)
    return cum_fit


def tail_threshold(vals, N=1000):
    """Determine a threshold above which there are N louder values"""
    vals = numpy.array(vals)
    if len(vals) < N:
        raise RuntimeError('Not enough input values to determine threshold')
    vals.sort()
    return min(vals[-N:])


def KS_test(distr, vals, alpha, thresh=None):
    """
    Perform Kolmogorov-Smirnov test for fitted distribution

    Compare the given set of discrete values above a given threshold to the
    fitted distribution function.
    If no threshold is specified, the minimum sample value will be used.
    Returns the KS test statistic and its p-value: lower p means less
    probable under the hypothesis of a perfect fit

    Parameters
    ----------
    distr : {'exponential', 'rayleigh', 'power'}
        Name of distribution
    vals : sequence of floats
        Values to compare to fit
    alpha :
        Fitted distribution parameter
    thresh : float
        Threshold to apply before fitting; if None, use min(vals)

    Returns
    -------
    D : float
        KS test statistic
    p-value : float
        p-value, assumed to be two-tailed
    """
    vals = numpy.array(vals)
    if thresh is None:
        thresh = min(vals)
    else:
        vals = vals[vals >= thresh]
    def cdf_fn(x):
        return 1 - cum_fndict[distr](x, alpha, thresh)
    return kstest(vals, cdf_fn)

