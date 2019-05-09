# Copyright (C) 2019 Steven Reyes
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""
This modules provides functions for estimating the marginal
likelihood or evidence of a model.
"""
import numpy
from scipy import integrate


def arithmetic_mean_estimator(log_likelihood):
    """Returns the log evidence via the prior arithmetic mean estimator (AME).

    The logarithm form of AME is used. This is the most basic
    evidence estimator, and often requires O(billions) of samples
    from the prior.

    Parameters
    ----------
    log_likelihood : 1d array of floats
        The log likelihood of the data sampled from the prior
        distribution.

    Returns
    -------
    float :
        Estimation of the log of the evidence.
    """
    num_samples = len(log_likelihood)
    logl_max = numpy.max(log_likelihood)

    log_evidence = 0.
    for i, _ in enumerate(log_likelihood):
        log_evidence += numpy.exp(log_likelihood[i] - logl_max)

    log_evidence = numpy.log(log_evidence)
    log_evidence += logl_max - numpy.log(num_samples)

    return log_evidence


def harmonic_mean_estimator(log_likelihood):
    """Returns the log evidence via posterior harmonic mean estimator (HME).

    The logarithm form of HME is used. This method is not
    recommended for general use. It is very slow to converge,
    formally, has infinite variance, and very error prone.

    Not recommended for general use.

    Parameters
    ----------
    log_likelihood : 1d array of floats
        The log likelihood of the data sampled from the posterior
        distribution.

    Returns
    -------
    float :
        Estimation of the log of the evidence.
    """
    num_samples = len(log_likelihood)
    logl_max = numpy.max(-1.0*log_likelihood)

    log_evidence = 0.
    for i, _ in enumerate(log_likelihood):
        log_evidence += numpy.exp(-1.0*log_likelihood[i] + logl_max)

    log_evidence = -1.0*numpy.log(log_evidence)
    log_evidence += logl_max
    log_evidence += numpy.log(num_samples)

    return log_evidence


def thermodynamic_integration(log_likelihood, betas,
                              method="simpsons"):
    """Returns the log evidence of the model via thermodynamic integration.
    Also returns an estimated standard deviation for the log evidence.

    Current options are integration through the trapezoid rule, a
    first-order corrected trapezoid rule, and Simpson's rule.

    Parameters
    ----------
    log_likelihood : 3d array of shape (betas, walker, iteration)
        The log likelihood for each temperature separated by
        temperature, walker, and iteration.

    betas : 1d array
        The inverse temperatures used in the MCMC.

    method : {"trapzoid", "trapezoid_corrected", "simpsons"},
             optional.
        The numerical integration method to use for the
        thermodynamic integration. Choices include: "trapezoid",
        "trapezoid_corrected", "simpsons", for the trapezoid rule,
        the first-order correction to the trapezoid rule, and
        Simpson's rule. [Default = "simpsons"]

    Returns
    -------
    log_evidence : float
        Estimation of the log of the evidence.

    mcmc_std : float
        The standard deviation of the log evidence estimate from
        Monte-Carlo spread.
    """
    # Check if the method of integration is in the list of choices
    method_list = ["trapezoid", "trapezoid_corrected", "simpsons"]

    if method not in method_list:
        raise ValueError("Method %s not supported. Expected %s"
                         % (method, method_list))
    # Read in the data and ensure ordering of data.
    # Ascending order sort
    order = numpy.argsort(betas)
    betas = betas[order]
    log_likelihood = log_likelihood[order]

    # Assume log likelihood is given in shape of beta, walker,
    # and iteration.
    log_likelihood = numpy.reshape(log_likelihood,
                                   (len(betas),
                                    len(log_likelihood[0].flatten())))

    average_logl = numpy.average(log_likelihood, axis=1)

    if method in ("trapezoid", "trapezoid_corrected"):
        log_evidence = numpy.trapz(average_logl, betas)

    if method == "trapezoid_corrected":
        # var_correction holds the derivative correction terms
        # See Friel et al. 2014 for expression and derivation.
        # https://link.springer.com/article/10.1007/s11222-013-9397-1
        var_correction = 0
        for i in range(len(betas) - 1):
            delta_beta = betas[i+1] - betas[i]
            pre_fac_var = (1. / 12.) * (delta_beta ** 2.0)
            var_diff = numpy.var(log_likelihood[i+1])
            var_diff -= numpy.var(log_likelihood[i])
            var_correction -= pre_fac_var * var_diff

        # Add the derivative correction term back to the log_evidence
        # from the first if statement.
        log_evidence += var_correction

    elif method == "simpsons":
        # beta -> 0 tends to contribute the least to the integral
        # so we can sacrifice precision there, rather than near
        # beta -> 1. Option even="last" puts trapezoid rule at
        # first few points.
        log_evidence = integrate.simps(average_logl, betas,
                                       even="last")

    # Estimate the Monte Carlo variance of the evidence calculation
    # See (Evans, Annis, 2019.)
    # https://www.sciencedirect.com/science/article/pii/S0022249617302651
    ti_vec = numpy.zeros(len(log_likelihood[0]))

    # Get log likelihood chains by sample and not by temperature.
    logl_per_samp = []
    for i, _ in enumerate(log_likelihood[0]):
        logl_per_samp.append([log_likelihood[x][i] for x in range(len(betas))])

    if method in ("trapezoid", "trapezoid_corrected"):
        for i, _ in enumerate(log_likelihood[0]):
            ti_vec[i] = numpy.trapz(logl_per_samp[i], betas)

    elif method == "simpsons":
        for i, _ in enumerate(log_likelihood[0]):
            ti_vec[i] = integrate.simps(logl_per_samp[i], betas,
                                        even="last")

    # Standard error is sample std / sqrt(number of samples)
    mcmc_std = numpy.std(ti_vec) / numpy.sqrt(float(len(log_likelihood[0])))

    return log_evidence, mcmc_std


def stepping_stone_algorithm(log_likelihood, betas):
    """Returns the log evidence of the model via stepping stone algorithm.
    Also returns an estimated standard deviation for the log evidence.

    Parameters
    ----------
    log_likelihood : 3d array of shape (betas, walker, iteration)
        The log likelihood for each temperature separated by
        temperature, walker, and iteration.

    betas          : 1d array
                     The inverse temperatures used in the MCMC.

    Returns
    -------
    log_evidence : float
        Estimation of the log of the evidence.
    mcmc_std : float
        The standard deviation of the log evidence estimate from
        Monte-Carlo spread.
    """
    # Reverse order sort
    order = numpy.argsort(betas)[::-1]
    betas = betas[order]
    log_likelihood = log_likelihood[order]

    # Assume log likelihood is given in shape of beta,
    # walker, iteration.
    log_likelihood = numpy.reshape(log_likelihood,
                                   (len(betas),
                                    len(log_likelihood[0].flatten())))

    log_rk_pb = numpy.zeros(len(betas) - 1)
    for i in range(len(betas) - 1):
        delta_beta = betas[i] - betas[i+1]
        # Max log likelihood for beta [i+1]
        max_logl_pb = numpy.max(log_likelihood[i+1])
        val_1 = delta_beta * max_logl_pb
        val_2 = delta_beta * (log_likelihood[i+1] - max_logl_pb)
        val_2 = numpy.log(numpy.average(numpy.exp(val_2)))
        log_rk_pb[i] = val_1 + val_2

    log_rk = numpy.sum(log_rk_pb)
    log_evidence = log_rk

    # Calculate the Monte Carlo variation
    mcmc_std = 0
    for i in range(len(betas) - 1):
        delta_beta = betas[i] - betas[i+1]
        pre_fact = (delta_beta * log_likelihood[i+1]) - log_rk_pb[i]
        pre_fact = numpy.exp(pre_fact) - 1.0
        val = numpy.sum(pre_fact ** 2)

        mcmc_std += val

    mcmc_std /= float(len(log_likelihood[0])) ** 2.0
    mcmc_std = numpy.sqrt(mcmc_std)

    return log_evidence, mcmc_std
