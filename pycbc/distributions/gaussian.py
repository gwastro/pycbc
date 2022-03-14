# Copyright (C) 2016  Christopher M. Biwer, Collin Capano
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
This modules provides classes for evaluating Gaussian distributions.
"""

import numpy
from scipy.special import erf, erfinv
import scipy.stats
from pycbc.distributions import bounded

class Gaussian(bounded.BoundedDist):
    r"""A Gaussian distribution on the given parameters; the parameters are
    independent of each other.

    Bounds can be provided on each parameter, in which case the distribution
    will be a truncated Gaussian distribution.  The PDF of a truncated
    Gaussian distribution is given by:

    .. math::
        p(x|a, b, \mu,\sigma) = \frac{1}{\sqrt{2 \pi \sigma^2}}\frac{e^{- \frac{\left( x - \mu \right)^2}{2 \sigma^2}}}{\Phi(b|\mu, \sigma) - \Phi(a|\mu, \sigma)},

    where :math:`\mu` is the mean, :math:`\sigma^2` is the variance,
    :math:`a,b` are the bounds, and :math:`\Phi` is the cumulative distribution
    of an unbounded normal distribution, given by:

    .. math::
        \Phi(x|\mu, \sigma) = \frac{1}{2}\left[1 + \mathrm{erf}\left(\frac{x-\mu}{\sigma \sqrt{2}}\right)\right].

    Note that if :math:`[a,b) = [-\infty, \infty)`, this reduces to a standard
    Gaussian distribution.


    Instances of this class can be called like a function. By default, logpdf
    will be called, but this can be changed by setting the class's __call__
    method to its pdf method.

    Parameters
    ----------
    \**params :
        The keyword arguments should provide the names of parameters and
        (optionally) some bounds, as either a tuple or a
        `boundaries.Bounds` instance. The mean and variance of each
        parameter can be provided by additional keyword arguments that have
        `_mean` and `_var` adding to the parameter name. For example,
        `foo=(-2,10), foo_mean=3, foo_var=2` would create a truncated Gaussian
        with mean 3 and variance 2, bounded between :math:`[-2, 10)`. If no
        mean or variance is provided, the distribution will have 0 mean and
        unit variance. If None is provided for the bounds, the distribution
        will be a normal, unbounded Gaussian (equivalent to setting the bounds
        to `[-inf, inf)`).

    Examples
    --------
    Create an unbounded Gaussian distribution with zero mean and unit variance:
    >>> dist = distributions.Gaussian(mass1=None)

    Create a bounded Gaussian distribution on :math:`[1,10)` with a mean of 3
    and a variance of 2:
    >>> dist = distributions.Gaussian(mass1=(1,10), mass1_mean=3, mass1_var=2)

    Create a bounded Gaussian distribution with the same parameters, but with
    cyclic boundary conditions:
    >>> dist = distributions.Gaussian(mass1=Bounds(1,10, cyclic=True), mass1_mean=3, mass1_var=2)
    """
    name = "gaussian"

    def __init__(self, **params):

        # save distribution parameters as dict
        # calculate the norm and exponential norm ahead of time
        # and save to self._norm, self._lognorm, and self._expnorm
        self._bounds = {}
        self._mean = {}
        self._var = {}
        self._norm = {}
        self._lognorm = {}
        self._expnorm = {}
        # pull out specified means, variance
        mean_args = [p for p in params if p.endswith('_mean')]
        var_args = [p for p in params if p.endswith('_var')]
        self._mean = dict([[p[:-5], params.pop(p)] for p in mean_args])
        self._var = dict([[p[:-4], params.pop(p)] for p in var_args])
        # initialize the bounds
        super(Gaussian, self).__init__(**params)

        # check that there are no params in mean/var that are not in params
        missing = set(self._mean.keys()) - set(params.keys())
        if any(missing):
            raise ValueError("means provided for unknow params {}".format(
                ', '.join(missing)))
        missing = set(self._var.keys()) - set(params.keys())
        if any(missing):
            raise ValueError("vars provided for unknow params {}".format(
                ', '.join(missing)))
        # set default mean/var for params not specified
        self._mean.update(dict([[p, 0.]
            for p in params if p not in self._mean]))
        self._var.update(dict([[p, 1.]
            for p in params if p not in self._var]))

        # compute norms
        for p,bnds in self._bounds.items():
            sigmasq = self._var[p]
            mu = self._mean[p]
            a,b = bnds
            invnorm = scipy.stats.norm.cdf(b, loc=mu, scale=sigmasq**0.5) \
                    - scipy.stats.norm.cdf(a, loc=mu, scale=sigmasq**0.5)
            invnorm *= numpy.sqrt(2*numpy.pi*sigmasq)
            self._norm[p] = 1./invnorm
            self._lognorm[p] = numpy.log(self._norm[p])
            self._expnorm[p] = -1./(2*sigmasq)


    @property
    def mean(self):
        return self._mean


    @property
    def var(self):
        return self._var

    def _normalcdf(self, param, value):
        """The CDF of the normal distribution, without bounds."""
        mu = self._mean[param]
        var = self._var[param]
        return 0.5*(1. + erf((value - mu)/(2*var)**0.5))

    def cdf(self, param, value):
        """Returns the CDF of the given parameter value."""
        a, b = self._bounds[param]
        if a != -numpy.inf:
            phi_a = self._normalcdf(param, a)
        else:
            phi_a = 0.
        if b != numpy.inf:
            phi_b = self._normalcdf(param, b)
        else:
            phi_b = 1.
        phi_x = self._normalcdf(param, value)
        return (phi_x - phi_a)/(phi_b - phi_a)

    def _normalcdfinv(self, param, p):
        """The inverse CDF of the normal distribution, without bounds."""
        mu = self._mean[param]
        var = self._var[param]
        return mu + (2*var)**0.5 * erfinv(2*p - 1.)

    def _cdfinv_param(self, param, p):
        """Return inverse of the CDF.
        """
        a, b = self._bounds[param]
        if a != -numpy.inf:
            phi_a = self._normalcdf(param, a)
        else:
            phi_a = 0.
        if b != numpy.inf:
            phi_b = self._normalcdf(param, b)
        else:
            phi_b = 1.
        adjusted_p = phi_a + p * (phi_b - phi_a)
        return self._normalcdfinv(param, adjusted_p)

    def _pdf(self, **kwargs):
        """Returns the pdf at the given values. The keyword arguments must
        contain all of parameters in self's params. Unrecognized arguments are
        ignored.
        """
        return numpy.exp(self._logpdf(**kwargs))


    def _logpdf(self, **kwargs):
        """Returns the log of the pdf at the given values. The keyword
        arguments must contain all of parameters in self's params. Unrecognized
        arguments are ignored.
        """
        if kwargs in self:
            return sum([self._lognorm[p] +
                        self._expnorm[p]*(kwargs[p]-self._mean[p])**2.
                        for p in self._params])
        else:
            return -numpy.inf

    @classmethod
    def from_config(cls, cp, section, variable_args):
        """Returns a Gaussian distribution based on a configuration file. The
        parameters for the distribution are retrieved from the section titled
        "[`section`-`variable_args`]" in the config file.

        Boundary arguments should be provided in the same way as described in
        `get_param_bounds_from_config`. In addition, the mean and variance of
        each parameter can be specified by setting `{param}_mean` and
        `{param}_var`, respectively. For example, the following would create a
        truncated Gaussian distribution between 0 and 6.28 for a parameter
        called `phi` with mean 3.14 and variance 0.5 that is cyclic:

        .. code-block:: ini

            [{section}-{tag}]
            min-phi = 0
            max-phi = 6.28
            phi_mean = 3.14
            phi_var = 0.5
            cyclic =

        Parameters
        ----------
        cp : pycbc.workflow.WorkflowConfigParser
            A parsed configuration file that contains the distribution
            options.
        section : str
            Name of the section in the configuration file.
        variable_args : str
            The names of the parameters for this distribution, separated by
            `prior.VARARGS_DELIM`. These must appear in the "tag" part
            of the section header.

        Returns
        -------
        Gaussian
            A distribution instance from the pycbc.inference.prior module.
        """
        return bounded.bounded_from_config(cls, cp, section, variable_args,
                                                  bounds_required=False)


__all__ = ['Gaussian']
