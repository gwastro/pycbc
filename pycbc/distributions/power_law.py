# Copyright (C) 2016 Christopher M. Biwer
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
This modules provides classes for evaluating distributions where the
probability density function is a power law.
"""

import numpy
from pycbc.distributions import bounded

class UniformPowerLaw(bounded.BoundedDist):
    r"""
    For a uniform distribution in power law. The parameters are
    independent of each other. Instances of this class can be called like
    a function. By default, logpdf will be called, but this can be changed
    by setting the class's __call__ method to its pdf method.

    The cumulative distribution function (CDF) will be the ratio of volumes:

    .. math::

        F(r) = \frac{V(r)}{V(R)}

    Where :math:`R` is the radius of the sphere. So we can write our
    probability density function (PDF) as:

    .. math::

        f(r) = c r^n

    For generality we use :math:`n` for the dimension of the volume element,
    eg. :math:`n=2` for a 3-dimensional sphere. And use
    :math:`c` as a general constant.

    So now we calculate the CDF in general for this type of PDF:

    .. math::

        F(r) = \int f(r) dr = \int c r^n dr = \frac{1}{n + 1} c r^{n + 1} + k

    Now with the definition of the CDF at radius :math:`r_{l}` is equal to 0
    and at radius :math:`r_{h}` is equal to 1 we find that the constant from
    integration from this system of equations:

    .. math::

        1 = \frac{1}{n + 1} c ((r_{h})^{n + 1} - (r_{l})^{n + 1}) + k

    Can see that :math:`c = (n + 1) / ((r_{h})^{n + 1} - (r_{l})^{n + 1}))`.
    And :math:`k` is:

    .. math::

        k = - \frac{r_{l}^{n + 1}}{(r_{h})^{n + 1} - (r_{l})^{n + 1}}

    Can see that :math:`c= \frac{n + 1}{R^{n + 1}}`. So can see that the CDF is:

    .. math::

        F(r) = \frac{1}{(r_{h})^{n + 1} - (r_{l})^{n + 1}} r^{n + 1} - \frac{r_{l}^{n + 1}}{(r_{h})^{n + 1} - (r_{l})^{n + 1}}

    And the PDF is the derivative of the CDF:

    .. math::

        f(r) = \frac{(n + 1)}{(r_{h})^{n + 1} - (r_{l})^{n + 1}} (r)^n

    Now we use the probabilty integral transform method to get sampling on
    uniform numbers from a continuous random variable. To do this we find
    the inverse of the CDF evaluated for uniform numbers:

    .. math::

        F(r) = u = \frac{1}{(r_{h})^{n + 1} - (r_{l})^{n + 1}} r^{n + 1} - \frac{r_{l}^{n + 1}}{(r_{h})^{n + 1} - (r_{l})^{n + 1}}

    And find :math:`F^{-1}(u)` gives:

    .. math::

        u = \frac{1}{n + 1} \frac{(r_{h})^{n + 1} - (r_{l})^{n + 1}} - \frac{r_{l}^{n + 1}}{(r_{h})^{n + 1} - (r_{l})^{n + 1}}

    And solving for :math:`r` gives:

    .. math::

        r = ( ((r_{h})^{n + 1} - (r_{l})^{n + 1}) u + (r_{l})^{n + 1})^{\frac{1}{n + 1}}

    Therefore the radius can be sampled by taking the n-th root of uniform
    numbers and multiplying by the radius offset by the lower bound radius.

    \**params :
        The keyword arguments should provide the names of parameters and their
        corresponding bounds, as either tuples or a `boundaries.Bounds`
        instance.

    Attributes
    ----------
    name : 'uniform_radius'
        The name of this distribution.
    dim : int
        The dimension of volume space. In the notation above `dim`
        is :math:`n+1`. For a 3-dimensional sphere this is 3.

    Attributes
    ----------
    params : list of strings
        The list of parameter names.
    bounds : dict
        A dictionary of the parameter names and their bounds.
    norm : float
        The normalization of the multi-dimensional pdf.
    lognorm : float
        The log of the normalization.
    """
    name = "uniform_power_law"
    def __init__(self, dim=None, **params):
        super(UniformPowerLaw, self).__init__(**params)
        self.dim = dim
        self._norm = 1.0
        self._lognorm = 0.0
        for p in self._params:
            self._norm *= self.dim  / \
                                   (self._bounds[p][1]**(self.dim) -
                                    self._bounds[p][0]**(self.dim))
            self._lognorm = numpy.log(self._norm)

    @property
    def norm(self):
        return self._norm

    @property
    def lognorm(self):
        return self._lognorm

    def rvs(self, size=1, param=None):
        """Gives a set of random values drawn from this distribution.

        Parameters
        ----------
        size : {1, int}
            The number of values to generate; default is 1.
        param : {None, string}
            If provided, will just return values for the given parameter.
            Otherwise, returns random values for each parameter.

        Returns
        -------
        structured array
            The random values in a numpy structured array. If a param was
            specified, the array will only have an element corresponding to the
            given parameter. Otherwise, the array will have an element for each
            parameter in self's params.
        """
        if param is not None:
            dtype = [(param, float)]
        else:
            dtype = [(p, float) for p in self.params]
        arr = numpy.zeros(size, dtype=dtype)
        for (p,_) in dtype:
            offset = numpy.power(self._bounds[p][0], self.dim)
            factor = numpy.power(self._bounds[p][1], self.dim) - \
                                      numpy.power(self._bounds[p][0], self.dim)
            arr[p] = numpy.random.uniform(0.0, 1.0, size=size)
            arr[p] = numpy.power(factor * arr[p] + offset, 1.0 / self.dim)
        return arr

    def _pdf(self, **kwargs):
        """Returns the pdf at the given values. The keyword arguments must
        contain all of parameters in self's params. Unrecognized arguments are
        ignored.
        """
        for p in self._params:
            if p not in kwargs.keys():
                raise ValueError(
                            'Missing parameter {} to construct pdf.'.format(p))
        if kwargs in self:
            pdf = self._norm * \
                  numpy.prod([(kwargs[p])**(self.dim - 1)
                              for p in self._params])
            return float(pdf)
        else:
            return 0.0

    def _logpdf(self, **kwargs):
        """Returns the log of the pdf at the given values. The keyword
        arguments must contain all of parameters in self's params. Unrecognized
        arguments are ignored.
        """
        for p in self._params:
            if p not in kwargs.keys():
                raise ValueError(
                            'Missing parameter {} to construct pdf.'.format(p))
        if kwargs in self:
            log_pdf = self._lognorm + \
                      (self.dim - 1) * \
                      numpy.log([kwargs[p] for p in self._params]).sum()
            return log_pdf
        else:
            return -numpy.inf

    @classmethod
    def from_config(cls, cp, section, variable_args):
        """Returns a distribution based on a configuration file. The parameters
        for the distribution are retrieved from the section titled
        "[`section`-`variable_args`]" in the config file.

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
        Uniform
            A distribution instance from the pycbc.inference.prior module.
        """
        return super(UniformPowerLaw, cls).from_config(cp, section,
                                                       variable_args,
                                                       bounds_required=True)


class UniformRadius(UniformPowerLaw):
    """ For a uniform distribution in volume using spherical coordinates, this
-   is the distriubtion to use for the radius.

    For more details see UniformPowerLaw.
    """
    name = "uniform_radius"
    def __init__(self, dim=3, **params):
        super(UniformRadius, self).__init__(dim=3, **params)

__all__ = ["UniformPowerLaw", "UniformRadius"]
