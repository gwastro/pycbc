# Copyright (C) 2017  Christopher M. Biwer
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
""" This modules provides classes for evaluating distributions whose logarithm
are uniform.
"""

import numpy
from pycbc.distributions import uniform

class UniformLog10(uniform.Uniform):
    """ A uniform distribution on the log base 10 of the given parameters.
    The parameters are independent of each other. Instances of this class can
    be called like a function. By default, logpdf will be called.

    Parameters
    ----------
    \**params :
        The keyword arguments should provide the names of parameters and their
        corresponding bounds, as either tuples or a `boundaries.Bounds`
        instance.

    Attributes
    ----------
    name : "uniform_log10"
        The name of this distribution.
    """
    name = "uniform_log10"

    def __init__(self, **params):
        super(UniformLog10, self).__init__(**params)
        self._norm = numpy.prod([numpy.log10(bnd[1]) - numpy.log10(bnd[0])
                                   for bnd in self._bounds.values()])
        self._lognorm = numpy.log(self._norm)

    def cdfinv(self, param, value):
        """Return the inverse cdf to map the unit interval to parameter bounds.
        """
        lower_bound = numpy.log10(self._bounds[param][0])
        upper_bound = numpy.log10(self._bounds[param][1])
        new_value = 10. ** ((upper_bound - lower_bound) * value + lower_bound)
        return new_value

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
            log_high = numpy.log10(self._bounds[p][0])
            log_low = numpy.log10(self._bounds[p][1])
            arr[p] = 10.0**(numpy.random.uniform(log_low, log_high, size=size))
        return arr


    def _pdf(self, **kwargs):
        """Returns the pdf at the given values. The keyword arguments must
        contain all of parameters in self's params. Unrecognized arguments are
        ignored.
        """
        if kwargs in self:
            vals = numpy.array([numpy.log(10) * self._norm * kwargs[param]
                                for param in kwargs.keys()])
            return 1.0 / numpy.prod(vals)
        else:
            return 0.

    def _logpdf(self, **kwargs):
        """Returns the log of the pdf at the given values. The keyword
        arguments must contain all of parameters in self's params. Unrecognized
        arguments are ignored.
        """
        if kwargs in self:
            return numpy.log(self._pdf(**kwargs))
        else:
            return -numpy.inf

__all__ = ["UniformLog10"]
