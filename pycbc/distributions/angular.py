# Copyright (C) 2016  Collin Capano
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
This modules provides classes for evaluating non-uniform angular distributions.
"""

import numpy
from pycbc.distributions import boundaries
from pycbc.distributions import uniform_angle

class SinAngle(uniform_angle.UniformAngle):
    r"""A sine distribution; the pdf of each parameter `\theta` is given by:

    ..math::
        p(\theta) = \frac{\sin \theta}{\cos\theta_0 - \cos\theta_1}, \theta_0 \leq \theta < \theta_1,

    and 0 otherwise. Here, :math:`\theta_0, \theta_1` are the bounds of the
    parameter.

    The domain of this distribution is `[0, pi]`. This is accomplished by
    putting hard boundaries at `[0, pi]`. Bounds may be provided to further
    limit the range for which the pdf has support.  As with `UniformAngle`,
    these are initizliaed as multiples of pi, while the stored bounds are in
    radians.

    Parameters
    ----------
    \**params :
        The keyword arguments should provide the names of parameters and
        (optionally) their corresponding bounds, as either
        `boundaries.Bounds` instances or tuples. The bounds must be
        in [0,1]. These are converted to radians for storage. None may also
        be passed; in that case, the domain bounds will be used.

    Attributes
    ----------------
    name : 'sin_angle'
        The name of this distribution.
    params : list of strings
        The list of parameter names.
    bounds : dict
        A dictionary of the parameter names and their bounds, in radians.
    """
    name = 'sin_angle'
    _func = numpy.cos
    _dfunc = numpy.sin
    _arcfunc = numpy.arccos
    # _domain applies the reflection off of 0, pi
    _domain = boundaries.Bounds(0, numpy.pi,
        btype_min='closed', btype_max='closed', cyclic=False)


    def _pdf(self, **kwargs):
        """Returns the pdf at the given values. The keyword arguments must
        contain all of parameters in self's params. Unrecognized arguments are
        ignored.
        """
        if kwargs not in self:
            return 0.
        return self._norm * \
            self._dfunc(numpy.array([kwargs[p] for p in self._params])).prod()


    def _logpdf(self, **kwargs):
        """Returns the log of the pdf at the given values. The keyword
        arguments must contain all of parameters in self's params. Unrecognized
        arguments are ignored.
        """
        if kwargs not in self:
            return -numpy.inf
        return self._lognorm + \
            numpy.log(self._dfunc(
                numpy.array([kwargs[p] for p in self._params]))).sum()


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
            arr[p] = self._arcfunc(numpy.random.uniform(
                                    self._func(self._bounds[p][0]),
                                    self._func(self._bounds[p][1]),
                                    size=size))
        return arr


class CosAngle(SinAngle):
    r"""A cosine distribution. This is the same thing as a sine distribution,
    but with the domain shifted to `[-pi/2, pi/2]`. See SinAngle for more
    details.

    Parameters
    ----------
    \**params :
        The keyword arguments should provide the names of parameters and
        (optionally) their corresponding bounds, as either
        `boundaries.Bounds` instances or tuples. The bounds must be
        in [-0.5, 0.5]. These are converted to radians for storage.
        None may also be passed; in that case, the domain bounds will be used.

    Attributes
    ----------------
    name : 'cos_angle'
        The name of this distribution.
    params : list of strings
        The list of parameter names.
    bounds : dict
        A dictionary of the parameter names and their bounds, in radians.
    """
    name = 'cos_angle'
    _func = numpy.sin
    _dfunc = numpy.cos
    _arcfunc = numpy.arcsin
    _domain = boundaries.Bounds(-numpy.pi/2., numpy.pi/2.,
        btype_min='closed', btype_max='closed', cyclic=False)


