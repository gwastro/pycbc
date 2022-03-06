# Copyright (C) 2021 Yifan Wang
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

"""This modules provides classes for evaluating distributions for mchirp and
q (i.e., mass ratio) from uniform component mass.
"""

import numpy
from scipy.interpolate import interp1d
from scipy.special import hyp2f1
from pycbc.distributions import power_law
from pycbc.distributions import bounded


class MchirpfromUniformMass1Mass2(power_law.UniformPowerLaw):
    r"""A distribution for chirp mass from uniform component mass +
    constraints given by chirp mass. This is a special case for UniformPowerLaw
    with index 1. For more details see UniformPowerLaw.

    The parameters (i.e. `**params`) are independent of each other. Instances
    of this class can be called like a function. By default, `logpdf` will be
    called, but this can be changed by setting the class's `__call__` method
    to its pdf method.

    Derivation for the probability density function:

    .. math::

        P(m_1,m_2)dm_1dm_2 = P(\mathcal{M}_c,q)d\mathcal{M}_cdq

    Where :math:`\mathcal{M}_c` is chirp mass and :math:`q` is mass ratio,
    :math:`m_1` and :math:`m_2` are component masses. The jacobian to transform
    chirp mass and mass ratio to component masses is

    .. math::

        \frac{\partial(m_1,m_2)}{\partial(\mathcal{M}_c,q)} = \
        \mathcal{M}_c \left(\frac{1+q}{q^3}\right)^{2/5}

    (https://github.com/gwastro/pycbc/blob/master/pycbc/transforms.py#L416.)

    Because :math:`P(m_1,m_2) = const`, then

    .. math::

        P(\mathcal{M}_c,q) = P(\mathcal{M}_c)P(q)\propto
        \mathcal{M}_c \left(\frac{1+q}{q^3}\right)^{2/5}`.

    Therefore,

    .. math::
        P(\mathcal{M}_c) \propto \mathcal{M}_c

    and

    .. math::
        P(q) \propto \left(\frac{1+q}{q^3}\right)^{2/5}

    Examples
    --------

    Generate 10000 random numbers from this distribution in [5,100]

    >>> from pycbc import distributions as dist
    >>> minmc = 5, maxmc = 100, size = 10000
    >>> mc = dist.MchirpfromUniformMass1Mass2(value=(minmc,maxmc)).rvs(size)

    The settings in the configuration file for pycbc_inference should be

    .. code-block:: ini

        [variable_params]
        mchirp =
        [prior-mchirp]
        name = mchirp_from_uniform_mass1_mass2
        min-mchirp = 10
        max-mchirp = 80

    Parameters
    ----------
    \**params :
        The keyword arguments should provide the names of parameters and their
        corresponding bounds, as either tuples or a `boundaries.Bounds`
        instance.
    """

    name = "mchirp_from_uniform_mass1_mass2"

    def __init__(self, dim=2, **params):
        super(MchirpfromUniformMass1Mass2, self).__init__(dim=2, **params)


class QfromUniformMass1Mass2(bounded.BoundedDist):
    r"""A distribution for mass ratio (i.e., q) from uniform component mass
    + constraints given by q.

    The parameters (i.e. `**params`) are independent of each other. Instances
    of this class can be called like a function. By default, `logpdf` will be
    called, but this can be changed by setting the class's `__call__` method
    to its pdf method.

    For mathematical derivation see the documentation above in the class
    `MchirpfromUniformMass1Mass2`.

    Parameters
    ----------
    \**params :
        The keyword arguments should provide the names of parameters and their
        corresponding bounds, as either tuples or a `boundaries.Bounds`
        instance.

    Examples
    --------

    Generate 10000 random numbers from this distribution in [1,8]

    >>> from pycbc import distributions as dist
    >>> minq = 1, maxq = 8, size = 10000
    >>> q = dist.QfromUniformMass1Mass2(value=(minq,maxq)).rvs(size)

    """

    name = 'q_from_uniform_mass1_mass2'

    def __init__(self, **params):
        super(QfromUniformMass1Mass2, self).__init__(**params)
        self._norm = 1.0
        self._lognorm = 0.0
        for p in self._params:
            self._norm /= self._cdf_param(p, self._bounds[p][1]) - \
                self._cdf_param(p, self._bounds[p][0])
        self._lognorm = numpy.log(self._norm)

    @property
    def norm(self):
        """float: The normalization of the multi-dimensional pdf."""
        return self._norm

    @property
    def lognorm(self):
        """float: The log of the normalization."""
        return self._lognorm

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
                numpy.prod([(1.+kwargs[p])**(2./5)/kwargs[p]**(6./5)
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
                    'Missing parameter {} to construct logpdf.'.format(p))
        if kwargs in self:
            return numpy.log(self._pdf(**kwargs))
        else:
            return -numpy.inf

    def _cdf_param(self, param, value):
        r""">>> from sympy import *
           >>> x = Symbol('x')
           >>> integrate((1+x)**(2/5)/x**(6/5))
           Output:
                             _
                      -0.2  |_  /-0.4, -0.2 |    I*pi\
                -5.0*x    * |   |           | x*e    |
                           2  1 \   0.8     |        /
        """
        if param in self._params:
            return -5. * value**(-1./5) * hyp2f1(-2./5, -1./5, 4./5, -value)
        else:
            raise ValueError('{} is not contructed yet.'.format(param))

    def _cdfinv_param(self, param, value):
        """Return the inverse cdf to map the unit interval to parameter bounds.
        Note that value should be uniform in [0,1]."""
        if (numpy.array(value) < 0).any() or (numpy.array(value) > 1).any():
            raise ValueError(
                'q_from_uniform_m1_m2 cdfinv requires input in [0,1].')
        if param in self._params:
            lower_bound = self._bounds[param][0]
            upper_bound = self._bounds[param][1]
            q_array = numpy.linspace(
                lower_bound, upper_bound, num=1000, endpoint=True)
            q_invcdf_interp = interp1d(self._cdf_param(param, q_array),
                                       q_array, kind='cubic',
                                       bounds_error=True)

            return q_invcdf_interp(
                (self._cdf_param(param, upper_bound) -
                 self._cdf_param(param, lower_bound)) * value +
                self._cdf_param(param, lower_bound))
        else:
            raise ValueError('{} is not contructed yet.'.format(param))

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
        for (p, _) in dtype:
            uniformcdfvalue = numpy.random.uniform(0, 1, size=size)
            arr[p] = self._cdfinv_param(p, uniformcdfvalue)
        return arr

    @classmethod
    def from_config(cls, cp, section, variable_args):
        """Returns a distribution based on a configuration file. The parameters
        for the distribution are retrieved from the section titled
        "[`section`-`variable_args`]" in the config file.

        Example:

        .. code-block:: ini

            [variable_params]
            q =
            [prior-q]
            name = q_from_uniform_mass1_mass2
            min-q = 1
            max-q = 8

        Parameters
        ----------
        cp : pycbc.workflow.WorkflowConfigParser
            A parsed configuration file that contains the distribution
            options.
        section : str
            Name of the section in the configuration file.
        variable_args : str
            The names of the parameters for this distribution, separated by
            ``VARARGS_DELIM``. These must appear in the "tag" part
            of the section header.

        Returns
        -------
        QfromUniformMass1Mass2
            A distribution instance from the pycbc.distributions.bounded
        module.
        """
        return super(QfromUniformMass1Mass2, cls).from_config(
            cp, section, variable_args, bounds_required=True)


__all__ = ["MchirpfromUniformMass1Mass2", "QfromUniformMass1Mass2"]
