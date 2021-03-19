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

""" This modules provides classes for evaluating distributions for mchirp and 
mass ratio from uniform component mass.
"""

import numpy
from pycbc.distributions import uniform
from scipy.interpolate import interp1d
from scipy.special import hyp2f1

class MchirpfromUniformMass1Mass2(UniformPowerLaw):
    """ For a uniform distribution in volume using spherical coordinates, this
    is the distriubtion to use for the radius.

    For more details see UniformPowerLaw.
    """
    name = "mchirp_from_uniform_mass1_mass2"
    def __init__(self, dim=2, **params):
        super(MchirpfromUniformMass1Mass2, self).__init__(dim=2, **params)

class QfromUniformMass1Mass2(bounded.BoundedDist):
    """
    A distribution for mass ratio (i.e., q) from uniform component mass. The parameters are
    independent of each other. Instances of this class can be called like
    a function. By default, logpdf will be called, but this can be changed
    by setting the class's __call__ method to its pdf method.

    Parameters
    ----------
    \**params :
        The keyword arguments should provide the names of parameters and their
        corresponding bounds, as either tuples or a `boundaries.Bounds`
        instance.

    Attributes
    ----------
    name : 'q_from_uniform_componentmass'
        The name of this distribution.

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

    Examples
    --------
    
    """
    name = 'q_from_uniform_componentmass'
    def __init__(self, **params):
        super(QfromUniformMass1Mass2, self).__init__(**params)
        self._norm = 1.0
        self._lognorm = 0.0
        for p in self._params:
            self._norm /= self._cdf_param(p,self._bounds[p][1]) - \
                          self._cdf_param(p,self._bounds[p][0]) 
            self._lognorm += numpy.log(self._norm)

    @property
    def norm(self):
        return self._norm

    @property
    def lognorm(self):
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
                  numpy.prod([(1.+kwargs[p])**(2./5)/kwargs[p]**(6./5) \
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
        if param in self._params:
            return -5. * value**(-1. / 5.) \
            * hyp2f1(-2. / 5., -1. / 5., 4. / 5., -value)
        else:
            raise ValueError('{} is not contructed yet.'.format(param))

    def _cdfinv_param(self, param, value):
        """Return the inverse cdf to map the unit interval to parameter bounds.
        """
        if param in self._params:
            lower_bound = self._bounds[param][0]
            upper_bound = self._bounds[param][1]
            q_array = np.linspace(lower_bound, upper_bound, 1000)
            q_invcdf_interp = interp1d(self._cdf_param(param,q_array), q_array, kind='cubic',
            bounds_error=False, fill_value=(minimum, maximum))  
            return q_invcdf_interp(value)
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
        for (p,_) in dtype:
            uniformcdf = numpy.random.uniform(self._cdf_param(p,self._bounds[p][0]),
                                        self._cdf_param(p,self._bounds[p][1]),
                                        size=size)
            arr[p] = self._cdfinv_param(p,uniformcdf)
        return arr

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
            ``VARARGS_DELIM``. These must appear in the "tag" part
            of the section header.

        Returns
        -------
        Uniform
            A distribution instance from the pycbc.inference.prior module.
        """
        return super(QfromUniformMass1Mass2, cls).from_config(cp, section, variable_args,
                     bounds_required=True)
