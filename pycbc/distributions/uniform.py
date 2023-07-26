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
This modules provides classes for evaluating uniform distributions.
"""

import numpy
from pycbc.distributions import bounded

class Uniform(bounded.BoundedDist):
    """
    A uniform distribution on the given parameters. The parameters are
    independent of each other. Instances of this class can be called like
    a function. By default, logpdf will be called, but this can be changed
    by setting the class's __call__ method to its pdf method.

    Parameters
    ----------
    \**params :
        The keyword arguments should provide the names of parameters and their
        corresponding bounds, as either tuples or a `boundaries.Bounds`
        instance.

    Examples
    --------
    Create a 2 dimensional uniform distribution:

    >>> from pycbc import distributions
    >>> dist = distributions.Uniform(mass1=(10.,50.), mass2=(10.,50.))

    Get the log of the pdf at a particular value:

    >>> dist.logpdf(mass1=25., mass2=10.)
        -7.3777589082278725

    Do the same by calling the distribution:

    >>> dist(mass1=25., mass2=10.)
        -7.3777589082278725

    Generate some random values:

    >>> dist.rvs(size=3)
        array([(36.90885758394699, 51.294212757995254),
               (39.109058546060346, 13.36220145743631),
               (34.49594465315212, 47.531953033719454)],
              dtype=[('mass1', '<f8'), ('mass2', '<f8')])

    Initialize a uniform distribution using a boundaries.Bounds instance,
    with cyclic bounds:

    >>> dist = distributions.Uniform(phi=Bounds(10, 50, cyclic=True))

    Apply boundary conditions to a value:

    >>> dist.apply_boundary_conditions(phi=60.)
        {'mass1': array(20.0)}

    The boundary conditions are applied to the value before evaluating the pdf;
    note that the following returns a non-zero pdf. If the bounds were not
    cyclic, the following would return 0:

    >>> dist.pdf(phi=60.)
        0.025
    """
    name = 'uniform'
    def __init__(self, **params):
        super(Uniform, self).__init__(**params)
        # compute the norm and save
        # temporarily suppress numpy divide by 0 warning
        with numpy.errstate(divide="ignore"):
            self._lognorm = -sum([numpy.log(abs(bnd[1]-bnd[0]))
                                  for bnd in self._bounds.values()])
            self._norm = numpy.exp(self._lognorm)

    @property
    def norm(self):
        """float: The normalization of the multi-dimensional pdf."""
        return self._norm

    @property
    def lognorm(self):
        """float: The log of the normalization"""
        return self._lognorm

    def _cdfinv_param(self, param, value):
        """Return the inverse cdf to map the unit interval to parameter bounds.
        """
        lower_bound = self._bounds[param][0]
        upper_bound = self._bounds[param][1]
        return (upper_bound - lower_bound) * value + lower_bound

    def _pdf(self, **kwargs):
        """Returns the pdf at the given values. The keyword arguments must
        contain all of parameters in self's params. Unrecognized arguments are
        ignored.
        """
        if kwargs in self:
            return self._norm
        else:
            return 0.

    def _logpdf(self, **kwargs):
        """Returns the log of the pdf at the given values. The keyword
        arguments must contain all of parameters in self's params. Unrecognized
        arguments are ignored.
        """
        if kwargs in self:
            return self._lognorm
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
            ``VARARGS_DELIM``. These must appear in the "tag" part
            of the section header.

        Returns
        -------
        Uniform
            A distribution instance from the pycbc.inference.prior module.
        """
        return super(Uniform, cls).from_config(cp, section, variable_args,
                     bounds_required=True)


__all__ = ['Uniform']
