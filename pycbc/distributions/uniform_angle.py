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
This modules provides classes for evaluating uniform angle distributions.
"""

import numpy
from ConfigParser import Error
from pycbc.distributions import boundaries
from pycbc.distributions import bounded
from pycbc.distributions import uniform

class UniformAngle(uniform.Uniform):
    """A uniform distribution in which the dependent variable is cyclic between
    `[0,2pi)`.
    
    Bounds may be provided to limit the range for which the pdf has support.
    If provided, the parameter bounds are initialized as multiples of pi,
    while the stored bounds are in radians.

    Parameters
    ----------
    \**params :
        The keyword arguments should provide the names of parameters and
        (optionally) their corresponding bounds, as either
        `boundaries.Bounds` instances or tuples. The bounds must be
        in [0,2). These are converted to radians for storage. None may also
        be passed; in that case, the domain bounds will be used.

    Attributes
    ----------------
    name : 'uniform_angle'
        The name of this distribution.
    params : list of strings
        The list of parameter names.
    bounds : dict
        A dictionary of the parameter names and their bounds, in radians.

    Notes
    ------
    For more information, see Uniform.
    """
    name = 'uniform_angle'
    # _domain is a bounds instance used apply the cyclic conditions; this is
    # applied first, before any bounds specified in the initialization are used
    _domain = boundaries.Bounds(0., 2*numpy.pi, cyclic=True)

    def __init__(self, **params):
        for p,bnds in params.items():
            if bnds is None:
                bnds = self._domain
            elif isinstance(bnds, boundaries.Bounds):
                # convert to radians
                bnds._min = bnds._min.__class__(bnds._min * numpy.pi)
                bnds._max = bnds._max.__class__(bnds._max * numpy.pi)
            else:
                # create a Bounds instance from the given tuple
                bnds = boundaries.Bounds(
                    bnds[0]*numpy.pi, bnds[1]*numpy.pi)
            # check that the bounds are in the domain
            if bnds.min < self._domain.min or bnds.max > self._domain.max:
                raise ValueError("bounds must be in [{x},{y}); "
                    "got [{a},{b})".format(x=self._domain.min/numpy.pi,
                    y=self._domain.max/numpy.pi, a=bnds.min/numpy.pi,
                    b=bnds.max/numpy.pi))
            # update
            params[p] = bnds
        super(UniformAngle, self).__init__(**params)

    def apply_boundary_conditions(self, **kwargs):
        """Maps values to be in [0, 2pi) (the domain) first, before applying
        any additional boundary conditions.

        Parameters
        ----------
        \**kwargs :
            The keyword args should be the name of a parameter and value to
            apply its boundary conditions to. The arguments need not include
            all of the parameters in self.

        Returns
        -------
        dict
            A dictionary of the parameter names and the conditioned values.
        """
        # map values to be within the domain
        kwargs = dict([[p, self._domain.apply_conditions(val)]
                      for p,val in kwargs.items()])
        # now apply additional conditions
        return super(UniformAngle, self).apply_boundary_conditions(**kwargs)

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
        UniformAngle
            A distribution instance from the pycbc.inference.prior module.
        """
        return bounded.bounded_from_config(cls, cp, section, variable_args,
                                           bounds_required=False)


