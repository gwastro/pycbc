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


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
"""
This modules provides classes and functions for evaluating the prior
for parameter estimation.
"""

import numpy

class PriorEvaluator(object):
    """
    Callable class that calculates the prior.

    Parameters
    ----------
    variable_args : list
        A list of strings that contain the names of the variable parameters and
        the order they are expected when the class is called.
    \*distributions :
        The rest of the arguments must be instances of distributions describing
        the priors on the variable parameters. A single distribution may contain
        multiple parameters. The set of all params across the distributions
        (retrieved from the distributions' params attribute) must be the same
        as the set of variable_args provided.
    \*\*kwargs :
        Valid keyword arguments include `constraints`. `constraints` is a list
        functions that accept a dict of parameters with the parameter name as
        the key. If the constraint is satisfied the function should return
        True, if the constraint is violated, then the function should return
        False.

    Attributes
    ----------
    variable_args : tuple
        The parameters expected when the evaluator is called.
    distributions : list
        The distributions for the parameters.
    constraints : list
        A list of functions to test if parameter values obey multi-dimensional
        constraints.

    Examples
    --------
    An example of creating a prior with constraint that total mass must
    be below 30.

    >>> from pycbc.inference import distributions
    >>> from pycbc.inference import prior
    >>> def mtotal_lt_30(params):
        ...    return True if params["mass1"] + params["mass2"] < 30 else False
    >>> mass_lim = (2, 50)
    >>> uniform_prior = distributions.Uniform(mass1=mass_lim, mass2=mass_lim)
    >>> prior_eval = prior.PriorEvaluator(["mass1", "mass2"], uniform_prior,
        ...                               constraints=[mtotal_lt_30])
    >>> print prior_eval([20, 1])

    """

    def __init__(self, variable_args, *distributions, **kwargs):

        # store the names of the variable params
        self.variable_args = tuple(variable_args)
        # store the distributions
        self.distributions = distributions
        # store the constraints
        self.constraints = kwargs["constraints"] \
                                  if "constraints" in kwargs.keys() else []

        # check that all of the variable args are described by the given
        # distributions
        distparams = set()
        [distparams.update(set(dist.params)) for dist in distributions]
        varset = set(self.variable_args)
        missing_params = distparams - varset
        if missing_params:
            raise ValueError("provided variable_args do not include "
                "parameters %s" %(','.join(missing_params)) + " which are "
                "required by the provided distributions")
        extra_params = varset - distparams
        if extra_params:
            raise ValueError("variable_args %s " %(','.join(extra_params)) +
                "are not in any of the provided distributions")

    def apply_boundary_conditions(self, params):
        """Applies each distributions' boundary conditions to the given list
        of parameters, returning a new list with the conditions applied.

        Parameters
        ----------
        params : list
            List of parameters to apply conditions to. The order of the
            parameters is assumed to be the same as self.variable_args.

        Returns
        -------
        list
            List of the parameters after each distribution's
            `apply_boundary_conditions` function has been applied.
        """
        params = dict(zip(self.variable_args, params))
        for dist in self.distributions:
            params.update(dist.apply_boundary_conditions(**params))
        return [params[p] for p in self.variable_args]

    def __call__(self, params):
        """ Evalualate prior for parameters.
        """
        params = dict(zip(self.variable_args, params))
        for constraint in self.constraints:
            if not constraint(params):
                return -numpy.inf
        return sum([d(**params) for d in self.distributions])


