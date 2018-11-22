# Copyright (C)  2017 Collin Capano, Christopher M. Biwer, Alex Nitz
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
""" This module provides classes to describe joint distributions
"""
import logging
import numpy
from pycbc.io import record

class JointDistribution(object):
    """
    Callable class that calculates the joint distribution built from a set of
    distributions.

    Parameters
    ----------
    variable_args : list
        A list of strings that contain the names of the variable parameters and
        the order they are expected when the class is called.
    \*distributions :
        The rest of the arguments must be instances of distributions describing
        the individual distributions on the variable parameters.
        A single distribution may contain
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
    An example of creating a joint distribution with constraint that total mass must
    be below 30.

    >>> from pycbc.distributions import uniform, JointDistribution
    >>> def mtotal_lt_30(params):
        ...    return True if params["mass1"] + params["mass2"] < 30 else False
    >>> mass_lim = (2, 50)
    >>> uniform_prior = Uniform(mass1=mass_lim, mass2=mass_lim)
    >>> prior_eval = JointDistribution(["mass1", "mass2"], uniform_prior,
        ...                               constraints=[mtotal_lt_30])
    >>> print(prior_eval(mass1=20, mass2=1))

    """
    name = 'joint'

    def __init__(self, variable_args, *distributions, **kwargs):

        # store the names of the parameters defined in the distributions
        self.variable_args = tuple(variable_args)

        # store the distributions
        self.distributions = distributions

        # store the constraints on the parameters defined inside the
        # distributions list
        self._constraints = kwargs["constraints"] \
                                  if "constraints" in kwargs.keys() else []

        # check that all of the supplied parameters are described by the given
        # distributions
        distparams = set()
        for dist in distributions:
            distparams.update(set(dist.params))

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

        # if there are constraints then find the renormalization factor
        # since a constraint will cut out part of the space
        # do this by random sampling the full space and find the percent
        # of samples rejected
        n_test_samples = kwargs["n_test_samples"] \
                             if "n_test_samples" in kwargs else int(1e6)
        if self._constraints:
            logging.info("Renormalizing distribution for constraints")

            # draw samples
            samples = {}
            for dist in self.distributions:
                draw = dist.rvs(n_test_samples)
                for param in dist.params:
                    samples[param] = draw[param][:]

            # evaluate constraints
            result = numpy.ones(len(samples.values()[0]), dtype=bool)
            for constraint in self._constraints:
                result = constraint(samples) & result

            # set new scaling factor for prior to be
            # the fraction of acceptances in random sampling of entire space
            self._pdf_scale = sum(result) / float(n_test_samples)

        else:
            self._pdf_scale = 1.0

        # since Distributions will return logpdf we keep the scale factor
        # in log scale as well for self.__call__
        self._logpdf_scale = numpy.log(self._pdf_scale)

    def apply_boundary_conditions(self, **params):
        """Applies each distributions' boundary conditions to the given list
        of parameters, returning a new list with the conditions applied.

        Parameters
        ----------
        **params :
            Keyword arguments should give the parameters to apply the
            conditions to.

        Returns
        -------
        dict
            A dictionary of the parameters after each distribution's
            `apply_boundary_conditions` function has been applied.
        """
        for dist in self.distributions:
            params.update(dist.apply_boundary_conditions(**params))
        return params

    def __call__(self, **params):
        """Evalualate joint distribution for parameters.
        """
        for constraint in self._constraints:
            if not constraint(params):
                return -numpy.inf
        return sum([d(**params)
                    for d in self.distributions]) - self._logpdf_scale

    def rvs(self, size=1):
        """ Rejection samples the parameter space.
        """

        # create output FieldArray
        out = record.FieldArray(size, dtype=[(arg, float)
                                    for arg in self.variable_args])

        # loop until enough samples accepted
        n = 0
        while n < size:

            # draw samples
            samples = {}
            for dist in self.distributions:
                draw = dist.rvs(1)
                for param in dist.params:
                    samples[param] = draw[param][0]
            vals = numpy.array([samples[arg] for arg in self.variable_args])

            # determine if all parameter values are in prior space
            # if they are then add to output
            if self(**dict(zip(self.variable_args, vals))) > -numpy.inf:
                out[n] = vals
                n += 1

        return out

