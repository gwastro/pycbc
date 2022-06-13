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
from pycbc.io.record import FieldArray

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
        Valid keyword arguments include:
        `constraints` : a list of functions that accept a dict of parameters
        with the parameter name as the key. If the constraint is satisfied the
        function should return True, if the constraint is violated, then the
        function should return False.
        `n_test_samples` : number of random draws used to fix pdf normalization
        factor after applying constraints.

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

    >>> from pycbc.distributions import Uniform, JointDistribution
    >>> def mtotal_lt_30(params):
        ...    return params["mass1"] + params["mass2"] < 30
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

        # store kwargs
        self.kwargs = kwargs

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
                    samples[param] = draw[param]
            samples = FieldArray.from_kwargs(**samples)

            # evaluate constraints
            result = self.within_constraints(samples)

            # set new scaling factor for prior to be
            # the fraction of acceptances in random sampling of entire space
            self._pdf_scale = result.sum() / float(n_test_samples)
            if self._pdf_scale == 0.0:
                raise ValueError("None of the random draws for pdf "
                    "renormalization satisfied the constraints. "
                    " You can try increasing the 'n_test_samples' keyword.")

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

    @staticmethod
    def _return_atomic(params):
        """Determines if an array or atomic value should be returned given a
        set of input params.

        Parameters
        ----------
        params : dict, numpy.record, array, or FieldArray
            The input to evaluate.

        Returns
        -------
        bool :
            Whether or not functions run on the parameters should be returned
            as atomic types or not.
        """
        if isinstance(params, dict):
            return not any(isinstance(val, numpy.ndarray)
                           for val in params.values())
        elif isinstance(params, numpy.record):
            return True
        elif isinstance(params, numpy.ndarray):
            return False
            params = params.view(type=FieldArray)
        elif isinstance(params, FieldArray):
            return False
        else:
            raise ValueError("params must be either dict, FieldArray, "
                             "record, or structured array")

    @staticmethod
    def _ensure_fieldarray(params):
        """Ensures the given params are a ``FieldArray``.

        Parameters
        ----------
        params : dict, FieldArray, numpy.record, or numpy.ndarray
            If the given object is a dict, it will be converted to a
            FieldArray.

        Returns
        -------
        FieldArray
            The given values as a FieldArray.
        """
        if isinstance(params, dict):
            return FieldArray.from_kwargs(**params)
        elif isinstance(params, numpy.record):
            return FieldArray.from_records(tuple(params),
                                           names=params.dtype.names)
        elif isinstance(params, numpy.ndarray):
            return params.view(type=FieldArray)
        elif isinstance(params, FieldArray):
            return params
        else:
            raise ValueError("params must be either dict, FieldArray, "
                             "record, or structured array")

    def within_constraints(self, params):
        """Evaluates whether the given parameters satisfy the constraints.

        Parameters
        ----------
        params : dict, FieldArray, numpy.record, or numpy.ndarray
            The parameter values to evaluate.

        Returns
        -------
        (array of) bool :
            If params was an array, or if params a dictionary and one or more
            of the parameters are arrays, will return an array of booleans.
            Otherwise, a boolean.
        """
        params = self._ensure_fieldarray(params)
        return_atomic = self._return_atomic(params)
        # convert params to a field array if it isn't one
        result = numpy.ones(params.shape, dtype=bool)
        for constraint in self._constraints:
            result &= constraint(params)
        if return_atomic:
            result = result.item()
        return result

    def contains(self, params):
        """Evaluates whether the given parameters satisfy the boundary
            conditions, boundaries, and constraints. This method is different
            from `within_constraints`, that method only check the constraints.

        Parameters
        ----------
        params : dict, FieldArray, numpy.record, or numpy.ndarray
            The parameter values to evaluate.

        Returns
        -------
        (array of) bool :
            If params was an array, or if params a dictionary and one or more
            of the parameters are arrays, will return an array of booleans.
            Otherwise, a boolean.
        """
        params = self.apply_boundary_conditions(**params)
        result = True
        for dist in self.distributions:
            param_name = dist.params[0]
            contain_array = numpy.ones(len(params[param_name]), dtype=bool)
            # note: enable `__contains__` in `pycbc.distributions.bounded`
            # to handle array-like input, it doesn't work now.
            for k in params[param_name]:
                index = numpy.where(params[param_name] == k)[0][0]
                contain_array[index] = {param_name: k} in dist
            result &= numpy.array(contain_array)
        result &= self.within_constraints(params)
        return result

    def __call__(self, **params):
        """Evaluate joint distribution for parameters.
        """
        return_atomic = self._return_atomic(params)
        # check if statisfies constraints
        if len(self._constraints) != 0:
            parray = self._ensure_fieldarray(params)
            isin = self.within_constraints(parray)
            if not isin.any():
                if return_atomic:
                    out = -numpy.inf
                else:
                    out = numpy.full(parray.shape, -numpy.inf)
                return out

        # evaluate
        # note: this step may fail if arrays of values were provided, as
        # not all distributions are vectorized currently
        logps = numpy.array([d(**params) for d in self.distributions])
        logp = logps.sum(axis=0)

        if len(self._constraints) != 0:
            logp += numpy.log(isin.astype(float))

        if return_atomic:
            logp = logp.item()

        return logp - self._logpdf_scale

    def rvs(self, size=1):
        """ Rejection samples the parameter space.
        """
        # create output FieldArray
        dtype = [(arg, float) for arg in self.variable_args]
        out = FieldArray(size, dtype=dtype)
        # loop until enough samples accepted
        remaining = size
        ndraw = size
        while remaining:
            # scratch space for evaluating constraints
            scratch = FieldArray(ndraw, dtype=dtype)
            for dist in self.distributions:
                # drawing samples from the distributions is generally faster
                # then evaluating constrants, so we'll always draw the full
                # size, even if that gives us more points than we need
                draw = dist.rvs(size=ndraw)
                for param in dist.params:
                    scratch[param] = draw[param]
            # apply any constraints
            keep = self.within_constraints(scratch)
            nkeep = keep.sum()
            kmin = size - remaining
            kmax = min(nkeep, remaining)
            out[kmin:kmin+kmax] = scratch[keep][:kmax]
            remaining = max(0, remaining - nkeep)
            # to try to speed up next go around, we'll increase the draw
            # size by the fraction of values that were kept, but cap at 1e6
            ndraw = int(min(1e6, ndraw * numpy.ceil(ndraw / (nkeep + 1.))))
        return out

    @property
    def well_reflected(self):
        """ Get list of which parameters are well reflected
        """
        reflect = []
        bounds = self.bounds
        for param in bounds:
            if bounds[param].reflected == 'well':
                reflect.append(param)
        return reflect

    @property
    def cyclic(self):
        """ Get list of which parameters are cyclic
        """
        cyclic = []
        bounds = self.bounds
        for param in bounds:
            if bounds[param].cyclic:
                cyclic.append(param)
        return cyclic

    @property
    def bounds(self):
        """ Get the dict of boundaries
        """
        bnds = {}
        for dist in self.distributions:
            if hasattr(dist, 'bounds'):
                bnds.update(dist.bounds)
        return bnds

    def cdfinv(self, **original):
        """ Apply the inverse cdf to the array of values [0, 1]. Every
        variable parameter must be given as a keyword argument.
        """
        updated = {}
        for dist in self.distributions:
            updated.update(dist.cdfinv(**original))
        return updated
