# Copyright (C) 2016  Christopher M. Biwer
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

#
#   Distributions for priors
#
class Uniform(object):
    """
    A uniform distribution on the given parameters. The parameters are
    independent of each other. Instances of this class can be called like
    a function. By default, logpdf will be called, but this can be changed
    by setting the class's __call__ method to its pdf method.

    Parameters
    ----------
    \**params :
        The keyword arguments should provide the names of parameters and their
        corresponding bounds, as tuples.

    Class Attributes
    ----------------
    name : 'uniform'
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

    Example
    -------
    Create a 2 dimensional uniform distribution:
    >>> dist = prior.Uniform(mass1=(10.,50.), mass2=(10.,50.))

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
    """
    name = 'uniform'
    def __init__(self, **params):
        self._bounds = params
        self._params = sorted(params.keys())
        # temporarily suppress numpy divide by 0 warning
        numpy.seterr(divide='ignore')
        self._lognorm = -sum([numpy.log(abs(bnd[1]-bnd[0]))
                                    for bnd in self._bounds.values()])
        self._norm = numpy.exp(self._lognorm)
        numpy.seterr(divide='warn')

    @property
    def params(self):
        return self._params

    @property
    def bounds(self):
        return self._bounds

    @property
    def norm(self):
        return self._norm

    @property
    def lognorm(self):
        return self._lognorm

    def __contains__(self, params):
        try:
            return all([(params[p] >= self._bounds[p][0]) &
                        (params[p] < self._bounds[p][1])
                       for p in self._params])
        except KeyError:
            raise ValueError("must provide all parameters [%s]" %(
                ', '.join(self._params)))

    def pdf(self, **kwargs):
        """
        Returns the pdf at the given values. The keyword arguments must contain
        all of parameters in self's params. Unrecognized arguments are ignored.
        """
        if kwargs in self:
            return self._norm
        else:
            return 0.

    def logpdf(self, **kwargs):
        """
        Returns the log of the pdf at the given values. The keyword arguments
        must contain all of parameters in self's params. Unrecognized arguments
        are ignored.
        """
        if kwargs in self:
            return self._lognorm
        else:
            return -numpy.inf

    __call__ = logpdf

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
            arr[p] = numpy.random.uniform(self._bounds[p][0],
                                        self._bounds[p][1],
                                        size=size)
        return arr

    @classmethod
    def from_config(cls, cp, section, tag):
        """ Returns a distribution based on a configuration file.

        Parameters
         ----------
         cp : pycbc.workflow.WorkflowConfigParser
             A parsed configuration file that contains the distribution
             options.
         section : str
             Name of the section in the configuration file.
         tag : str
             Name of the tag in the configuration file to use, eg. the
             configuration file has a [section-tag] section.

         Returns
         -------
         Uniform
             A distribution instance from the pycbc.inference.prior module.
        """

        # seperate out tags from section tag
        variable_args = tag.split("+")

        # list of args that are used to construct distribution
        special_args = ["name"] + ["min-%s"%param for param in variable_args] \
                                + ["max-%s"%param for param in variable_args]

        # get a dict with bounds as value
        dist_args = {}
        for param in variable_args:
            low = float(cp.get_opt_tag(section, "min-%s"%param, tag))
            high = float(cp.get_opt_tag(section, "max-%s"%param, tag))
            dist_args[param] = (low,high)

        # add any additional options that user put in that section
        for key in cp.options( "-".join([section,tag]) ):

            # ignore options that are already included
            if key in special_args:
                continue

            # check if option can be cast as a float
            val = cp.get_opt_tag("prior", key, tag)
            try:
                val = float(val)
            except:
                pass

            # add option
            dist_args.update({key:val})

        # construction distribution and add to list
        return cls(**dist_args)

priors = {Uniform.name: Uniform}

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

    Attributes
    ----------
    variable_args : tuple
        The parameters expected when the evaluator is called.
    distributions : list
        The distributions for the parameters.
    """

    def __init__(self, variable_args, *distributions):

        # store the names of the variable params
        self.variable_args = tuple(variable_args)
        # store the distributions
        self.distributions = distributions

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

    def __call__(self, *params):
        """ Evalualate prior for parameters.
        """
        params = dict(zip(self.variable_args, params))
        return sum([d(**params) for d in self.distributions])

