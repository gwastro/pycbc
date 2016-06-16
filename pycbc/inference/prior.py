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
import scipy.stats

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

class Gaussian(Uniform):
    """
    A gaussian distribution on the given parameters. The parameters are
    independent of each other. Instances of this class can be called like
    a function. By default, logpdf will be called, but this can be changed
    by setting the class's __call__ method to its pdf method.

    Parameters
    ----------
    variable_args : {list, str}
        A list of str for each parameter name. Other parameters are list types
        and the index for a particular param is used as a map.
    low : {list, float}
        A list of lower bounds as float.
    high : {list, float}
        A list of higher bounds as float.
    mean : {list, float}
        A list of means as float.
    var : {list, float}
        A list of variances as float.

    Class Attributes
    ----------------
    name : 'guassian'
        The name of this distribution.
    """
    name = "gaussian"
    def __init__(self, variable_args, low, high, mean, var, **kwargs):

        # save distribution parameters as dict
        # calculate the norm and exponential norm ahead of time
        # and save to self._norm, self._lognorm, and self._expnorm
        self._bounds = {}
        self._mean = {}
        self._var = {}
        self._norm = {}
        self._lognorm = {}
        self._expnorm = {}
        for i,param in enumerate(variable_args):
            self._bounds[param] = (low[i],high[i])
            self._mean[param] = mean[i]
            self._var[param] = var[i]
            norm = numpy.sqrt( 2 * self._var[param] * numpy.pi )
            self._norm[param] = 1.0 / norm
            self._lognorm[param] = numpy.log(self._norm[param])
            self._expnorm[param] = 2 * self._var[param]

        # save variable parameters
        self._params = sorted(variable_args)

    @property
    def params(self):
        return self._params

    @property
    def bounds(self):
        return self._bounds

    @property
    def mean(self):
        return self._mean

    @property
    def var(self):
        return self._var

    def pdf(self, **kwargs):
        """ Returns the probability density function (PDF) at the given values.
        The keyword arguments must contain all of parameters in self's params.
        Unrecognized arguments are ignored.

        The PDF is calculated using

            p(x) = \frac{1}{\sqrt{2 \pi \sigma^2}} \exp{- \frac{\left( x - \mu \right^2}{2 \sigma^2} }

        Parameters
        ----------
        kwargs : **dict
            An unpacked dict with variable parameter as key and parameter value as value.

        Returns
        -------
        _pdf : float
            The PDF.
        """
        if kwargs in self:
            _pdf = 1.0
            for param in kwargs.keys():
                if param in self._params:
                    _pdf *= self._norm[param]
                    _pdf *= numpy.exp( -1 * (kwargs[param] - self._mean[param])**2 / self._expnorm[param] )
            return _pdf
        else:
            return 0.0

    def logpdf(self, **kwargs):
       """ Returns the natural logarithm of the probability density function
        (PDF) at the given values. For PDF formula see self._pdf docstring.
        The keyword arguments must contain all of parameters in self's params.
        Unrecognized arguments are ignored.

        Parameters
        ----------
        kwargs : **dict
            An unpacked dict with variable parameter as key and parameter value as value.

        Returns
        -------
        logpdf : float
            The natural logarithm of the PDF.
        """
        if kwargs in self:
            logpdf = 0.0
            for param in kwargs.keys():
                if param in self._params:
                    logpdf += self._lognorm[param]
                    logpdf += -1 * (kwargs[param] - self._mean[param])**2 / self._expnorm[param]
            return logpdf
        else:
            return -numpy.inf

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
            The random values in a numpy array with shape (niterations,ndim).
            If a param was specified, the array will only have an element
            corresponding to the given parameter. Otherwise, the array will
            have an element for each parameter in self's params.
        """
        params = [param] if param != None else self._params
        vals = numpy.zeros(shape=(size,len(params)))
        for i,param in enumerate(params):
            sigma = numpy.sqrt(self._var[param])
            vals[:,i] = scipy.stats.truncnorm.rvs(
                              (self._bounds[param][0]-self._mean[param])/sigma,
                              (self._bounds[param][1]-self._mean[param])/sigma,
                              loc=self._mean[param], scale=sigma, size=size)
        return vals

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
         Gaussian
             A distribution instance from the pycbc.inference.prior module.
        """

        # seperate out tags from section tag
        variable_args = tag.split("+")

        # list of args that are used to construct distribution
        special_args = ["name"] + ["min-%s"%param for param in variable_args] \
                                + ["max-%s"%param for param in variable_args] \
                                + ["mean-%s"%param for param in variable_args] \
                                + ["var-%s"%param for param in variable_args]

        # get input kwargs
        low = []
        high = []
        mean = []
        var = []
        for param in variable_args:
            low.append( float(cp.get_opt_tag(section, "min-%s"%param, tag)) )
            high.append( float(cp.get_opt_tag(section, "max-%s"%param, tag)) )
            mean.append( float(cp.get_opt_tag(section, "mean-%s"%param, tag)) )
            var.append( float(cp.get_opt_tag(section, "var-%s"%param, tag)) )

        # add any additional options that user put in that section
        dist_args = {}
        for key in cp.options( "-".join([section,tag]) ):

            # ignore options that are already included
            if key in special_args:
                continue

            # check if option can be cast as a float then add option
            val = cp.get_opt_tag("prior", key, tag)
            try:
                val = float(val)
            except:
                pass
            dist_args.update({key:val})

        # construction distribution and add to list
        return cls(variable_args, low, high, mean, var, **dist_args)

priors = {
    Uniform.name : Uniform,
    Gaussian.name : Gaussian,
}

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

