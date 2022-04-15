# Copyright (C) 2016  Collin Capano, Christopher M. Biwer
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
This modules provides classes for evaluating distributions with bounds.
"""

import warnings
from configparser import Error

import numpy
from pycbc import boundaries
from pycbc import VARARGS_DELIM

#
#   Distributions for priors
#
def get_param_bounds_from_config(cp, section, tag, param):
    """Gets bounds for the given parameter from a section in a config file.

    Minimum and maximum values for bounds are specified by adding
    `min-{param}` and `max-{param}` options, where `{param}` is the name of
    the parameter. The types of boundary (open, closed, or reflected) to create
    may also be specified by adding options `btype-min-{param}` and
    `btype-max-{param}`. Cyclic conditions can be adding option
    `cyclic-{param}`. If no `btype` arguments are provided, the
    left bound will be closed and the right open.

    For example, the following will create right-open bounds for parameter
    `foo`:

    .. code-block:: ini

        [{section}-{tag}]
        min-foo = -1
        max-foo = 1

    This would make the boundaries cyclic:

    .. code-block:: ini

        [{section}-{tag}]
        min-foo = -1
        max-foo = 1
        cyclic-foo =

    For more details on boundary types and their meaning, see
    `boundaries.Bounds`.

    If the parameter is not found in the section will just return None (in
    this case, all `btype` and `cyclic` arguments are ignored for that
    parameter).  If bounds are specified, both a minimum and maximum must be
    provided, else a Value or Type Error will be raised.

    Parameters
    ----------
    cp : ConfigParser instance
        The config file.
    section : str
        The name of the section.
    tag : str
        Any tag in the section name. The full section name searched for in
        the config file is `{section}(-{tag})`.
    param : str
        The name of the parameter to retrieve bounds for.

    Returns
    -------
    bounds : {Bounds instance | None}
        If bounds were provided, a `boundaries.Bounds` instance
        representing the bounds. Otherwise, `None`.
    """
    try:
        minbnd = float(cp.get_opt_tag(section, 'min-'+param, tag))
    except Error:
        minbnd = None
    try:
        maxbnd = float(cp.get_opt_tag(section, 'max-'+param, tag))
    except Error:
        maxbnd = None
    if minbnd is None and maxbnd is None:
        bnds = None
    elif minbnd is None or maxbnd is None:
        raise ValueError("if specifying bounds for %s, " %(param) +
            "you must provide both a minimum and a maximum")
    else:
        bndargs = {'min_bound': minbnd, 'max_bound': maxbnd}
        # try to get  any other conditions, if provided
        try:
            minbtype = cp.get_opt_tag(section, 'btype-min-{}'.format(param),
                                      tag)
        except Error:
            minbtype = 'closed'
        try:
            maxbtype = cp.get_opt_tag(section, 'btype-max-{}'.format(param),
                                      tag)
        except Error:
            maxbtype = 'open'
        bndargs.update({'btype_min': minbtype, 'btype_max': maxbtype})
        cyclic = cp.has_option_tag(section, 'cyclic-{}'.format(param), tag)
        bndargs.update({'cyclic': cyclic})
        bnds = boundaries.Bounds(**bndargs)
    return bnds


def bounded_from_config(cls, cp, section, variable_args,
        bounds_required=False, additional_opts=None):
    """Returns a bounded distribution based on a configuration file. The
    parameters for the distribution are retrieved from the section titled
    "[`section`-`variable_args`]" in the config file.

    Parameters
    ----------
    cls : pycbc.prior class
        The class to initialize with.
    cp : pycbc.workflow.WorkflowConfigParser
        A parsed configuration file that contains the distribution
        options.
    section : str
        Name of the section in the configuration file.
    variable_args : str
        The names of the parameters for this distribution, separated by
        `prior.VARARGS_DELIM`. These must appear in the "tag" part
        of the section header.
    bounds_required : {False, bool}
       If True, raise a ValueError if a min and max are not provided for
       every parameter. Otherwise, the prior will be initialized with the
       parameter set to None. Even if bounds are not required, a
       ValueError will be raised if only one bound is provided; i.e.,
       either both bounds need to provided or no bounds.
    additional_opts : {None, dict}
        Provide additional options to be passed to the distribution class;
        should be a dictionary specifying option -> value. If an option is
        provided that also exists in the config file, the value provided will
        be used instead of being read from the file.

    Returns
    -------
    cls
        An instance of the given class.
    """
    tag = variable_args
    variable_args = variable_args.split(VARARGS_DELIM)

    if additional_opts is None:
        additional_opts = {}

    # list of args that are used to construct distribution
    special_args = ["name"] + \
        ['min-{}'.format(arg) for arg in variable_args] + \
        ['max-{}'.format(arg) for arg in variable_args] + \
        ['btype-min-{}'.format(arg) for arg in variable_args] + \
        ['btype-max-{}'.format(arg) for arg in variable_args] + \
        ['cyclic-{}'.format(arg) for arg in variable_args] + \
        list(additional_opts.keys())

    # get a dict with bounds as value
    dist_args = {}
    for param in variable_args:
        bounds = get_param_bounds_from_config(cp, section, tag, param)
        if bounds_required and bounds is None:
            raise ValueError("min and/or max missing for parameter %s"%(
                param))
        dist_args[param] = bounds

    # add any additional options that user put in that section
    for key in cp.options("-".join([section, tag])):

        # ignore options that are already included
        if key in special_args:
            continue

        # check if option can be cast as a float
        val = cp.get_opt_tag(section, key, tag)
        try:
            val = float(val)
        except ValueError:
            pass

        # add option
        dist_args.update({key:val})

    dist_args.update(additional_opts)

    # construction distribution and add to list
    return cls(**dist_args)


class BoundedDist(object):
    """
    A generic class for storing common properties of distributions in which
    each parameter has a minimum and maximum value.

    Parameters
    ----------
    \**params :
        The keyword arguments should provide the names of parameters and their
        corresponding bounds, as either tuples or a `boundaries.Bounds`
        instance.
    """
    def __init__(self, **params):
        # convert input bounds to Bounds class, if necessary
        for param,bnds in params.items():
            if bnds is None:
                params[param] = boundaries.Bounds()
            elif not isinstance(bnds, boundaries.Bounds):
                params[param] = boundaries.Bounds(bnds[0], bnds[1])
            # warn the user about reflected boundaries
            if isinstance(bnds, boundaries.Bounds) and (
                    bnds.min.name == 'reflected' or
                    bnds.max.name == 'reflected'):
                warnings.warn("Param {} has one or more ".format(param) +
                              "reflected boundaries. Reflected boundaries "
                              "can cause issues when used in an MCMC.")
        self._bounds = params
        self._params = sorted(list(params.keys()))

    @property
    def params(self):
        """list of strings: The list of parameter names."""
        return self._params

    @property
    def bounds(self):
        """dict: A dictionary of the parameter names and their bounds."""
        return self._bounds

    def __contains__(self, params):
        try:
            return all(self._bounds[p].contains_conditioned(params[p])
                       for p in self._params)
        except KeyError:
            raise ValueError("must provide all parameters [%s]" %(
                ', '.join(self._params)))

    def apply_boundary_conditions(self, **kwargs):
        """Applies any boundary conditions to the given values (e.g., applying
        cyclic conditions, and/or reflecting values off of boundaries). This
        is done by running `apply_conditions` of each bounds in self on the
        corresponding value. See `boundaries.Bounds.apply_conditions` for
        details.

        Parameters
        ----------
        \**kwargs :
            The keyword args should be the name of a parameter and value to
            apply its boundary conditions to. The arguments need not include
            all of the parameters in self. Any unrecognized arguments are
            ignored.

        Returns
        -------
        dict
            A dictionary of the parameter names and the conditioned values.
        """
        return dict([[p, self._bounds[p].apply_conditions(val)]
                     for p,val in kwargs.items() if p in self._bounds])

    def pdf(self, **kwargs):
        """Returns the pdf at the given values. The keyword arguments must
        contain all of parameters in self's params. Unrecognized arguments are
        ignored. Any boundary conditions are applied to the values before the
        pdf is evaluated.
        """
        return self._pdf(**self.apply_boundary_conditions(**kwargs))

    def _pdf(self, **kwargs):
        """The underlying pdf function called by `self.pdf`. This must be set
        by any class that inherits from this class. Otherwise, a
        `NotImplementedError` is raised.
        """
        raise NotImplementedError("pdf function not set")

    def logpdf(self, **kwargs):
        """Returns the log of the pdf at the given values. The keyword
        arguments must contain all of parameters in self's params.
        Unrecognized arguments are ignored. Any boundary conditions are
        applied to the values before the pdf is evaluated.
        """
        return self._logpdf(**self.apply_boundary_conditions(**kwargs))

    def _logpdf(self, **kwargs):
        """The underlying log pdf function called by `self.logpdf`. This must
        be set by any class that inherits from this class. Otherwise, a
        `NotImplementedError` is raised.
        """
        raise NotImplementedError("pdf function not set")

    __call__ = logpdf

    def _cdfinv_param(self, param, value):
        """Return the cdfinv for a single given parameter """
        raise NotImplementedError("inverse cdf not set")

    def cdfinv(self, **kwds):
        """Return the inverse cdf to map the unit interval to parameter bounds.
        You must provide a keyword for every parameter.
        """
        updated = {}
        for param in self.params:
            updated[param] = self._cdfinv_param(param, kwds[param])
        return updated

    def rvs(self, size=1, **kwds):
        "Draw random value"
        dtype = [(p, float) for p in self.params]
        arr = numpy.zeros(size, dtype=dtype)
        draw = {}
        for param in self.params:
            draw[param] = numpy.random.uniform(0, 1, size=size)
        exp = self.cdfinv(**draw)
        for param in self.params:
            arr[param] = exp[param]
        return arr

    @classmethod
    def from_config(cls, cp, section, variable_args, bounds_required=False):
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
        bounds_required : {False, bool}
           If True, raise a ValueError if a min and max are not provided for
           every parameter. Otherwise, the prior will be initialized with the
           parameter set to None. Even if bounds are not required, a
           ValueError will be raised if only one bound is provided; i.e.,
           either both bounds need to provided or no bounds.

        Returns
        -------
        BoundedDist
            A distribution instance from the pycbc.distribution subpackage.
        """
        return bounded_from_config(cls, cp, section, variable_args,
                                    bounds_required=bounds_required)
