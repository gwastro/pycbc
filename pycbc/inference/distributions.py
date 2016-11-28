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
from ConfigParser import Error
import warnings
from pycbc.inference import boundaries

VARARGS_DELIM = '+'

#
#   Distributions for priors
#
def get_param_bounds_from_config(cp, section, tag, param):
    """Gets bounds for the given parameter from a section in a config file.

    Minimum and maximum values for bounds are specified by adding
    `min-{param}` and `max-{param}` options, where `{param}` is the name of
    the paramter. The types of boundary (open, closed, or reflected) to create
    may also be specified by adding options `bytime-min-{param}` and
    `btype-max-{param}`. Cyclic conditions can be adding option
    `cyclic-{param}`. If no `btype` arguments are provided, the
    left bound will be closed and the right open.

    For example, the following will create right-open bounds for parameter
    `foo`:

    .. code::
        [{section}-{tag}]
        min-foo = -1
        max-foo = 1

    This would make the boundaries cyclic:

    .. code::
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


def _bounded_from_config(cls, cp, section, variable_args,
        bounds_required=False):
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

    Returns
    -------
    cls
        An instance of the given class.
    """
    tag = variable_args
    variable_args = variable_args.split(VARARGS_DELIM)

    # list of args that are used to construct distribution
    special_args = ["name"] + \
        ['min-{}'.format(arg) for arg in variable_args] + \
        ['max-{}'.format(arg) for arg in variable_args] + \
        ['btype-min-{}'.format(arg) for arg in variable_args] + \
        ['btype-max-{}'.format(arg) for arg in variable_args] + \
        ['cyclic-{}'.format(arg) for arg in variable_args]

    # get a dict with bounds as value
    dist_args = {}
    for param in variable_args:
        bounds = get_param_bounds_from_config(cp, section, tag, param)
        if bounds_required and bounds is None:
            raise ValueError("min and/or max missing for parameter %s"%(
                param))
        dist_args[param] = bounds

    # add any additional options that user put in that section
    for key in cp.options( "-".join([section,tag]) ):
        # ignore options that are already included
        if key in special_args:
            continue
        # check if option can be cast as a float
        val = cp.get_opt_tag("prior", key, tag)
        try:
            val = float(val)
        except ValueError:
            pass
        # add option
        dist_args.update({key:val})

    # construction distribution and add to list
    return cls(**dist_args)


class _BoundedDist(object):
    """
    A generic class for storing common properties of distributions in which
    each parameter has a minimum and maximum value.

    Parameters
    ----------
    \**params :
        The keyword arguments should provide the names of parameters and their
        corresponding bounds, as either tuples or a `boundaries.Bounds`
        instance.

    Attributes
    ----------
    params : list of strings
        The list of parameter names.
    bounds : dict
        A dictionary of the parameter names and their bounds.
    """
    def __init__(self, **params):
        # convert input bounds to Bounds class, if necessary
        for param,bnds in params.items():
            if not isinstance(bnds, boundaries.Bounds):
                params[param] = boundaries.Bounds(bnds[0], bnds[1])
            # warn the user about reflected boundaries
            if isinstance(bnds, boundaries.Bounds) and (
                    bnds.min.name == 'reflected' or
                    bnds.max.name == 'reflected'):
                warnings.warn("Param {} has one or more ".format(param) +
                              "reflected boundaries. Reflected boundaries "
                              "can cause issues when used in an MCMC.")
        self._bounds = params
        self._params = sorted(params.keys())

    @property
    def params(self):
        return self._params

    @property
    def bounds(self):
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
        _BoundedDist
            A distribution instance from the pycbc.inference.prior module.
        """
        return _bounded_from_config(cls, cp, section, variable_args,
            bounds_required=bounds_required)


class Uniform(_BoundedDist):
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
        numpy.seterr(divide='ignore')
        self._lognorm = -sum([numpy.log(abs(bnd[1]-bnd[0]))
                                    for bnd in self._bounds.values()])
        self._norm = numpy.exp(self._lognorm)
        numpy.seterr(divide='warn')

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
        Uniform
            A distribution instance from the pycbc.inference.prior module.
        """
        return super(Uniform, cls).from_config(cp, section, variable_args,
            bounds_required=True)


class UniformAngle(Uniform):
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

    Class Attributes
    ----------------
    name : 'uniform_angle'
        The name of this distribution.

    Attributes
    ----------
    params : list of strings
        The list of parameter names.
    bounds : dict
        A dictionary of the parameter names and their bounds, in radians.

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
        return _bounded_from_config(cls, cp, section, variable_args,
            bounds_required=False)


class SinAngle(UniformAngle):
    r"""A sine distribution; the pdf of each parameter `\theta` is given by:

    ..math::
        p(\theta) = \frac{\sin \theta}{\cos\theta_0 - \cos\theta_1}, \theta_0 \leq \theta < \theta_1,

    and 0 otherwise. Here, :math:`\theta_0, \theta_1` are the bounds of the
    parameter.

    The domain of this distribution is `[0, pi]`. This is accomplished by
    putting hard boundaries at `[0, pi]`. Bounds may be provided to further
    limit the range for which the pdf has support.  As with `UniformAngle`,
    these are initizliaed as multiples of pi, while the stored bounds are in
    radians.

    Parameters
    ----------
    \**params :
        The keyword arguments should provide the names of parameters and
        (optionally) their corresponding bounds, as either
        `boundaries.Bounds` instances or tuples. The bounds must be
        in [0,1]. These are converted to radians for storage. None may also
        be passed; in that case, the domain bounds will be used.

    Class Attributes
    ----------------
    name : 'sin_angle'
        The name of this distribution.

    Attributes
    ----------
    params : list of strings
        The list of parameter names.
    bounds : dict
        A dictionary of the parameter names and their bounds, in radians.
    """
    name = 'sin_angle'
    _func = numpy.cos
    _dfunc = numpy.sin
    _arcfunc = numpy.arccos
    # _domain applies the reflection off of 0, pi
    _domain = boundaries.Bounds(0, numpy.pi,
        btype_min='closed', btype_max='closed', cyclic=False)


    def _pdf(self, **kwargs):
        """Returns the pdf at the given values. The keyword arguments must
        contain all of parameters in self's params. Unrecognized arguments are
        ignored.
        """
        if kwargs not in self:
            return 0.
        return self._norm * \
            self._dfunc(numpy.array([kwargs[p] for p in self._params])).prod()


    def _logpdf(self, **kwargs):
        """Returns the log of the pdf at the given values. The keyword
        arguments must contain all of parameters in self's params. Unrecognized
        arguments are ignored.
        """
        if kwargs not in self:
            return -numpy.inf
        return self._lognorm + \
            numpy.log(self._dfunc(
                numpy.array([kwargs[p] for p in self._params]))).sum()


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
            arr[p] = self._arcfunc(numpy.random.uniform(
                                    self._func(self._bounds[p][0]),
                                    self._func(self._bounds[p][1]),
                                    size=size))
        return arr


class CosAngle(SinAngle):
    r"""A cosine distribution. This is the same thing as a sine distribution,
    but with the domain shifted to `[-pi/2, pi/2]`. See SinAngle for more
    details.

    Parameters
    ----------
    \**params :
        The keyword arguments should provide the names of parameters and
        (optionally) their corresponding bounds, as either
        `boundaries.Bounds` instances or tuples. The bounds must be
        in [-0.5, 0.5]. These are converted to radians for storage.
        None may also be passed; in that case, the domain bounds will be used.

    Class Attributes
    ----------------
    name : 'cos_angle'
        The name of this distribution.

    Attributes
    ----------
    params : list of strings
        The list of parameter names.
    bounds : dict
        A dictionary of the parameter names and their bounds, in radians.
    """
    name = 'cos_angle'
    _func = numpy.sin
    _dfunc = numpy.cos
    _arcfunc = numpy.arcsin
    _domain = boundaries.Bounds(-numpy.pi/2., numpy.pi/2.,
        btype_min='closed', btype_max='closed', cyclic=False)


class UniformSolidAngle(_BoundedDist):
    """A distribution that is uniform in the solid angle of a sphere. The names
    of the two angluar parameters can be specified on initalization.

    Parameters
    ----------
    polar_angle : {'theta', str}
        The name of the polar angle.
    azimuthal_angle : {'phi', str}
        The name of the azimuthal angle.
    polar_bounds : {None, (min, max)}
        Limit the polar angle to the given bounds. If None provided, the polar
        angle will vary from 0 (the north pole) to pi (the south pole). The
        bounds should be specified as factors of pi. For example, to limit
        the distribution to the northern hemisphere, set
        `polar_bounds=(0,0.5)`.
    azimuthal_bounds : {None, (min, max)}
        Limit the azimuthal angle to the given bounds. If None provided, the
        azimuthal angle will vary from 0 to 2pi. The
        bounds should be specified as factors of pi. For example, to limit
        the distribution to the one hemisphere, set `azimuthal_bounds=(0,1)`.

    Class Attributes
    ----------------
    name : 'uniform_solidangle'
        The name of the distribution.

    Attributes
    ----------
    bounds : dict
        The bounds on each angle. The keys are the names of the polar and
        azimuthal angles, the values are the minimum and maximum of each, in
        radians. For example, if the distribution was initialized with
        `polar_angle='theta', polar_bounds=(0,0.5)` then the bounds will have
        `'theta': 0, 1.5707963267948966` as an entry.
    params : list
        The names of the polar and azimuthal angles.
    polar_angle : str
        The name of the polar angle.
    azimuthal_angle : str
        The name of the azimuthal angle.
    """
    name = 'uniform_solidangle'
    _polardistcls = SinAngle
    _azimuthaldistcls = UniformAngle
    _default_polar_angle = 'theta'
    _default_azimuthal_angle = 'phi'

    def __init__(self, polar_angle=_default_polar_angle,
                 azimuthal_angle=_default_azimuthal_angle,
                 polar_bounds=None, azimuthal_bounds=None):
        self._polardist = self._polardistcls(**{
            polar_angle: polar_bounds}) 
        self._azimuthaldist = self._azimuthaldistcls(**{
            azimuthal_angle: azimuthal_bounds})
        self._polar_angle = polar_angle
        self._azimuthal_angle = azimuthal_angle
        self._bounds = dict(self._polardist.bounds.items() +
                            self._azimuthaldist.bounds.items())
        self._params = sorted(self._bounds.keys())


    @property
    def polar_angle(self):
        return self._polar_angle


    @property
    def azimuthal_angle(self):
        return self._azimuthal_angle


    def apply_boundary_conditions(self, **kwargs):
        """Maps the given values to be within the domain of the azimuthal and
        polar angles, before applying any other boundary conditions.
        
        Parameters
        ----------
        \**kwargs :
            The keyword args must include values for both the azimuthal and
            polar angle, using the names they were initilialized with. For
            example, if `polar_angle='theta'` and `azimuthal_angle=`phi`, then
            the keyword args must be `theta={val1}, phi={val2}`.

        Returns
        -------
        dict
            A dictionary of the parameter names and the conditioned values.
        """
        polarval = kwargs[self._polar_angle]
        azval = kwargs[self._azimuthal_angle]
        # constrain each angle to its domain
        polarval = self._polardist._domain.apply_conditions(polarval)
        azval = self._azimuthaldist._domain.apply_conditions(azval)
        # apply any other boundary conditions
        polarval = self._bounds[self._polar_angle].apply_conditions(polarval)
        azval = self._bounds[self._azimuthal_angle].apply_conditions(azval)
        return {self._polar_angle: polarval, self._azimuthal_angle: azval}


    def _pdf(self, **kwargs):
        """
        Returns the pdf at the given angles.

        Parameters
        ----------
        \**kwargs:
            The keyword arguments should specify the value for each angle,
            using the names of the polar and azimuthal angles as the keywords.
            Unrecognized arguments are ignored.

        Returns
        -------
        float
            The value of the pdf at the given values.
        """
        return self._polardist._pdf(**kwargs) * \
            self._azimuthaldist._pdf(**kwargs)
        

    def _logpdf(self, **kwargs):
        """
        Returns the logpdf at the given angles.

        Parameters
        ----------
        \**kwargs:
            The keyword arguments should specify the value for each angle,
            using the names of the polar and azimuthal angles as the keywords.
            Unrecognized arguments are ignored.

        Returns
        -------
        float
            The value of the pdf at the given values.
        """
        return self._polardist._logpdf(**kwargs) +\
            self._azimuthaldist._logpdf(**kwargs)


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
            if p == self._polar_angle:
                arr[p] = self._polardist.rvs(size=size)
            elif p == self._azimuthal_angle:
                arr[p] = self._azimuthaldist.rvs(size=size)
            else:
                raise ValueError("unrecognized parameter %s" %(p))
        return arr

    @classmethod
    def from_config(cls, cp, section, variable_args):
        """Returns a distribution based on a configuration file. The section
        must have the names of the polar and azimuthal angles in the tag part
        of the section header. For example:

        .. code-block::
            [prior-theta+phi]
            name = uniform_solidangle

        If nothing else is provided, the default names and bounds of the polar
        and azimuthal angles will be used. To specify a different name for
        each angle, set the `polar-angle` and `azimuthal-angle` attributes. For
        example: 

        ..code-block::
            [prior-foo+bar]
            name = uniform_solidangle
            polar-angle = foo
            azimuthal-angle = bar
        
        Note that the names of the variable args in the tag part of the section
        name must match the names of the polar and azimuthal angles.

        Bounds may also be specified for each angle, as factors of pi. For
        example:

        .. code-block::
            [prior-theta+phi]
            polar-angle = theta
            azimuthal-angle = phi
            min-theta = 0
            max-theta = 0.5

        This will return a distribution that is uniform in the upper
        hemisphere.

        Parameters
        ----------
        cp : ConfigParser instance
            The config file.
        section : str
            The name of the section.
        variable_args : str
            The names of the parameters for this distribution, separated by
            `prior.VARARGS_DELIM`. These must appear in the "tag" part
            of the section header.

        Returns
        -------
        UniformSolidAngle
            A distribution instance from the pycbc.inference.prior module.
        """
        tag = variable_args
        variable_args = variable_args.split(VARARGS_DELIM)

        # get the variables that correspond to the polar/azimuthal angles
        try:
            polar_angle = cp.get_opt_tag(section, 'polar-angle', tag)
        except Error:
            polar_angle = cls._default_polar_angle
        try:
            azimuthal_angle = cp.get_opt_tag(section, 'azimuthal-angle', tag)
        except Error:
            azimuthal_angle = cls._default_azimuthal_angle

        if polar_angle not in variable_args:
            raise Error("polar-angle %s is not one of the variable args (%s)"%(
                polar_angle, ', '.join(variable_args)))
        if azimuthal_angle not in variable_args:
            raise Error("azimuthal-angle %s is not one of the variable args "%(
                azimuthal_angle) + "(%s)"%(', '.join(variable_args)))

        # get the bounds, if provided
        polar_bounds = get_param_bounds_from_config(cp, section, tag,
            polar_angle)
        azimuthal_bounds = get_param_bounds_from_config(cp, section, tag,
            azimuthal_angle)

        return cls(polar_angle=polar_angle, azimuthal_angle=azimuthal_angle,
            polar_bounds=polar_bounds, azimuthal_bounds=azimuthal_bounds)


class UniformSky(UniformSolidAngle):
    """A distribution that is uniform on the sky. This is the same as
    UniformSolidAngle, except that the polar angle varies from pi/2 (the north
    pole) to -pi/2 (the south pole) instead of 0 to pi. Also, the default
    names are "dec" (declination) for the polar angle and "ra" (right
    ascension) for the azimuthal angle, instead of "theta" and "phi".
    """
    name = 'uniform_sky'
    _polardistcls = CosAngle
    _default_polar_angle = 'dec'
    _default_azimuthal_angle = 'ra'


class Gaussian(_BoundedDist):
    r"""A Gaussian distribution on the given parameters; the parameters are
    independent of each other.
    
    Bounds can be provided on each parameter, in which case the distribution
    will be a truncated Gaussian distribution.  The PDF of a truncated
    Gaussian distribution is given by:

    .. math::
        p(x|a, b, \mu,\sigma) = \frac{1}{\sqrt{2 \pi \sigma^2}}\frac{e^{- \frac{\left( x - \mu \right)^2}{2 \sigma^2}}}{\Phi(b|\mu, \sigma) - \Phi(a|\mu, \sigma)},

    where :math:`\mu` is the mean, :math:`\sigma^2` is the variance,
    :math:`a,b` are the bounds, and :math:`\Phi` is the cumulative distribution
    of an unbounded normal distribution, given by:

    .. math::
        \Phi(x|\mu, \sigma) = \frac{1}{2}\left[1 + \mathrm{erf}\left(\frac{x-\mu}{\sigma \sqrt{2}}\right)\right].

    Note that if :math:`[a,b) = [-\infty, \infty)`, this reduces to a standard
    Gaussian distribution.

    
    Instances of this class can be called like a function. By default, logpdf
    will be called, but this can be changed by setting the class's __call__
    method to its pdf method.

    Parameters
    ----------
    \**params :
        The keyword arguments should provide the names of parameters and
        (optionally) some bounds, as either a tuple or a
        `boundaries.Bounds` instance. The mean and variance of each
        parameter can be provided by additional keyword arguments that have
        `_mean` and `_var` adding to the parameter name. For example,
        `foo=(-2,10), foo_mean=3, foo_var=2` would create a truncated Gaussian
        with mean 3 and variance 2, bounded between :math:`[-2, 10)`. If no
        mean or variance is provided, the distribution will have 0 mean and
        unit variance. If None is provided for the bounds, the distribution
        will be a normal, unbounded Gaussian (equivalent to setting the bounds
        to `[-inf, inf)`).

    Class Attributes
    ----------------
    name : 'guassian'
        The name of this distribution.

    Example
    -------
    Create an unbounded Gaussian distribution with zero mean and unit variance:
    >>> dist = distributions.Gaussian(mass1=None)

    Create a bounded Gaussian distribution on :math:`[1,10)` with a mean of 3
    and a variance of 2:
    >>> dist = distributions.Gaussian(mass1=(1,10), mass1_mean=3, mass1_var=2)
   
    Create a bounded Gaussian distribution with the same parameters, but with
    cyclic boundary conditions:
    >>> dist = distributions.Gaussian(mass1=Bounds(1,10, cyclic=True), mass1_mean=3, mass1_var=2)
    """
    name = "gaussian"

    def __init__(self, **params):

        # save distribution parameters as dict
        # calculate the norm and exponential norm ahead of time
        # and save to self._norm, self._lognorm, and self._expnorm
        self._bounds = {}
        self._mean = {}
        self._var = {}
        self._norm = {}
        self._lognorm = {}
        self._expnorm = {}
        # pull out specified means, variance
        mean_args = [p for p in params if p.endswith('_mean')]
        var_args = [p for p in params if p.endswith('_var')]
        self._mean = dict([[p[:-5], params.pop(p)] for p in mean_args])
        self._var = dict([[p[:-4], params.pop(p)] for p in var_args])
        # if any param is set to None, make its bounds -inf, inf
        for param in params:
            if params[param] is None:
                # Bounds defaults to -inf, inf
                params[param] = boundaries.Bounds()
        # initialize the bounds
        super(Gaussian, self).__init__(**params)

        # check that there are no params in mean/var that are not in params
        missing = set(self._mean.keys()) - set(params.keys())
        if any(missing):
            raise ValueError("means provided for unknow params {}".format(
                ', '.join(missing)))
        missing = set(self._var.keys()) - set(params.keys())
        if any(missing):
            raise ValueError("vars provided for unknow params {}".format(
                ', '.join(missing)))
        # set default mean/var for params not specified
        self._mean.update(dict([[p, 0.]
            for p in params if p not in self._mean]))
        self._var.update(dict([[p, 1.]
            for p in params if p not in self._var]))

        # compute norms
        for p,bnds in self._bounds.items():
            sigmasq = self._var[p]
            mu = self._mean[p]
            a,b = bnds
            invnorm = scipy.stats.norm.cdf(b, loc=mu, scale=sigmasq**0.5) \
                    - scipy.stats.norm.cdf(a, loc=mu, scale=sigmasq**0.5)
            invnorm *= numpy.sqrt(2*numpy.pi*sigmasq)
            self._norm[p] = 1./invnorm
            self._lognorm[p] = numpy.log(self._norm[p])
            self._expnorm[p] = -1./(2*sigmasq)


    @property
    def mean(self):
        return self._mean


    @property
    def var(self):
        return self._var


    def _pdf(self, **kwargs):
        """Returns the pdf at the given values. The keyword arguments must
        contain all of parameters in self's params. Unrecognized arguments are
        ignored.
        """
        return numpy.exp(self._logpdf(**kwargs))


    def _logpdf(self, **kwargs):
        """Returns the log of the pdf at the given values. The keyword
        arguments must contain all of parameters in self's params. Unrecognized
        arguments are ignored.
        """
        if kwargs in self:
            return sum([self._lognorm[p] +
                        self._expnorm[p]*(kwargs[p]-self._mean[p])**2.
                        for p in self._params])
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
            sigma = numpy.sqrt(self._var[p])
            mu = self._mean[p]
            a,b = self._bounds[p]
            arr[p][:] = scipy.stats.truncnorm.rvs((a-mu)/sigma, (b-mu)/sigma,
                loc=self._mean[p], scale=sigma, size=size)
        return arr


    @classmethod
    def from_config(cls, cp, section, variable_args):
        """Returns a Gaussian distribution based on a configuration file. The
        parameters for the distribution are retrieved from the section titled
        "[`section`-`variable_args`]" in the config file.

        Boundary arguments should be provided in the same way as described in
        `get_param_bounds_from_config`. In addition, the mean and variance of
        each parameter can be specified by setting `{param}_mean` and
        `{param}_var`, respectively. For example, the following would create a
        truncated Gaussian distribution between 0 and 6.28 for a parameter
        called `phi` with mean 3.14 and variance 0.5 that is cyclic:

        .. code::
            [{section}-{tag}]
            min-phi = 0
            max-phi = 6.28
            phi_mean = 3.14
            phi_var = 0.5
            cyclic =

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
        Gaussain
            A distribution instance from the pycbc.inference.prior module.
        """
        return _bounded_from_config(cls, cp, section, variable_args,
            bounds_required=False)


distribs = {
    Uniform.name : Uniform,
    UniformAngle.name : UniformAngle,
    CosAngle.name : CosAngle,
    SinAngle.name : SinAngle,
    UniformSolidAngle.name : UniformSolidAngle,
    UniformSky.name : UniformSky,
    Gaussian.name : Gaussian,
}

def read_distributions_from_config(cp, section="prior"):
    """Returns a list of PyCBC distribution instances for a section in the
    given configuration file.

    Parameters
    ----------
    cp : WorflowConfigParser
        An open config file to read.
    section : {"prior", string}
        Prefix on section names from which to retrieve the distributions.

    Returns
    -------
    list
        A list of the parsed distributions.
    """
    dists = []
    variable_args = []
    for subsection in cp.get_subsections(section):
        name = cp.get_opt_tag(section, "name", subsection)
        dist = distribs[name].from_config(cp, section, subsection)
        if set(dist.params).isdisjoint(variable_args):
            dists.append(dist)
            variable_args += dist.params
        else:
            raise ValueError("Same parameter in more than one distribution.")
    return dists
