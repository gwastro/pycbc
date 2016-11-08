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
from pycbc.inference import boundary_utils

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

    This would result in two reflected boundaries:

    .. code::
        [{section}-{tag}]
        min-foo = -1
        max-foo = 1
        btype-min-foo = reflected
        btype-max-foo = reflected
    
    For more details on boundary types and their meaning, see
    `boundary_utils.Bounds`.

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
        If bounds were provided, a `boundary_utils.Bounds` instance
        representing the bounds. Otherwise, `None`.
    """
    try:
        minbnd = cp.get_opt_tag(section, 'min-'+param, tag)
    except Error:
        minbnd = None
    try:
        maxbnd = cp.get_opt_tag(section, 'max-'+param, tag)
    except Error:
        maxbnd = None
    if minbnd is None and maxbnd is None:
        bnds = None
    elif minbnd is None or maxbnd is None:
        raise ValueError("if specifying bounds for %s, " %(param) +
            "you must provide both a minimum and a maximum")
    else:
        bndargs = {'min-bound': minbnd, 'max-bnd': maxbnd}
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
        bndargs.update({'btype-min': minbtype, 'btype-max': maxbtype})
        cyclic = cp.has_option_tag(section, 'cyclic-{}'.format(param), tag)
        bndargs.update({'cyclic': cyclic})
        bnds = boundary_utils.Bounds(**bndargs)
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
        corresponding bounds, as tuples.

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
            if not isinstance(bnds, boundary_utils.Bounds):
                params[param] = boundary_utils.Bounds(bnds[0], bnds[1])
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
        corresponding value. See `boundary_utils.Bounds.apply_conditions` for
        details.

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
        return dict([[p, self._bounds[p].apply_conditions(val)]
                     for p,val in kwargs.items()])

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
        `boundary_utils.Bounds` instances or tuples. The bounds must be
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
    _domain = boundary_utils.Bounds(0., 2*numpy.pi, cyclic=True)

    def __init__(self, **params):
        for p,bnds in params.items():
            if bnds is None:
                bnds = self._domain
            elif isinstance(bnds, boundary_utils.Bounds):
                # convert to radians
                bnds._min *= numpy.pi
                bnds._max *= numpy.pi
            else:
                # create a Bounds instance from the given tuple
                bnds = boundary_utils.Bounds(
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
    reflecting a value between `0` and `pi` until it falls between `[0, pi]`.
    For example, if `pdf` or `logpdf` are evaluated at `3.25 pi`, the value is
    mapped to `0.75 pi` (`3.25 pi -> -1.25 pi -> 0.75 pi`) before being
    evaluated.

    Bounds may be provided to limit the range for which the pdf has support.
    As with `UniformAngle`, these are initizliaed as multiples of pi, while
    the stored bounds are in radians.

    Parameters
    ----------
    \**params :
        The keyword arguments should provide the names of parameters and
        (optionally) their corresponding bounds, as either
        `boundary_utils.Bounds` instances or tuples. The bounds must be
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
    _domain = boundary_utils.Bounds(0, numpy.pi,
        btype_min='reflected', btype_max='reflected', cyclic=False)


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
        `boundary_utils.Bounds` instances or tuples. The bounds must be
        in [-0.5, 0.5]. These are converted to radians for storage.
        None may also be passed; in that case, the domain bounds will be used.

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
    name = 'cos_angle'
    _func = numpy.sin
    _dfunc = numpy.cos
    _arcfunc = numpy.arcsin
    _domain = boundary_utils.Bounds(-numpy.pi/2., numpy.pi/2.,
        btype_min='reflected', btype_max='reflected', cyclic=False)


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
    """
    A gaussian distribution on the given parameters. The parameters are
    independent of each other. Instances of this class can be called like
    a function. By default, logpdf will be called, but this can be changed
    by setting the class's __call__ method to its pdf method.

    Parameters
    ----------
    \**params :
        The keyword arguments should provide the names of parameters and their
        corresponding bounds, mean, and variance, as tuples,
        eg. parameter=(low,high,mean,var).

    Class Attributes
    ----------------
    name : 'guassian'
        The name of this distribution.

    Example
    -------
    >> dist = Gaussian(mass1=(1.3,1.5,1.4,0.01))
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
        for param in params.keys():
            self._bounds[param] = (params[param][0], params[param][1])
            self._mean[param] = params[param][2]
            self._var[param] = params[param][3]
            norm = numpy.sqrt( 2 * self._var[param] * numpy.pi )
            self._norm[param] = 1.0 / norm
            self._lognorm[param] = numpy.log(norm)
            self._expnorm[param] = 2 * self._var[param]

        # save variable parameters
        self._params = sorted(params.keys())

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
                    logpdf -= (kwargs[param] - self._mean[param])**2 / self._expnorm[param]
            return logpdf
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
        vals : numpy array
            The random values in a numpy.array with shape (niterations,ndim).
            If a param was specified, the array will only have an element
            corresponding to the given parameter. Otherwise, the array will
            have an element for each parameter in self's params.
        """
        params = [param] if param is not None else self._params
        vals = numpy.zeros(shape=(size,len(params)))
        for i,param in enumerate(params):
            sigma = numpy.sqrt(self._var[param])
            vals[:,i] = scipy.stats.truncnorm.rvs(
                              (self._bounds[param][0]-self._mean[param])/sigma,
                              (self._bounds[param][1]-self._mean[param])/sigma,
                              loc=self._mean[param], scale=sigma, size=size)
        return vals

    @classmethod
    def from_config(cls, cp, section, variable_args):
        """ Returns a distribution based on a configuration file.

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
        Gaussian
            A distribution instance from the pycbc.inference.prior module.
        """
        tag = variable_args
        variable_args = variable_args.split(VARARGS_DELIM)

        # list of args that are used to construct distribution
        special_args = ["name"] + ["min-%s"%param for param in variable_args] \
                                + ["max-%s"%param for param in variable_args] \
                                + ["mean-%s"%param for param in variable_args] \
                                + ["var-%s"%param for param in variable_args]

        # get input kwargs
        params = {}
        for param in variable_args:
            params[param] = [None, None, None, None]
            params[param][0] = float(cp.get_opt_tag(section,
                                                    "min-%s" % param, tag))
            params[param][1] = float(cp.get_opt_tag(section,
                                                    "max-%s" % param, tag))
            params[param][2] = float(cp.get_opt_tag(section,
                                                    "mean-%s" % param, tag))
            params[param][3] = float(cp.get_opt_tag(section,
                                                    "var-%s" % param, tag))

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
            except ValueError:
                pass
            dist_args.update({key:val})

        # construction distribution and add to list
        return cls(**params)

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
