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

    Attributes
    ----------------
    name : 'sin_angle'
        The name of this distribution.
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

    Attributes
    ----------------
    name : 'cos_angle'
        The name of this distribution.
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

class UniformSolidAngle(bounded.BoundedDist):
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

    Attributes
    ----------------
    name : 'uniform_solidangle'
        The name of the distribution.
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

        .. code-block:: ini

            [prior-theta+phi]
            name = uniform_solidangle

        If nothing else is provided, the default names and bounds of the polar
        and azimuthal angles will be used. To specify a different name for
        each angle, set the `polar-angle` and `azimuthal-angle` attributes. For
        example: 

        .. code-block:: ini

            [prior-foo+bar]
            name = uniform_solidangle
            polar-angle = foo
            azimuthal-angle = bar
        
        Note that the names of the variable args in the tag part of the section
        name must match the names of the polar and azimuthal angles.

        Bounds may also be specified for each angle, as factors of pi. For
        example:

        .. code-block:: ini

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
        variable_args = variable_args.split(bounded.VARARGS_DELIM)

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

