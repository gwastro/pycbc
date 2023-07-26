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


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#

"""Base class for models.
"""

import numpy
import logging
from abc import (ABCMeta, abstractmethod)
from configparser import NoSectionError
from pycbc import (transforms, distributions)
from pycbc.io import FieldArray


#
# =============================================================================
#
#                               Support classes
#
# =============================================================================
#


class _NoPrior(object):
    """Dummy class to just return 0 if no prior is given to a model.
    """
    @staticmethod
    def apply_boundary_conditions(**params):
        return params

    def __call__(self, **params):
        return 0.


class ModelStats(object):
    """Class to hold model's current stat values."""

    @property
    def statnames(self):
        """Returns the names of the stats that have been stored."""
        return list(self.__dict__.keys())

    def getstats(self, names, default=numpy.nan):
        """Get the requested stats as a tuple.

        If a requested stat is not an attribute (implying it hasn't been
        stored), then the default value is returned for that stat.

        Parameters
        ----------
        names : list of str
            The names of the stats to get.
        default : float, optional
            What to return if a requested stat is not an attribute of self.
            Default is ``numpy.nan``.

        Returns
        -------
        tuple
            A tuple of the requested stats.
        """
        return tuple(getattr(self, n, default) for n in names)

    def getstatsdict(self, names, default=numpy.nan):
        """Get the requested stats as a dictionary.

        If a requested stat is not an attribute (implying it hasn't been
        stored), then the default value is returned for that stat.

        Parameters
        ----------
        names : list of str
            The names of the stats to get.
        default : float, optional
            What to return if a requested stat is not an attribute of self.
            Default is ``numpy.nan``.

        Returns
        -------
        dict
            A dictionary of the requested stats.
        """
        return dict(zip(names, self.getstats(names, default=default)))


class SamplingTransforms(object):
    """Provides methods for transforming between sampling parameter space and
    model parameter space.
    """

    def __init__(self, variable_params, sampling_params,
                 replace_parameters, sampling_transforms):
        assert len(replace_parameters) == len(sampling_params), (
            "number of sampling parameters must be the "
            "same as the number of replace parameters")
        # pull out the replaced parameters
        self.sampling_params = [arg for arg in variable_params
                                if arg not in replace_parameters]
        # add the sampling parameters
        self.sampling_params += sampling_params
        # sort to make sure we have a consistent order
        self.sampling_params.sort()
        self.sampling_transforms = sampling_transforms

    def logjacobian(self, **params):
        r"""Returns the log of the jacobian needed to transform pdfs in the
        ``variable_params`` parameter space to the ``sampling_params``
        parameter space.

        Let :math:`\mathbf{x}` be the set of variable parameters,
        :math:`\mathbf{y} = f(\mathbf{x})` the set of sampling parameters, and
        :math:`p_x(\mathbf{x})` a probability density function defined over
        :math:`\mathbf{x}`.
        The corresponding pdf in :math:`\mathbf{y}` is then:

        .. math::

            p_y(\mathbf{y}) =
                p_x(\mathbf{x})\left|\mathrm{det}\,\mathbf{J}_{ij}\right|,

        where :math:`\mathbf{J}_{ij}` is the Jacobian of the inverse transform
        :math:`\mathbf{x} = g(\mathbf{y})`. This has elements:

        .. math::

            \mathbf{J}_{ij} = \frac{\partial g_i}{\partial{y_j}}

        This function returns
        :math:`\log \left|\mathrm{det}\,\mathbf{J}_{ij}\right|`.


        Parameters
        ----------
        \**params :
            The keyword arguments should specify values for all of the variable
            args and all of the sampling args.

        Returns
        -------
        float :
            The value of the jacobian.
        """
        return numpy.log(abs(transforms.compute_jacobian(
            params, self.sampling_transforms, inverse=True)))

    def apply(self, samples, inverse=False):
        """Applies the sampling transforms to the given samples.

        Parameters
        ----------
        samples : dict or FieldArray
            The samples to apply the transforms to.
        inverse : bool, optional
            Whether to apply the inverse transforms (i.e., go from the sampling
            args to the ``variable_params``). Default is False.

        Returns
        -------
        dict or FieldArray
            The transformed samples, along with the original samples.
        """
        return transforms.apply_transforms(samples, self.sampling_transforms,
                                           inverse=inverse)

    @classmethod
    def from_config(cls, cp, variable_params):
        """Gets sampling transforms specified in a config file.

        Sampling parameters and the parameters they replace are read from the
        ``sampling_params`` section, if it exists. Sampling transforms are
        read from the ``sampling_transforms`` section(s), using
        ``transforms.read_transforms_from_config``.

        An ``AssertionError`` is raised if no ``sampling_params`` section
        exists in the config file.

        Parameters
        ----------
        cp : WorkflowConfigParser
            Config file parser to read.
        variable_params : list
            List of parameter names of the original variable params.

        Returns
        -------
        SamplingTransforms
            A sampling transforms class.
        """
        # Check if a sampling_params section is provided
        try:
            sampling_params, replace_parameters = \
                read_sampling_params_from_config(cp)
        except NoSectionError as e:
            logging.warning("No sampling_params section read from config file")
            raise e
        # get sampling transformations
        sampling_transforms = transforms.read_transforms_from_config(
            cp, 'sampling_transforms')
        logging.info("Sampling in {} in place of {}".format(
            ', '.join(sampling_params), ', '.join(replace_parameters)))
        return cls(variable_params, sampling_params,
                   replace_parameters, sampling_transforms)


def read_sampling_params_from_config(cp, section_group=None,
                                     section='sampling_params'):
    """Reads sampling parameters from the given config file.

    Parameters are read from the `[({section_group}_){section}]` section.
    The options should list the variable args to transform; the parameters they
    point to should list the parameters they are to be transformed to for
    sampling. If a multiple parameters are transformed together, they should
    be comma separated. Example:

    .. code-block:: ini

        [sampling_params]
        mass1, mass2 = mchirp, logitq
        spin1_a = logitspin1_a

    Note that only the final sampling parameters should be listed, even if
    multiple intermediate transforms are needed. (In the above example, a
    transform is needed to go from mass1, mass2 to mchirp, q, then another one
    needed to go from q to logitq.) These transforms should be specified
    in separate sections; see ``transforms.read_transforms_from_config`` for
    details.

    Parameters
    ----------
    cp : WorkflowConfigParser
        An open config parser to read from.
    section_group : str, optional
        Append `{section_group}_` to the section name. Default is None.
    section : str, optional
        The name of the section. Default is 'sampling_params'.

    Returns
    -------
    sampling_params : list
        The list of sampling parameters to use instead.
    replaced_params : list
        The list of variable args to replace in the sampler.
    """
    if section_group is not None:
        section_prefix = '{}_'.format(section_group)
    else:
        section_prefix = ''
    section = section_prefix + section
    replaced_params = set()
    sampling_params = set()
    for args in cp.options(section):
        map_args = cp.get(section, args)
        sampling_params.update(set(map(str.strip, map_args.split(','))))
        replaced_params.update(set(map(str.strip, args.split(','))))
    return sorted(sampling_params), sorted(replaced_params)


#
# =============================================================================
#
#                               Base model definition
#
# =============================================================================
#


class BaseModel(metaclass=ABCMeta):
    r"""Base class for all models.

    Given some model :math:`h` with parameters :math:`\Theta`, Bayes Theorem
    states that the probability of observing parameter values :math:`\vartheta`
    is:

    .. math::

        p(\vartheta|h) = \frac{p(h|\vartheta) p(\vartheta)}{p(h)}.

    Here:

     * :math:`p(\vartheta|h)` is the **posterior** probability;

     * :math:`p(h|\vartheta)` is the **likelihood**;

     * :math:`p(\vartheta)` is the **prior**;

     * :math:`p(h)` is the **evidence**.

    This class defines properties and methods for evaluating the log
    likelihood, log prior, and log posteror. A set of parameter values is set
    using the ``update`` method. Calling the class's
    ``log(likelihood|prior|posterior)`` properties will then evaluate the model
    at those parameter values.

    Classes that inherit from this class must implement a ``_loglikelihood``
    function that can be called by ``loglikelihood``.

    Parameters
    ----------
    variable_params : (tuple of) string(s)
        A tuple of parameter names that will be varied.
    static_params : dict, optional
        A dictionary of parameter names -> values to keep fixed.
    prior : callable, optional
        A callable class or function that computes the log of the prior. If
        None provided, will use ``_noprior``, which returns 0 for all parameter
        values.
    sampling_params : list, optional
        Replace one or more of the ``variable_params`` with the given
        parameters for sampling.
    replace_parameters : list, optional
        The ``variable_params`` to replace with sampling parameters. Must be
        the same length as ``sampling_params``.
    sampling_transforms : list, optional
        List of transforms to use to go between the ``variable_params`` and the
        sampling parameters. Required if ``sampling_params`` is not None.
    waveform_transforms : list, optional
        A list of transforms to convert the ``variable_params`` into something
        understood by the likelihood model. This is useful if the prior is
        more easily parameterized in parameters that are different than what
        the likelihood is most easily defined in. Since these are used solely
        for converting parameters, and not for rescaling the parameter space,
        a Jacobian is not required for these transforms.
    """
    name = None

    def __init__(self, variable_params, static_params=None, prior=None,
                 sampling_transforms=None, waveform_transforms=None, **kwargs):
        # store variable and static args
        self.variable_params = variable_params
        self.static_params = static_params
        # store prior
        if prior is None:
            self.prior_distribution = _NoPrior()
        elif set(prior.variable_args) != set(variable_params):
            raise ValueError("variable params of prior and model must be the "
                             "same")
        else:
            self.prior_distribution = prior
        # store transforms
        self.sampling_transforms = sampling_transforms
        self.waveform_transforms = waveform_transforms
        # initialize current params to None
        self._current_params = None
        # initialize a model stats
        self._current_stats = ModelStats()

    @property
    def variable_params(self):
        """Returns the model parameters."""
        return self._variable_params

    @variable_params.setter
    def variable_params(self, variable_params):
        if isinstance(variable_params, str):
            variable_params = (variable_params,)
        if not isinstance(variable_params, tuple):
            variable_params = tuple(variable_params)
        self._variable_params = variable_params

    @property
    def static_params(self):
        """Returns the model's static arguments."""
        return self._static_params

    @static_params.setter
    def static_params(self, static_params):
        if static_params is None:
            static_params = {}
        self._static_params = static_params

    @property
    def sampling_params(self):
        """Returns the sampling parameters.

        If ``sampling_transforms`` is None, this is the same as the
        ``variable_params``.
        """
        if self.sampling_transforms is None:
            sampling_params = self.variable_params
        else:
            sampling_params = self.sampling_transforms.sampling_params
        return sampling_params

    def update(self, **params):
        """Updates the current parameter positions and resets stats.

        If any sampling transforms are specified, they are applied to the
        params before being stored.
        """
        # add the static params
        values = self.static_params.copy()
        values.update(params)
        self._current_params = self._transform_params(**values)
        self._current_stats = ModelStats()

    @property
    def current_params(self):
        if self._current_params is None:
            raise ValueError("no parameters values currently stored; "
                             "run update to add some")
        return self._current_params

    @property
    def default_stats(self):
        """The stats that ``get_current_stats`` returns by default."""
        return ['logjacobian', 'logprior', 'loglikelihood'] + self._extra_stats

    @property
    def _extra_stats(self):
        """Allows child classes to add more stats to the default stats.

        This returns an empty list; classes that inherit should override this
        property if they want to add extra stats.
        """
        return []

    def get_current_stats(self, names=None):
        """Return one or more of the current stats as a tuple.

        This function does no computation. It only returns what has already
        been calculated. If a stat hasn't been calculated, it will be returned
        as ``numpy.nan``.

        Parameters
        ----------
        names : list of str, optional
            Specify the names of the stats to retrieve. If ``None`` (the
            default), will return ``default_stats``.

        Returns
        -------
        tuple :
            The current values of the requested stats, as a tuple. The order
            of the stats is the same as the names.
        """
        if names is None:
            names = self.default_stats
        return self._current_stats.getstats(names)

    @property
    def current_stats(self):
        """Return the ``default_stats`` as a dict.

        This does no computation. It only returns what has already been
        calculated. If a stat hasn't been calculated, it will be returned
        as ``numpy.nan``.

        Returns
        -------
        dict :
            Dictionary of stat names -> current stat values.
        """
        return self._current_stats.getstatsdict(self.default_stats)

    def _trytoget(self, statname, fallback, apply_transforms=False, **kwargs):
        r"""Helper function to get a stat from ``_current_stats``.

        If the statistic hasn't been calculated, ``_current_stats`` will raise
        an ``AttributeError``. In that case, the ``fallback`` function will
        be called. If that call is successful, the ``statname`` will be added
        to ``_current_stats`` with the returned value.

        Parameters
        ----------
        statname : str
            The stat to get from ``current_stats``.
        fallback : method of self
            The function to call if the property call fails.
        apply_transforms : bool, optional
            Apply waveform transforms to the current parameters before calling
            the fallback function. Default is False.
        \**kwargs :
            Any other keyword arguments are passed through to the function.

        Returns
        -------
        float :
            The value of the property.
        """
        try:
            return getattr(self._current_stats, statname)
        except AttributeError:
            # apply waveform transforms if requested
            if apply_transforms and self.waveform_transforms is not None:
                self._current_params = transforms.apply_transforms(
                    self._current_params, self.waveform_transforms,
                    inverse=False)
            val = fallback(**kwargs)
            setattr(self._current_stats, statname, val)
            return val

    @property
    def loglikelihood(self):
        """The log likelihood at the current parameters.

        This will initially try to return the ``current_stats.loglikelihood``.
        If that raises an ``AttributeError``, will call `_loglikelihood`` to
        calculate it and store it to ``current_stats``.
        """
        return self._trytoget('loglikelihood', self._loglikelihood,
                              apply_transforms=True)

    @abstractmethod
    def _loglikelihood(self):
        """Low-level function that calculates the log likelihood of the current
        params."""
        pass

    @property
    def logjacobian(self):
        """The log jacobian of the sampling transforms at the current postion.

        If no sampling transforms were provided, will just return 0.

        Parameters
        ----------
        \**params :
            The keyword arguments should specify values for all of the variable
            args and all of the sampling args.

        Returns
        -------
        float :
            The value of the jacobian.
        """
        return self._trytoget('logjacobian', self._logjacobian)

    def _logjacobian(self):
        """Calculates the logjacobian of the current parameters."""
        if self.sampling_transforms is None:
            logj = 0.
        else:
            logj = self.sampling_transforms.logjacobian(
                **self.current_params)
        return logj

    @property
    def logprior(self):
        """Returns the log prior at the current parameters."""
        return self._trytoget('logprior', self._logprior)

    def _logprior(self):
        """Calculates the log prior at the current parameters."""
        logj = self.logjacobian
        logp = self.prior_distribution(**self.current_params) + logj
        if numpy.isnan(logp):
            logp = -numpy.inf
        return logp

    @property
    def logposterior(self):
        """Returns the log of the posterior of the current parameter values.

        The logprior is calculated first. If the logprior returns ``-inf``
        (possibly indicating a non-physical point), then the ``loglikelihood``
        is not called.
        """
        logp = self.logprior
        if logp == -numpy.inf:
            return logp
        else:
            return logp + self.loglikelihood

    def prior_rvs(self, size=1, prior=None):
        """Returns random variates drawn from the prior.

        If the ``sampling_params`` are different from the ``variable_params``,
        the variates are transformed to the `sampling_params` parameter space
        before being returned.

        Parameters
        ----------
        size : int, optional
            Number of random values to return for each parameter. Default is 1.
        prior : JointDistribution, optional
            Use the given prior to draw values rather than the saved prior.

        Returns
        -------
        FieldArray
            A field array of the random values.
        """
        # draw values from the prior
        if prior is None:
            prior = self.prior_distribution
        p0 = prior.rvs(size=size)
        # transform if necessary
        if self.sampling_transforms is not None:
            ptrans = self.sampling_transforms.apply(p0)
            # pull out the sampling args
            p0 = FieldArray.from_arrays([ptrans[arg]
                                         for arg in self.sampling_params],
                                        names=self.sampling_params)
        return p0

    def _transform_params(self, **params):
        """Applies sampling transforms and boundary conditions to parameters.

        Parameters
        ----------
        \**params :
            Key, value pairs of parameters to apply the transforms to.

        Returns
        -------
        dict
            A dictionary of the transformed parameters.
        """
        # apply inverse transforms to go from sampling parameters to
        # variable args
        if self.sampling_transforms is not None:
            params = self.sampling_transforms.apply(params, inverse=True)
        # apply boundary conditions
        params = self.prior_distribution.apply_boundary_conditions(**params)
        return params

    #
    # Methods for initiating from a config file.
    #
    @staticmethod
    def extra_args_from_config(cp, section, skip_args=None, dtypes=None):
        """Gets any additional keyword in the given config file.

        Parameters
        ----------
        cp : WorkflowConfigParser
            Config file parser to read.
        section : str
            The name of the section to read.
        skip_args : list of str, optional
            Names of arguments to skip.
        dtypes : dict, optional
            A dictionary of arguments -> data types. If an argument is found
            in the dict, it will be cast to the given datatype. Otherwise, the
            argument's value will just be read from the config file (and thus
            be a string).

        Returns
        -------
        dict
            Dictionary of keyword arguments read from the config file.
        """
        kwargs = {}
        if dtypes is None:
            dtypes = {}
        if skip_args is None:
            skip_args = []
        read_args = [opt for opt in cp.options(section)
                     if opt not in skip_args]
        for opt in read_args:
            val = cp.get(section, opt)
            # try to cast the value if a datatype was specified for this opt
            try:
                val = dtypes[opt](val)
            except KeyError:
                pass
            kwargs[opt] = val
        return kwargs

    @staticmethod
    def prior_from_config(cp, variable_params, static_params, prior_section,
                          constraint_section):
        """Gets arguments and keyword arguments from a config file.

        Parameters
        ----------
        cp : WorkflowConfigParser
            Config file parser to read.
        variable_params : list
            List of variable model parameter names.
        static_params : dict
            Dictionary of static model parameters and their values.
        prior_section : str
            Section to read prior(s) from.
        constraint_section : str
            Section to read constraint(s) from.

        Returns
        -------
        pycbc.distributions.JointDistribution
            The prior.
        """
        # get prior distribution for each variable parameter
        logging.info("Setting up priors for each parameter")
        dists = distributions.read_distributions_from_config(cp, prior_section)
        constraints = distributions.read_constraints_from_config(
            cp, constraint_section, static_args=static_params)
        return distributions.JointDistribution(variable_params, *dists,
                                               constraints=constraints)

    @classmethod
    def _init_args_from_config(cls, cp):
        """Helper function for loading parameters.

        This retrieves the prior, variable parameters, static parameterss,
        constraints, sampling transforms, and waveform transforms
        (if provided).

        Parameters
        ----------
        cp : ConfigParser
            Config parser to read.

        Returns
        -------
        dict :
            Dictionary of the arguments. Has keys ``variable_params``,
            ``static_params``, ``prior``, and ``sampling_transforms``. If
            waveform transforms are in the config file, will also have
            ``waveform_transforms``.
        """
        section = "model"
        prior_section = "prior"
        vparams_section = 'variable_params'
        sparams_section = 'static_params'
        constraint_section = 'constraint'
        # check that the name exists and matches
        name = cp.get(section, 'name')
        if name != cls.name:
            raise ValueError("section's {} name does not match mine {}".format(
                             name, cls.name))
        # get model parameters
        variable_params, static_params = distributions.read_params_from_config(
            cp, prior_section=prior_section, vargs_section=vparams_section,
            sargs_section=sparams_section)
        # get prior
        prior = cls.prior_from_config(
            cp, variable_params, static_params, prior_section,
            constraint_section)
        args = {'variable_params': variable_params,
                'static_params': static_params,
                'prior': prior}
        # try to load sampling transforms
        try:
            sampling_transforms = SamplingTransforms.from_config(
                cp, variable_params)
        except NoSectionError:
            sampling_transforms = None
        args['sampling_transforms'] = sampling_transforms
        # get any waveform transforms
        if any(cp.get_subsections('waveform_transforms')):
            logging.info("Loading waveform transforms")
            waveform_transforms = transforms.read_transforms_from_config(
                cp, 'waveform_transforms')
            args['waveform_transforms'] = waveform_transforms
        else:
            waveform_transforms = []
        # safety check for spins
        # we won't do this if the following exists in the config file
        ignore = "no_err_on_missing_cartesian_spins"
        check_for_cartesian_spins(1, variable_params, static_params,
                                  waveform_transforms, cp, ignore)
        check_for_cartesian_spins(2, variable_params, static_params,
                                  waveform_transforms, cp, ignore)
        return args

    @classmethod
    def from_config(cls, cp, **kwargs):
        """Initializes an instance of this class from the given config file.

        Parameters
        ----------
        cp : WorkflowConfigParser
            Config file parser to read.
        \**kwargs :
            All additional keyword arguments are passed to the class. Any
            provided keyword will over ride what is in the config file.
        """
        args = cls._init_args_from_config(cp)
        # get any other keyword arguments provided in the model section
        args.update(cls.extra_args_from_config(cp, "model",
                                               skip_args=['name']))
        args.update(kwargs)
        return cls(**args)

    def write_metadata(self, fp, group=None):
        """Writes metadata to the given file handler.

        Parameters
        ----------
        fp : pycbc.inference.io.BaseInferenceFile instance
            The inference file to write to.
        group : str, optional
            If provided, the metadata will be written to the attrs specified
            by group, i.e., to ``fp[group].attrs``. Otherwise, metadata is
            written to the top-level attrs (``fp.attrs``).
        """
        attrs = fp.getattrs(group=group)
        attrs['model'] = self.name
        attrs['variable_params'] = list(map(str, self.variable_params))
        attrs['sampling_params'] = list(map(str, self.sampling_params))
        fp.write_kwargs_to_attrs(attrs, static_params=self.static_params)


def check_for_cartesian_spins(which, variable_params, static_params,
                              waveform_transforms, cp, ignore):
    """Checks that if any spin parameters exist, cartesian spins also exist.

    This looks for parameters starting with ``spinN`` in the variable and
    static params, where ``N`` is either  1 or 2 (specified by the ``which``
    argument). If any parameters are found with those names, the params and
    the output of the waveform transforms are checked to see that there is
    at least one of ``spinN(x|y|z)``. If not, a ``ValueError`` is raised.

    This check will not be done if the config file has an section given by
    the ignore argument.

    Parameters
    ----------
    which : {1, 2}
        Which component to check for. Must be either 1 or 2.
    variable_params : list
        List of the variable parameters.
    static_params : dict
        The dictionary of static params.
    waveform_transforms : list
        List of the transforms that will be applied to the variable and
        static params before being passed to the waveform generator.
    cp : ConfigParser
        The config file.
    ignore : str
        The section to check for in the config file. If the section is
        present in the config file, the check will not be done.
    """
    # don't do this check if the config file has the ignore section
    if cp.has_section(ignore):
        logging.info("[{}] found in config file; not performing check for "
                     "cartesian spin{} parameters".format(ignore, which))
        return
    errmsg = (
        "Spin parameters {sp} found in variable/static "
        "params for component {n}, but no Cartesian spin parameters ({cp}) "
        "found in either the variable/static params or "
        "the waveform transform outputs. Most waveform "
        "generators only recognize Cartesian spin "
        "parameters; without them, all spins are set to "
        "zero. If you are using spherical spin coordinates, add "
        "the following waveform_transform to your config file:\n\n"
        "[waveform_transforms-spin{n}x+spin{n}y+spin{n}z]\n"
        "name = spherical_to_cartesian\n"
        "x = spin{n}x\n"
        "y = spin{n}y\n"
        "z = spin{n}z\n"
        "radial = spin{n}_a\n"
        "azimuthal = spin{n}_azimuthal\n"
        "polar = spin{n}_polar\n\n"
        "Here, spin{n}_a, spin{n}_azimuthal, and spin{n}_polar are the names "
        "of your radial, azimuthal, and polar coordinates, respectively. "
        "If you intentionally did not include Cartesian spin parameters, "
        "(e.g., you are using a custom waveform or model) add\n\n"
        "[{ignore}]\n\n"
        "to your config file as an empty section and rerun. This check will "
        "not be performed in that case.")
    allparams = set(variable_params) | set(static_params.keys())
    spinparams = set(p for p in allparams
                     if p.startswith('spin{}'.format(which)))
    if any(spinparams):
        cartspins = set('spin{}{}'.format(which, coord)
                        for coord in ['x', 'y', 'z'])
        # add any parameters to all params that will be output by waveform
        # transforms
        allparams = allparams.union(*[t.outputs for t in waveform_transforms])
        if not any(allparams & cartspins):
            raise ValueError(errmsg.format(sp=', '.join(spinparams),
                                           cp=', '.join(cartspins),
                                           n=which, ignore=ignore))
