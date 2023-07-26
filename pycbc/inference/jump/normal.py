# Copyright (C) 2019  Collin Capano
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

"""Jump proposals that use a normal distribution."""


import numpy

from epsie import proposals as epsie_proposals
from epsie.proposals import Boundaries

from pycbc import VARARGS_DELIM


class EpsieNormal(epsie_proposals.Normal):
    """Adds ``from_config`` method to epsie's normal proposal."""

    @classmethod
    def from_config(cls, cp, section, tag):
        """Loads a proposal from a config file.

        This calls :py:func:`epsie_from_config` with ``cls`` set to
        :py:class:`epsie.proposals.Normal` and ``with_boundaries`` set to
        False. See that function for details on options that can be read.

        Example::

            [jump_proposal-mchrip+q]
            name = normal
            var-q = 0.1

        Parameters
        ----------
        cp : WorkflowConfigParser instance
            Config file to read from.
        section : str
            The name of the section to look in.
        tag : str
            :py:const:`pycbc.VARARGS_DELIM` separated list of parameter names
            to create proposals for.

        Returns
        -------
        :py:class:`epsie.proposals.Normal`:
            A normal proposal for use with ``epsie`` samplers.
        """
        return epsie_from_config(cls, cp, section, tag, with_boundaries=False)


class EpsieAdaptiveNormal(epsie_proposals.AdaptiveNormal):
    """Adds ``from_config`` method to epsie's adaptive normal proposal."""

    @classmethod
    def from_config(cls, cp, section, tag):
        """Loads a proposal from a config file.

        This calls :py:func:`epsie_adaptive_from_config` with ``cls`` set to
        :py:class:`epsie.proposals.AdaptiveNormal`. See that function for
        details on options that can be read.

        Example::

            [jump_proposal-mchirp+q]
            name = adaptive_normal
            adaptation-duration = 1000
            min-q = 1
            max-q = 8
            min-mchirp = 20
            max-mchirp = 80

        Parameters
        ----------
        cp : WorkflowConfigParser instance
            Config file to read from.
        section : str
            The name of the section to look in.
        tag : str
            :py:const:`pycbc.VARARGS_DELIM` separated list of parameter names
            to create proposals for.

        Returns
        -------
        :py:class:`epsie.proposals.AdaptiveNormal`:
            An adaptive normal proposal for use with ``epsie`` samplers.
        """
        return epsie_adaptive_from_config(cls, cp, section, tag,
                                          boundary_arg_name='prior_widths')


class EpsieATAdaptiveNormal(epsie_proposals.ATAdaptiveNormal):
    """Adds ``from_config`` method to epsie's ATAdaptiveProposal."""

    @classmethod
    def from_config(cls, cp, section, tag):
        """Loads a proposal from a config file.

        This calls :py:func:`epsie_from_config` with ``cls`` set to
        :py:class:`epsie.proposals.AdaptiveProposal` and ``with_boundaries``
        set to False. See that function for details on options that can be
        read.

        Example::

            [jump_proposal-mchrip+q]
            name = adaptive_proposal
            diagonal =

        Parameters
        ----------
        cp : WorkflowConfigParser instance
            Config file to read from.
        section : str
            The name of the section to look in.
        tag : str
            :py:const:`pycbc.VARARGS_DELIM` separated list of parameter names
            to create proposals for.

        Returns
        -------
        :py:class:`epsie.proposals.AdaptiveProposal`:
            An adaptive proposal for use with ``epsie`` samplers.
        """
        return epsie_at_adaptive_from_config(cls, cp, section, tag,
                                             with_boundaries=False)


def epsie_from_config(cls, cp, section, tag, with_boundaries=False):
    r"""Generic function for loading epsie proposals from a config file.

    This should be used for proposals that are not adaptive.

    The section that is read should have the format ``[{section}-{tag}]``,
    where ``{tag}`` is a :py:const:`pycbc.VARARGS_DELIM` separated list
    of the parameters to create the jump proposal for.

    Options that are read:

    * name : str
        Required. Must match the name of the proposal.
    * var-{param} : float
        Optional. Variance to use for parameter {param}. If ``with_boundaries``
        is True, then any parameter not specified will use a default variance
        of :math:`(\Delta p/10)^2`, where :math:`\Delta p` is the boundary
        width for that parameter. If ``with_boundaries`` is False, will use
        a default value of 1.
    * min-{param} : float
    * max-{param} : float
        The bounds on each parameter. Required if ``with_boundaries`` is set to
        True, in which case bounds must be provided for every parameter.

    Parameters
    ----------
    cls : epsie.proposals.BaseProposal
        The epsie proposal class to initialize.
    cp : WorkflowConfigParser instance
        Config file to read from.
    section : str
        The name of the section to look in.
    tag : str
        :py:const:`pycbc.VARARGS_DELIM` separated list of parameter names
        to create proposals for.
    with_boundaries : bool, optional
        Try to load boundaries from the section and pass a ``boundaries``
        argument to the class's initialization. This should be set to true
        for bounded proposals. Default is False.

    Returns
    -------
    cls :
        The class initialized with the options read from the config file.
    """
    # check that the name matches
    assert cp.get_opt_tag(section, "name", tag) == cls.name, (
        "name in specified section must match mine")
    params, opts = load_opts(cp, section, tag, skip=['name'])
    args = {'parameters': params}
    if with_boundaries:
        boundaries = get_param_boundaries(params, opts)
        args['boundaries'] = boundaries
    if 'discrete' in cls.name.split('_'):
        args.update({'successive':
                     get_epsie_discrete_successive_settings(params, opts)})
    # if there are any options left, assume they are for setting the variance
    if opts:
        cov = get_variance(params, opts)
    elif with_boundaries:
        cov = numpy.array([abs(boundaries[p])/10. for p in params])**2.
    else:
        cov = None
    args['cov'] = cov
    # no other options should remain
    if opts:
        raise ValueError("unrecognized options {}"
                         .format(', '.join(opts.keys())))
    return cls(**args)


def epsie_adaptive_from_config(cls, cp, section, tag, with_boundaries=True,
                               boundary_arg_name='boundaries'):
    """Generic function for loading adaptive epsie proposals from a config
    file.

    The section that is read should have the format ``[{section}-{tag}]``,
    where ``{tag}`` is a :py:const:`pycbc.VARARGS_DELIM` separated list
    of the parameters to create the jump proposal for.

    Options that are read:

    * name : str
        Required. Must match the name of the proposal.
    * adaptation-duration : int
        Required. Sets the ``adaptation_duration``.
    * min-{param} : float
    * max-{param} : float
        The bounds on each parameter. Required if ``with_boundaries`` is set to
        True, in which case bounds must be provided for every parameter.
    * var-{param} : float
        Optional. Initial variance to use. If not provided, will use a
        default based on the bounds (see
        :py:class:`epsie.proposals.AdaptiveSupport` for details).
    * adaptation-decay : int
        Optional. Sets the ``adaptation_decay``. If not provided, will use
        the class's default.
    * start-iteration : int
        Optional. Sets the ``start_iteration``.If not provided, will use
        the class's default.
    * target-rate : float
        Optional. Sets the ``target_rate``. If not provided, will use
        the class's default.

    Parameters
    ----------
    cls : epsie.proposals.BaseProposal
        The epsie proposal class to initialize. The class should have
        :py:class:`epsie.proposals.normal.AdaptiveSupport`.
    cp : WorkflowConfigParser instance
        Config file to read from.
    section : str
        The name of the section to look in.
    tag : str
        :py:const:`pycbc.VARARGS_DELIM` separated list of parameter names
        to create proposals for.
    with_boundaries : bool, optional
        Try to load boundaries from the section and pass a ``boundaries``
        argument to the class's initialization. Default is True.
    boundary_arg_name : str, optional
        The name of the argument for the boundaries (only used if
        ``with_boundaries`` is True). Provided because some adaptive proposals
        that only need the boundary widths call this ``prior_widths``. Default
        is ``'boundaries'``.

    Returns
    -------
    cls :
        The class initialized with the options read from the config file.
    """
    # check that the name matches
    assert cp.get_opt_tag(section, "name", tag) == cls.name, (
        "name in specified section must match mine")
    params, opts = load_opts(cp, section, tag, skip=['name'])
    args = {'parameters': params}
    # get the bounds
    if with_boundaries:
        args[boundary_arg_name] = get_param_boundaries(params, opts)
    if 'discrete' in cls.name.split('_'):
        args.update({'successive':
                     get_epsie_discrete_successive_settings(params, opts)})
    # get the adaptation parameters
    args.update(get_epsie_adaptation_settings(opts))
    # if there are any other options, assume they are for setting the
    # initial standard deviation
    if opts:
        var = get_variance(params, opts)
        args['initial_std'] = var**0.5
        # at this point, there should be no options left
        if opts:
            raise ValueError('unrecognized options {} in section {}'
                             .format(', '.join(opts.keys()),
                                     '-'.join([section, tag])))
    return cls(**args)


def epsie_at_adaptive_from_config(cls, cp, section, tag,
                                  with_boundaries=False):
    """Generic function for loading AT Adaptive Normal proposals from a config
    file.

    The section that is read should have the format ``[{section}-{tag}]``,
    where ``{tag}`` is a :py:const:`pycbc.VARARGS_DELIM` separated list
    of the parameters to create the jump proposal for.

    Options that are read:

    * name : str
        Required. Must match the name of the proposal.
    * adaptation-duration : int
        Sets the ``adaptation_duration``. If not provided will use the class's
        default.
    * diagonal : bool, optional
        Determines whether only to adapt the variance. If True will only train
        the diagonal elements.
    * componentwise : bool, optional
        Whether to include a componentwise scaling of the parameters.
        By default set to False. Componentwise scaling `ndim` times more
        expensive than global scaling.
    * min-{param} : float
    * max-{param} : float
        The bounds on each parameter. Required if ``with_boundaries`` is set to
        True, in which case bounds must be provided for every parameter.
    * start-iteration : int
        Optional. Sets the ``start_iteration``. If not provided, will use
        the class's default.
    * target-rate : float
        Optional. Sets the ``target_rate``. If not provided, will use
        the class's default.

    Parameters
    ----------
    cls : epsie.proposals.BaseProposal
        The epsie proposal class to initialize. The class should have
        :py:class:`epsie.proposals.normal.AdaptiveSupport`.
    cp : WorkflowConfigParser instance
        Config file to read from.
    section : str
        The name of the section to look in.
    tag : str
        :py:const:`pycbc.VARARGS_DELIM` separated list of parameter names
        to create proposals for.
    with_boundaries : bool, optional
        Try to load boundaries from the section and pass a ``boundaries``
        argument to the class's initialization. Default is True.

    Returns
    -------
    cls :
        The class initialized with the options read from the config file.
    """
    # check that the name matches
    assert cp.get_opt_tag(section, "name", tag) == cls.name, (
        "name in specified section must match mine")
    params, opts = load_opts(cp, section, tag, skip=['name'])
    args = {'parameters': params}
    # get the bounds
    if with_boundaries:
        args['boundaries'] = get_param_boundaries(params, opts)
    if 'discrete' in cls.name.split('_'):
        args.update({'successive':
                     get_epsie_discrete_successive_settings(params, opts)})
    # get the adaptation parameters
    args.update(get_epsie_adaptation_settings(opts, cls.name))
    # bounded and angular adaptive proposals support diagonal-only
    diagonal = opts.pop('diagonal', None)
    if not any(p in cls.name.split('_') for p in ['bounded', 'angular']):
        args.update({'diagonal': diagonal is not None})
    componentwise = opts.pop('componentwise', None)
    if componentwise is not None:
        args.update({'componentwise': True})
    if opts:
        raise ValueError("unrecognized options {}"
                         .format(', '.join(opts.keys())))
    return cls(**args)


def load_opts(cp, section, tag, skip=None):
    """Loads config options for jump proposals.

    All `-` in option names are converted to `_` before returning.

    Parameters
    ----------
    cp : WorkflowConfigParser instance
        Config file to read from.
    section : str
        The name of the section to look in.
    tag : str
        :py:const:`pycbc.VARARGS_DELIM` separated list of parameter names
        to create proposals for.
    skip : list, optional
        List of option names to skip loading.

    Returns
    -------
    params : list
        List of parameter names the jump proposal is for.
    opts : dict
        Dictionary of option names -> values, where all values are strings.
    """
    if skip is None:
        skip = []
    params = tag.split(VARARGS_DELIM)
    # get options
    readsection = '-'.join([section, tag])
    opts = {opt.replace('-', '_'): cp.get(readsection, opt)
            for opt in cp.options(readsection) if opt not in skip}
    return params, opts


def get_variance(params, opts, default=1.):
    """Gets variance for jump proposals from the dictionary of options.

    This looks for ``var_{param}`` for every parameter listed in ``params``.
    If found, the argument is popped from the given ``opts`` dictionary. If not
    found, ``default`` will be used.

    Parameters
    ----------
    params : list of str
        List of parameter names to look for.
    opts : dict
        Dictionary of option -> value that was loaded from a config file
        section.
    default : float, optional
        Default value to use for parameters that do not have variances
        provided. Default is 1.

    Returns
    -------
    numpy.array
        Array of variances to use. Order is the same as the parameter names
        given in ``params``.
    """
    varfmt = 'var_{}'
    cov = numpy.array([float(opts.pop(varfmt.format(param), default))
                       for param in params])
    return cov


def get_param_boundaries(params, opts):
    """Gets parameter boundaries for jump proposals.

    The syntax for the options should be ``(min|max)_{param} = value``. Both
    a minimum and maximum should be provided for every parameter in ``params``.
    If the opts are created using ``load_opts``, then the options can be
    formatted as ``(min|max)-{param}``, since that function will turn all ``-``
    to ``_`` in option names.

    Arguments will be popped from the given ``opts`` dictionary.

    Parameters
    ----------
    params : list of str
        List of parameter names to get boundaries for.
    opts : dict
        Dictionary of option -> value that was loaded from a config file
        section.

    Returns
    -------
    dict :
        Dictionary of parameter names -> :py:class:`epsie.proposals.Boundaries`
    """
    boundaries = {}
    for param in params:
        minbound = opts.pop('min_{}'.format(param), None)
        if minbound is None:
            raise ValueError("Must provide a minimum bound for {p}."
                             "Syntax is min_{p} = val".format(p=param))
        maxbound = opts.pop('max_{}'.format(param), None)
        if maxbound is None:
            raise ValueError("Must provide a maximum bound for {p}."
                             "Syntax is max_{p} = val".format(p=param))
        boundaries[param] = Boundaries((float(minbound), float(maxbound)))
    return boundaries


def get_epsie_adaptation_settings(opts, name=None):
    """Get settings for Epsie adaptive proposals from a config file.

    This requires that ``adaptation_duration`` is in the given dictionary.
    It will also look for ``adaptation_decay``, ``start_iteration``, and
    ``target_rate``, but these are optional. Arguments will be popped from the
    given dictionary.

    Parameters
    ----------
    opts : dict
        Dictionary of option -> value that was loaded from a config file
        section.
    name : str (optional)
        Proposal name

    Returns
    -------
    dict :
        Dictionary of argument name -> values.
    """
    args = {}
    adaptation_duration = opts.pop('adaptation_duration', None)
    if adaptation_duration is None:
        if name is not None:
            if all(p in name.split('_') for p in ['at', 'adaptive']):
                args.update({'adaptation_duration': None})
        else:
            raise ValueError("No adaptation_duration specified")
    else:
        args.update({'adaptation_duration': int(adaptation_duration)})
    # optional args
    adaptation_decay = opts.pop('adaptation_decay', None)
    if adaptation_decay is not None:
        args.update({'adaptation_decay': int(adaptation_decay)})
    start_iteration = opts.pop('start_iteration', None)
    if start_iteration is not None:
        args.update({'start_iteration': int(start_iteration)})
    target_rate = opts.pop('target_rate', None)
    if target_rate is not None:
        args.update({'target_rate': float(target_rate)})
    return args


def get_epsie_discrete_successive_settings(params, opts):
    """Get settings for Epsie successive discrete proposal successive jumps
    from a config file.

    If ``successive`` is not defined for a parameter then assumes successive
    jumps are not allowed (i.e. jumps from an integer to the same integer).
    Arguments will be popped from the given dictionary.

    Example::
        [jump_proposal-k+n]
        name = discrete
        successive-k =

    This example sets successive jumps for ``k`` but does not do so for ``n``.

    Parameters
    ----------
    params : list of str
        List of parameter names to get the successive option for.
    opts : dict
        Dictionary of option -> value that was loaded from a config file
        section.

    Returns
    -------
    dict :
        Dictionary of parameter names -> bools
    """
    successive = {}
    for param in params:
        successive.update(
            {param: opts.pop('successive_{}'.format(param), None) is not None})
    return successive
