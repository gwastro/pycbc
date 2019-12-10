# Copyright (C) 2016 Collin Capano, Duncan Brown
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Generals
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

"""This module contains standard options used for inference-related programs.
"""

import logging
import argparse

from six import string_types

from pycbc.psd import from_cli_multi_ifos as psd_from_cli_multi_ifos
from pycbc.strain import from_cli_multi_ifos as strain_from_cli_multi_ifos
from pycbc.strain import (gates_from_cli, psd_gates_from_cli,
                          apply_gates_to_td, apply_gates_to_fd)
from pycbc.types import MultiDetOptionAction
from pycbc import waveform
from pycbc import distributions


# -----------------------------------------------------------------------------
#
#                       Utilities for loading data
#
# -----------------------------------------------------------------------------


def add_low_frequency_cutoff_opt(parser):
    """Adds the low-frequency-cutoff option to the given parser."""
    parser.add_argument("--data-conditioning-low-freq", type=float,
                        nargs="+", action=MultiDetOptionAction,
                        metavar='IFO:FLOW', dest="low_frequency_cutoff",
                        help="Low frequency cutoff of the data. Needed for "
                             "PSD estimation and when creating fake strain.")


def data_from_cli(opts):
    """Loads the data needed for a model from the given
    command-line options. Gates specifed on the command line are also applied.

    Parameters
    ----------
    opts : ArgumentParser parsed args
        Argument options parsed from a command line string (the sort of thing
        returned by `parser.parse_args`).

    Returns
    -------
    strain_dict : dict
        Dictionary of instruments -> `TimeSeries` strain.
    stilde_dict : dict
        Dictionary of instruments -> `FrequencySeries` strain.
    psd_dict : dict
        Dictionary of instruments -> `FrequencySeries` psds.
    """
    # get gates to apply
    gates = gates_from_cli(opts)
    psd_gates = psd_gates_from_cli(opts)

    # get strain time series
    instruments = opts.instruments if opts.instruments is not None else []
    strain_dict = strain_from_cli_multi_ifos(opts, instruments,
                                             precision="double")
    # apply gates if not waiting to overwhiten
    if not opts.gate_overwhitened:
        strain_dict = apply_gates_to_td(strain_dict, gates)

    # get strain time series to use for PSD estimation
    # if user has not given the PSD time options then use same data as analysis
    if opts.psd_start_time and opts.psd_end_time:
        logging.info("Will generate a different time series for PSD "
                     "estimation")
        psd_opts = opts
        psd_opts.gps_start_time = psd_opts.psd_start_time
        psd_opts.gps_end_time = psd_opts.psd_end_time
        psd_strain_dict = strain_from_cli_multi_ifos(psd_opts,
                                                     instruments,
                                                     precision="double")
        # apply any gates
        logging.info("Applying gates to PSD data")
        psd_strain_dict = apply_gates_to_td(psd_strain_dict, psd_gates)

    elif opts.psd_start_time or opts.psd_end_time:
        raise ValueError("Must give --psd-start-time and --psd-end-time")
    else:
        psd_strain_dict = strain_dict

    # FFT strain and save each of the length of the FFT, delta_f, and
    # low frequency cutoff to a dict
    stilde_dict = {}
    length_dict = {}
    delta_f_dict = {}
    for ifo in instruments:
        stilde_dict[ifo] = strain_dict[ifo].to_frequencyseries()
        length_dict[ifo] = len(stilde_dict[ifo])
        delta_f_dict[ifo] = stilde_dict[ifo].delta_f

    # get PSD as frequency series
    psd_dict = psd_from_cli_multi_ifos(
        opts, length_dict, delta_f_dict, opts.low_frequency_cutoff,
        instruments, strain_dict=psd_strain_dict, precision="double")

    # apply any gates to overwhitened data, if desired
    if opts.gate_overwhitened and opts.gate is not None:
        logging.info("Applying gates to overwhitened data")
        # overwhiten the data
        for ifo in gates:
            stilde_dict[ifo] /= psd_dict[ifo]
        stilde_dict = apply_gates_to_fd(stilde_dict, gates)
        # unwhiten the data for the model
        for ifo in gates:
            stilde_dict[ifo] *= psd_dict[ifo]

    return strain_dict, stilde_dict, psd_dict


# -----------------------------------------------------------------------------
#
#                Utilities for plotting results
#
# -----------------------------------------------------------------------------


class ParseLabelArg(argparse.Action):
    """Argparse action that will parse arguments that can accept labels.

    This assumes that the values set on the command line for its assigned
    argument are strings formatted like ``PARAM[:LABEL]``. When the arguments
    are parsed, the ``LABEL`` bit is stripped off and added to a dictionary
    mapping ``PARAM -> LABEL``. This dictionary is stored to the parsed
    namespace called ``{dest}_labels``, where ``{dest}`` is the argument's
    ``dest`` setting (by default, this is the same as the option string).
    Likewise, the argument's ``dest`` in the parsed namespace is updated so
    that it is just ``PARAM``.

    If no ``LABEL`` is provided, then ``PARAM`` will be used for ``LABEL``.

    This action can work on arguments that have ``nargs != 0`` and ``type`` set
    to ``str``.
    """
    def __init__(self, type=str, nargs=None,
                 **kwargs):  # pylint: disable=redefined-builtin
        # check that type is string
        if type != str:
            raise ValueError("the type for this action must be a string")
        if nargs == 0:
            raise ValueError("nargs must not be 0 for this action")
        super(ParseLabelArg, self).__init__(type=type, nargs=nargs,
                                            **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        singlearg = isinstance(values, string_types)
        if singlearg:
            values = [values]
        params = []
        labels = {}
        for param in values:
            psplit = param.split(':')
            if len(psplit) == 2:
                param, label = psplit
            else:
                label = param
            labels[param] = label
            params.append(param)
        # update the namespace
        if singlearg:
            params = params[0]
        setattr(namespace, self.dest, params)
        setattr(namespace, '{}_labels'.format(self.dest), labels)


class ParseParametersArg(ParseLabelArg):
    """Argparse action that will parse parameters and labels from an opton.

    Does the same as ``ParseLabelArg``, with the additional functionality that
    if ``LABEL`` is a known parameter in ``pycbc.waveform.parameters``, then
    the label attribute there will be used in the labels dictionary.
    Otherwise, ``LABEL`` will be used.

    Examples
    --------
    Create a parser and add two arguments that use this action (note that the
    first argument accepts multiple inputs while the second only accepts a
    single input):

    >>> import argparse
    >>> parser = argparse.ArgumentParser()
    >>> parser.add_argument('--parameters', type=str, nargs="+",
                            action=ParseParametersArg)
    >>> parser.add_argument('--z-arg', type=str, action=ParseParametersArg)

    Parse a command line that uses these options:

    >>> import shlex
    >>> cli = "--parameters 'mass1+mass2:mtotal' ra ni --z-arg foo:bar"
    >>> opts = parser.parse_args(shlex.split(cli))
    >>> opts.parameters
    ['mass1+mass2', 'ra', 'ni']
    >>> opts.parameters_labels
    {'mass1+mass2': '$M~(\\mathrm{M}_\\odot)$', 'ni': 'ni', 'ra': '$\\alpha$'}
    >>> opts.z_arg
    'foo'
    >>> opts.z_arg_labels
    {'foo': 'bar'}

    In the above, the first argument to ``--parameters`` was ``mtotal``. Since
    this is a recognized parameter in ``pycbc.waveform.parameters``, the label
    dictionary contains the latex string associated with the ``mtotal``
    parameter. A label was not provided for the second argument, and so ``ra``
    was used. Since ``ra`` is also a recognized parameter, its associated latex
    string was used in the labels dictionary. Since ``ni`` and ``bar`` (the
    label for ``z-arg``) are not recognized parameters, they were just used
    as-is in the labels dictionaries.
    """
    def __call__(self, parser, namespace, values, option_string=None):
        super(ParseParametersArg, self).__call__(parser, namespace, values,
                                                 option_string=option_string)
        # try to replace the labels with a label from waveform.parameters
        labels = getattr(namespace, '{}_labels'.format(self.dest))
        for param, label in labels.items():
            try:
                label = getattr(waveform.parameters, label).label
                labels[param] = label
            except AttributeError:
                pass


def add_plot_posterior_option_group(parser):
    """Adds the options needed to configure plots of posterior results.

    Parameters
    ----------
    parser : object
        ArgumentParser instance.
    """
    pgroup = parser.add_argument_group("Options for what plots to create and "
                                       "their formats.")
    pgroup.add_argument('--plot-marginal', action='store_true', default=False,
                        help="Plot 1D marginalized distributions on the "
                             "diagonal axes.")
    pgroup.add_argument('--marginal-percentiles', nargs='+', default=None,
                        type=float,
                        help="Percentiles to draw lines at on the 1D "
                             "histograms.")
    pgroup.add_argument("--plot-scatter", action='store_true', default=False,
                        help="Plot each sample point as a scatter plot.")
    pgroup.add_argument("--plot-density", action="store_true", default=False,
                        help="Plot the posterior density as a color map.")
    pgroup.add_argument("--plot-contours", action="store_true", default=False,
                        help="Draw contours showing the 50th and 90th "
                             "percentile confidence regions.")
    pgroup.add_argument('--contour-percentiles', nargs='+', default=None,
                        type=float,
                        help="Percentiles to draw contours if different "
                             "than 50th and 90th.")
    # add mins, maxs options
    pgroup.add_argument('--mins', nargs='+', metavar='PARAM:VAL', default=[],
                        help="Specify minimum parameter values to plot. This "
                             "should be done by specifying the parameter name "
                             "followed by the value. Parameter names must be "
                             "the same as the PARAM argument in --parameters "
                             "(or, if no parameters are provided, the same as "
                             "the parameter name specified in the variable "
                             "args in the input file. If none provided, "
                             "the smallest parameter value in the posterior "
                             "will be used.")
    pgroup.add_argument('--maxs', nargs='+', metavar='PARAM:VAL', default=[],
                        help="Same as mins, but for the maximum values to "
                             "plot.")
    # add expected parameters options
    pgroup.add_argument('--expected-parameters', nargs='+',
                        metavar='PARAM:VAL',
                        default=[],
                        help="Specify expected parameter values to plot. If "
                             "provided, a cross will be plotted in each axis "
                             "that an expected parameter is provided. "
                             "Parameter names must be "
                             "the same as the PARAM argument in --parameters "
                             "(or, if no parameters are provided, the same as "
                             "the parameter name specified in the variable "
                             "args in the input file.")
    pgroup.add_argument('--expected-parameters-color', default='r',
                        help="What to color the expected-parameters cross. "
                             "Default is red.")
    pgroup.add_argument('--plot-injection-parameters', action='store_true',
                        default=False,
                        help="Get the expected parameters from the injection "
                             "in the input file. There must be only a single "
                             "injection in the file to work. Any values "
                             "specified by expected-parameters will override "
                             "the values obtained for the injection.")
    return pgroup


def plot_ranges_from_cli(opts):
    """Parses the mins and maxs arguments from the `plot_posterior` option
    group.

    Parameters
    ----------
    opts : ArgumentParser
        The parsed arguments from the command line.

    Returns
    -------
    mins : dict
        Dictionary of parameter name -> specified mins. Only parameters that
        were specified in the --mins option will be included; if no parameters
        were provided, will return an empty dictionary.
    maxs : dict
        Dictionary of parameter name -> specified maxs. Only parameters that
        were specified in the --mins option will be included; if no parameters
        were provided, will return an empty dictionary.
    """
    mins = {}
    for x in opts.mins:
        x = x.split(':')
        if len(x) != 2:
            raise ValueError("option --mins not specified correctly; see help")
        mins[x[0]] = float(x[1])
    maxs = {}
    for x in opts.maxs:
        x = x.split(':')
        if len(x) != 2:
            raise ValueError("option --maxs not specified correctly; see help")
        maxs[x[0]] = float(x[1])
    return mins, maxs


def expected_parameters_from_cli(opts):
    """Parses the --expected-parameters arguments from the `plot_posterior`
    option group.

    Parameters
    ----------
    opts : ArgumentParser
        The parsed arguments from the command line.

    Returns
    -------
    dict
        Dictionary of parameter name -> expected value. Only parameters that
        were specified in the --expected-parameters option will be included; if
        no parameters were provided, will return an empty dictionary.
    """
    expected = {}
    for x in opts.expected_parameters:
        x = x.split(':')
        if len(x) != 2:
            raise ValueError("option --expected-paramters not specified "
                             "correctly; see help")
        expected[x[0]] = float(x[1])
    return expected


def add_scatter_option_group(parser):
    """Adds the options needed to configure scatter plots.

    Parameters
    ----------
    parser : object
        ArgumentParser instance.
    """
    scatter_group = parser.add_argument_group("Options for configuring the "
                                              "scatter plot.")

    scatter_group.add_argument(
        '--z-arg', type=str, default=None, action=ParseParametersArg,
        help='What to color the scatter points by. Syntax is the same as the '
             'parameters option.')
    scatter_group.add_argument(
        "--vmin", type=float, help="Minimum value for the colorbar.")
    scatter_group.add_argument(
        "--vmax", type=float, help="Maximum value for the colorbar.")
    scatter_group.add_argument(
        "--scatter-cmap", type=str, default='plasma',
        help="Specify the colormap to use for points. Default is plasma.")

    return scatter_group


def add_density_option_group(parser):
    """Adds the options needed to configure contours and density colour map.

    Parameters
    ----------
    parser : object
        ArgumentParser instance.
    """
    density_group = parser.add_argument_group("Options for configuring the "
                                              "contours and density color map")

    density_group.add_argument(
        "--density-cmap", type=str, default='viridis',
        help="Specify the colormap to use for the density. "
             "Default is viridis.")
    density_group.add_argument(
        "--contour-color", type=str, default=None,
        help="Specify the color to use for the contour lines. Default is "
             "white for density plots and black for scatter plots.")
    density_group.add_argument(
        '--use-kombine-kde', default=False, action="store_true",
        help="Use kombine's KDE for determining contours. "
             "Default is to use scipy's gaussian_kde.")

    return density_group


def prior_from_config(cp, prior_section='prior'):
    """Loads a prior distribution from the given config file.

    Parameters
    ----------
    cp : pycbc.workflow.WorkflowConfigParser
        The config file to read.
    sections : list of str, optional
        The sections to retrieve the prior from. If ``None`` (the default),
        will look in sections starting with 'prior'.

    Returns
    -------
    distributions.JointDistribution
        The prior distribution.
    """
    # Read variable and static parameters from the config file
    variable_params, _ = distributions.read_params_from_config(
        cp, prior_section=prior_section, vargs_section='variable_params',
        sargs_section='static_params')
    # Read constraints to apply to priors from the config file
    constraints = distributions.read_constraints_from_config(cp)
    # Get PyCBC distribution instances for each variable parameter in the
    # config file
    dists = distributions.read_distributions_from_config(cp, prior_section)
    # construct class that will return draws from the prior
    return distributions.JointDistribution(variable_params, *dists,
                                           **{"constraints": constraints})
