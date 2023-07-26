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

import argparse
from pycbc import waveform

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
        singlearg = isinstance(values, str)
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


def add_injsamples_map_opt(parser):
    """Adds option to parser to specify a mapping between injection parameters
    an sample parameters.
    """
    parser.add_argument('--injection-samples-map', nargs='+',
                        metavar='INJECTION_PARAM:SAMPLES_PARAM',
                        help='Rename/apply functions to the injection '
                             'parameters and name them the same as one of the '
                             'parameters in samples. This can be used if the '
                             'injection parameters are not the same as the '
                             'samples parameters. INJECTION_PARAM may be a '
                             'function of the injection parameters; '
                             'SAMPLES_PARAM must a name of one of the '
                             'parameters in the samples group.')


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
    pgroup.add_argument('--no-marginal-lines', action='store_true',
                        default=False,
                        help="Do not add vertical lines in the 1D marginal "
                             "plots showing the marginal percentiles.")
    pgroup.add_argument('--no-marginal-titles', action='store_true',
                        default=False,
                        help="Do not add titles giving the 1D credible range "
                             "over the 1D marginal plots.")
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
    pgroup.add_argument('--pick-injection-by-time', action='store_true',
                        default=False,
                        help="In the case of multiple injections, pick one"
                             " for plotting based on its proximity in time.")
    add_injsamples_map_opt(pgroup)
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
        "--contour-linestyles", type=str, default=None, nargs="+",
        help="Specify the linestyles to use for the contour lines. Defaut "
             "is solid for all.")
    density_group.add_argument(
        "--no-contour-labels", action="store_true", default=False,
        help="Don't put labels on the contours.")
    density_group.add_argument(
        '--use-kombine-kde', default=False, action="store_true",
        help="Use kombine's clustered KDE for determining 2D marginal "
             "contours and density instead of scipy's gaussian_kde (the "
             "default). This is better at distinguishing bimodal "
             "distributions, but is much slower than the default. For speed, "
             "suggest setting --kde-args 'max_samples:20000' or smaller if "
             "using this. Requires kombine to be installed.")
    density_group.add_argument(
        '--max-kde-samples', type=int, default=None,
        help="Limit the number of samples used for KDE construction to the "
             "given value. This can substantially speed up plot generation "
             "(particularly when plotting multiple parameters). Suggested "
             "values: 5000 to 10000.")
    density_group.add_argument(
        '--kde-args', metavar="ARG:VALUE", nargs='+', default=None,
        help="Pass the given argrument, value pairs to the KDE function "
             "(either scipy's or kombine's) when setting it up.")
    return density_group
