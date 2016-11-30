# Copyright (C) 2016 Collin Capano
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

"""
This module contains standard options used for inference-related programs.
"""

import logging
from pycbc.io import InferenceFile
import pycbc.inference.sampler
from pycbc.inference import likelihood


def add_sampler_option_group(parser):
    """Adds the options needed to set up an inference sampler.

    Parameters
    ----------
    parser : object
        ArgumentParser instance.
    """
    sampler_group = parser.add_argument_group("Arguments for setting up "
        "a sampler")

    # required options
    sampler_group.add_argument("--sampler", required=True,
        choices=pycbc.inference.sampler.samplers.keys(),
        help="Sampler class to use for finding posterior.")
    sampler_group.add_argument("--niterations", type=int, required=True,
        help="Number of iterations to perform after burn in.")
    sampler_group.add_argument("--nprocesses", type=int, default=None,
        help="Number of processes to use. If not given then use maximum.")
    # sampler-specific options
    sampler_group.add_argument("--nwalkers", type=int, default=None,
        help="Number of walkers to use in sampler. Required for MCMC "
             "samplers.")
    sampler_group.add_argument("--ntemps", type=int, default=None,
        help="Number of temperatures to use in sampler. Required for parallel "
             "tempered MCMC samplers.")
    sampler_group.add_argument("--min-burn-in", type=int, default=None,
        help="Force the burn-in to be at least the given number of "
             "iterations. If a sampler has an internal algorithm for "
             "determining the burn-in size (e.g., kombine), and it returns "
             "a value < this, the burn-in will be repeated until the "
             "number of iterations is at least this value.")
    sampler_group.add_argument("--skip-burn-in", action="store_true",
        default=False,
        help="Do not burn in with sampler. An error will be raised if "
             "min-burn-in is also provided.")

    return sampler_group


def sampler_from_cli(opts, likelihood_evaluator):
    """Parses the given command-line options to set up a sampler.

    Parameters
    ----------
    opts : object
        ArgumentParser options.
    likelihood_evaluator : LikelihoodEvaluator
        The likelihood evaluator to use with the sampler.

    Returns
    -------
    pycbc.inference.sampler
        A sampler initialized based on the given arguments.
    """
    sclass = pycbc.inference.sampler.samplers[opts.sampler]
    # check for consistency
    if opts.skip_burn_in and opts.min_burn_in is not None:
        raise ValueError("both skip-burn-in and min-burn-in specified")
    return sclass.from_cli(opts, likelihood_evaluator)


def add_inference_results_option_group(parser):
    """Adds the options used to call pycbc.inference.results_from_cli function
    to an argument parser. These are options releated to loading the results
    from a run of pycbc_inference, for purposes of plotting and/or creating
    tables.

    Parameters
    ----------
    parser : object
        ArgumentParser instance.
    """

    results_reading_group = parser.add_argument_group("Arguments for loading "
        "inference results")

    # required options
    results_reading_group.add_argument("--input-file", type=str, required=True,
        help="Path to input HDF file.")
    results_reading_group.add_argument("--parameters", type=str, nargs="+",
        metavar="PARAM[:LABEL]",
        help="Name of parameters to plot. If none provided will load all of "
             "the variable args in the input-file. If provided, the "
             "parameters can be any of the variable args or posteriors in "
             "the input file, derived parameters from them, or any function "
             "of them. Syntax for functions is python; any math functions in "
             "the numpy libary may be used. Can optionally also specify a "
             "label for each parameter. If no label is provided, will try to "
             "retrieve a label from the input-file. If no label can be found "
             "in the input-file, will try to get a label from "
             "pycbc.waveform.parameters. If no label can be found in either "
             "place, will just use the parameter.")

    # optionals
    results_reading_group.add_argument("--thin-start", type=int, default=None,
        help="Sample number to start collecting samples to plot. If none "
             "provided, will start at the end of the burn-in.")
    results_reading_group.add_argument("--thin-interval", type=int,
        default=None,
        help="Interval to use for thinning samples. If none provided, will "
             "use the auto-correlation length found in the file.")
    results_reading_group.add_argument("--thin-end", type=int, default=None,
        help="Sample number to stop collecting samples to plot. If none "
             "provided, will stop at the last sample from the sampler.")
    results_reading_group.add_argument("--iteration", type=int, default=None,
        help="Only retrieve the given iteration. To load the last n-th sampe "
             "use -n, e.g., -1 will load the last iteration. This overrides "
             "the thin-start/interval/end options.")

    return results_reading_group


def results_from_cli(opts, load_samples=True, walkers=None):
    """
    Loads an inference result file along with any labels associated with it
    from the command line options.

    Parameters
    ----------
    opts : ArgumentParser options
        The options from the command line.
    load_samples : {True, bool}
        Load samples from the results file using the parameters, thin_start,
        and thin_interval specified in the options. The samples are returned
        as a WaveformArray instance.
    walkers : {None, (list of) int}
        If loading samples, the walkers to load from. If None, will load from
        all walkers.

    Returns
    -------
    result_file : pycbc.io.InferenceFile
        The result file as an InferenceFile.
    parameters : list
        List of the parameters to use, parsed from the parameters option.
    labels : list
        List of labels to associate with the parameters.
    samples : {None, WaveformArray}
        If load_samples, the samples as a WaveformArray; otherwise, None.
    """
    logging.info("Reading input file")
    fp = InferenceFile(opts.input_file, "r")
    parameters = fp.variable_args if opts.parameters is None \
                 else opts.parameters
    # load the labels
    labels = []
    for ii,p in enumerate(parameters):
        if len(p.split(':')) == 2:
            p, label = p.split(':')
            parameters[ii] = p
        else:
            label = fp.read_label(p)
        labels.append(label)
    if load_samples:
        logging.info("Loading samples")
        samples = fp.read_samples(parameters, walkers=walkers,
            thin_start=opts.thin_start, thin_interval=opts.thin_interval,
            thin_end=opts.thin_end, iteration=opts.iteration)
    else:
        samples = None
    return fp, parameters, labels, samples


def get_zvalues(fp, arg, likelihood_stats):
    """Reads the data for the z-value of the plots from the inference file.

    Parameters
    ----------
    fp : InferenceFile
        An open inference file; needed to get the value of the log noise
        likelihood.
    arg : str
        The argument to plot; must be one of `loglr`, `snr`, `logplr`,
        `logposterior`, or `prior`. If not one of these, a ValueError is
        raised.
    likelihood_stats : FieldArray
        The likelihood stats; the sort of thing returned by
        `fp.read_likelihood_stats`.

    Returns
    -------
    zvals : numpy.array
        An array of the desired likelihood values to plot.
    zlbl : str
        The label to use for the values on a plot.
    """
    if arg == 'loglr':
        zvals = likelihood_stats.loglr
        zlbl = r'$\log\mathcal{L}(\vec{\vartheta})$'
    elif arg == 'snr':
        zvals = likelihood.snr_from_loglr(likelihood_stats.loglr)
        zlbl = r'$\rho(\vec{\vartheta})$'
    elif arg == 'logplr':
        zvals = likelihood_stats.loglr + likelihood_stats.prior
        zlbl = r'$\log[\mathcal{L}(\vec{\vartheta})p(\vec{\vartheta})]$'
    elif arg == 'logposterior':
        zvals = likelihood_stats.loglr + likelihood_stats.prior + fp.lognl
        zlbl = r'$\log[p(d|\vec{\vartheta})p(\vec{\vartheta})]$'
    elif arg == 'prior':
        zvals = likelihood_stats.prior
        zlbl = r'$\log p(\vec{\vartheta})$'
    else:
        raise ValueError("Unrecognized arg {}".format(arg))
    return zvals, zlbl


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
    pgroup.add_argument("--plot-scatter", action='store_true', default=False,
                        help="Plot each sample point as a scatter plot.")
    pgroup.add_argument("--plot-density", action="store_true", default=False,
                        help="Plot the posterior density as a color map.")
    pgroup.add_argument("--plot-contours", action="store_true", default=False,
                        help="Draw contours showing the 50th and 90th "
                             "percentile confidence regions.")
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


def add_scatter_option_group(parser):
    """Adds the options needed to configure scatter plots.

    Parameters
    ----------
    parser : object
        ArgumentParser instance.
    """
    scatter_group = parser.add_argument_group("Options for configuring the "
                                              "scatter plot.")

    scatter_group.add_argument('--z-arg', type=str, default=None,
                    choices=['loglr', 'snr', 'logplr', 'logposterior',
                             'prior'],
                    help='What to color the scatter points by. If not set, '
                         'all points will be the same color. Choices are: '
                         'loglr: the log likelihood ratio; snr: SNR; '
                         'logplr: loglr + log of the prior; '
                         'logposterior: log likelihood function + log prior; '
                         'prior: the log of the prior.')
    scatter_group.add_argument("--vmin", type=float,
                    help="Minimum value for the colorbar.")
    scatter_group.add_argument("--vmax", type=float,
                    help="Maximum value for the colorbar.")
    scatter_group.add_argument("--scatter-cmap", type=str, default='plasma',
                    help="Specify the colormap to use for points. Default is "
                         "plasma.")

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

    density_group.add_argument("--density-cmap", type=str, default='viridis',
                    help="Specify the colormap to use for the density. "
                         "Default is viridis.")
    density_group.add_argument("--contour-color", type=str,
                    help="Specify the color to use for the contour lines. "
                         "Default is white for density plots and black "
                         "for scatter plots.")
    density_group.add_argument('--use-kombine-kde', default=False,
                    action="store_true",
                    help="Use kombine's KDE for determining contours. Default "
                         "is to use scipy's gaussian_kde.")

    return density_group
