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

"""This module contains standard options used for inference-related programs.
"""

import logging
import numpy
import pycbc.inference.sampler
from pycbc import conversions
from pycbc import transforms
from pycbc.distributions import bounded
from pycbc.distributions import constraints
from pycbc.io import InferenceFile
from pycbc.inference import likelihood
from pycbc.workflow import WorkflowConfigParser
from pycbc.pool import choose_pool
from pycbc.psd import from_cli_multi_ifos as psd_from_cli_multi_ifos
from pycbc.strain import from_cli_multi_ifos as strain_from_cli_multi_ifos
from pycbc.gate import gates_from_cli, psd_gates_from_cli, apply_gates_to_td, \
                       apply_gates_to_fd


#-----------------------------------------------------------------------------
#
#                   Utilities for loading config files
#
#-----------------------------------------------------------------------------

def convert_liststring_to_list(lstring):
    """Checks if an argument of the configuration file is a string of a list
    and returns the corresponding list (of strings).

    The argument is considered to be a list if it starts with '[' and ends
    with ']'. List elements should be comma separated. For example, passing
    `'[foo bar, cat]'` will result in `['foo bar', 'cat']` being returned. If
    the argument does not start and end with '[' and ']', the argument will
    just be returned as is.
    """
    if lstring[0]=='[' and lstring[-1]==']':
        lstring = [str(lstring[1:-1].split(',')[n].strip().strip("'"))
                      for n in range(len(lstring[1:-1].split(',')))]
    return lstring


def add_config_opts_to_parser(parser):
    """Adds options for the configuration files to the given parser.
    """
    parser.add_argument("--config-files", type=str, nargs="+", required=True,
                        help="A file parsable by "
                             "pycbc.workflow.WorkflowConfigParser.")
    parser.add_argument("--config-overrides", type=str, nargs="+",
                        default=None, metavar="SECTION:OPTION:VALUE",
                        help="List of section:option:value combinations to "
                             "add into the configuration file.")


def config_parser_from_cli(opts):
    """Loads a config file from the given options, applying any overrides
    specified. Specifically, config files are loaded from the `--config-files`
    options while overrides are loaded from `--config-overrides`.
    """
    # read configuration file
    logging.info("Reading configuration file")
    if opts.config_overrides is not None:
        overrides = [override.split(":") for override in opts.config_overrides]
    else:
        overrides = None
    return WorkflowConfigParser(opts.config_files, overrides)


def read_args_from_config(cp, section_group=None):
    """Given an open config file, loads the static and variable arguments to
    use in the parameter estmation run.

    Parameters
    ----------
    cp : WorkflowConfigParser
        An open config parser to read from.
    section_group : {None, str}
        When reading the config file, only read from sections that begin with
        `{section_group}_`. For example, if `section_group='foo'`, the
        variable arguments will be retrieved from section
        `[foo_variable_args]`. If None, no prefix will be appended to section
        names.

    Returns
    -------
    variable_args : list
        The names of the parameters to vary in the PE run.
    static_args : dict
        Dictionary of names -> values giving the parameters to keep fixed.
    """
    logging.info("Loading arguments")
    if section_group is not None:
        section_prefix = '{}_'.format(section_group)
    else:
        section_prefix = ''

    # sanity check that each parameter in [variable_args] has a priors section
    variable_args = cp.options("{}variable_args".format(section_prefix))
    subsections = cp.get_subsections("{}prior".format(section_prefix))
    tags = numpy.concatenate([tag.split("+") for tag in subsections])
    if not any(param in tags for param in variable_args):
        raise KeyError("You are missing a priors section in the config file.")

    # get parameters that do not change in sampler
    static_args = dict([(key,cp.get_opt_tags(
        "{}static_args".format(section_prefix), key, []))
        for key in cp.options("{}static_args".format(section_prefix))])
    # try converting values to float
    for key,val in static_args.iteritems():
        try:
            # the following will raise a ValueError if it cannot be cast to
            # float (as we would expect for string arguments)
            static_args[key] = float(val)
        except ValueError:
            # try converting to a list of strings; this function will just
            # return val if it does not begin (end) with [ (])
            static_args[key] = convert_liststring_to_list(val) 

    # get additional constraints to apply in prior
    cons = []
    section = "{}constraint".format(section_prefix)
    for subsection in cp.get_subsections(section):
        name = cp.get_opt_tag(section, "name", subsection)
        constraint_arg = cp.get_opt_tag(section, "constraint_arg", subsection)
        kwargs = {}
        for key in cp.options(section + "-" + subsection):
            if key in ["name", "constraint_arg"]:
                continue
            val = cp.get_opt_tag(section, key, subsection)
            if key == "required_parameters":
                kwargs["required_parameters"] = val.split(
                                                        bounded.VARARGS_DELIM)
                continue
            try:
                val = float(val)
            except ValueError:
                pass
            kwargs[key] = val
        cons.append(constraints.constraints[name](variable_args,
                                                  constraint_arg, **kwargs))

    return variable_args, static_args, cons

def read_sampling_args_from_config(cp, section_group=None):
    """Reads sampling parameters from the given config file.

    Parameters are read from the `[({section_group}_)sampling_params]` section.
    The options should list the variable args to transform; the parameters they
    point to should list the parameters they are to be transformed to for
    sampling. If a multiple parameters are transformed together, they should
    be comma separated. Example:

    .. code-block:: ini

        [sampling_parameters]
        mass1, mass2 = mchirp, logitq
        spin1_a = logitspin1_a

    Note that only the final sampling parameters should be listed, even if
    multiple intermediate transforms are needed. (In the above example, a
    transform is needed to go from mass1, mass2 to mchirp, q, then another one
    needed to go from q to logitq.) These transforms should be specified
    in separate sections; see `transforms.read_transforms_from_config` for
    details.

    Parameters
    ----------
    cp : WorkflowConfigParser
        An open config parser to read from.
    section_group : str, optional
        Append `{section_group}_` to the section name. Default is None.

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
    section = section_prefix + 'sampling_parameters'
    replaced_params = set()
    sampling_params = set()
    for args in cp.options(section):
        map_args = cp.get(section, args)
        sampling_params.update(set(map(str.strip, map_args.split(','))))
        replaced_params.update(set(map(str.strip, args.split(',')))) 
    return list(sampling_params), list(replaced_params)


#-----------------------------------------------------------------------------
#
#                    Utilities for setting up a sampler
#
#-----------------------------------------------------------------------------

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
    sampler_group.add_argument("--update-interval", type=int, default=None,
        help="If using kombine, specify the number of steps to take between "
             "proposal updates. Note: for purposes of updating, kombine "
             "counts iterations since the last checkpoint. This interval "
             "should therefore be less than the checkpoint interval, else "
             "no updates will occur. To ensure that updates happen at equal "
             "intervals, make checkpoint-interval a multiple of "
             "update-interval.")
    sampler_group.add_argument("--nprocesses", type=int, default=None,
        help="Number of processes to use. If not given then use maximum.")
    sampler_group.add_argument("--use-mpi", action='store_true', default=False,
        help="Use MPI to parallelize the sampler")
    return sampler_group


def sampler_from_cli(opts, likelihood_evaluator, pool=None):
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
    # Used to help paralleize over multiple cores / MPI
    if opts.nprocesses > 1:
        likelihood._global_instance = likelihood_evaluator
        likelihood_call = likelihood._call_global_likelihood
    else:
        likelihood_call = None

    sclass = pycbc.inference.sampler.samplers[opts.sampler]
    # check for consistency
    if opts.skip_burn_in and opts.min_burn_in is not None:
        raise ValueError("both skip-burn-in and min-burn-in specified")

    pool = choose_pool(mpi=opts.use_mpi, processes=opts.nprocesses)

    if pool is not None:
        pool.count = opts.nprocesses

    return sclass.from_cli(opts, likelihood_evaluator,
                           pool=pool, likelihood_call=likelihood_call)


#-----------------------------------------------------------------------------
#
#                       Utilities for loading data
#
#-----------------------------------------------------------------------------

def add_low_frequency_cutoff_opt(parser):
    """Adds the low-frequency-cutoff option to the given parser."""
    # FIXME: this just uses the same frequency cutoff for every instrument for
    # now. We should allow for different frequency cutoffs to be used; that
    # will require (minor) changes to the Likelihood class
    parser.add_argument("--low-frequency-cutoff", type=float, required=True,
                        help="Low frequency cutoff for each IFO.")


def low_frequency_cutoff_from_cli(opts):
    """Parses the low frequency cutoff from the given options.

    Returns
    -------
    dict
        Dictionary of instruments -> low frequency cutoff.
    """
    # FIXME: this just uses the same frequency cutoff for every instrument for
    # now. We should allow for different frequency cutoffs to be used; that
    # will require (minor) changes to the Likelihood class
    return {ifo: opts.low_frequency_cutoff for ifo in opts.instruments}


def data_from_cli(opts):
    """Loads the data needed for a likelihood evaluator from the given
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
    strain_dict = strain_from_cli_multi_ifos(opts, opts.instruments,
                                             precision="double")
    # apply gates if not waiting to overwhiten
    if not opts.gate_overwhitened:
        logging.info("Applying gates to strain data")
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
                                                    opts.instruments,
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
    logging.info("FFT strain")
    stilde_dict = {}
    length_dict = {}
    delta_f_dict = {}
    low_frequency_cutoff_dict = low_frequency_cutoff_from_cli(opts)
    for ifo in opts.instruments:
        stilde_dict[ifo] = strain_dict[ifo].to_frequencyseries()
        length_dict[ifo] = len(stilde_dict[ifo])
        delta_f_dict[ifo] = stilde_dict[ifo].delta_f

    # get PSD as frequency series
    psd_dict = psd_from_cli_multi_ifos(opts, length_dict, delta_f_dict,
                               low_frequency_cutoff_dict, opts.instruments,
                               strain_dict=psd_strain_dict, precision="double")

    # apply any gates to overwhitened data, if desired
    if opts.gate_overwhitened and opts.gate is not None:
        logging.info("Applying gates to overwhitened data")
        # overwhiten the data
        for ifo in gates:
            stilde_dict[ifo] /= psd_dict[ifo]
        stilde_dict = apply_gates_to_fd(stilde_dict, gates)
        # unwhiten the data for the likelihood generator
        for ifo in gates:
            stilde_dict[ifo] *= psd_dict[ifo]

    return strain_dict, stilde_dict, psd_dict



#-----------------------------------------------------------------------------
#
#                Utilities for loading and plotting results
#
#-----------------------------------------------------------------------------

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
    results_reading_group.add_argument(
        "--parameters-group", type=str, default=InferenceFile.samples_group,
        choices=[InferenceFile.samples_group, InferenceFile.stats_group],
        help="Group in the HDF InferenceFile to look for parameters.")
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
        as a FieldArray instance.
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
    samples : {None, FieldArray}
        If load_samples, the samples as a FieldArray; otherwise, None.
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

    # load the samples
    if load_samples:
        logging.info("Loading samples")
        # check if need extra parameters for a non-sampling parameter
        file_parameters, ts = transforms.get_common_cbc_transforms(
                                                 parameters, fp.variable_args)
        # read samples from file
        samples = fp.read_samples(
            file_parameters, walkers=walkers,
            thin_start=opts.thin_start, thin_interval=opts.thin_interval,
            thin_end=opts.thin_end, iteration=opts.iteration,
            samples_group=opts.parameters_group)
        # add parameters not included in file
        samples = transforms.apply_transforms(samples, ts)
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
        zvals = conversions.snr_from_loglr(likelihood_stats.loglr)
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
    # add expected parameters options
    pgroup.add_argument('--expected-parameters', nargs='+', metavar='PARAM:VAL',
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
