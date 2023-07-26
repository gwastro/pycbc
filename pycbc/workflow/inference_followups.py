# Copyright (C) 2016 Christopher M. Biwer, Alexander Harvey Nitz, Collin Capano
#
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
Module that contains functions for setting up the inference workflow.
"""
from pycbc.workflow.core import (Executable, makedir)
from pycbc.workflow.plotting import PlotExecutable
from pycbc.results import layout


def make_inference_plot(workflow, input_file, output_dir,
                        name, analysis_seg=None,
                        tags=None, input_file_opt='input-file',
                        output_file_extension='.png',
                        add_to_workflow=False):
    """Boiler-plate function for creating a standard plotting job.

    Parameters
    ----------
    workflow: pycbc.workflow.Workflow
        The core workflow instance we are populating
    input_file: (list of) pycbc.workflow.File
        The file used for the input. May provide either a single file or a
        list of files.
    output_dir: str
        The directory to store result plots.
    name: str
        The name in the [executables] section of the configuration file
        to use.
    analysis_segs: ligo.segments.Segment, optional
       The segment this job encompasses. If None then use the total analysis
       time from the workflow.
    tags: list, optional
        Tags to add to the inference executables.
    input_file_opt : str, optional
        The name of the input-file option used by the executable. Default
        is ``input-file``.
    output_file_extension : str, optional
        What file type to create. Default is ``.png``.
    add_to_workflow : bool, optional
        If True, the node will be added to the workflow before being returned.
        **This means that no options may be added to the node afterward.**
        Default is ``False``.

    Returns
    -------
    pycbc.workflow.plotting.PlotExecutable
        The job node for creating the plot.
    """
    # default values
    if tags is None:
        tags = []
    if analysis_seg is None:
        analysis_seg = workflow.analysis_time
    # make the directory that will contain the output files
    makedir(output_dir)
    # Catch if a parameters option was specified:
    # we need to do this because PlotExecutable will automatically add any
    # option in the section to the node. However, we need to add the
    # appropriate escapes to the parameters option so pegasus will render it
    # properly (see _params_for_pegasus for details).
    parameters = None
    if workflow.cp.has_option(name, 'parameters'):
        parameters = workflow.cp.get(name, 'parameters')
        workflow.cp.remove_option(name, 'parameters')
    # make a node for plotting the posterior as a corner plot
    node = PlotExecutable(workflow.cp, name, ifos=workflow.ifos,
                          out_dir=output_dir,
                          tags=tags).create_node()
    # add back the parameters option if it was specified
    if parameters is not None:
        node.add_opt("--parameters", _params_for_pegasus(parameters))
        # and put the opt back in the config file in memory
        workflow.cp.set(name, 'parameters', parameters)
    # add input and output options
    if isinstance(input_file, list):
        # list of input files are given, use input_list_opt
        node.add_input_list_opt("--{}".format(input_file_opt), input_file)
    else:
        # assume just a single file
        node.add_input_opt("--{}".format(input_file_opt), input_file)
    node.new_output_file_opt(analysis_seg, output_file_extension,
                             "--output-file")
    # add node to workflow
    if add_to_workflow:
        workflow += node
    return node


def make_inference_prior_plot(workflow, config_file, output_dir,
                              name="plot_prior",
                              analysis_seg=None, tags=None):
    """Sets up the corner plot of the priors in the workflow.

    Parameters
    ----------
    workflow: pycbc.workflow.Workflow
        The core workflow instance we are populating
    config_file: pycbc.workflow.File
        The WorkflowConfigParser parasable inference configuration file..
    output_dir: str
        The directory to store result plots and files.
    name: str
        The name in the [executables] section of the configuration file
        to use, and the section to read for additional arguments to pass to
        the executable. Default is ``plot_prior``.
    analysis_segs: ligo.segments.Segment, optional
       The segment this job encompasses. If None then use the total analysis
       time from the workflow.
    tags: list, optional
        Tags to add to the inference executables.

    Returns
    -------
    pycbc.workflow.FileList
        A list of the output files.
    """
    node = make_inference_plot(workflow, config_file, output_dir,
                               name, analysis_seg=analysis_seg, tags=tags,
                               input_file_opt='config-file',
                               add_to_workflow=True)
    return node.output_files


def create_posterior_files(workflow, samples_files, output_dir,
                           parameters=None, name="extract_posterior",
                           analysis_seg=None, tags=None):
    """Sets up job to create posterior files from some given samples files.

    Parameters
    ----------
    workflow: pycbc.workflow.Workflow
        The workflow instance we are populating
    samples_files : str or list of str
        One or more files to extract the posterior samples from.
    output_dir: str
        The directory to store result plots and files.
    name: str, optional
        The name in the [executables] section of the configuration file
        to use, and the section to read for additional arguments to pass to
        the executable. Default is ``extract_posterior``.
    analysis_segs: ligo.segments.Segment, optional
       The segment this job encompasses. If None then use the total analysis
       time from the workflow.
    tags: list, optional
        Tags to add to the inference executables.

    Returns
    -------
    pycbc.workflow.FileList
        A list of output files.
    """
    if analysis_seg is None:
        analysis_seg = workflow.analysis_time
    if tags is None:
        tags = []
    # Catch if a parameters option was specified:
    # we need to do this because Executable will automatically add any
    # option in the section to the node. However, we need to add the
    # appropriate escapes to the parameters option so pegasus will render it
    # properly (see _params_for_pegasus for details).
    parameters = None
    if workflow.cp.has_option(name, 'parameters'):
        parameters = workflow.cp.get(name, 'parameters')
        workflow.cp.remove_option(name, 'parameters')
    extract_posterior_exe = Executable(workflow.cp, name,
                                       ifos=workflow.ifos,
                                       out_dir=output_dir)
    node = extract_posterior_exe.create_node()
    # add back the parameters option if it was specified
    if parameters is not None:
        node.add_opt("--parameters", _params_for_pegasus(parameters))
        # and put the opt back in the config file in memory
        workflow.cp.set(name, 'parameters', parameters)
    if not isinstance(samples_files, list):
        samples_files = [samples_files]
    node.add_input_list_opt("--input-file", samples_files)
    node.new_output_file_opt(analysis_seg, ".hdf", "--output-file", tags=tags)
    # add node to workflow
    workflow += node
    return node.output_files


def create_fits_file(workflow, inference_file, output_dir,
                     name="create_fits_file",
                     analysis_seg=None, tags=None):
    """Sets up job to create fits files from some given samples files.

    Parameters
    ----------
    workflow: pycbc.workflow.Workflow
        The workflow instance we are populating
    inference_file: pycbc.workflow.File
        The file with posterior samples.
    output_dir: str
        The directory to store result plots and files.
    name: str, optional
        The name in the [executables] section of the configuration file
        to use, and the section to read for additional arguments to pass to
        the executable. Default is ``create_fits_file``.
    analysis_segs: ligo.segments.Segment, optional
       The segment this job encompasses. If None then use the total analysis
       time from the workflow.
    tags: list, optional
        Tags to add to the inference executables.

    Returns
    -------
    pycbc.workflow.FileList
        A list of output files.
    """
    if analysis_seg is None:
        analysis_seg = workflow.analysis_time
    if tags is None:
        tags = []
    create_fits_exe = Executable(workflow.cp, name,
                                 ifos=workflow.ifos,
                                 out_dir=output_dir)
    node = create_fits_exe.create_node()
    node.add_input_opt("--input-file", inference_file)
    node.new_output_file_opt(analysis_seg, ".fits", "--output-file", tags=tags)
    # add node to workflow
    workflow += node
    return node.output_files


def make_inference_skymap(workflow, fits_file, output_dir,
                          name="plot_skymap", analysis_seg=None,
                          tags=None):
    """Sets up the skymap plot.

    Parameters
    ----------
    workflow: pycbc.workflow.Workflow
        The core workflow instance we are populating
    fits_file: pycbc.workflow.File
        The fits file with the sky location.
    output_dir: str
        The directory to store result plots and files.
    name: str, optional
        The name in the [executables] section of the configuration file
        to use, and the section to read for additional arguments to pass to
        the executable. Default is ``plot_skymap``.
    analysis_segs: ligo.segments.Segment, optional
       The segment this job encompasses. If None then use the total analysis
       time from the workflow.
    tags: list, optional
        Tags to add to the inference executables.

    Returns
    -------
    pycbc.workflow.FileList
        A list of result and output files.
    """
    node = make_inference_plot(workflow, fits_file, output_dir,
                               name, analysis_seg=analysis_seg, tags=tags,
                               add_to_workflow=True)
    return node.output_files


def make_inference_summary_table(workflow, inference_file, output_dir,
                                 parameters=None, print_metadata=None,
                                 name="table_summary",
                                 analysis_seg=None, tags=None):
    """Sets up the html table summarizing parameter estimates.

    Parameters
    ----------
    workflow: pycbc.workflow.Workflow
        The core workflow instance we are populating
    inference_file: pycbc.workflow.File
        The file with posterior samples.
    output_dir: str
        The directory to store result plots and files.
    parameters : list or str
        A list or string of parameters to generate the table for. If a string
        is provided, separate parameters should be space or new-line separated.
    print_metadata : list or str
        A list or string of metadata parameters to print. Syntax is the same
        as for ``parameters``.
    name: str, optional
        The name in the [executables] section of the configuration file
        to use, and the section to read for additional arguments to pass to
        the executable. Default is ``table_summary``.
    analysis_segs: ligo.segments.Segment, optional
       The segment this job encompasses. If None then use the total analysis
       time from the workflow.
    tags: list, optional
        Tags to add to the inference executables.

    Returns
    -------
    pycbc.workflow.FileList
        A list of output files.
    """
    # we'll use make_inference_plot even though this isn't a plot; the
    # setup is the same, we just change the file extension
    node = make_inference_plot(workflow, inference_file, output_dir,
                               name, analysis_seg=analysis_seg, tags=tags,
                               output_file_extension='.html',
                               add_to_workflow=False)
    # now add the parameters and print metadata options; these are pulled
    # from separate sections in the workflow config file, which is why we
    # add them separately here
    if parameters is not None:
        node.add_opt("--parameters", _params_for_pegasus(parameters))
    if print_metadata is not None:
        node.add_opt("--print-metadata", _params_for_pegasus(print_metadata))
    workflow += node
    return node.output_files


def make_inference_posterior_plot(workflow, inference_file, output_dir,
                                  parameters=None, plot_prior_from_file=None,
                                  name="plot_posterior",
                                  analysis_seg=None, tags=None):
    """Sets up the corner plot of the posteriors in the workflow.

    Parameters
    ----------
    workflow: pycbc.workflow.Workflow
        The core workflow instance we are populating
    inference_file: pycbc.workflow.File
        The file with posterior samples.
    output_dir: str
        The directory to store result plots and files.
    parameters : list or str
        The parameters to plot.
    plot_prior_from_file : str, optional
        Plot the prior from the given config file on the 1D marginal plots.
    name: str, optional
        The name in the [executables] section of the configuration file
        to use, and the section to read for additional arguments to pass to
        the executable. Default is ``plot_posterior``.
    analysis_segs: ligo.segments.Segment, optional
       The segment this job encompasses. If None then use the total analysis
       time from the workflow.
    tags: list, optional
        Tags to add to the inference executables.

    Returns
    -------
    pycbc.workflow.FileList
        A list of output files.
    """
    # create the node, but delay adding it to the workflow so we can add
    # the prior file if it is requested
    node = make_inference_plot(workflow, inference_file, output_dir,
                               name, analysis_seg=analysis_seg, tags=tags,
                               add_to_workflow=False)
    if parameters is not None:
        node.add_opt("--parameters", _params_for_pegasus(parameters))
    if plot_prior_from_file is not None:
        node.add_input_opt('--plot-prior', plot_prior_from_file)
    # now add the node to workflow
    workflow += node
    return node.output_files


def make_inference_samples_plot(workflow, inference_file, output_dir,
                                name="plot_samples",
                                analysis_seg=None, tags=None):
    """Sets up a plot of the samples versus iteration (for MCMC samplers).

    Parameters
    ----------
    workflow: pycbc.workflow.Workflow
        The core workflow instance we are populating
    inference_file: pycbc.workflow.File
        The file with posterior samples.
    output_dir: str
        The directory to store result plots and files.
    name: str, optional
        The name in the [executables] section of the configuration file
        to use, and the section to read for additional arguments to pass to
        the executable. Default is ``plot_samples``.
    analysis_segs: ligo.segments.Segment, optional
       The segment this job encompasses. If None then use the total analysis
       time from the workflow.
    tags: list, optional
        Tags to add to the inference executables.

    Returns
    -------
    pycbc.workflow.FileList
        A list of output files.
    """
    node = make_inference_plot(workflow, inference_file, output_dir,
                               name, analysis_seg=analysis_seg, tags=tags,
                               add_to_workflow=True)
    return node.output_files


def make_inference_acceptance_rate_plot(workflow, inference_file, output_dir,
                                        name="plot_acceptance_rate",
                                        analysis_seg=None, tags=None):
    """Sets up a plot of the acceptance rate (for MCMC samplers).

    Parameters
    ----------
    workflow: pycbc.workflow.Workflow
        The core workflow instance we are populating
    inference_file: pycbc.workflow.File
        The file with posterior samples.
    output_dir: str
        The directory to store result plots and files.
    name: str, optional
        The name in the [executables] section of the configuration file
        to use, and the section to read for additional arguments to pass to
        the executable. Default is ``plot_acceptance_rate``.
    analysis_segs: ligo.segments.Segment, optional
       The segment this job encompasses. If None then use the total analysis
       time from the workflow.
    tags: list, optional
        Tags to add to the inference executables.

    Returns
    -------
    pycbc.workflow.FileList
        A list of output files.
    """
    node = make_inference_plot(workflow, inference_file, output_dir,
                               name, analysis_seg=analysis_seg, tags=tags,
                               add_to_workflow=True)
    return node.output_files


def make_inference_plot_mcmc_history(workflow, inference_file, output_dir,
                                     name="plot_mcmc_history",
                                     analysis_seg=None, tags=None):
    """Sets up a plot showing the checkpoint history of an MCMC sampler.

    Parameters
    ----------
    workflow: pycbc.workflow.Workflow
        The core workflow instance we are populating
    inference_file: pycbc.workflow.File
        The file with posterior samples.
    output_dir: str
        The directory to store result plots and files.
    name: str, optional
        The name in the [executables] section of the configuration file
        to use, and the section to read for additional arguments to pass to
        the executable. Default is ``plot_mcmc_history``.
    analysis_segs: ligo.segments.Segment, optional
       The segment this job encompasses. If None then use the total analysis
       time from the workflow.
    tags: list, optional
        Tags to add to the inference executables.

    Returns
    -------
    pycbc.workflow.FileList
        A list of output files.
    """
    node = make_inference_plot(workflow, inference_file, output_dir,
                               name, analysis_seg=analysis_seg, tags=tags,
                               add_to_workflow=True)
    return node.output_files


def make_inference_dynesty_run_plot(workflow, inference_file, output_dir,
                                    name="plot_dynesty_run",
                                    analysis_seg=None, tags=None):
    """Sets up a debugging plot for the dynesty run (for Dynesty sampler).

    Parameters
    ----------
    workflow: pycbc.workflow.Workflow
        The core workflow instance we are populating
    inference_file: pycbc.workflow.File
        The file with posterior samples.
    output_dir: str
        The directory to store result plots and files.
    name: str, optional
        The name in the [executables] section of the configuration file
        to use, and the section to read for additional arguments to pass to
        the executable. Default is ``plot_dynesty_run``.
    analysis_segs: ligo.segments.Segment, optional
       The segment this job encompasses. If None then use the total analysis
       time from the workflow.
    tags: list, optional
        Tags to add to the inference executables.

    Returns
    -------
    pycbc.workflow.FileList
        A list of output files.
    """
    node = make_inference_plot(workflow, inference_file, output_dir,
                               name, analysis_seg=analysis_seg, tags=tags,
                               add_to_workflow=True)
    return node.output_files


def make_inference_dynesty_trace_plot(workflow, inference_file, output_dir,
                                      name="plot_dynesty_traceplot",
                                      analysis_seg=None, tags=None):
    """Sets up a trace plot for the dynesty run (for Dynesty sampler).

    Parameters
    ----------
    workflow: pycbc.workflow.Workflow
        The core workflow instance we are populating
    inference_file: pycbc.workflow.File
        The file with posterior samples.
    output_dir: str
        The directory to store result plots and files.
    name: str, optional
        The name in the [executables] section of the configuration file
        to use, and the section to read for additional arguments to pass to
        the executable. Default is ``plot_dynesty_traceplot``.
    analysis_segs: ligo.segments.Segment, optional
       The segment this job encompasses. If None then use the total analysis
       time from the workflow.
    tags: list, optional
        Tags to add to the inference executables.

    Returns
    -------
    pycbc.workflow.FileList
        A list of output files.
    """
    node = make_inference_plot(workflow, inference_file, output_dir,
                               name, analysis_seg=analysis_seg, tags=tags,
                               add_to_workflow=True)
    return node.output_files


def make_inference_pp_table(workflow, posterior_files, output_dir,
                            parameters=None, injection_samples_map=None,
                            name="pp_table_summary",
                            analysis_seg=None, tags=None):
    """Performs a PP, writing results to an html table.

    Parameters
    ----------
    workflow : pycbc.workflow.Workflow
        The core workflow instance we are populating
    posterior_files : pycbc.workflow.core.FileList
        List of files with posteriors of injections.
    output_dir : str
        The directory to store result plots and files.
    parameters : list or str, optional
        A list or string of parameters to generate the table for. If a string
        is provided, separate parameters should be space or new-line separated.
    injection_samples_map : (list of) str, optional
        Map between injection parameters and parameters in the posterior file.
        Format is ``INJECTION_PARAM:SAMPLES_PARAM``.
    name : str, optional
        The name in the [executables] section of the configuration file
        to use, and the section to read for additional arguments to pass to
        the executable. Default is ``table_summary``.
    analysis_segs : ligo.segments.Segment, optional
       The segment this job encompasses. If None then use the total analysis
       time from the workflow.
    tags : list, optional
        Tags to add to the inference executables.

    Returns
    -------
    pycbc.workflow.FileList
        A list of output files.
    """
    # we'll use make_inference_plot even though this isn't a plot; the
    # setup is the same, we just change the file extension
    node = make_inference_plot(workflow, posterior_files, output_dir,
                               name, analysis_seg=analysis_seg, tags=tags,
                               output_file_extension='.html',
                               add_to_workflow=False)
    # add the parameters and inj/samples map
    if parameters is not None:
        node.add_opt("--parameters", _params_for_pegasus(parameters))
    if injection_samples_map is not None:
        node.add_opt("--injection-samples-map",
                     _params_for_pegasus(injection_samples_map))
    workflow += node
    return node.output_files


def make_inference_pp_plot(workflow, posterior_files, output_dir,
                           parameters=None, injection_samples_map=None,
                           name="plot_pp",
                           analysis_seg=None, tags=None):
    """Sets up a pp plot in the workflow.

    Parameters
    ----------
    workflow: pycbc.workflow.Workflow
        The core workflow instance we are populating
    posterior_files: pycbc.workflow.core.FileList
        List of files with posteriors of injections.
    output_dir: str
        The directory to store result plots and files.
    parameters : list or str, optional
        The parameters to plot.
    injection_samples_map : (list of) str, optional
        Map between injection parameters and parameters in the posterior file.
        Format is ``INJECTION_PARAM:SAMPLES_PARAM``.
    name: str, optional
        The name in the [executables] section of the configuration file
        to use, and the section to read for additional arguments to pass to
        the executable. Default is ``plot_pp``.
    analysis_segs: ligo.segments.Segment, optional
       The segment this job encompasses. If None then use the total analysis
       time from the workflow.
    tags: list, optional
        Tags to add to the inference executables.

    Returns
    -------
    pycbc.workflow.FileList
        A list of output files.
    """
    node = make_inference_plot(workflow, posterior_files, output_dir,
                               name, analysis_seg=analysis_seg, tags=tags,
                               add_to_workflow=False)
    # add the parameters and inj/samples map
    if parameters is not None:
        node.add_opt("--parameters", _params_for_pegasus(parameters))
    if injection_samples_map is not None:
        node.add_opt("--injection-samples-map",
                     _params_for_pegasus(injection_samples_map))
    # now add the node to workflow
    workflow += node
    return node.output_files


def make_inference_inj_recovery_plot(workflow, posterior_files, output_dir,
                                     parameter, injection_samples_map=None,
                                     name="inj_recovery",
                                     analysis_seg=None, tags=None):
    """Sets up the recovered versus injected parameter plot in the workflow.

    Parameters
    ----------
    workflow: pycbc.workflow.Workflow
        The core workflow instance we are populating
    inference_files: pycbc.workflow.core.FileList
        List of files with posteriors of injections.
    output_dir: str
        The directory to store result plots and files.
    parameter : str
        The parameter to plot.
    injection_samples_map : (list of) str, optional
        Map between injection parameters and parameters in the posterior file.
        Format is ``INJECTION_PARAM:SAMPLES_PARAM``.
    name: str, optional
        The name in the [executables] section of the configuration file
        to use, and the section to read for additional arguments to pass to
        the executable. Default is ``inj_recovery``.
    analysis_segs: ligo.segments.Segment, optional
       The segment this job encompasses. If None then use the total analysis
       time from the workflow.
    tags: list, optional
        Tags to add to the inference executables.

    Returns
    -------
    pycbc.workflow.FileList
        A list of output files.
    """
    # arguments are the same as plot_pp, so just call that with the
    # different executable name
    return make_inference_pp_plot(
        workflow, posterior_files, output_dir, parameters=parameter,
        injection_samples_map=injection_samples_map,
        name=name, analysis_seg=analysis_seg, tags=tags)


def get_plot_group(cp, section_tag):
    """Gets plotting groups from ``[workflow-section_tag]``."""
    group_prefix = "plot-group-"
    # parameters for the summary plots
    plot_groups = {}
    opts = [opt for opt in cp.options("workflow-{}".format(section_tag))
            if opt.startswith(group_prefix)]
    for opt in opts:
        group = opt.replace(group_prefix, "").replace("-", "_")
        plot_groups[group] = cp.get_opt_tag("workflow", opt, section_tag)
    return plot_groups


def get_diagnostic_plots(workflow):
    """Determines what diagnostic plots to create based on workflow.

    The plots to create are based on what executable's are specified in the
    workflow's config file. A list of strings is returned giving the diagnostic
    plots to create. This list may contain:

    * ``samples``: For MCMC samplers, a plot of the sample chains as a function
      of iteration. This will be created if ``plot_samples`` is in the
      executables section.
    * ``acceptance_rate``: For MCMC samplers, a plot of the acceptance rate.
      This will be created if ``plot_acceptance_rate`` is in the executables
      section.

    Returns
    -------
    list :
        List of names of diagnostic plots.
    """
    diagnostics = []
    if "plot_samples" in workflow.cp.options("executables"):
        diagnostics.append('samples')
    if "plot_acceptance_rate" in workflow.cp.options("executables"):
        diagnostics.append('acceptance_rate')
    if "plot_mcmc_history" in workflow.cp.options("executables"):
        diagnostics.append('mcmc_history')
    if "plot_dynesty_run" in workflow.cp.options("executables"):
        diagnostics.append('dynesty_run')
    if "plot_dynesty_traceplot" in workflow.cp.options("executables"):
        diagnostics.append('dynesty_traceplot')
    return diagnostics


def make_diagnostic_plots(workflow, diagnostics, samples_file, label, rdir,
                          tags=None):
    """Makes diagnostic plots.

    Diagnostic plots are sampler-specific plots the provide information on
    how the sampler performed. All diagnostic plots use the output file
    produced by ``pycbc_inference`` as their input. Diagnostic plots are added
    to the results directory ``rdir/NAME`` where ``NAME`` is the name of the
    diagnostic given in ``diagnostics``.

    Parameters
    ----------
    workflow : pycbc.workflow.core.Workflow
        The workflow to add the plotting jobs to.
    diagnostics : list of str
        The names of the diagnostic plots to create. See
        :py:func:`get_diagnostic_plots` for recognized names.
    samples_file : (list of) pycbc.workflow.File
        One or more samples files with which to create the diagnostic plots.
        If a list of files is provided, a diagnostic plot for each file will
        be created.
    label : str
        Event label for the diagnostic plots.
    rdir : pycbc.results.layout.SectionNumber
        Results directory layout.
    tags : list of str, optional
        Additional tags to add to the file names.

    Returns
    -------
    dict :
        Dictionary of diagnostic name -> list of files giving the plots that
        will be created.
    """
    if tags is None:
        tags = []
    out = {}
    if not isinstance(samples_file, list):
        samples_file = [samples_file]
    if 'samples' in diagnostics:
        # files for samples summary subsection
        base = "samples/{}".format(label)
        samples_plots = []
        for kk, sf in enumerate(samples_file):
            samples_plots += make_inference_samples_plot(
                workflow, sf, rdir[base],
                analysis_seg=workflow.analysis_time,
                tags=tags+[label, str(kk)])
        out['samples'] = samples_plots
        layout.group_layout(rdir[base], samples_plots)

    if 'acceptance_rate' in diagnostics:
        # files for samples acceptance_rate subsection
        base = "acceptance_rate/{}".format(label)
        acceptance_plots = []
        for kk, sf in enumerate(samples_file):
            acceptance_plots += make_inference_acceptance_rate_plot(
                workflow, sf, rdir[base],
                analysis_seg=workflow.analysis_time,
                tags=tags+[label, str(kk)])
        out['acceptance_rate'] = acceptance_plots
        layout.single_layout(rdir[base], acceptance_plots)

    if 'mcmc_history' in diagnostics:
        # files for samples mcmc history subsection
        base = "mcmc_history/{}".format(label)
        history_plots = []
        for kk, sf in enumerate(samples_file):
            history_plots += make_inference_plot_mcmc_history(
                workflow, sf, rdir[base],
                analysis_seg=workflow.analysis_time,
                tags=tags+[label, str(kk)])
        out['mcmc_history'] = history_plots
        layout.single_layout(rdir[base], history_plots)

    if 'dynesty_run' in diagnostics:
        # files for dynesty run subsection
        base = "dynesty_run/{}".format(label)
        dynesty_run_plots = []
        for kk, sf in enumerate(samples_file):
            dynesty_run_plots += make_inference_dynesty_run_plot(
                workflow, sf, rdir[base],
                analysis_seg=workflow.analysis_time,
                tags=tags+[label, str(kk)])
        out['dynesty_run'] = dynesty_run_plots
        layout.single_layout(rdir[base], dynesty_run_plots)

    if 'dynesty_traceplot' in diagnostics:
        # files for samples dynesty tyrace plots subsection
        base = "dynesty_traceplot/{}".format(label)
        dynesty_trace_plots = []
        for kk, sf in enumerate(samples_file):
            dynesty_trace_plots += make_inference_dynesty_trace_plot(
                workflow, sf, rdir[base],
                analysis_seg=workflow.analysis_time,
                tags=tags+[label, str(kk)])
        out['dynesty_traceplot'] = dynesty_trace_plots
        layout.single_layout(rdir[base], dynesty_trace_plots)

    return out


def make_posterior_workflow(workflow, samples_files, config_file, label,
                            rdir, posterior_file_dir='posterior_files',
                            tags=None):
    """Adds jobs to a workflow that make a posterior file and subsequent plots.

    A posterior file is first created from the given samples file(s). The
    settings for extracting the posterior are set by the
    ``[extract_posterior]`` section. If that section has a ``parameters``
    argument, then the parameters in the posterior file (and for use in all
    subsequent plotting) will be whatever that option is set to. Otherwise,
    the parameters in the posterior file will be whatever is common to
    all of the given samples file(s).

    Except for prior plots (which use the given inference config file), all
    subsequent jobs use the posterior file. The following are created:

    * **Summary table**: an html table created using the ``table_summary``
      executable. The parameters to print in the table are retrieved from the
      ``table-params`` option in the ``[workflow-summary_table]`` section.
      Metadata may also be printed by adding a ``print-metadata`` option to
      that section.
    * **Summary posterior plots**: a collection of posterior plots to include
      in the summary page, after the summary table. The parameters to plot
      are read from ``[workflow-summary_plots]``. Parameters should be grouped
      together by providing
      ``plot-group-NAME = PARAM1[:LABEL1] PARAM2[:LABEL2]`` in that section,
      where ``NAME`` is a unique name for each group. One posterior plot will
      be created for each plot group. For clarity, only one or two parameters
      should be plotted in each summary group, but this is not enforced.
      Settings for the plotting executable are read from the
      ``plot_posterior_summary`` section; likewise, the executable used
      is read from ``plot_posterior_summary`` in the
      ``[executables]`` section.
    * **Sky maps**: if *both* ``create_fits_file`` and ``plot_skymap``
      are listed in the ``[executables]`` section, then a ``.fits`` file and
      sky map plot will be produced. The sky map plot will be included in
      the summary plots. You must be running in a python 3 environment to
      create these.
    * **Prior plots**: plots of the prior will be created using the
      ``plot_prior`` executable. By default, all of the variable
      parameters will be plotted. The prior plots are added to
      ``priors/LALBEL/`` in the results directory, where ``LABEL`` is the
      given ``label``.
    * **Posterior plots**: additional posterior plots are created using the
      ``plot_posterior`` executable. The parameters to plot are
      read from ``[workflow-plot_params]`` section. As with the summary
      posterior plots, parameters are grouped together by providing
      ``plot-group-NAME`` options in that section. A posterior plot will be
      created for each group, and added to the ``posteriors/LABEL/`` directory.
      Plot settings are read from the ``[plot_posterior]`` section; this
      is kept separate from the posterior summary so that different settings
      can be used. For example, you may want to make a density plot for the
      summary plots, but a scatter plot colored by SNR for the posterior plots.


    Parameters
    ----------
    samples_file : pycbc.workflow.core.FileList
        List of samples files to combine into a single posterior file.
    config_file : pycbc.worfkow.File
        The inference configuration file used to generate the samples file(s).
        This is needed to make plots of the prior.
    label : str
        Unique label for the plots. Used in file names.
    rdir : pycbc.results.layout.SectionNumber
        The results directory to save the plots to.
    posterior_file_dir : str, optional
        The name of the directory to save the posterior file to. Default is
        ``posterior_files``.
    tags : list of str, optional
        Additional tags to add to the file names.

    Returns
    -------
    posterior_file : pycbc.workflow.File
        The posterior file that was created.
    summary_files : list
        List of files to go on the summary results page.
    prior_plots : list
        List of prior plots that will be created. These will be saved to
        ``priors/LABEL/`` in the resuls directory, where ``LABEL`` is the
        provided label.
    posterior_plots : list
        List of posterior plots that will be created. These will be saved to
        ``posteriors/LABEL/`` in the results directory.
    """
    # the list of plots to go in the summary
    summary_files = []

    if tags is None:
        tags = []

    analysis_seg = workflow.analysis_time

    # figure out what parameters user wants to plot from workflow configuration
    # parameters for the summary plots
    summary_plot_params = get_plot_group(workflow.cp, 'summary_plots')
    # parameters to plot in large corner plots
    plot_params = get_plot_group(workflow.cp, 'plot_params')
    # get parameters for the summary tables
    table_params = workflow.cp.get_opt_tag('workflow', 'table-params',
                                           'summary_table')
    # get any metadata that should be printed
    if workflow.cp.has_option('workflow-summary_table', 'print-metadata'):
        table_metadata = workflow.cp.get_opt_tag('workflow', 'print-metadata',
                                                 'summary_table')
    else:
        table_metadata = None

    # figure out if we are making a skymap
    make_skymap = ("create_fits_file" in workflow.cp.options("executables") and
                   "plot_skymap" in workflow.cp.options("executables"))

    make_prior = ("plot_prior" in workflow.cp.options("executables"))
    _config = None
    if make_prior:
        _config = config_file

    # make node for running extract samples
    posterior_file = create_posterior_files(
        workflow, samples_files, posterior_file_dir,
        analysis_seg=analysis_seg, tags=tags+[label])[0]

    # summary table
    summary_files += (make_inference_summary_table(
        workflow, posterior_file, rdir.base,
        parameters=table_params, print_metadata=table_metadata,
        analysis_seg=analysis_seg,
        tags=tags+[label]),)

    # summary posteriors
    summary_plots = []
    for group, params in summary_plot_params.items():
        summary_plots += make_inference_posterior_plot(
            workflow, posterior_file, rdir.base,
            name='plot_posterior_summary',
            parameters=params, plot_prior_from_file=_config,
            analysis_seg=analysis_seg,
            tags=tags+[label, group])

    # sky map
    if make_skymap:
        # create the fits file
        fits_file = create_fits_file(
            workflow, posterior_file, rdir.base, analysis_seg=analysis_seg,
            tags=tags+[label])[0]
        # now plot the skymap
        skymap_plot = make_inference_skymap(
            workflow, fits_file, rdir.base, analysis_seg=analysis_seg,
            tags=tags+[label])
        summary_plots += skymap_plot

    summary_files += list(layout.grouper(summary_plots, 2))

    # files for posteriors summary subsection
    base = "posteriors/{}".format(label)
    posterior_plots = []
    for group, params in plot_params.items():
        posterior_plots += make_inference_posterior_plot(
            workflow, posterior_file, rdir[base],
            parameters=params, plot_prior_from_file=_config,
            analysis_seg=analysis_seg,
            tags=tags+[label, group])
    layout.single_layout(rdir[base], posterior_plots)

    prior_plots = []
    # files for priors summary section
    if make_prior:
        base = "priors/{}".format(label)
        prior_plots += make_inference_prior_plot(
            workflow, config_file, rdir[base],
            analysis_seg=workflow.analysis_time, tags=tags+[label])
        layout.single_layout(rdir[base], prior_plots)
    return posterior_file, summary_files, prior_plots, posterior_plots


def _params_for_pegasus(parameters):
    """Escapes $ and escapes in parameters string for pegasus.

    Pegaus kickstart tries to do variable substitution if it sees a ``$``, and
    it will strip away back slashes. This can be problematic when trying to use
    LaTeX in parameter labels. This function adds escapes to all ``$`` and
    backslashes in a parameters argument, so the argument can be safely passed
    through pegasus-kickstart.

    Parameters
    ----------
    parameters : list or str
        The parameters argument to modify. If a list, the output will be
        converted to a space-separated string.
    """
    if isinstance(parameters, list):
        parameters = " ".join(parameters)
    return parameters.replace('\\', '\\\\').replace('$', '\$')
