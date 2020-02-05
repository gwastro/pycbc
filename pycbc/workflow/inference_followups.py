# Copyright (C) 2016 Christopher M. Biwer, Alexander Harvey Nitz
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

import logging, os.path
from Pegasus import DAX3 as dax
from pycbc.workflow.core import (Executable, FileList, Node, makedir, File,
                                 Workflow)
from pycbc.workflow.plotting import PlotExecutable, requirestr, excludestr
from pycbc.workflow import WorkflowConfigParser
from pycbc.results import layout
from pycbc.workflow import pegasus_workflow as wdax

def setup_foreground_inference(workflow, coinc_file, single_triggers,
                       tmpltbank_file, insp_segs, insp_data_name,
                       insp_anal_name, dax_output, out_dir, tags=None):
    """ Creates workflow node that will run the inference workflow.

    Parameters
    ----------
    workflow: pycbc.workflow.Workflow
        The core workflow instance we are populating
    coinc_file: pycbc.workflow.File
        The file associated with coincident triggers.
    single_triggers: list of pycbc.workflow.File
        A list cointaining the file objects associated with the merged
        single detector trigger files for each ifo.
    tmpltbank_file: pycbc.workflow.File
        The file object pointing to the HDF format template bank
    insp_segs: SegFile
       The segment file containing the data read and analyzed by each inspiral
       job.
    insp_data_name: str
        The name of the segmentlist storing data read.
    insp_anal_name: str
        The name of the segmentlist storing data analyzed.
    dax_output : str
        The name of the output DAX file.
    out_dir: path
        The directory to store inference result plots and files
    tags: {None, optional}
        Tags to add to the inference executables
    """

    logging.info("Entering inference module")

    # check if configuration file has inference section
    if not workflow.cp.has_section("workflow-inference"):
        logging.info("There is no [workflow-inference] section in "
                     "configuration file")
        logging.info("Leaving inference module")
        return

    # default tags is a list
    tags = [] if tags is None else tags

    # make the directory that will contain the dax file
    makedir(dax_output)

    # turn the config file into a File class
    config_path = os.path.abspath(dax_output + "/" + "_".join(tags) \
                                        + "foreground_inference.ini")
    workflow.cp.write(open(config_path, "w"))
    config_file = wdax.File(os.path.basename(config_path))
    config_file.PFN(config_path, "local")

    # create an Executable for the inference workflow generator
    exe = Executable(workflow.cp, "foreground_inference", ifos=workflow.ifos,
                     out_dir=dax_output)

    # create the node that will run in the workflow
    node = exe.create_node()
    node.add_input_opt("--config-files", config_file)
    node.add_input_opt("--bank-file", tmpltbank_file)
    node.add_input_opt("--statmap-file", coinc_file)
    node.add_multiifo_input_list_opt("--single-detector-triggers",
                                     single_triggers)
    node.new_output_file_opt(workflow.analysis_time, ".dax", "--output-file",
                                     tags=tags)
    node.new_output_file_opt(workflow.analysis_time, ".dax.map",
                                     "--output-map", tags=tags)
    node.new_output_file_opt(workflow.analysis_time, ".tc.txt",
                                     "--transformation-catalog", tags=tags)

    # get dax name and use it for the workflow name
    name = node.output_files[0].name
    node.add_opt("--workflow-name", name)

    # get output map name and use it for the output dir name
    map_file = node.output_files[1]
    node.add_opt("--output-dir", out_dir)

    # get the transformation catalog name
    tc_file = node.output_files[2]

    # add this node to the workflow
    workflow += node

    # create job for dax that will run a sub-workflow
    # and add it to the workflow
    fil = node.output_files[0]
    job = dax.DAX(fil)
    job.addArguments("--basename {}".format(
        os.path.splitext(os.path.basename(name))[0]))
    Workflow.set_job_properties(job, map_file, tc_file)
    workflow._adag.addJob(job)

    # make dax a child of the inference workflow generator node
    dep = dax.Dependency(parent=node._dax_node, child=job)
    workflow._adag.addDependency(dep)

    logging.info("Leaving inference module")

def make_inference_prior_plot(workflow, config_file, output_dir,
                              sections=None, name="inference_prior",
                              parameters=None, analysis_seg=None, tags=None):
    """ Sets up the corner plot of the priors in the workflow.

    Parameters
    ----------
    workflow: pycbc.workflow.Workflow
        The core workflow instance we are populating
    config_file: pycbc.workflow.File
        The WorkflowConfigParser parasable inference configuration file..
    output_dir: str
        The directory to store result plots and files.
    sections : list
        A list of subsections to use.
    name: str
        The name in the [executables] section of the configuration file
        to use.
    analysis_segs: {None, ligo.segments.Segment}
       The segment this job encompasses. If None then use the total analysis
       time from the workflow.
    tags: {None, optional}
        Tags to add to the inference executables.

    Returns
    -------
    pycbc.workflow.FileList
        A list of result and output files.
    """

    # default values
    tags = [] if tags is None else tags
    analysis_seg = workflow.analysis_time \
                       if analysis_seg is None else analysis_seg

    # make the directory that will contain the output files
    makedir(output_dir)

    # make a node for plotting the posterior as a corner plot
    node = PlotExecutable(workflow.cp, name, ifos=workflow.ifos,
                      out_dir=output_dir,
                      tags=tags).create_node()

    # add command line options
    node.add_input_opt("--config-file", config_file)
    node.new_output_file_opt(analysis_seg, ".png", "--output-file")
    if sections is not None:
        node.add_opt("--sections", " ".join(sections))
    if parameters is not None:
        node.add_opt("--parameters", _params_for_pegasus(parameters))

    # add node to workflow
    workflow += node

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
    parameters : list, optional
        A list of the parameters to extract, and (optionally) a name for them
        to be mapped to. This is passed to the program's ``--parameters``
        argument.
    name: str, optional
        The name in the [executables] section of the configuration file
        to use.
    analysis_segs: {None, ligo.segments.Segment}
       The segment this job encompasses. If None then use the total analysis
       time from the workflow.
    tags: list, optional
        Tags to add to the inference executables.

    Returns
    -------
    pycbc.workflow.FileList
        A list of result and output files.
    """
    if analysis_seg is None:
        analysis_seg = workflow.analysis_time
    if tags is None:
        tags = []
    extract_posterior_exe = Executable(workflow.cp, name,
                                       ifos=workflow.ifos,
                                       out_dir=output_dir)
    node = extract_posterior_exe.create_node()
    if not isinstance(samples_files, list):
        samples_files = [samples_files]
    node.add_input_list_opt("--input-file", samples_files)
    if parameters is not None:
        node.add_opt("--parameters", _params_for_pegasus(parameters))
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
    name: str
        The name in the [executables] section of the configuration file
        to use.
    analysis_segs: {None, ligo.segments.Segment}
       The segment this job encompasses. If None then use the total analysis
       time from the workflow.
    tags: list, optional
        Tags to add to the inference executables.

    Returns
    -------
    pycbc.workflow.FileList
        A list of result and output files.
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
                          name="inference_skymap", analysis_seg=None,
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
    name: str
        The name in the [executables] section of the configuration file
        to use.
    analysis_segs: {None, ligo.segments.Segment}
       The segment this job encompasses. If None then use the total analysis
       time from the workflow.
    tags: {None, optional}
        Tags to add to the inference executables.

    Returns
    -------
    pycbc.workflow.FileList
        A list of result and output files.
    """
    # default values
    tags = [] if tags is None else tags
    analysis_seg = workflow.analysis_time \
                       if analysis_seg is None else analysis_seg
    # make the directory that will contain the output files
    makedir(output_dir)
    # make a node for plotting the posterior as a corner plot
    node = PlotExecutable(workflow.cp, name, ifos=workflow.ifos,
                          out_dir=output_dir,
                          tags=tags).create_node()
    # add command line options
    node.add_input_opt("--input-file", fits_file)
    node.new_output_file_opt(analysis_seg, ".png", "--output-file")
    # add node to workflow
    workflow += node
    return node.output_files


def make_inference_summary_table(workflow, inference_file, output_dir,
                    parameters=None, print_metadata=None,
                    name="inference_table",
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
    parameters : list
        A list of parameters to generate the table for.
    name: str
        The name in the [executables] section of the configuration file
        to use.
    analysis_segs: {None, ligo.segments.Segment}
       The segment this job encompasses. If None then use the total analysis
       time from the workflow.
    tags: {None, optional}
        Tags to add to the inference executables.

    Returns
    -------
    pycbc.workflow.FileList
        A list of result and output files.
    """

    # default values
    tags = [] if tags is None else tags
    analysis_seg = workflow.analysis_time \
                       if analysis_seg is None else analysis_seg

    # make the directory that will contain the output files
    makedir(output_dir)

    # make a node for plotting the posterior as a corner plot
    node = PlotExecutable(workflow.cp, name, ifos=workflow.ifos,
                      out_dir=output_dir, tags=tags).create_node()

    # add command line options
    node.add_input_opt("--input-file", inference_file)
    node.new_output_file_opt(analysis_seg, ".html", "--output-file")
    if parameters is not None:
        node.add_opt("--parameters", _params_for_pegasus(parameters))
    if print_metadata is not None:
        node.add_opt("--print-metadata", _params_for_pegasus(print_metadata))

    # add node to workflow
    workflow += node

    return node.output_files

def make_inference_posterior_plot(
                    workflow, inference_file, output_dir, parameters=None,
                    plot_prior_from_file=None,
                    name="inference_posterior", analysis_seg=None, tags=None):
    """ Sets up the corner plot of the posteriors in the workflow.

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
    name: str
        The name in the [executables] section of the configuration file
        to use.
    analysis_segs: {None, ligo.segments.Segment}
       The segment this job encompasses. If None then use the total analysis
       time from the workflow.
    tags: {None, optional}
        Tags to add to the inference executables.

    Returns
    -------
    pycbc.workflow.FileList
        A list of result and output files.
    """

    # default values
    tags = [] if tags is None else tags
    analysis_seg = workflow.analysis_time \
                       if analysis_seg is None else analysis_seg

    # make the directory that will contain the output files
    makedir(output_dir)

    # make a node for plotting the posterior as a corner plot
    node = PlotExecutable(workflow.cp, name, ifos=workflow.ifos,
                      out_dir=output_dir,
                      tags=tags).create_node()

    # add command line options
    node.add_input_opt("--input-file", inference_file)
    node.new_output_file_opt(analysis_seg, ".png", "--output-file")
    if parameters is not None:
        node.add_opt("--parameters", _params_for_pegasus(parameters))
    if plot_prior_from_file is not None:
        node.add_input_opt('--plot-prior', plot_prior_from_file)

    # add node to workflow
    workflow += node

    return node.output_files

def make_inference_1d_posterior_plots(
                    workflow, inference_file, output_dir, parameters=None,
                    analysis_seg=None, tags=None):
    parameters = [] if parameters is None else parameters
    files = FileList([])
    for (ii, parameter) in enumerate(parameters):
        files += make_inference_posterior_plot(
                    workflow, inference_file, output_dir,
                    parameters=[parameter], analysis_seg=analysis_seg,
                    tags=tags+['param{}'.format(ii)])
    return files

def make_inference_samples_plot(
                    workflow, inference_file, output_dir, parameters=None,
                    name="inference_samples", analysis_seg=None, tags=None):
    # default values
    tags = [] if tags is None else tags
    analysis_seg = workflow.analysis_time \
                       if analysis_seg is None else analysis_seg

    # make the directory that will contain the output files
    makedir(output_dir)

    # make a node for plotting the posterior as a corner plot
    node = PlotExecutable(workflow.cp, name, ifos=workflow.ifos,
                      out_dir=output_dir,
                      tags=tags).create_node()

    # add command line options
    node.add_input_opt("--input-file", inference_file)
    node.new_output_file_opt(analysis_seg, ".png", "--output-file")
    if parameters is not None:
        node.add_opt("--parameters", _params_for_pegasus(parameters))

    # add node to workflow
    workflow += node

    return node.output_files


def make_inference_acceptance_rate_plot(workflow, inference_file, output_dir,
                    name="inference_rate", analysis_seg=None, tags=None):
    """ Sets up the acceptance rate plot in the workflow.

    Parameters
    ----------
    workflow: pycbc.workflow.Workflow
        The core workflow instance we are populating
    inference_file: pycbc.workflow.File
        The file with posterior samples.
    output_dir: str
        The directory to store result plots and files.
    name: str
        The name in the [executables] section of the configuration file
        to use.
    analysis_segs: {None, ligo.segments.Segment}
       The segment this job encompasses. If None then use the total analysis
       time from the workflow.
    tags: {None, optional}
        Tags to add to the inference executables.

    Returns
    -------
    pycbc.workflow.FileList
        A list of result and output files.
    """

    # default values
    tags = [] if tags is None else tags
    analysis_seg = workflow.analysis_time \
                       if analysis_seg is None else analysis_seg

    # make the directory that will contain the output files
    makedir(output_dir)

    # make a node for plotting the acceptance rate
    node = PlotExecutable(workflow.cp, name, ifos=workflow.ifos,
                      out_dir=output_dir, tags=tags).create_node()

    # add command line options
    node.add_input_opt("--input-file", inference_file)
    node.new_output_file_opt(analysis_seg, ".png", "--output-file")

    # add node to workflow
    workflow += node

    return node.output_files

def make_inference_inj_plots(workflow, inference_files, output_dir,
                             parameters, name="inference_recovery",
                             analysis_seg=None, tags=None):
    """ Sets up the recovered versus injected parameter plot in the workflow.

    Parameters
    ----------
    workflow: pycbc.workflow.Workflow
        The core workflow instance we are populating
    inference_files: pycbc.workflow.FileList
        The files with posterior samples.
    output_dir: str
        The directory to store result plots and files.
    parameters : list
        A ``list`` of parameters. Each parameter gets its own plot.
    name: str
        The name in the [executables] section of the configuration file
        to use.
    analysis_segs: {None, ligo.segments.Segment}
       The segment this job encompasses. If None then use the total analysis
       time from the workflow.
    tags: {None, optional}
        Tags to add to the inference executables.

    Returns
    -------
    pycbc.workflow.FileList
        A list of result and output files.
    """

    # default values
    tags = [] if tags is None else tags
    analysis_seg = workflow.analysis_time \
                       if analysis_seg is None else analysis_seg
    output_files = FileList([])

    # make the directory that will contain the output files
    makedir(output_dir)

    # add command line options
    for (ii, param) in enumerate(parameters):
        plot_exe = PlotExecutable(workflow.cp, name, ifos=workflow.ifos,
                                  out_dir=output_dir,
                                  tags=tags+['param{}'.format(ii)])
        node = plot_exe.create_node()
        node.add_input_list_opt("--input-file", inference_files)
        node.new_output_file_opt(analysis_seg, ".png", "--output-file")
        node.add_opt("--parameters", param)
        workflow += node
        output_files += node.output_files

    return output_files


def get_posterior_params(cp, section='workflow-posterior_params'):
    """Gets the posterior parameters from the given config file.

    The posterior parameters are read from the given ``section``. Parameters
    should be specified as ``OUTPUT = [INPUT]``, where ``OUTPUT`` is what
    the parameter should be named in the posterior file and ``INPUT`` is the
    (function of) parameter(s) to read from the samples file. If no ``INPUT``
    is provided, the ``INPUT`` name will assumed to be the same as the
    ``OUTPUT``. Example:

    .. code-block:: ini

       [workflow-posterior_params]
       mass1 = primary_mass(mass1, mass2)
       mass2 = secondary_mass(mass1, mass2)
       distance =

    Parameters
    ----------
    cp : pycbc.workflow.configuration.WorkflowConfigParser
        Config parser to read.
    section : str, optional
        The name of the section to load the parameters from. Default is
        ``workflow-posterior_params``.

    Returns
    -------
    list :
        List of strings giving ``INPUT:OUTPUT``. This can be passed as the
        ``parameters`` argument to :py:func:`create_posterior_files`.
    """
    params = []
    for opt in cp.options(section):
        val = cp.get(section, opt)
        if val == '':
            val = opt
        params.append('{}:{}'.format(val, opt))
    return params


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
      of iteration. This will be created if ``inference_samples`` is in the
      executables section.
    * ``acceptance_rate``: For MCMC samplers, a plot of the acceptance rate.
      This will be created if ``inference_rate`` is in the executables section.

    Returns
    -------
    list :
        List of names of diagnostic plots.
    """
    diagnostics = []
    if "inference_samples" in workflow.cp.options("executables"):
        diagnostics.append('samples')
    if "inference_rate" in workflow.cp.options("executables"):
        diagnostics.append('acceptance_rate')
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

    return out


def make_posterior_workflow(workflow, samples_files, config_file, label,
                            rdir, posterior_file_dir='posterior_files',
                            tags=None):
    """Adds jobs to a workflow that make a posterior file and subsequent plots.

    The parameters to be written to the posterior file are read from the
    ``[workflow-posterior_params]`` section of the workflow's config file; see
    :py:func:`get_posterior_params` for details.

    Except for prior plots (which use the given inference config file), all
    subsequent jobs use the posterior file, and so may use the parameters
    provided in ``[workflow-posterior_params]``. The following are created:

    * **Summary table**: an html table created using the ``inference_table``
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
      ``inference_posterior_summary`` section; likewise, the executable used
      is read from ``inference_posterior_summary`` in the
      ``[executables]`` section.
    * **Sky maps**: if *both* ``create_fits_file`` and ``inference_skymap``
      are listed in the ``[executables]`` section, then a ``.fits`` file and
      sky map plot will be produced. The sky map plot will be included in
      the summary plots. You must be running in a python 3 environment to
      create these.
    * **Prior plots**: plots of the prior will be created using the
      ``inference_prior`` executable. By default, all of the variable
      parameters will be plotted. The prior plots are added to
      ``priors/LALBEL/`` in the results directory, where ``LABEL`` is the
      given ``label``.
    * **Posterior plots**: additional posterior plots are created using the
      ``inference_posterior`` executable. The parameters to plot are
      read from ``[workflow-plot_params]`` section. As with the summary
      posterior plots, parameters are grouped together by providing
      ``plot-group-NAME`` options in that section. A posterior plot will be
      created for each group, and added to the ``posteriors/LABEL/`` directory.
      Plot settings are read from the ``[inference_posterior]`` section; this
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

    # get the parameters that will be used for posterior files and plots
    posterior_params = get_posterior_params(workflow.cp)

    # figure out what parameters user wants to plot from workflow configuration
    group_prefix = "plot-group-"
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
                   "inference_skymap" in workflow.cp.options("executables"))

    # make node for running extract samples
    posterior_file = create_posterior_files(
        workflow, samples_files, posterior_file_dir,
        parameters=posterior_params, analysis_seg=analysis_seg,
        tags=tags+[label])[0]

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
            name='inference_posterior_summary',
            parameters=params, plot_prior_from_file=config_file,
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

    # files for priors summary section
    base = "priors/{}".format(label)
    prior_plots = make_inference_prior_plot(
        workflow, config_file, rdir[base],
        analysis_seg=workflow.analysis_time,
        tags=tags+[label])
    layout.single_layout(rdir[base], prior_plots)

    # files for posteriors summary subsection
    base = "posteriors/{}".format(label)
    posterior_plots = []
    for group, params in plot_params.items():
        posterior_plots += make_inference_posterior_plot(
            workflow, posterior_file, rdir[base],
            parameters=params, plot_prior_from_file=config_file,
            analysis_seg=analysis_seg,
            tags=tags+[label, group])
    layout.single_layout(rdir[base], posterior_plots)

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
