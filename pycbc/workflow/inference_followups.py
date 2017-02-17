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
from pycbc.workflow.core import Executable, FileList, Node, makedir, File, Workflow
from pycbc.workflow.plotting import PlotExecutable, requirestr, excludestr
from pycbc.workflow import WorkflowConfigParser
from itertools import izip_longest
from Pegasus import DAX3 as dax
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
        logging.info("There is no [workflow-inference] section in configuration file")
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

    # get dax name and use it for the workflow name
    name = node.output_files[0].name
    node.add_opt("--workflow-name", name)

    # get output map name and use it for the output dir name
    map_loc = node.output_files[1].name
    node.add_opt("--output-dir", out_dir)

    # add this node to the workflow
    workflow += node

    # create job for dax that will run a sub-workflow
    # and add it to the workflow
    fil = node.output_files[0]
    job = dax.DAX(fil)
    job.addArguments("--basename %s" % os.path.splitext(os.path.basename(name))[0])
    Workflow.set_job_properties(job, map_loc)
    workflow._adag.addJob(job)

    # make dax a child of the inference workflow generator node
    dep = dax.Dependency(parent=node._dax_node, child=job)
    workflow._adag.addDependency(dep)

    logging.info("Leaving inference module")

def make_inference_prior_plot(workflow, config_file, output_dir,
                    sections=None, name="inference_prior",
                    analysis_seg=None, tags=None):
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
    analysis_segs: {None, glue.segments.Segment}
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
                      out_dir=output_dir, universe="local",
                      tags=tags).create_node()

    # add command line options
    node.add_input_opt("--config-file", config_file)
    node.new_output_file_opt(analysis_seg, ".png", "--output-file")
    if sections is not None:
        node.add_opt("--sections", " ".join(sections))

    # add node to workflow
    workflow += node

    return node.output_files

def make_inference_summary_table(workflow, inference_file, output_dir,
                    variable_args=None, name="inference_table",
                    analysis_seg=None, tags=None):
    """ Sets up the corner plot of the posteriors in the workflow.

    Parameters
    ----------
    workflow: pycbc.workflow.Workflow
        The core workflow instance we are populating
    inference_file: pycbc.workflow.File
        The file with posterior samples.
    output_dir: str
        The directory to store result plots and files.
    variable_args : list
        A list of parameters to use instead of [variable_args].
    name: str
        The name in the [executables] section of the configuration file
        to use.
    analysis_segs: {None, glue.segments.Segment}
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
    node.add_opt("--parameters", " ".join(variable_args))

    # add node to workflow
    workflow += node

    return node.output_files

def make_inference_posterior_plot(workflow, inference_file, output_dir,
                    variable_args=None,
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
    variable_args : list
        A list of parameters to use instead of [variable_args].
    name: str
        The name in the [executables] section of the configuration file
        to use.
    analysis_segs: {None, glue.segments.Segment}
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
                      out_dir=output_dir, universe="local",
                      tags=tags).create_node()

    # add command line options
    node.add_input_opt("--input-file", inference_file)
    node.new_output_file_opt(analysis_seg, ".png", "--output-file")
    node.add_opt("--parameters", " ".join(variable_args))

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
    analysis_segs: {None, glue.segments.Segment}
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

def make_inference_single_parameter_plots(workflow, inference_file, output_dir,
                    variable_args=None, samples_name="inference_samples",
                    acf_name="inference_acf", acl_name="inference_acl",
                    analysis_seg=None, tags=None):
    """ Sets up single-parameter plots for inference workflow.

    Parameters
    ----------
    workflow: pycbc.workflow.Workflow
        The core workflow instance we are populating
    inference_file: pycbc.workflow.File
        The file with posterior samples.
    output_dir: str
        The directory to store result plots and files.
    variable_args : list
        A list of parameters to use instead of [variable_args].
    samples_name: str
        The name in the [executables] section of the configuration file
        to use for the plot that shows all samples.
    acf_name: str
        The name in the [executables] section of the configuration file
        to use for the autocorrelation function plot.
    acl_name: str
        The name in the [executables] section of the configuration file
        to use for the autocorrelation length histogram.
    analysis_segs: {None, glue.segments.Segment}
       The segment this job encompasses. If None then use the total analysis
       time from the workflow.
    tags: {None, optional}
        Tags to add to the inference executables.

    Returns
    -------
    files: pycbc.workflow.FileList
        A list of result and output files. 
    """

    # default values
    tags = [] if tags is None else tags
    analysis_seg = workflow.analysis_time \
                       if analysis_seg is None else analysis_seg

    # make the directory that will contain the output files
    makedir(output_dir)

    # list of all output files
    files = FileList()

    # make a set of plots for each parameter
    for arg in variable_args:

        # plot posterior distribution
        corner_files = make_inference_posterior_plot(workflow, inference_file,
                          output_dir, variable_args=[arg],
                          analysis_seg=analysis_seg, tags=tags + [arg])

        # make a node for plotting all the samples for each walker
        samples_node = PlotExecutable(workflow.cp, samples_name,
                          ifos=workflow.ifos, out_dir=output_dir,
                          tags=tags + [arg]).create_node()
        samples_node.add_input_opt("--input-file", inference_file)
        samples_node.new_output_file_opt(analysis_seg, ".png", "--output-file")
        samples_node.add_opt("--parameters", arg)

        # make node for plotting the autocorrelation function for each walker
        acf_node = PlotExecutable(workflow.cp, acf_name, ifos=workflow.ifos,
                          out_dir=output_dir, tags=tags + [arg]).create_node()
        acf_node.add_input_opt("--input-file", inference_file)
        acf_node.new_output_file_opt(analysis_seg, ".png", "--output-file")
        acf_node.add_opt("--parameters", arg)

        # make node for plotting the autocorrelation function for each walker
        acl_node = PlotExecutable(workflow.cp, acl_name, ifos=workflow.ifos,
                          out_dir=output_dir, tags=tags + [arg]).create_node()
        acl_node.add_input_opt("--input-file", inference_file)
        acl_node.new_output_file_opt(analysis_seg, ".png", "--output-file")
        acl_node.add_opt("--parameters", arg)

        # add nodes to workflow
        workflow += samples_node
        workflow += acf_node
        workflow += acl_node

        # add files to output files list
        files += corner_files
        files += samples_node.output_files
        files += acf_node.output_files
        files += acl_node.output_files

    return files
