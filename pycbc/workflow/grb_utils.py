# Copyright (C) 2015  Andrew Williamson, Francesco Pannarale
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

#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#

"""
This library code contains functions and classes that are used in the
generation of pygrb workflows. For details about pycbc.workflow see here:
http://pycbc.org/pycbc/latest/html/workflow.html
"""

import logging
import os

import numpy as np
from gwdatafind.utils import filename_metadata
from scipy.stats import rayleigh

from pycbc import makedir
from pycbc.results import layout
from pycbc.workflow.core import Executable, File, FileList, Node, resolve_url_to_file
from pycbc.workflow.jobsetup import select_generic_executable
from pycbc.workflow.minifollowups import add_wellfile_opt, write_config_for_subworkflow
from pycbc.workflow.pegasus_workflow import SubWorkflow
from pycbc.workflow.plotting import PlotExecutable
from pycbc.workflow.versioning import make_versioning_page

logger = logging.getLogger(__name__)


def _select_grb_pp_class(wflow, curr_exe):
    """
    This function returns the class for PyGRB post-processing scripts.

    Parameters
    ----------
    curr_exe : string
        The name of the executable

    Returns
    -------
    exe_class : Sub-class of pycbc.workflow.core.Executable that holds utility
        functions appropriate for the given executable.  Instances of the class
        ('jobs') **must** have methods
        * job.create_node()
        and
        * job.get_valid_times(ifo, )

    """
    exe_path = wflow.cp.get("executables", curr_exe)
    exe_name = os.path.basename(exe_path)
    exe_to_class_map = {
        "pycbc_grb_trig_combiner": PycbcGrbTrigCombinerExecutable,
        "pycbc_grb_trig_cluster": PycbcGrbTrigClusterExecutable,
        "pycbc_grb_inj_finder": PycbcGrbInjFinderExecutable,
    }
    if exe_name not in exe_to_class_map:
        raise ValueError(f"No job class exists for executable {curr_exe}")

    return exe_to_class_map[exe_name]


def set_grb_start_end(cp, start, end):
    """
    Function to update analysis boundaries as workflow is generated

    Parameters
    ----------
    cp : pycbc.workflow.configuration.WorkflowConfigParser object
    The parsed configuration options of a pycbc.workflow.core.Workflow.

    start : int
    The start of the workflow analysis time.

    end : int
    The end of the workflow analysis time.

    Returns
    -------
    cp : pycbc.workflow.configuration.WorkflowConfigParser object
    The modified WorkflowConfigParser object.

    """
    cp.set("workflow", "start-time", str(start))
    cp.set("workflow", "end-time", str(end))

    return cp


def make_gating_node(workflow, datafind_files, outdir=None, tags=None):
    """
    Generate jobs for autogating the data for PyGRB runs.

    Parameters
    ----------
    workflow: pycbc.workflow.core.Workflow
        An instanced class that manages the constructed workflow.
    datafind_files : pycbc.workflow.core.FileList
        A FileList containing the frame files to be gated.
    outdir : string
        Path of the output directory
    tags : list of strings
        If given these tags are used to uniquely name and identify output files
        that would be produced in multiple calls to this function.

    Returns
    -------
    condition_strain_nodes : list
        List containing the pycbc.workflow.core.Node objects representing the
        autogating jobs.
    condition_strain_outs : pycbc.workflow.core.FileList
        FileList containing the pycbc.workflow.core.File objects representing
        the gated frame files.

    """
    cp = workflow.cp
    if tags is None:
        tags = []

    condition_strain_class = select_generic_executable(workflow, "condition_strain")
    condition_strain_nodes = []
    condition_strain_outs = FileList([])
    for ifo in workflow.ifos:
        input_files = FileList(
            [
                datafind_file
                for datafind_file in datafind_files
                if datafind_file.ifo == ifo
            ]
        )
        condition_strain_jobs = condition_strain_class(
            cp, "condition_strain", ifos=ifo, out_dir=outdir, tags=tags
        )
        condition_strain_node, condition_strain_out = condition_strain_jobs.create_node(
            input_files, tags=tags
        )
        condition_strain_nodes.append(condition_strain_node)
        condition_strain_outs.extend(FileList([condition_strain_out]))

    return condition_strain_nodes, condition_strain_outs


def fermi_core_tail_model(
    sky_err, rad, core_frac=0.98, core_sigma=3.6, tail_sigma=29.6
):
    """
    Fermi systematic error model following
    https://arxiv.org/abs/1909.03006, with default values valid
    before 11 September 2019.

    Parameters
    ----------
    core_frac : float
        Fraction of the systematic uncertainty contained within the core
        component.
    core_sigma : float
        Size of the GBM systematic core component.
    tail_sigma : float
        Size of the GBM systematic tail component.

    Returns
    -------
    _______
    tuple
        Tuple containing the core and tail probability distributions
        as a function of radius.

    """
    scaledsq = sky_err**2 / -2 / np.log(0.32)
    return (
        frac * (1 - np.exp(-0.5 * (rad / np.sqrt(scaledsq + sigma**2)) ** 2))
        for frac, sigma in zip(
            [core_frac, 1 - core_frac], [core_sigma, tail_sigma], strict=True
        )
    )


def get_sky_grid_scale(
    sky_error=0.0,
    containment=0.9,
    upscale=False,
    fermi_sys=False,
    precision=1e-3,
    **kwargs,
):
    """
    Calculate the angular radius corresponding to a desired
    localization uncertainty level. This is used to generate the search
    grid and involves scaling up the standard 1-sigma value provided to
    the workflow, assuming a normal probability profile. Fermi
    systematic errors can be included, following
    https://arxiv.org/abs/1909.03006, with default values valid before
    11 September 2019. The default probability coverage is 90%.

    Parameters
    ----------
    sky_error : float
        The reported statistical 1-sigma sky error of the trigger.
    containment : float
        The desired localization probability to be covered by the sky
        grid.
    upscale : bool, optional
        Whether to apply rescale to convert from 1 sigma -> containment
        for non-Fermi triggers. Default = True as Swift reports 90%
        radius directly.
    fermi_sys : bool, optional
        Whether to apply Fermi-GBM systematics via
        ``fermi_core_tail_model``. Default = False.
    precision : float, optional
        Precision (in degrees) for calculating the error radius via
        Fermi-GBM model.
    **kwargs
        Additional keyword arguments passed to `fermi_core_tail_model`.

    Returns
    -------
    _______

    float
        Sky error radius in degrees.

    """
    if fermi_sys:
        lims = (0.5, 4)
        radii = np.linspace(
            lims[0] * sky_error,
            lims[1] * sky_error,
            int((lims[1] - lims[0]) * sky_error / precision) + 1,
        )
        core, tail = fermi_core_tail_model(sky_error, radii, **kwargs)
        out = radii[(abs(core + tail - containment)).argmin()]
    else:
        # Use Rayleigh distribution to go from 1 sigma containment to
        # containment given by function variable. Interval method returns
        # bounds of equal probability about the median, but we want 1-sided
        # bound, hence use (2 * containment - 1)
        out = sky_error
        if upscale:
            out *= rayleigh.interval(2 * containment - 1)[-1]
    return out


def make_skygrid_node(workflow, out_dir, tags=None):
    """
    Adds a job to the workflow to produce the PyGRB search skygrid."""
    tags = [] if tags is None else tags

    # Initialize job node
    grb_name = workflow.cp.get("workflow", "trigger-name")
    extra_tags = ["GRB" + grb_name]
    node = Executable(
        workflow.cp,
        "make_sky_grid",
        ifos=workflow.ifos,
        out_dir=out_dir,
        tags=tags + extra_tags,
    ).create_node()
    node.add_opt("--instruments", " ".join(workflow.ifos))
    node.new_output_file_opt(
        workflow.analysis_time, ".h5", "--output", tags=extra_tags, store_file=True
    )

    # Add job node to the workflow
    workflow += node

    return node.output_files


def generate_tc_prior(wflow, tc_path, buffer_seg):
    """
    Generate the configuration file for the prior on the coalescence
    time of injections, ensuring that these times fall in the analysis
    time and avoid the onsource and its buffer.

    Parameters
    ----------
    tc_path : str
        Path where the configuration file for the prior needs to be written.
    buffer_seg : segmentlist
        Start and end times of the buffer segment encapsulating the onsource.

    """
    # Write the tc-prior configuration file if it does not exist
    if os.path.exists(tc_path):
        raise ValueError(f"Refusing to overwrite {tc_path}.")
    tc_file = open(tc_path, "w")
    tc_file.write("[prior-tc]\n")
    tc_file.write("name = uniform\n")
    tc_file.write(f"min-tc = {wflow.analysis_time[0]}\n")
    tc_file.write(f"max-tc = {wflow.analysis_time[1]}\n\n")
    tc_file.write("[constraint-tc]\n")
    tc_file.write("name = custom\n")
    tc_file.write(f"constraint_arg = (tc < {buffer_seg[0]}) | (tc > {buffer_seg[1]})\n")
    tc_file.close()

    # Add the tc-prior configuration file url to wflow.cp if necessary
    tc_file_path = "file://" + tc_path
    for inj_sec in wflow.cp.get_subsections("injections"):
        config_urls = wflow.cp.get("workflow-injections", inj_sec + "-config-files")
        config_urls = [url.strip() for url in config_urls.split(",")]
        if tc_file_path not in config_urls:
            config_urls += [tc_file_path]
        config_urls = ", ".join([str(item) for item in config_urls])
        wflow.cp.set("workflow-injections", inj_sec + "-config-files", config_urls)


def setup_pygrb_pp_workflow(
    wf,
    pp_dir,
    seg_files,
    segment,
    bank_file,
    insp_files,
    inj_files,
    inj_insp_files,
    inj_tags,
):
    """
    Generate post-processing section of PyGRB offline workflow

    Parameters
    ----------
    wf : The workflow object
    pp_dir : The directory where the post-processing files will be stored
    seg_files : FileList of the segment files
    segment : The segment to be analyzed
    bank_file : The full template bank file
    insp_files : The list of inspiral files
    inj_files : The list of injection files
    inj_insp_files : The list of inspiral files for injections
    inj_tags : The list of injection tags

    Returns
    -------
    trig_files : FileList
        The list of combined trigger files
        [ALL_TIMES, ONSOURCE, OFFSOURCE, OFFTRIAL_1, ..., OFFTRIAL_N]
        FileList (N can be set by the user and is 6 by default)
    clustered_files : FileList
        CLUSTERED FileList, same order as trig_files
        Contains triggers after clustering
    inj_find_files : FileList
        FOUNDMISSED FileList covering all injection sets

    """
    # Begin setting up trig combiner job(s)
    # Select executable class and initialize
    exe_class = _select_grb_pp_class(wf, "trig_combiner")
    job_instance = exe_class(wf.cp, "trig_combiner")
    # Create node for coherent no injections jobs
    node, trig_files = job_instance.create_node(
        wf.ifo_string, seg_files, segment, insp_files, pp_dir, bank_file
    )
    wf.add_node(node)

    # Trig clustering for each trig file
    exe_class = _select_grb_pp_class(wf, "trig_cluster")
    job_instance = exe_class(wf.cp, "trig_cluster")
    clustered_files = FileList([])
    for trig_file in trig_files:
        # Create and add nodes
        node, out_file = job_instance.create_node(trig_file, pp_dir)
        wf.add_node(node)
        clustered_files.append(out_file)

    # Find injections from triggers
    exe_class = _select_grb_pp_class(wf, "inj_finder")
    job_instance = exe_class(wf.cp, "inj_finder")
    inj_find_files = FileList([])
    for inj_tag in inj_tags:
        tag_inj_files = FileList([f for f in inj_files if inj_tag in f.tags])
        # The here stems from the injection group information
        # being stored in the second tag. This could be improved
        # depending on the final implementation of injections
        tag_insp_files = FileList([f for f in inj_insp_files if inj_tag in f.tags[1]])
        node, inj_find_file = job_instance.create_node(
            tag_inj_files, tag_insp_files, bank_file, pp_dir, tags=[inj_tag]
        )
        wf.add_node(node)
        inj_find_files.append(inj_find_file)

    return trig_files, clustered_files, inj_find_files


class PycbcGrbTrigCombinerExecutable(Executable):
    """
    The class responsible for creating jobs
    for ''pycbc_grb_trig_combiner''.
    """

    current_retention_level = Executable.ALL_TRIGGERS

    def __init__(self, cp, name):
        super().__init__(cp=cp, name=name)
        self.trigger_name = cp.get("workflow", "trigger-name")
        self.trig_start_time = cp.get("workflow", "start-time")
        self.num_trials = int(cp.get("trig_combiner", "num-trials"))

    def create_node(
        self, ifo_tag, seg_files, segment, insp_files, out_dir, bank_file, tags=None
    ):
        node = Node(self)
        node.add_opt("--verbose")
        node.add_opt("--ifo-tag", ifo_tag)
        node.add_opt("--grb-name", self.trigger_name)
        node.add_opt("--trig-start-time", self.trig_start_time)
        node.add_input_list_opt("--seg-files", seg_files)
        node.add_input_list_opt("--input-files", insp_files)
        node.add_opt("--user-tag", "PYGRB")
        node.add_input_opt("--bank-file", bank_file)
        # Prepare output file tag
        user_tag = f"PYGRB_GRB{self.trigger_name}"
        if tags:
            user_tag += f"_{tags}"
        # Add on/off source and off trial outputs
        output_files = FileList([])
        outfile_types = ["ALL_TIMES", "ONSOURCE", "OFFSOURCE"]
        for i in range(self.num_trials):
            outfile_types.append(f"OFFTRIAL_{i + 1}")
        for out_type in outfile_types:
            out_name = f"{ifo_tag}-{user_tag}_{out_type}-{segment[0]}-{segment[1] - segment[0]}.h5"
            out_file = File(
                ifo_tag,
                "trig_combiner",
                segment,
                file_url=os.path.join(out_dir, out_name),
            )
            node.add_output(out_file)
            output_files.append(out_file)

        return node, output_files


class PycbcGrbTrigClusterExecutable(Executable):
    """
    The class responsible for creating jobs
    for ''pycbc_grb_trig_cluster''.
    """

    current_retention_level = Executable.ALL_TRIGGERS

    def __init__(self, cp, name):
        super().__init__(cp=cp, name=name)

    def create_node(self, in_file, out_dir):
        node = Node(self)
        node.add_input_opt("--trig-file", in_file)
        # Determine output file name
        ifotag, filetag, segment = filename_metadata(in_file.name)
        start, end = segment
        out_name = f"{ifotag}-{filetag}_CLUSTERED-{start}-{end - start}.h5"
        out_file = File(
            ifotag, "trig_cluster", segment, file_url=os.path.join(out_dir, out_name)
        )
        node.add_output(out_file)

        return node, out_file


class PycbcGrbInjFinderExecutable(Executable):
    """The class responsible for creating jobs for ``pycbc_grb_inj_finder``"""

    current_retention_level = Executable.ALL_TRIGGERS

    def __init__(self, cp, exe_name):
        super().__init__(cp=cp, name=exe_name)

    def create_node(self, inj_files, inj_insp_files, bank_file, out_dir, tags=None):
        if tags is None:
            tags = []
        node = Node(self)
        node.add_input_list_opt("--input-files", inj_insp_files)
        node.add_input_list_opt("--inj-files", inj_files)
        node.add_input_opt("--bank-file", bank_file)
        ifo_tag, desc, segment = filename_metadata(inj_files[0].name)
        desc = "_".join(desc.split("_")[:-1])
        out_name = f"{ifo_tag}-{desc}_FOUNDMISSED-{segment[0]}-{abs(segment)}.h5"
        out_file = File(
            ifo_tag, "inj_finder", segment, os.path.join(out_dir, out_name), tags=tags
        )
        node.add_output(out_file)
        return node, out_file


def build_segment_filelist(seg_dir):
    """Construct a FileList instance containing all segments txt files"""
    # Needs to be in this order for consistency with _read_seg_files
    file_names = ["bufferSeg.txt", "offSourceSeg.txt", "onSourceSeg.txt"]
    seg_files = FileList([])
    for file_name in file_names:
        path = os.path.join(seg_dir, file_name)
        seg_file = File.from_path(path)
        seg_file.add_pfn(path, site="local")
        seg_files.append(seg_file)

    return seg_files


def make_pygrb_plot(
    workflow,
    exec_name,
    out_dir,
    ifo=None,
    inj_file=None,
    trig_file=None,
    onsource_file=None,
    bank_file=None,
    seg_files=None,
    sky_grid_file=None,
    veto_file=None,
    tags=None,
    **kwargs,
):
    """Adds a node for a plot of PyGRB results to the workflow"""
    tags = [] if tags is None else tags

    # Initialize job node with its tags
    grb_name = workflow.cp.get("workflow", "trigger-name")
    extra_tags = ["GRB" + grb_name]
    if ifo:
        extra_tags.append(ifo)
    node = PlotExecutable(
        workflow.cp,
        exec_name,
        ifos=workflow.ifos,
        out_dir=out_dir,
        tags=tags + extra_tags,
    ).create_node()
    if trig_file:
        node.add_input_opt("--trig-file", trig_file)
    # Pass the veto and segment files and options
    if seg_files:
        node.add_input_list_opt("--seg-files", seg_files)
    if sky_grid_file:
        node.add_input_opt("--sky-grid", sky_grid_file)
    if veto_file:
        node.add_input_opt("--veto-file", veto_file)
    # Option to show the onsource trial if this is a plot of all data
    if exec_name == "pygrb_plot_snr_timeseries" and "alltimes" in tags:
        node.add_opt("--onsource")
    if exec_name in ["pygrb_plot_injs_results", "pygrb_plot_snr_timeseries"]:
        trig_time = workflow.cp.get("workflow", "trigger-time")
        node.add_opt("--trigger-time", trig_time)
    # Pass the injection file as an input File instance
    if inj_file is not None and exec_name not in [
        "pygrb_plot_skygrid",
        "pygrb_plot_stats_distribution",
        "pycbc_plot_bank_corner",
    ]:
        node.add_input_opt("--found-missed-file", inj_file)
    if exec_name == "pycbc_plot_bank_corner":
        if inj_file is not None:
            node.add_input_opt("--bank-file", inj_file)
        else:
            node.add_input_opt("--bank-file", bank_file)
    # IFO option
    if ifo:
        node.add_opt("--ifo", ifo)
    # Output files and final input file (passed as a File instance)
    if exec_name == "pygrb_efficiency":
        # In this case tags[0] is the offtrial number
        node.add_input_opt("--bank-file", bank_file)
        node.add_opt("--trial-name", tags[0])
        node.add_opt("--injection-set-name", tags[1])
        # Output the sensitivity plot
        if kwargs["plot_bkgd"]:
            node.new_output_file_opt(
                workflow.analysis_time,
                ".png",
                "--background-output-file",
                tags=[*extra_tags, "max_background"],
            )
        # Output the exclusion distance plot and table
        else:
            node.add_input_opt("--onsource-file", onsource_file)
            node.new_output_file_opt(
                workflow.analysis_time,
                ".png",
                "--onsource-output-file",
                tags=["onsource", *extra_tags],
            )
            node.new_output_file_opt(
                workflow.analysis_time,
                ".json",
                "--exclusion-dist-output-file",
                tags=extra_tags,
            )
    elif exec_name == "pycbc_plot_bank_corner":
        node.new_output_file_opt(
            workflow.analysis_time, ".png", "--output-plot-file", tags=extra_tags + tags
        )
    else:
        node.new_output_file_opt(
            workflow.analysis_time, ".png", "--output-file", tags=extra_tags
        )
        if (
            exec_name in ["pygrb_plot_coh_ifosnr", "pygrb_plot_null_stats"]
            and "zoomin" in tags
        ):
            node.add_opt("--zoom-in")
    # Quantity to be displayed on the y-axis of the plot
    if exec_name in [
        "pygrb_plot_chisq_veto",
        "pygrb_plot_null_stats",
        "pygrb_plot_snr_timeseries",
    ]:
        node.add_opt("--y-variable", tags[0])
    # Quantity to be displayed on the x-axis of the plot
    elif exec_name == "pygrb_plot_stats_distribution":
        node.add_opt("--x-variable", tags[0])
    elif exec_name == "pygrb_plot_injs_results":
        # Variables to plot on x and y axes
        node.add_opt("--y-variable", tags[0])
        node.add_opt("--x-variable", tags[1])
        # Flag to plot found over missed or missed over found
        if tags[2] == "missed-on-top":
            node.add_opt("--" + tags[2])
        # Enable log axes
        subsection = "_".join(tags[0:2])
        for log_flag in ["x-log", "y-log"]:
            if workflow.cp.has_option_tags(exec_name, log_flag, tags=[subsection]):
                node.add_opt("--" + log_flag)
    elif exec_name == "pycbc_plot_bank_corner":
        node.add_opt("--no-suptitle")
        # if inj_file:
        #    node.add_opt('--title', f'\"{tags[0]} injections\"')
        #    params = workflow.cp.get_opt_tags(exec_name, 'parameters', ['injs'])
        # else:
        #    node.add_opt('--title', '\"Template bank\"')
        #    params = workflow.cp.get_opt_tags(exec_name, 'parameters', ['bank'])
        # node.add_opt('--parameters', params)

    # Add job node to workflow
    workflow += node

    return node, node.output_files


def make_pygrb_info_table(workflow, exec_name, out_dir, in_files=None, tags=None):
    """
    Setup a job to create an html snippet with the GRB trigger information
    or exlusion distances information.
    """
    # Organize tags
    tags = [] if tags is None else tags
    grb_name = workflow.cp.get("workflow", "trigger-name")
    extra_tags = ["GRB" + grb_name]

    # Initialize job node
    node = PlotExecutable(
        workflow.cp,
        exec_name,
        ifos=workflow.ifos,
        out_dir=out_dir,
        tags=tags + extra_tags,
    ).create_node()

    # Options
    if exec_name == "pygrb_grb_info_table":
        node.add_opt("--ifos", " ".join(workflow.ifos))
    elif exec_name == "pygrb_exclusion_dist_table":
        node.add_input_opt("--input-files", in_files)

    # Output
    node.new_output_file_opt(
        workflow.analysis_time, ".html", "--output-file", tags=extra_tags
    )

    # Add job node to workflow
    workflow += node

    return node, node.output_files


def make_pygrb_injs_tables(
    workflow,
    out_dir,
    bank_file,
    off_file,
    seg_files,
    inj_file=None,
    on_file=None,
    veto_file=None,
    tags=None,
):
    """
    Adds a job to make quiet-found and missed-found injection tables,
    or loudest trigger(s) table.
    """
    tags = [] if tags is None else tags

    # Executable
    exec_name = "pygrb_page_tables"
    # Initialize job node
    grb_name = workflow.cp.get("workflow", "trigger-name")
    extra_tags = ["GRB" + grb_name]
    node = PlotExecutable(
        workflow.cp,
        exec_name,
        ifos=workflow.ifos,
        out_dir=out_dir,
        tags=tags + extra_tags,
    ).create_node()
    # Pass the bank-file
    node.add_input_opt("--bank-file", bank_file)
    # Offsource input file (or equivalently trigger file for injections)
    offsource_file = off_file
    node.add_input_opt("--offsource-file", offsource_file)
    # Pass the veto and segment files (as File instances)
    if veto_file:
        node.add_input_opt("--veto-file", veto_file)
    node.add_input_list_opt("--seg-files", seg_files)
    # Handle input/output for injections
    if inj_file:
        # Found-missed injection file (passed as File instance)
        node.add_input_opt("--found-missed-file", inj_file)
        # Missed-found and quiet-found injections html output files
        for mf_or_qf in ["missed-found", "quiet-found"]:
            mf_or_qf_tags = [mf_or_qf.upper().replace("-", "_")]
            node.new_output_file_opt(
                workflow.analysis_time,
                ".html",
                "--" + mf_or_qf + "-injs-output-file",
                tags=extra_tags + mf_or_qf_tags,
            )
        # Quiet-found injections h5 output file
        node.new_output_file_opt(
            workflow.analysis_time,
            ".h5",
            "--quiet-found-injs-h5-output-file",
            tags=[*extra_tags, "QUIET_FOUND"],
        )
    # Handle input/output for onsource/offsource
    else:
        src_type = "offsource-trigs"
        if on_file:
            src_type = "onsource-trig"
            # Pass onsource input File instance
            node.add_input_opt("--onsource-file", on_file)
        # Loudest offsource/onsource triggers html and h5 output files
        src_type_tags = [src_type.upper().replace("-", "_")]
        node.new_output_file_opt(
            workflow.analysis_time,
            ".html",
            "--loudest-" + src_type + "-output-file",
            tags=extra_tags + src_type_tags,
        )
        node.new_output_file_opt(
            workflow.analysis_time,
            ".h5",
            "--loudest-" + src_type + "-h5-output-file",
            tags=extra_tags + src_type_tags,
        )

    # Add job node to the workflow
    workflow += node

    return node, node.output_files


# Based on setup_single_det_minifollowups
def setup_pygrb_minifollowups(
    workflow,
    followups_file,
    trigger_file,
    dax_output,
    out_dir,
    seg_files=None,
    veto_file=None,
    tags=None,
):
    """
    Create plots that followup the the loudest PyGRB triggers or
    missed injections from an HDF file.

    Parameters
    ----------
    workflow: pycbc.workflow.Workflow
        The core workflow instance we are populating
    followups_file: pycbc.workflow.File
        The File class holding the triggers/injections to follow up
    trigger_file: pycbc.workflow.File
        The File class holding the triggers
    dax_output: The directory that will contain the dax file
    out_dir: path
        The directory to store minifollowups result plots and files
    seg_files: {pycbc.workflow.FileList, optional}
        The list of segments Files
    veto_file: {pycbc.workflow.File, optional}
        The veto definer file
    tags: {None, optional}
        Tags to add to the minifollowups executables

    """
    logger.info("Entering minifollowups module")

    if not workflow.cp.has_section("workflow-minifollowups"):
        msg = "There is no [workflow-minifollowups] section in "
        msg += "the configuration file"
        logger.info(msg)
        logger.info("Leaving minifollowups")
        return None

    tags = [] if tags is None else tags

    # Turn the config file into a File instance
    config_path = os.path.abspath(
        dax_output + "/" + "_".join(tags) + "_minifollowup.ini"
    )
    write_config_for_subworkflow(workflow.cp, config_path)
    config_file = resolve_url_to_file(config_path)

    # Create the node
    exe = Executable(
        workflow.cp,
        "pygrb_minifollowups",
        ifos=workflow.ifos,
        out_dir=dax_output,
        tags=tags,
    )
    node = exe.create_node()

    node.add_input_opt("--trig-file", trigger_file)

    # Grab and pass all necessary files as File instances
    if seg_files:
        node.add_input_list_opt("--seg-files", seg_files)
    if veto_file:
        node.add_input_opt("--veto-file", veto_file)
    node.add_input_opt("--config-files", config_file)
    node.add_input_opt("--followups-file", followups_file)
    if tags:
        node.add_list_opt("--tags", tags)
    node.new_output_file_opt(workflow.analysis_time, ".dax", "--dax-file")
    node.new_output_file_opt(workflow.analysis_time, ".dax.map", "--output-map")

    name = node.output_files[0].name
    assert name.endswith(".dax")
    map_file = node.output_files[1]
    assert map_file.name.endswith(".map")

    wiki_file = File(
        workflow.ifos,
        "pygrb_minifollowups",
        workflow.analysis_time,
        file_url=os.path.join(out_dir, "_".join(tags) + "loudest_table.txt"),
    )
    node.add_output_opt("--wiki-file", wiki_file)
    add_wellfile_opt(node, workflow, "--html-layout-file", out_dir)

    node.add_opt("--workflow-name", name)
    node.add_opt("--output-dir", out_dir)
    node.add_opt("--dax-file-directory", ".")

    workflow += node

    # Execute this in a sub-workflow
    fil = node.output_files[0]
    job = SubWorkflow(fil.name, is_planned=False)
    # Declare the files the sub-workflow's jobs read, so that the planner is
    # given their locations from this workflow.
    input_files = [trigger_file, *(seg_files or [])]
    if veto_file:
        input_files.append(veto_file)
    job.add_inputs(*input_files)
    job.set_subworkflow_properties(
        map_file, staging_site=workflow.staging_site, cache_file=workflow.cache_file
    )
    job.add_into_workflow(workflow)
    logger.info("Leaving minifollowups module")

    return job


def _add_seg_plot_metadata(plot_path):
    """Add the metadata used by the results page to the segments plot."""
    from PIL import Image, PngImagePlugin

    im = Image.open(plot_path)
    meta = PngImagePlugin.PngInfo()
    meta.add_text("title", "Segments and analysis time")
    meta.add_text(
        "caption",
        "Segments (horizontal bands) available around the trigger time (orange "
        "vertical line) and analysis time used (orange box).",
    )
    meta.add_text(
        "cmd",
        "This plot is generated by the make_grb_segments_plot function during "
        "the pre-processing stage of the analysis.",
    )
    im.save(plot_path, "png", pnginfo=meta)


def setup_pygrb_results_workflow(
    workflow,
    rdir,
    trig_files,
    full_inj_files,
    inj_files,
    bank_file,
    seg_files,
    seg_plot,
    sky_grid_file,
    veto_file=None,
):
    """
    Add the jobs that produce the plots and tables of a PyGRB analysis, and
    write the layout of the results webpage that displays them.

    Parameters
    ----------
    workflow: pycbc.workflow.Workflow
        The core workflow instance we are populating
    rdir: pycbc.results.layout.SectionNumber
        The sections of the results webpage
    trig_files: FileList of clustered trigger files, in the order
        [ALL_TIMES, ONSOURCE, OFFSOURCE, OFFTRIAL_1, ..., OFFTRIAL_N]
    full_inj_files: FileList of generated injections
    inj_files: FileList of injection results
    bank_file: The template bank File object
    seg_files: FileList of the segment files
    seg_plot: The File object of the segments plot
    sky_grid_file: The sky grid File object
    veto_file: {None, optional}
        The veto File object

    """
    [all_times_file, onsource_file, offsource_file] = trig_files[0:3]
    offtrial_files = trig_files[3:]
    ifos = workflow.ifos

    # Sanity check and File instances of the input injection files
    inj_sets = workflow.cp.get_subsections("injections")
    eff_secs = workflow.cp.get_subsections("pygrb_efficiency")
    if not set(inj_sets).issubset(eff_secs):
        err_msg = "Each [injections] subsection requires at "
        err_msg += "least one dedicated [pygrb_efficiency] subsection. "
        err_msg += "[injections] subsections found: %s. "
        err_msg += "[pygrb_efficiency] subsections found: %s. "
        logger.error(err_msg, inj_sets, eff_secs)
    inj_sets = [i.upper() for i in inj_sets]
    out_dir = rdir.base

    #
    # Opening page
    #
    # Create the information table about the GRB trigger
    _, grb_info_table = make_pygrb_info_table(workflow, "pygrb_grb_info_table", out_dir)
    # Plot the search grid
    _, skygrid_plot = make_pygrb_plot(
        workflow, "pygrb_plot_skygrid", out_dir, sky_grid_file=sky_grid_file
    )
    # Plot the template bank
    _, bank_plot = make_pygrb_plot(
        workflow, "pycbc_plot_bank_corner", out_dir, bank_file=bank_file, tags=["bank"]
    )
    # The segments plot was produced in the preprocessing stage
    _add_seg_plot_metadata(seg_plot.storage_path)
    summary_layout = [
        (grb_info_table[0],),
        (seg_plot, skygrid_plot[0]),
        (bank_plot[0],),
    ]
    layout.two_column_layout(out_dir, summary_layout)

    #
    # Plot SNR timeseries
    #
    out_dir = rdir["offsource_triggers_vs_time"]

    # Retrieve the name of the preferred injection set for these plots
    tuning_inj_set = workflow.cp.get(
        "workflow-pygrb_results_workflow", "tuning-inj-set"
    )
    tuning_inj_set = tuning_inj_set.upper()
    logger.info("The tuning injections set is %s", tuning_inj_set)

    # Retrieve the injections result File corresponding to the tuning set
    tuning_inj_file = inj_files.find_output_with_tag(
        tuning_inj_set, fail_if_not_single_file=True
    )

    # Loop over timeseries request by the user
    timeseries = workflow.cp.get_subsections("pygrb_plot_snr_timeseries")
    out_dirs_dict = {
        "coherent": "offsource_triggers_vs_time/coh_snr_timeseries",
        "reweighted": "offsource_triggers_vs_time/reweighted_snr_timeseries",
        "single": "offsource_triggers_vs_time/single_ifo_snr_timeseries",
        "null": "offsource_triggers_vs_time/null_snr_timeseries",
    }
    for snr_type in timeseries:
        out_dir = rdir[out_dirs_dict[snr_type]]
        files = FileList([])
        # Only single SNR timeseries requires looping over IFOs
        ifos_to_loop = ifos if snr_type == "single" else [None]
        for ifo in ifos_to_loop:
            timeseries_plots = FileList([])
            # Plots without and with injections
            for inj_file in {None, tuning_inj_file}:
                # If the tuning injection file is used for the plot, include its
                # tags in assembling the plot name
                inj_tags = inj_file.tags if inj_file else []
                tags = [snr_type, *inj_tags]
                _, output_files = make_pygrb_plot(
                    workflow,
                    "pygrb_plot_snr_timeseries",
                    out_dir,
                    trig_file=offsource_file,
                    inj_file=inj_file,
                    ifo=ifo,
                    seg_files=seg_files,
                    veto_file=veto_file,
                    tags=tags,
                )
                # We want a File, not a 1-element list with a File
                # pycbc_pygrb_plot_snr_timeseries produces only one plot: take [0]
                timeseries_plots.append(output_files[0])
            files.append(timeseries_plots)
        layout.two_column_layout(out_dir, files)

    #
    # Signal consistency plots
    #
    out_dir = rdir["signal_consistency"]
    # Chisq veto vs SNR plots
    out_dir = rdir["signal_consistency/chi_squared_tests"]
    files = FileList([])
    # Loop over vetoes request by the user
    vetoes = workflow.cp.get_subsections("pygrb_plot_chisq_veto")
    for veto in vetoes:
        # Loop over IFOs only in the case of single detector chi-squares
        ifo_loop = [ifos[0]] if veto == "network" else ifos
        for ifo in ifo_loop:
            # Plot with and without injections
            for inj_file in {tuning_inj_file, None}:
                # If the tuning injection file is used for the plot, include its
                # tags in assembling the plot name
                inj_tags = inj_file.tags if inj_file else []
                tags = [veto, *inj_tags]
                ifo_arg = None if veto == "network" else ifo
                _, output_files = make_pygrb_plot(
                    workflow,
                    "pygrb_plot_chisq_veto",
                    out_dir,
                    trig_file=offsource_file,
                    inj_file=inj_file,
                    ifo=ifo_arg,
                    seg_files=seg_files,
                    veto_file=veto_file,
                    tags=tags,
                )
                # We want a File, not a 1-element list with a File
                # pycbc_pygrb_plot_chisq_veto produces only one plot: take [0]
                files.append(output_files[0])
    chisq_layout = list(layout.grouper(files, 2))
    layout.two_column_layout(out_dir, chisq_layout)

    # Single detector chi-square plots: zoomed in and zoomed out
    out_dir = rdir["signal_consistency/individual_detector_snrs"]
    files = FileList([])
    # Single IFO SNR vs Coherent SNR plots: zoomed in and zoomed out
    # Requires looping over IFOs
    if workflow.cp.has_section("pygrb_plot_coh_ifosnr"):
        for ifo in ifos:
            # Plots with and without injections
            for inj_file in {tuning_inj_file, None}:
                sngl_snr_plots = FileList([])
                for zoom_tag in [["zoomin"], []]:
                    # If the tuning injection file is used for the plot, include
                    # its tags in assembling the plot name
                    inj_tags = inj_file.tags if inj_file else []
                    tags = zoom_tag + inj_tags
                    # Single IFO SNR vs Coherent SNR
                    _, output_files = make_pygrb_plot(
                        workflow,
                        "pygrb_plot_coh_ifosnr",
                        out_dir,
                        trig_file=offsource_file,
                        inj_file=inj_file,
                        ifo=ifo,
                        seg_files=seg_files,
                        veto_file=veto_file,
                        tags=tags,
                    )
                    # We want a File, not a 1-element list with a File
                    # pycbc_pygrb_plot_cohifo_snr produces only one plot: take [0]
                    sngl_snr_plots.append(output_files[0])
                files.append(sngl_snr_plots)
        layout.two_column_layout(out_dir, files)
    else:
        msg = "No pygrb_plot_coh_ifosnr section found in the configuration file. "
        msg += "No coherent vs single detector SNR plots will be generated."
        logger.info(msg)

    # Null SNR/Overwhitened null stat vs Coherent SNR plots
    null_snr_out_dir = rdir["signal_consistency/null_snrs"]
    null_snr_files = FileList([])
    # Coincident SNR vs Coherent SNR plots
    coinc_out_dir = rdir["signal_consistency/coincident_snr"]
    coinc_files = FileList([])
    # Loop over null statistics requested by the user (including coincident SNR)
    nstats = workflow.cp.get_subsections("pygrb_plot_null_stats")
    for nstat in nstats:
        # Plots with and without injections, zoomed in and zoomed out
        for inj_file in {tuning_inj_file, None}:
            if nstat == "coincident":
                out_dir = coinc_out_dir
                files = coinc_files
            else:
                out_dir = null_snr_out_dir
                files = null_snr_files
            null_stats_plots = FileList([])
            for zoom_tag in [["zoomin"], []]:
                # If the tuning injection file is used for the plot, include its
                # tags in assembling the plot name
                inj_tags = inj_file.tags if inj_file else []
                tags = [nstat, *zoom_tag, *inj_tags]
                _, output_files = make_pygrb_plot(
                    workflow,
                    "pygrb_plot_null_stats",
                    out_dir,
                    trig_file=offsource_file,
                    inj_file=inj_file,
                    seg_files=seg_files,
                    veto_file=veto_file,
                    tags=tags,
                )
                # We want a File, not a 1-element list with a File
                # pycbc_pygrb_plot_null_stats produces only one plot: take [0]
                null_stats_plots.append(output_files[0])
            files.append(null_stats_plots)
    layout.two_column_layout(null_snr_out_dir, null_snr_files)
    layout.two_column_layout(coinc_out_dir, coinc_files)

    # layout.group_layout(rdir['coincident_triggers'],
    #                     closed_box_ifars + all_snrifar + [bank_plot[0][0]])

    #
    # Found/missed injections plots and tables
    #
    out_dir = rdir["injections"]
    ## Create corner plot
    inj_corner_plot_files = FileList([])
    for inj_set in inj_sets:
        # Retrieve the injections full files
        inj_file = full_inj_files.find_output_with_tag(
            inj_set, fail_if_not_single_file=True
        )
        # Generate a corner bank plot for each injection set
        tags = [inj_set]
        _, output_files = make_pygrb_plot(
            workflow,
            "pycbc_plot_bank_corner",
            out_dir,
            inj_file=inj_file,
            tags=tags,
        )
        inj_corner_plot_files.append(output_files[0])
    layout.single_layout(out_dir, inj_corner_plot_files)
    # Loop over injection plots requested by the user
    inj_plots = workflow.cp.get_subsections("pygrb_plot_injs_results")
    # The command above also picks up the injection set names so we remove them
    # from the set of requested injection plot types
    inj_plots = [inj_plot for inj_plot in inj_plots if inj_plot not in inj_sets]
    for inj_set in inj_sets:
        # Retrieve the injections result File corresponding to the inj_set label
        inj_file = inj_files.find_output_with_tag(inj_set, fail_if_not_single_file=True)
        out_dir = rdir["injections/" + inj_set]
        files = FileList([])
        # Generate plots: loop over inj_plots and found-missed/missed-found
        for inj_plot in inj_plots:
            y_qty, x_qty = inj_plot.split("_")
            ifos_to_loop = [None]
            if y_qty == "effsitedist":
                ifos_to_loop = ifos
            for ifo in ifos_to_loop:
                for fm_or_mf in ["missed-on-top", "found-on-top"]:
                    tags = [y_qty, x_qty, fm_or_mf, inj_set]
                    _, output_files = make_pygrb_plot(
                        workflow,
                        "pygrb_plot_injs_results",
                        out_dir,
                        trig_file=offsource_file,
                        ifo=ifo,
                        inj_file=inj_file,
                        seg_files=seg_files,
                        veto_file=veto_file,
                        tags=tags,
                    )
                    # A File is needed, not a 1-element list with a File
                    # pycbc_pygrb_plot_injs_results produces only one plot, so
                    # take the first output file
                    files.append(output_files[0])
        # Generate quiet-found and missed-found html tables
        _, [mf_table, qf_table, qf_h5] = make_pygrb_injs_tables(
            workflow,
            out_dir,
            bank_file,
            offsource_file,
            seg_files,
            inj_file=inj_file,
            veto_file=veto_file,
            tags=inj_file.tags,
        )

        inj_layout = [*list(layout.grouper(files, 2)), (mf_table,), (qf_table,)]
        layout.two_column_layout(out_dir, inj_layout)

        # Follow up of loudest N quiet/missed-found injections:
        # qf_h5 includes quiet and missed found injections
        out_dir = rdir["injections/" + inj_set + "loudest_quiet_found_followups"]
        if not workflow.cp.has_section("workflow-minifollowups"):
            msg = "No [workflow-minifollowups] section found in "
            msg += "configuration file"
            logger.info(msg)
        else:
            logger.info("Entering minifollowups module for injections")
            setup_pygrb_minifollowups(
                workflow,
                qf_h5,
                offsource_file,
                "daxes",
                out_dir,
                seg_files=seg_files,
                veto_file=veto_file,
                tags=[*inj_file.tags, "loudest_quiet_found_injs"],
            )
            logger.info("Leaving minifollowups")

    #
    # FAP distributions
    #
    out_dir = rdir["loudest_offsource_events"]
    files = []
    # Loop over statistics requested by the user
    stats = workflow.cp.get_subsections("pygrb_plot_stats_distribution")
    for stat in stats:
        _, output_file = make_pygrb_plot(
            workflow,
            "pygrb_plot_stats_distribution",
            out_dir,
            trig_file=offsource_file,
            inj_file=inj_file,
            seg_files=seg_files,
            veto_file=veto_file,
            tags=[stat],
        )
        # We want a File, not a 1-element list with a File
        # pycbc_pygrb_plot_stats_distribution produces only one plot: take [0]
        files.append(output_file[0])

    # Generate loudest off-source triggers table
    _, [lofft_table, lofft_h5] = make_pygrb_injs_tables(
        workflow, out_dir, bank_file, offsource_file, seg_files
    )

    lofft_layout = [*list(layout.grouper(files, 2)), (lofft_table,)]
    layout.two_column_layout(out_dir, lofft_layout)

    # Follow up N loudest offsource triggers (parent-child)
    out_dir = rdir["loudest_offsource_events/followups"]
    if not workflow.cp.has_section("workflow-minifollowups"):
        msg = "No [workflow-minifollowups] section found in "
        msg += "configuration file"
        logger.info(msg)
    else:
        logger.info("Entering minifollowups module for offsource")
        setup_pygrb_minifollowups(
            workflow,
            lofft_h5,
            offsource_file,
            "daxes",
            out_dir,
            seg_files=seg_files,
            veto_file=veto_file,
            tags=["loudest_offsource_events"],
        )
        logger.info("Leaving minifollowups")

    #
    # Exclusion distance and efficiency plots based on offtrials
    #
    out_dir = rdir["exclusion_distances"]

    # Offtrials and injection sets requested by the user
    num_trials = int(workflow.cp.get("trig_combiner", "num-trials"))
    offtrials = [f"offtrial_{i + 1}" for i in range(num_trials)]

    # Sensitivity plots of each injection set
    out_dir = rdir["exclusion_distances"]
    bkgd_plots = []
    for inj_set in inj_sets:
        # Retrieve the injections result File for the inj_set label
        inj_file = inj_files.find_output_with_tag(inj_set, fail_if_not_single_file=True)
        tags = [offtrials[0], *inj_file.tags]
        _, output_files = make_pygrb_plot(
            workflow,
            "pygrb_efficiency",
            out_dir,
            trig_file=offsource_file,
            inj_file=inj_file,
            bank_file=bank_file,
            seg_files=seg_files,
            veto_file=veto_file,
            tags=tags,
            plot_bkgd=True,
        )
        bkgd_plots += [output_files[0]]
    bkgd_plots = list(layout.grouper(bkgd_plots, 2))
    layout.two_column_layout(out_dir, bkgd_plots)

    # Exclusion distances of each injection set
    for i, offtrial in enumerate(offtrials):
        out_dir = rdir[f"exclusion_distances/{offtrial}"]
        eff_layout = []
        for inj_set in inj_sets:
            # Retrieve the injections result File for the inj_set label
            inj_file = inj_files.find_output_with_tag(
                inj_set, fail_if_not_single_file=True
            )
            tags = [offtrial, *inj_file.tags]
            _, output_files = make_pygrb_plot(
                workflow,
                "pygrb_efficiency",
                out_dir,
                trig_file=offsource_file,
                onsource_file=offtrial_files[i],
                inj_file=inj_file,
                bank_file=bank_file,
                seg_files=seg_files,
                veto_file=veto_file,
                tags=tags,
                plot_bkgd=False,
            )
            eff_plot = output_files[0]
            json_input = output_files[1]
            # Create information table about exclusion distances
            _, excl_dist_table = make_pygrb_info_table(
                workflow,
                "pygrb_exclusion_dist_table",
                out_dir,
                in_files=json_input,
                tags=tags,
            )
            eff_layout += [(eff_plot, excl_dist_table[0])]
        layout.two_column_layout(out_dir, eff_layout, offtrial)

    # Save global config file
    base = rdir["workflow/configuration"]
    makedir(base)
    ini_file_path = os.path.join(base, workflow.name + ".ini")
    with open(ini_file_path, "w") as ini_fh:
        workflow.cp.write(ini_fh)
    ini_file = FileList(
        [
            File(
                workflow.ifos,
                "",
                workflow.analysis_time,
                file_url="file://" + ini_file_path,
            )
        ]
    )
    layout.single_layout(base, ini_file)

    # Create versioning information
    make_versioning_page(
        workflow,
        workflow.cp,
        rdir["workflow/version"],
    )

    # Protect the open box results folder
    out_dir = rdir["open_box"]
    makedir(out_dir)
    os.chmod(out_dir, 0o0700)

    # Generate loudest on-source trigger table
    _, [lont_table, lont_h5] = make_pygrb_injs_tables(
        workflow, out_dir, bank_file, offsource_file, seg_files, on_file=onsource_file
    )

    # Reweighted and coherent SNR timeseries for the full data stretch
    timeseries_plots = FileList([])
    for snr_type in ["reweighted", "coherent"]:
        _, output_files = make_pygrb_plot(
            workflow,
            "pygrb_plot_snr_timeseries",
            out_dir,
            trig_file=all_times_file,
            seg_files=seg_files,
            veto_file=veto_file,
            tags=[snr_type, "alltimes"],
        )
        timeseries_plots.extend(output_files)
    openbox_layout = [(lont_table,), *list(layout.grouper(timeseries_plots, 2))]
    layout.two_column_layout(out_dir, openbox_layout)

    # Loudest on-source trigger follow up
    out_dir = rdir["open_box/loudest_event_followup"]
    if not workflow.cp.has_section("workflow-minifollowups"):
        msg = "No [workflow-minifollowups] section found in "
        msg += "configuration file"
        logger.info(msg)
    else:
        logger.info("Entering minifollowups module for loudest onsource")
        setup_pygrb_minifollowups(
            workflow,
            lont_h5,
            onsource_file,
            "daxes",
            out_dir,
            seg_files=seg_files,
            veto_file=veto_file,
            tags=["loudest_onsource_event"],
        )
        logger.info("Leaving onsource minifollowups")

    # Exclusion distances and efficiency plots based on the on-source
    out_dir = rdir["open_box/exclusion_distances"]
    eff_layout = []
    for inj_set in inj_sets:
        # Retrieve the injections result File corresponding to the inj_set label
        inj_file = inj_files.find_output_with_tag(inj_set, fail_if_not_single_file=True)
        tags = ["onsource", *inj_file.tags]
        _, output_files = make_pygrb_plot(
            workflow,
            "pygrb_efficiency",
            out_dir,
            trig_file=offsource_file,
            onsource_file=onsource_file,
            inj_file=inj_file,
            bank_file=bank_file,
            seg_files=seg_files,
            veto_file=veto_file,
            tags=tags,
            plot_bkgd=False,
        )
        eff_plot = output_files[0]
        json_input = output_files[1]
        # Create information table about exclusion distances
        _, excl_dist_table = make_pygrb_info_table(
            workflow,
            "pygrb_exclusion_dist_table",
            out_dir,
            in_files=json_input,
            tags=tags,
        )
        eff_layout += [(eff_plot, excl_dist_table[0])]
    layout.two_column_layout(out_dir, eff_layout)
