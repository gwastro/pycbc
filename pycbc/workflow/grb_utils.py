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

import os
import logging
import numpy as np
from scipy.stats import rayleigh
from gwdatafind.utils import filename_metadata

from pycbc import makedir
from pycbc.workflow.core import \
    File, FileList, resolve_url_to_file, \
    Executable, Node
from pycbc.workflow.jobsetup import select_generic_executable
from pycbc.workflow.pegasus_workflow import SubWorkflow
from pycbc.workflow.plotting import PlotExecutable

logger = logging.getLogger('pycbc.workflow.grb_utils')


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
    exe_path = wflow.cp.get('executables', curr_exe)
    exe_name = os.path.basename(exe_path)
    exe_to_class_map = {
        'pycbc_grb_trig_combiner': PycbcGrbTrigCombinerExecutable,
        'pycbc_grb_trig_cluster': PycbcGrbTrigClusterExecutable,
        'pycbc_grb_inj_finder': PycbcGrbInjFinderExecutable
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
    --------
    cp : pycbc.workflow.configuration.WorkflowConfigParser object
    The modified WorkflowConfigParser object.

    """
    cp.set("workflow", "start-time", str(start))
    cp.set("workflow", "end-time", str(end))

    return cp


def make_gating_node(workflow, datafind_files, outdir=None, tags=None):
    '''
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
    --------
    condition_strain_nodes : list
        List containing the pycbc.workflow.core.Node objects representing the
        autogating jobs.
    condition_strain_outs : pycbc.workflow.core.FileList
        FileList containing the pycbc.workflow.core.File objects representing
        the gated frame files.
    '''

    cp = workflow.cp
    if tags is None:
        tags = []

    condition_strain_class = select_generic_executable(workflow,
                                                       "condition_strain")
    condition_strain_nodes = []
    condition_strain_outs = FileList([])
    for ifo in workflow.ifos:
        input_files = FileList([datafind_file for datafind_file in
                                datafind_files if datafind_file.ifo == ifo])
        condition_strain_jobs = condition_strain_class(cp, "condition_strain",
                                                       ifos=ifo,
                                                       out_dir=outdir,
                                                       tags=tags)
        condition_strain_node, condition_strain_out = \
            condition_strain_jobs.create_node(input_files, tags=tags)
        condition_strain_nodes.append(condition_strain_node)
        condition_strain_outs.extend(FileList([condition_strain_out]))

    return condition_strain_nodes, condition_strain_outs


def fermi_core_tail_model(
        sky_err, rad, core_frac=0.98, core_sigma=3.6, tail_sigma=29.6):
    """Fermi systematic error model following
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
    _______
    tuple
        Tuple containing the core and tail probability distributions
        as a function of radius.
    """
    scaledsq = sky_err**2 / -2 / np.log(0.32)
    return (
        frac * (1 - np.exp(-0.5 * (rad / np.sqrt(scaledsq + sigma**2))**2))
        for frac, sigma
        in zip([core_frac, 1 - core_frac], [core_sigma, tail_sigma]))


def get_sky_grid_scale(
        sky_error=0.0, containment=0.9, upscale=False, fermi_sys=False,
        precision=1e-3, **kwargs):
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
    _______

    float
        Sky error radius in degrees.
    """
    if fermi_sys:
        lims = (0.5, 4)
        radii = np.linspace(
            lims[0] * sky_error, lims[1] * sky_error,
            int((lims[1] - lims[0]) * sky_error / precision) + 1)
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
    grb_name = workflow.cp.get('workflow', 'trigger-name')
    extra_tags = ['GRB'+grb_name]
    node = Executable(workflow.cp, 'make_sky_grid',
                      ifos=workflow.ifos, out_dir=out_dir,
                      tags=tags+extra_tags).create_node()
    node.add_opt('--instruments', ' '.join(workflow.ifos))
    node.new_output_file_opt(workflow.analysis_time, '.h5', '--output',
                             tags=extra_tags, store_file=True)

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
        raise ValueError("Refusing to overwrite %s." % tc_path)
    tc_file = open(tc_path, "w")
    tc_file.write("[prior-tc]\n")
    tc_file.write("name = uniform\n")
    tc_file.write("min-tc = %s\n" % wflow.analysis_time[0])
    tc_file.write("max-tc = %s\n\n" % wflow.analysis_time[1])
    tc_file.write("[constraint-tc]\n")
    tc_file.write("name = custom\n")
    tc_file.write("constraint_arg = (tc < %s) | (tc > %s)\n" %
                  (buffer_seg[0], buffer_seg[1]))
    tc_file.close()

    # Add the tc-prior configuration file url to wflow.cp if necessary
    tc_file_path = "file://"+tc_path
    for inj_sec in wflow.cp.get_subsections("injections"):
        config_urls = wflow.cp.get("workflow-injections",
                                   inj_sec+"-config-files")
        config_urls = [url.strip() for url in config_urls.split(",")]
        if tc_file_path not in config_urls:
            config_urls += [tc_file_path]
        config_urls = ', '.join([str(item) for item in config_urls])
        wflow.cp.set("workflow-injections",
                     inj_sec+"-config-files",
                     config_urls)


def setup_pygrb_pp_workflow(wf, pp_dir, seg_dir, segment, bank_file,
                            insp_files, inj_files, inj_insp_files, inj_tags):
    """
    Generate post-processing section of PyGRB offline workflow

    Parameters
    ----------
    wf : The workflow object
    pp_dir : The directory where the post-processing files will be stored
    seg_dir : The directory where the segment files are stored
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
    node, trig_files = job_instance.create_node(wf.ifo_string,
                                                seg_dir,
                                                segment,
                                                insp_files,
                                                pp_dir,
                                                bank_file)
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
        tag_inj_files = FileList([f for f in inj_files
                                  if inj_tag in f.tags])
        # The here stems from the injection group information
        # being stored in the second tag. This could be improved
        # depending on the final implementation of injections
        tag_insp_files = FileList([f for f in inj_insp_files
                                   if inj_tag in f.tags[1]])
        node, inj_find_file = job_instance.create_node(
                                           tag_inj_files, tag_insp_files,
                                           bank_file, pp_dir)
        wf.add_node(node)
        inj_find_files.append(inj_find_file)

    return trig_files, clustered_files, inj_find_files


class PycbcGrbTrigCombinerExecutable(Executable):
    """ The class responsible for creating jobs
    for ''pycbc_grb_trig_combiner''.
    """

    current_retention_level = Executable.ALL_TRIGGERS

    def __init__(self, cp, name):
        super().__init__(cp=cp, name=name)
        self.trigger_name = cp.get('workflow', 'trigger-name')
        self.trig_start_time = cp.get('workflow', 'start-time')
        self.num_trials = int(cp.get('trig_combiner', 'num-trials'))

    def create_node(self, ifo_tag, seg_dir, segment, insp_files,
                    out_dir, bank_file, tags=None):
        node = Node(self)
        node.add_opt('--verbose')
        node.add_opt("--ifo-tag", ifo_tag)
        node.add_opt("--grb-name", self.trigger_name)
        node.add_opt("--trig-start-time", self.trig_start_time)
        node.add_opt("--segment-dir", seg_dir)
        node.add_input_list_opt("--input-files", insp_files)
        node.add_opt("--user-tag", "PYGRB")
        node.add_input_opt("--bank-file", bank_file)
        # Prepare output file tag
        user_tag = f"PYGRB_GRB{self.trigger_name}"
        if tags:
            user_tag += "_{}".format(tags)
        # Add on/off source and off trial outputs
        output_files = FileList([])
        outfile_types = ['ALL_TIMES', 'ONSOURCE', 'OFFSOURCE']
        for i in range(self.num_trials):
            outfile_types.append("OFFTRIAL_{}".format(i+1))
        for out_type in outfile_types:
            out_name = "{}-{}_{}-{}-{}.h5".format(
                       ifo_tag, user_tag, out_type,
                       segment[0], segment[1]-segment[0])
            out_file = File(ifo_tag, 'trig_combiner', segment,
                            file_url=os.path.join(out_dir, out_name))
            node.add_output(out_file)
            output_files.append(out_file)

        return node, output_files


class PycbcGrbTrigClusterExecutable(Executable):
    """ The class responsible for creating jobs
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
        out_name = "{}-{}_CLUSTERED-{}-{}.h5".format(ifotag, filetag,
                                                     start, end-start)
        out_file = File(ifotag, 'trig_cluster', segment,
                        file_url=os.path.join(out_dir, out_name))
        node.add_output(out_file)

        return node, out_file


class PycbcGrbInjFinderExecutable(Executable):
    """The class responsible for creating jobs for ``pycbc_grb_inj_finder``
    """
    current_retention_level = Executable.ALL_TRIGGERS

    def __init__(self, cp, exe_name):
        super().__init__(cp=cp, name=exe_name)

    def create_node(self, inj_files, inj_insp_files, bank_file,
                    out_dir, tags=None):
        if tags is None:
            tags = []
        node = Node(self)
        node.add_input_list_opt('--input-files', inj_insp_files)
        node.add_input_list_opt('--inj-files', inj_files)
        node.add_input_opt('--bank-file', bank_file)
        ifo_tag, desc, segment = filename_metadata(inj_files[0].name)
        desc = '_'.join(desc.split('_')[:-1])
        out_name = "{}-{}_FOUNDMISSED-{}-{}.h5".format(
            ifo_tag, desc, segment[0], abs(segment))
        out_file = File(ifo_tag, 'inj_finder', segment,
                        os.path.join(out_dir, out_name), tags=tags)
        node.add_output(out_file)
        return node, out_file


def build_segment_filelist(seg_dir):
    """Construct a FileList instance containing all segments txt files"""

    # Needs to be in this order for consistency with _read_seg_files
    file_names = ["bufferSeg.txt", "offSourceSeg.txt", "onSourceSeg.txt"]
    seg_files = [os.path.join(seg_dir, fn) for fn in file_names]
    seg_files = [resolve_url_to_file(sf) for sf in seg_files]
    seg_files = FileList(seg_files)

    return seg_files


def make_pygrb_plot(workflow, exec_name, out_dir,
                    ifo=None, inj_file=None, trig_file=None,
                    onsource_file=None, bank_file=None,
                    seg_files=None, veto_file=None, tags=None, **kwargs):
    """Adds a node for a plot of PyGRB results to the workflow"""

    tags = [] if tags is None else tags

    # Initialize job node with its tags
    grb_name = workflow.cp.get('workflow', 'trigger-name')
    extra_tags = ['GRB'+grb_name]
    if ifo:
        extra_tags.append(ifo)
    node = PlotExecutable(workflow.cp, exec_name, ifos=workflow.ifos,
                          out_dir=out_dir,
                          tags=tags+extra_tags).create_node()
    if trig_file:
        node.add_input_opt('--trig-file', trig_file)
    # Pass the veto and segment files and options
    if seg_files:
        node.add_input_list_opt('--seg-files', seg_files)
    if veto_file:
        node.add_input_opt('--veto-file', veto_file)
    # Option to show the onsource trial if this is a plot of all data
    if exec_name == 'pygrb_plot_snr_timeseries' and 'alltimes' in tags:
        node.add_opt('--onsource')
    if exec_name in ['pygrb_plot_injs_results',
                     'pygrb_plot_snr_timeseries']:
        trig_time = workflow.cp.get('workflow', 'trigger-time')
        node.add_opt('--trigger-time', trig_time)
    # Pass the injection file as an input File instance
    if inj_file is not None and exec_name not in \
            ['pygrb_plot_skygrid', 'pygrb_plot_stats_distribution']:
        node.add_input_opt('--found-missed-file', inj_file)
    # IFO option
    if ifo:
        node.add_opt('--ifo', ifo)
    # Output files and final input file (passed as a File instance)
    if exec_name == 'pygrb_efficiency':
        # In this case tags[0] is the offtrial number
        node.add_input_opt('--bank-file', bank_file)
        node.add_opt('--trial-name', tags[0])
        node.add_opt('--injection-set-name', tags[1])
        # Output the sensitivity plot
        if kwargs['plot_bkgd']:
            node.new_output_file_opt(workflow.analysis_time, '.png',
                                     '--background-output-file',
                                     tags=extra_tags+['max_background'])
        # Output the exclusion distance plot and table
        else:
            node.add_input_opt('--onsource-file',
                               onsource_file)
            node.new_output_file_opt(workflow.analysis_time, '.png',
                                     '--onsource-output-file',
                                     tags=['onsource']+extra_tags)
            node.new_output_file_opt(workflow.analysis_time, '.json',
                                     '--exclusion-dist-output-file',
                                     tags=extra_tags)
    else:
        node.new_output_file_opt(workflow.analysis_time, '.png',
                                 '--output-file', tags=extra_tags)
        if exec_name in ['pygrb_plot_coh_ifosnr', 'pygrb_plot_null_stats'] \
                and 'zoomin' in tags:
            node.add_opt('--zoom-in')
    # Quantity to be displayed on the y-axis of the plot
    if exec_name in ['pygrb_plot_chisq_veto', 'pygrb_plot_null_stats',
                     'pygrb_plot_snr_timeseries']:
        node.add_opt('--y-variable', tags[0])
    # Quantity to be displayed on the x-axis of the plot
    elif exec_name == 'pygrb_plot_stats_distribution':
        node.add_opt('--x-variable', tags[0])
    elif exec_name == 'pygrb_plot_injs_results':
        # Variables to plot on x and y axes
        node.add_opt('--y-variable', tags[0])
        node.add_opt('--x-variable', tags[1])
        # Flag to plot found over missed or missed over found
        if tags[2] == 'missed-on-top':
            node.add_opt('--'+tags[2])
        # Enable log axes
        subsection = '_'.join(tags[0:2])
        for log_flag in ['x-log', 'y-log']:
            if workflow.cp.has_option_tags(exec_name, log_flag,
                                           tags=[subsection]):
                node.add_opt('--'+log_flag)

    # Add job node to workflow
    workflow += node

    return node, node.output_files


def make_pygrb_info_table(workflow, exec_name, out_dir, in_files=None,
                          tags=None):
    """
    Setup a job to create an html snippet with the GRB trigger information
    or exlusion distances information.
    """

    # Organize tags
    tags = [] if tags is None else tags
    grb_name = workflow.cp.get('workflow', 'trigger-name')
    extra_tags = ['GRB'+grb_name]

    # Initialize job node
    node = PlotExecutable(workflow.cp, exec_name,
                          ifos=workflow.ifos, out_dir=out_dir,
                          tags=tags+extra_tags).create_node()

    # Options
    if exec_name == 'pygrb_grb_info_table':
        node.add_opt('--ifos', ' '.join(workflow.ifos))
    elif exec_name == 'pygrb_exclusion_dist_table':
        node.add_input_opt('--input-files', in_files)

    # Output
    node.new_output_file_opt(workflow.analysis_time, '.html',
                             '--output-file', tags=extra_tags)

    # Add job node to workflow
    workflow += node

    return node, node.output_files


def make_pygrb_injs_tables(workflow, out_dir, bank_file, off_file, seg_files,
                           inj_file=None, on_file=None, veto_file=None,
                           tags=None):
    """
    Adds a job to make quiet-found and missed-found injection tables,
    or loudest trigger(s) table."""

    tags = [] if tags is None else tags

    # Executable
    exec_name = 'pygrb_page_tables'
    # Initialize job node
    grb_name = workflow.cp.get('workflow', 'trigger-name')
    extra_tags = ['GRB'+grb_name]
    node = PlotExecutable(workflow.cp, exec_name,
                          ifos=workflow.ifos, out_dir=out_dir,
                          tags=tags+extra_tags).create_node()
    # Pass the bank-file
    node.add_input_opt('--bank-file', bank_file)
    # Offsource input file (or equivalently trigger file for injections)
    offsource_file = off_file
    node.add_input_opt('--offsource-file', offsource_file)
    # Pass the veto and segment files (as File instances)
    if veto_file:
        node.add_input_opt('--veto-file', veto_file)
    node.add_input_list_opt('--seg-files', seg_files)
    # Handle input/output for injections
    if inj_file:
        # Found-missed injection file (passed as File instance)
        node.add_input_opt('--found-missed-file', inj_file)
        # Missed-found and quiet-found injections html output files
        for mf_or_qf in ['missed-found', 'quiet-found']:
            mf_or_qf_tags = [mf_or_qf.upper().replace('-', '_')]
            node.new_output_file_opt(workflow.analysis_time, '.html',
                                     '--'+mf_or_qf+'-injs-output-file',
                                     tags=extra_tags+mf_or_qf_tags)
        # Quiet-found injections h5 output file
        node.new_output_file_opt(workflow.analysis_time, '.h5',
                                 '--quiet-found-injs-h5-output-file',
                                 tags=extra_tags+['QUIET_FOUND'])
    # Handle input/output for onsource/offsource
    else:
        src_type = 'offsource-trigs'
        if on_file:
            src_type = 'onsource-trig'
            # Pass onsource input File instance
            node.add_input_opt('--onsource-file', on_file)
        # Loudest offsource/onsource triggers html and h5 output files
        src_type_tags = [src_type.upper().replace('-', '_')]
        node.new_output_file_opt(workflow.analysis_time, '.html',
                                 '--loudest-'+src_type+'-output-file',
                                 tags=extra_tags+src_type_tags)
        node.new_output_file_opt(workflow.analysis_time, '.h5',
                                 '--loudest-'+src_type+'-h5-output-file',
                                 tags=extra_tags+src_type_tags)

    # Add job node to the workflow
    workflow += node

    return node, node.output_files


# Based on setup_single_det_minifollowups
def setup_pygrb_minifollowups(workflow, followups_file, trigger_file,
                              dax_output, out_dir, seg_files=None,
                              veto_file=None, tags=None):
    """ Create plots that followup the the loudest PyGRB triggers or
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

    logging.info('Entering minifollowups module')

    if not workflow.cp.has_section('workflow-minifollowups'):
        msg = 'There is no [workflow-minifollowups] section in '
        msg += 'the configuration file'
        logging.info(msg)
        logging.info('Leaving minifollowups')
        return

    tags = [] if tags is None else tags
    makedir(dax_output)

    # Turn the config file into a File instance
    config_path = os.path.abspath(dax_output + '/' +
                                  '_'.join(tags) + '_minifollowup.ini')
    workflow.cp.write(open(config_path, 'w'))
    config_file = resolve_url_to_file(config_path)

    # wikifile = curr_ifo + '_'.join(tags) + 'loudest_table.txt'
    wikifile = '_'.join(tags) + 'loudest_table.txt'

    # Create the node
    exe = Executable(workflow.cp, 'pygrb_minifollowups',
                     ifos=workflow.ifos, out_dir=dax_output,
                     tags=tags)
    node = exe.create_node()

    node.add_input_opt('--trig-file', trigger_file)

    # Grab and pass all necessary files as File instances
    if seg_files:
        node.add_input_list_opt('--seg-files', seg_files)
    if veto_file:
        node.add_input_opt('--veto-file', veto_file)
    node.add_input_opt('--config-files', config_file)
    node.add_input_opt('--followups-file', followups_file)
    node.add_opt('--wiki-file', wikifile)
    if tags:
        node.add_list_opt('--tags', tags)
    node.new_output_file_opt(workflow.analysis_time, '.dax', '--dax-file')
    node.new_output_file_opt(workflow.analysis_time, '.dax.map',
                             '--output-map')

    name = node.output_files[0].name
    assert name.endswith('.dax')
    map_file = node.output_files[1]
    assert map_file.name.endswith('.map')

    node.add_opt('--workflow-name', name)
    node.add_opt('--output-dir', out_dir)
    node.add_opt('--dax-file-directory', '.')

    workflow += node

    # Execute this in a sub-workflow
    fil = node.output_files[0]
    job = SubWorkflow(fil.name, is_planned=False)
    job.set_subworkflow_properties(map_file,
                                   staging_site=workflow.staging_site,
                                   cache_file=workflow.cache_file)
    job.add_into_workflow(workflow)
    logging.info('Leaving minifollowups module')


def setup_pygrb_results_workflow(workflow, res_dir, trig_files,
                                 inj_files, bank_file, seg_dir,
                                 veto_file=None, tags=None,
                                 explicit_dependencies=None):
    """Create subworkflow to produce plots, tables,
    and results webpage for a PyGRB analysis.

    Parameters
    ----------
    workflow: pycbc.workflow.Workflow
        The core workflow instance we are populating
    res_dir: The post-processing directory where
        results (plots, etc.) will be stored
    trig_files: FileList of trigger files
    inj_files: FileList of injection results
    bank_file: The template bank File object
    seg_dir: The directory path with the segments files
    veto_file: {None, optional}
        The veto File object
    tags: {None, optional}
        Tags to add to the executables
    explicit_dependencies: {None, optional}
        nodes that must precede this
    """

    tags = [] if tags is None else tags
    dax_output = res_dir+'/webpage_daxes'
    # _workflow.makedir(dax_output)
    makedir(dax_output)

    # Create the node
    exe = Executable(workflow.cp, 'pygrb_results_workflow',
                     ifos=workflow.ifo_string, out_dir=dax_output,
                     tags=tags)
    node = exe.create_node()
    # Grab and pass all necessary files
    node.add_input_list_opt('--trig-files', trig_files)
    # node.add_input_opt('--config-files', config_file)
    node.add_input_list_opt('--inj-files', inj_files)
    node.add_input_opt('--bank-file', bank_file)
    node.add_opt('--segment-dir', seg_dir)
    if veto_file:
        node.add_input_opt('--veto-file', veto_file)

    if tags:
        node.add_list_opt('--tags', tags)

    node.new_output_file_opt(workflow.analysis_time, '.dax',
                             '--dax-file', tags=tags)
    node.new_output_file_opt(workflow.analysis_time, '.map',
                             '--output-map', tags=tags)
    # + ['MAP'], use_tmp_subdirs=True)
    name = node.output_files[0].name
    assert name.endswith('.dax')
    map_file = node.output_files[1]
    assert map_file.name.endswith('.map')
    node.add_opt('--workflow-name', name)
    # This is the output dir for the products of this node, namely dax and map
    node.add_opt('--output-dir', res_dir)
    node.add_opt('--dax-file-directory', '.')

    # Turn the config file into a File instance
    config_path = os.path.abspath(dax_output + '/' +
                                  '_'.join(tags) + 'webpage.ini')
    workflow.cp.write(open(config_path, 'w'))
    config_file = resolve_url_to_file(config_path)
    node.add_input_opt('--config-files', config_file)

    # Track additional ini file produced by pycbc_pygrb_results_workflow
    out_file = File(workflow.ifos,
                    'pygrb_results_workflow',
                    workflow.analysis_time,
                    file_url=os.path.join(dax_output, name+'.ini'))
    node.add_output(out_file)

    # Add node to the workflow
    workflow += node
    if explicit_dependencies is not None:
        for dep in explicit_dependencies:
            workflow.add_explicit_dependancy(dep, node)

    # Execute this in a sub-workflow
    job = SubWorkflow(name, is_planned=False)  # , _id='results')
    job.set_subworkflow_properties(map_file,
                                   staging_site=workflow.staging_site,
                                   cache_file=workflow.cache_file)
    job.add_into_workflow(workflow)

    return node.output_files
