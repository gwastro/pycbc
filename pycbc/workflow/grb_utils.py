# Copyright (C) 2015  Andrew Williamson
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
import shutil
from urllib.request import pathname2url
from urllib.parse import urljoin
import numpy as np
from scipy.stats import rayleigh
from ligo import segments
from ligo.lw import ligolw, lsctables, utils
from gwdatafind.utils import filename_metadata
from pycbc.workflow.core import \
    File, FileList, resolve_url_to_file, Executable, Node
from pycbc.workflow.jobsetup import select_generic_executable


def select_grb_pp_class(wflow, curr_exe):
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
        'pycbc_grb_inj_finder': PycbcGrbInjFinderExecutable,
        'pycbc_grb_inj_combiner': PycbcGrbInjCombinerExecutable
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


def get_coh_PTF_files(cp, ifos, run_dir, bank_veto=False, summary_files=False):
    """
    Retrieve files needed to run coh_PTF jobs within a PyGRB workflow

    Parameters
    ----------
    cp : pycbc.workflow.configuration.WorkflowConfigParser object
    The parsed configuration options of a pycbc.workflow.core.Workflow.

    ifos : str
    String containing the analysis interferometer IDs.

    run_dir : str
    The run directory, destination for retrieved files.

    bank_veto : Boolean
    If true, will retrieve the bank_veto_bank.xml file.

    summary_files : Boolean
    If true, will retrieve the summary page style files.

    Returns
    -------
    file_list : pycbc.workflow.FileList object
    A FileList containing the retrieved files.
    """
    if os.getenv("LAL_SRC") is None:
        raise ValueError("The environment variable LAL_SRC must be set to a "
                         "location containing the file lalsuite.git")
    else:
        lalDir = os.getenv("LAL_SRC")
        sci_seg = segments.segment(int(cp.get("workflow", "start-time")),
                                   int(cp.get("workflow", "end-time")))
        file_list = FileList([])

        # Bank veto
        if bank_veto:
            shutil.copy("%s/lalapps/src/ring/coh_PTF_config_files/" \
                        "bank_veto_bank.xml" % lalDir, "%s" % run_dir)
            bank_veto_url = "file://localhost%s/bank_veto_bank.xml" % run_dir
            bank_veto = File(ifos, "bank_veto_bank", sci_seg,
                             file_url=bank_veto_url)
            # FIXME: Is this an input file? If so use the from_path classmethod
            bank_veto.add_pfn(bank_veto.cache_entry.path, site="local")
            file_list.extend(FileList([bank_veto]))

        if summary_files:
            # summary.js file
            shutil.copy("%s/lalapps/src/ring/coh_PTF_config_files/" \
                        "coh_PTF_html_summary.js" % lalDir, "%s" % run_dir)
            summary_js_url = "file://localhost%s/coh_PTF_html_summary.js" \
                             % run_dir
            summary_js = File(ifos, "coh_PTF_html_summary_js", sci_seg,
                              file_url=summary_js_url)
            summary_js.add_pfn(summary_js.cache_entry.path, site="local")
            file_list.extend(FileList([summary_js]))

            # summary.css file
            shutil.copy("%s/lalapps/src/ring/coh_PTF_config_files/" \
                        "coh_PTF_html_summary.css" % lalDir, "%s" % run_dir)
            summary_css_url = "file://localhost%s/coh_PTF_html_summary.css" \
                              % run_dir
            summary_css = File(ifos, "coh_PTF_html_summary_css", sci_seg,
                               file_url=summary_css_url)
            summary_css.add_pfn(summary_css.cache_entry.path, site="local")
            file_list.extend(FileList([summary_css]))

        return file_list


def make_exttrig_file(cp, ifos, sci_seg, out_dir):
    '''
    Make an ExtTrig xml file containing information on the external trigger

    Parameters
    ----------
    cp : pycbc.workflow.configuration.WorkflowConfigParser object
    The parsed configuration options of a pycbc.workflow.core.Workflow.

    ifos : str
    String containing the analysis interferometer IDs.

    sci_seg : ligo.segments.segment
    The science segment for the analysis run.

    out_dir : str
    The output directory, destination for xml file.

    Returns
    -------
    xml_file : pycbc.workflow.File object
    The xml file with external trigger information.

    '''
    # Initialise objects
    xmldoc = ligolw.Document()
    xmldoc.appendChild(ligolw.LIGO_LW())
    tbl = lsctables.New(lsctables.ExtTriggersTable)
    cols = tbl.validcolumns
    xmldoc.childNodes[-1].appendChild(tbl)
    row = tbl.appendRow()

    # Add known attributes for this GRB
    setattr(row, "event_ra", float(cp.get("workflow", "ra")))
    setattr(row, "event_dec", float(cp.get("workflow", "dec")))
    setattr(row, "start_time", int(cp.get("workflow", "trigger-time")))
    setattr(row, "event_number_grb", str(cp.get("workflow", "trigger-name")))

    # Fill in all empty rows
    for entry in cols.keys():
        if hasattr(row, entry):
            continue
        if cols[entry] in ['real_4', 'real_8']:
            setattr(row, entry, 0.)
        elif cols[entry] in ['int_4s', 'int_8s']:
            setattr(row, entry, 0)
        elif cols[entry] == 'lstring':
            setattr(row, entry, '')
        elif entry == 'process_id':
            row.process_id = 0
        elif entry == 'event_id':
            row.event_id = 0
        else:
            raise ValueError("Column %s not recognized" % entry)

    # Save file
    xml_file_name = "triggerGRB%s.xml" % str(cp.get("workflow",
                                                    "trigger-name"))
    xml_file_path = os.path.join(out_dir, xml_file_name)
    utils.write_filename(xmldoc, xml_file_path)
    xml_file_url = urljoin("file:", pathname2url(xml_file_path))
    xml_file = File(ifos, xml_file_name, sci_seg, file_url=xml_file_url)
    xml_file.add_pfn(xml_file_url, site="local")

    return xml_file


def get_ipn_sky_files(workflow, file_url, tags=None):
    '''
    Retreive the sky point files for searching over the IPN error box and
    populating it with injections.

    Parameters
    ----------
    workflow: pycbc.workflow.core.Workflow
        An instanced class that manages the constructed workflow.
    file_url : string
        The URL of the IPN sky points file.
    tags : list of strings
        If given these tags are used to uniquely name and identify output files
        that would be produced in multiple calls to this function.

    Returns
    --------
    sky_points_file : pycbc.workflow.core.File
        File object representing the IPN sky points file.
    '''
    tags = tags or []
    file_attrs = {
        'ifos': workflow.ifos,
        'segs': workflow.analysis_time,
        'exe_name': "IPN_SKY_POINTS",
        'tags': tags
    }
    sky_points_file = resolve_url_to_file(file_url, attrs=file_attrs)

    return sky_points_file

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
        input_files = FileList([datafind_file for datafind_file in \
                                datafind_files if datafind_file.ifo == ifo])
        condition_strain_jobs = condition_strain_class(cp, "condition_strain",
                ifos=ifo, out_dir=outdir, tags=tags)
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


def setup_pygrb_pp_workflow(wf, pp_dir, seg_dir, segment, insp_files,
                            inj_files, inj_insp_files, inj_tags):
    """
    Generate post-processing section of PyGRB offline workflow
    """
    pp_outs = FileList([])
    # Begin setting up trig combiner job(s)
    # Select executable class and initialize
    exe_class = select_grb_pp_class(wf, "trig_combiner")
    job_instance = exe_class(wf.cp, "trig_combiner")
    # Create node for coherent no injections jobs
    node, trig_files = job_instance.create_node(wf.ifos, seg_dir, segment,
                                    insp_files, pp_dir)
    wf.add_node(node)
    pp_outs.append(trig_files)

    # Trig clustering for each trig file
    exe_class = select_grb_pp_class(wf, "trig_cluster")
    job_instance = exe_class(wf.cp, "trig_cluster")
    for trig_file in trig_files:
        # Create and add nodes
        node, out_file = job_instance.create_node(trig_file, pp_dir)
        wf.add_node(node)
        pp_outs.append(out_file)

    # Find injections from triggers
    exe_class = select_grb_pp_class(wf, "inj_finder")
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
                                           pp_dir)
        wf.add_node(node)
        inj_find_files.append(inj_find_file)
    pp_outs.append(inj_find_files)

    # Combine injections
    exe_class = select_grb_pp_class(wf, "inj_combiner")
    job_instance = exe_class(wf.cp, "inj_combiner")
    inj_comb_files = FileList([])
    for in_file in inj_find_files:
        if 'DETECTION' not in in_file.tags:
            node, inj_comb_file = job_instance.create_node(in_file,
                                                           pp_dir,
                                                           in_file.tags,
                                                           segment)
            wf.add_node(node)
            inj_comb_files.append(inj_comb_file)
    pp_outs.append(inj_comb_files)

    return pp_outs


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
                    out_dir, tags=None):
        node = Node(self)
        node.add_opt('--verbose')
        node.add_opt("--ifo-tag", ifo_tag)
        node.add_opt("--grb-name", self.trigger_name)
        node.add_opt("--trig-start-time", self.trig_start_time)
        node.add_opt("--segment-dir", seg_dir)
        node.add_input_list_opt("--input-files", insp_files)
        node.add_opt("--user-tag", "PYGRB")
        node.add_opt("--num-trials", self.num_trials)
        # Prepare output file tag
        user_tag = f"PYGRB_GRB{self.trigger_name}"
        if tags:
            user_tag += "_{}".format(tags)
        # Add on/off source and off trial outputs
        output_files = FileList([])
        outfile_types = ['ALL_TIMES', 'OFFSOURCE', 'ONSOURCE']
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

    def create_node(self, inj_files, inj_insp_files,
                    out_dir, tags=None):
        if tags is None:
            tags = []
        node = Node(self)
        node.add_input_list_opt('--input-files', inj_insp_files)
        node.add_input_list_opt('--inj-files', inj_files)
        ifo_tag, desc, segment = filename_metadata(inj_files[0].name)
        desc = '_'.join(desc.split('_')[:-1])
        out_name = "{}-{}_FOUNDMISSED-{}-{}.h5".format(
            ifo_tag, desc, segment[0], abs(segment))
        out_file = File(ifo_tag, 'inj_finder', segment,
                        os.path.join(out_dir, out_name), tags=tags)
        node.add_output(out_file)
        return node, out_file


class PycbcGrbInjCombinerExecutable(Executable):
    """The class responsible for creating jobs ``pycbc_grb_inj_combiner``
    """
    current_retention_level = Executable.ALL_TRIGGERS

    def __init__(self, cp, exe_name):
        super().__init__(cp=cp, name=exe_name)

    def create_node(self, input_file, out_dir, ifo_tag, segment, tags=None):
        if tags is None:
            tags = []
        node = Node(self)
        node.add_input_opt('--input-files', input_file)
        out_name = input_file.name.replace('.h5', '-FILTERED.h5')
        out_file = File(ifo_tag, 'inj_combiner', segment,
                        os.path.join(out_dir, out_name), tags=tags)
        node.add_output_opt('--output-file', out_file)
        return node, out_file
