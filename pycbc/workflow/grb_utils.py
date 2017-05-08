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
http://ligo-cbc.github.io/pycbc/latest/html/workflow.html
"""

from __future__ import print_function
import sys
import os
import shutil
import urlparse, urllib
from pycbc_glue import segments
from pycbc_glue.ligolw import ligolw, lsctables, utils, ilwd
from pycbc.workflow.core import File, FileList, resolve_url
from pycbc.workflow.jobsetup import select_generic_executable


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
            bank_veto.PFN(bank_veto.cache_entry.path, site="local")
            file_list.extend(FileList([bank_veto]))

        if summary_files:
            # summary.js file
            shutil.copy("%s/lalapps/src/ring/coh_PTF_config_files/" \
                        "coh_PTF_html_summary.js" % lalDir, "%s" % run_dir)
            summary_js_url = "file://localhost%s/coh_PTF_html_summary.js" \
                             % run_dir
            summary_js = File(ifos, "coh_PTF_html_summary_js", sci_seg,
                              file_url=summary_js_url)
            summary_js.PFN(summary_js.cache_entry.path, site="local")
            file_list.extend(FileList([summary_js]))

            # summary.css file
            shutil.copy("%s/lalapps/src/ring/coh_PTF_config_files/" \
                        "coh_PTF_html_summary.css" % lalDir, "%s" % run_dir)
            summary_css_url = "file://localhost%s/coh_PTF_html_summary.css" \
                              % run_dir
            summary_css = File(ifos, "coh_PTF_html_summary_css", sci_seg,
                               file_url=summary_css_url)
            summary_css.PFN(summary_css.cache_entry.path, site="local")
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

    sci_seg : glue.segments.segment
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
        if not hasattr(row, entry):
            if cols[entry] in ['real_4','real_8']:
                setattr(row,entry,0.)
            elif cols[entry] == 'int_4s':
                setattr(row,entry,0)
            elif cols[entry] == 'lstring':
                setattr(row,entry,'')
            elif entry == 'process_id':
                row.process_id = ilwd.ilwdchar("external_trigger:process_id:0")
            elif entry == 'event_id':
                row.event_id = ilwd.ilwdchar("external_trigger:event_id:0")
            else:
                print("Column %s not recognized" %(entry), file=sys.stderr)
                raise ValueError

    # Save file
    xml_file_name = "triggerGRB%s.xml" % str(cp.get("workflow",
                                                    "trigger-name"))
    xml_file_path = os.path.join(out_dir, xml_file_name)
    utils.write_filename(xmldoc, xml_file_path)
    xml_file_url = urlparse.urljoin("file:", urllib.pathname2url(xml_file_path))
    xml_file = File(ifos, xml_file_name, sci_seg, file_url=xml_file_url)
    xml_file.PFN(xml_file_url, site="local")
    
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
    ipn_sky_points = resolve_url(file_url)
    sky_points_url = urlparse.urljoin("file:",
            urllib.pathname2url(ipn_sky_points))
    sky_points_file = File(workflow.ifos, "IPN_SKY_POINTS",
            workflow.analysis_time, file_url=sky_points_url, tags=tags)
    sky_points_file.PFN(sky_points_url, site="local")

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
                ifo=ifo, out_dir=outdir, tags=tags)
        condition_strain_node, condition_strain_out = \
                condition_strain_jobs.create_node(input_files, tags=tags)
        condition_strain_nodes.append(condition_strain_node)
        condition_strain_outs.extend(FileList([condition_strain_out]))

    return condition_strain_nodes, condition_strain_outs
