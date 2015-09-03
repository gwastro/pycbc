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
https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/ahope.html
"""

import os
import shutil
from glue import segments
from pycbc.workflow.core import File, FileList


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

