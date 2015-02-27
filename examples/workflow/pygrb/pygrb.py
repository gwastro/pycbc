#!/usr/bin/env python
# coding=utf-8
"""
Coherent inspiral pipeline development script
"""

import pycbc.version

__author__ = "Andrew Williamson <andrew.williamson@ligo.org>"
__version__ = pycbc.version.git_verbose_msg
__date__ = pycbc.version.date
__program__ = "pyGRB"

import shutil
import os
import argparse
import logging
import pycbc.workflow as _workflow
from glue.segments import segment
###############################################################################


def set_start_end(cp, a, b):
    """
    Function to update analysis boundaries as workflow is generated
    """
    cp.set("workflow", "start-time", str(a))
    cp.set("workflow", "end-time", str(b))
    return cp


def get_bank_veto(cp, run_dir):
    """
    Retrieve the bank_veto_bank.xml file needed to run coh_PTF_inspiral jobs.
    """
    if os.getenv("LAL_SRC") is None:
        raise ValueError("The environment variable LAL_SRC must be set to a "
                         "location containing the file lalsuite.git")
    else:
        lalDir = os.getenv("LAL_SRC")
        shutil.copy("%s/lalapps/src/ring/coh_PTF_config_files/" \
                    "bank_veto_bank.xml" % lalDir, "%s" % run_dir)
        sci_seg = segment(int(cp.get("workflow", "start-time")),
                          int(cp.get("workflow", "end-time")))
        bank_veto_url = "file://localhost%s/bank_veto_bank.xml" % run_dir
        bank_veto = _workflow.File("H1L1", "bank_veto_bank",
                                   sci_seg, file_url=bank_veto_url)
        bank_veto.PFN(bank_veto.cache_entry.path, site="local")
        return bank_veto


def make_cache(input_files, ifos, name, seg, cacheDir, tags=[]):
    """
    Place input_files into a cache file.
    """
    cache_file = _workflow.File(ifos, name, seg, extension="lcf",
                                directory=cacheDir, tags=tags)
    cache_file.PFN(cache_file.cache_entry.path, site="local")

    # Dump output to file
    fP = open(cache_file.storage_path, "w")
    for entry in input_files:
        start = str(int(entry.segment[0]))
        duration = str(int(abs(entry.segment)))
        print >> fP, "%s %s %s %s file://localhost%s" \
            %(ifos, entry.description.upper(), start, duration,
              entry.storage_path)
    fP.close()
    
    return cache_file


logging.basicConfig(format="%(asctime)s:%(levelname)s : %(message)s",
                    level=logging.INFO)

# Parse command line options and instantiate pycbc workflow object
parser = argparse.ArgumentParser()
parser.add_argument("--version", action="version", version=__version__)
parser.add_argument("-d", "--output-dir", default=None,
                    help="Path to output directory.")
_workflow.add_workflow_command_line_group(parser)
args = parser.parse_args()
wflow = _workflow.Workflow(args, "pygrb")
all_files = _workflow.FileList([])
tags = []

# Setup run directory
if args.output_dir:
    baseDir = args.output_dir
else:
    baseDir = os.getcwd()
runDir = os.path.join(baseDir, "GRB%s" % str(wflow.cp.get("workflow",
                                                          "trigger-name")))
if not os.path.exists(runDir):
    os.makedirs(runDir)
os.chdir(runDir)

# Hack to set start and end times based on maximum allowed duration
start = int(wflow.cp.get("workflow", "trigger-time")) - int(wflow.cp.get(
            "workflow-exttrig_segments", "max-duration"))
end = int(wflow.cp.get("workflow", "trigger-time")) + int(wflow.cp.get(
    "workflow-exttrig_segments", "max-duration"))
wflow.cp = set_start_end(wflow.cp, start, end)


# Retrieve segments ahope-style
currDir = os.getcwd()
segDir = os.path.join(currDir, "segments")
sciSegs, segsFileList = _workflow.setup_segment_generation(wflow, segDir)

# Make coherent network segments
onSrc, sciSegs = _workflow.get_triggered_coherent_segment(wflow, segDir,
                                                          sciSegs)
# FIXME: The following two lines are/were crude hacks.
ifo = sciSegs.keys()[0]
wflow.analysis_time = sciSegs[ifo][0]
wflow.cp = set_start_end(wflow.cp, sciSegs[ifo][0][0], sciSegs[ifo][0][1])

# Datafind
dfDir = os.path.join(currDir, "datafind")
datafind_files, sciSegs = _workflow.setup_datafind_workflow(wflow, sciSegs,
                                                            dfDir,
                                                            segsFileList)

# If using coh_PTF_inspiral we need bank_veto_bank.xml
if os.path.basename(wflow.cp.get("executables", "inspiral")) \
                    == "lalapps_coh_PTF_inspiral":
    datafind_veto_files = _workflow.FileList([])
    bank_veto_file = _workflow.FileList([get_bank_veto(wflow.cp, runDir)])
    datafind_veto_files.extend(datafind_files)
    datafind_veto_files.extend(bank_veto_file)

all_files.extend(datafind_files)

# Template bank and splitting the bank
# TODO: Move from pregenerated to generated coherent network bank
bank_files = _workflow.setup_tmpltbank_workflow(wflow, sciSegs,
                                                datafind_files, dfDir)
splitbank_files = _workflow.setup_splittable_workflow(wflow, bank_files, dfDir)
all_files.extend(bank_files)
all_files.extend(splitbank_files)

# Injections
if wflow.cp.has_section("inspiral-inj"):
    injDir = os.path.join(currDir, "inj_files")
    inj_files, inj_tags = _workflow.setup_injection_workflow(wflow,
                                                             output_dir=injDir)
    all_files.extend(inj_files)
else:
    inj_files = None

# Matched-filtering
# TODO: Write coherent matched filtering code
inspDir = os.path.join(currDir, "inspiral")
inspiral_files = _workflow.setup_matchedfltr_workflow(wflow, sciSegs,
                                                      datafind_veto_files,
                                                      splitbank_files, inspDir)
all_files.extend(inspiral_files)

# Post processing
ppDir = os.path.join(currDir, "post_processing")
ifos = ''.join(sciSegs.keys())
post_proc_method = wflow.cp.get_opt_tags("workflow-postproc",
                                         "postproc-method", tags)
if post_proc_method == "COH_PTF_PP_DAG_WORKFLOW":
    pp_config_file_name = wflow.cp.get("workflow-postproc", "config-file")
    pp_config_file_url = "file://localhost%s/%s" % (baseDir,
                                                    pp_config_file_name)
    pp_config_file = _workflow.File(ifos, pp_config_file_name,
                                    wflow.analysis_time,
                                    file_url=pp_config_file_url)
    pp_config_file.PFN(pp_config_file.cache_entry.path, site="local")
    inspiral_files.extend(_workflow.FileList([pp_config_file]))
    pp_files = _workflow.setup_coh_PTF_post_processing(wflow, inspiral_files,
                                                       inj_files, ppDir,
                                                       segDir, run_dir=runDir,
                                                       ifos=ifos)

elif post_proc_method == "COH_PTF_WORKFLOW":
    inspiral_cache = make_cache(inspiral_files, ifos, "INSPIRAL",
                                sciSegs[ifo][0], inspDir)
    inspiral_files.extend(_workflow.FileList([inspiral_cache]))
    pp_files = _workflow.setup_coh_PTF_post_processing(wflow, inspiral_files,
                                                       inj_files, ppDir,
                                                       segDir, ifos=ifos)

all_files.extend(pp_files)

# Compile workflow and write DAX
wflow.save()
logging.info("Written dax.")
