#!/usr/bin/env python
# coding=utf-8
"""
Coherent inspiral pipeline development script
"""

import pycbc.version

__author__ = 'Andrew Williamson <andrew.williamson@ligo.org>'
__version__ = pycbc.version.git_verbose_msg
__date__ = pycbc.version.date
__program__ = "pyGRB"

import os
import argparse
import logging
import pycbc.workflow as _workflow
###############################################################################

logging.basicConfig(format='%(asctime)s:%(levelname)s : %(message)s',
                    level=logging.INFO)

# Parse command line options and instantiate pycbc workflow object
parser = argparse.ArgumentParser()
parser.add_argument('--version', action='version', version=__version__)
parser.add_argument("-d", "--output-dir", default=None,
                    help="Path to output directory.")
_workflow.add_workflow_command_line_group(parser)
args = parser.parse_args()
wflow = _workflow.Workflow(args, 'pygrb')
all_files = _workflow.FileList([])
tags = []

if args.output_dir:
    baseDir = args.output_dir
else:
    baseDir = os.getcwd()
runDir = os.path.join(baseDir, 'GRB%s' % str(wflow.cp.get('workflow',
                                                          'trigger-name')))
if not os.path.exists(runDir):
    os.makedirs(runDir)
os.chdir(runDir)

# Hack to set start and end times based on maximum allowed duration
start = int(wflow.cp.get('workflow', 'trigger-time')) - int(wflow.cp.get(
            'workflow-exttrig_segments', 'max-duration'))
end = int(wflow.cp.get('workflow', 'trigger-time')) + int(wflow.cp.get(
    'workflow-exttrig_segments', 'max-duration'))
wflow.cp.set('workflow', 'start-time', str(start))
wflow.cp.set('workflow', 'end-time', str(end))

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

# Datafind
dfDir = os.path.join(currDir, "datafind")
datafind_files, sciSegs = _workflow.setup_datafind_workflow(wflow, sciSegs,
                                                            dfDir,
                                                            segsFileList)
all_files.extend(datafind_files)

# Template bank and splitting the bank
# TODO: Move from pregenerated to generated coherent network bank
bank_files = _workflow.setup_tmpltbank_workflow(wflow, sciSegs,
                                                datafind_files, dfDir)
splitbank_files = _workflow.setup_splittable_workflow(wflow, bank_files, dfDir)
all_files.extend(bank_files)
all_files.extend(splitbank_files)

# Injections
"""
injDir = os.path.join(currDir, "inj_files")
inj_files, inj_tags = ahope.setup_injection_workflow(wflow,
                                                     output_dir=injDir)
all_files.extend(inj_files)
"""

# Matched-filtering
# TODO: Write coherent matched filtering code
inspDir = os.path.join(currDir, "inspiral")
inspiral_files = _workflow.setup_matchedfltr_workflow(wflow, sciSegs,
                                                      datafind_files,
                                                      splitbank_files, inspDir)
all_files.extend(inspiral_files)

# Compile workflow and write DAX
wflow.save()
logging.info("Written dax.")
