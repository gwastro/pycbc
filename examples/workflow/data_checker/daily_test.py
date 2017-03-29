#!/usr/bin/env python

# Copyright (C) 2013 Ian W. Harry
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
Description of daily_test
"""

import pycbc
import pycbc.version
__author__  = "Ian Harry <ian.harry@astro.cf.ac.uk>"
__version__ = pycbc.version.git_verbose_msg
__date__    = pycbc.version.date
__program__ = "daily_test"

import os
import copy
import logging
import argparse
from pycbc_glue import segments
import pycbc.workflow as _workflow

logging.basicConfig(format='%(asctime)s:%(levelname)s : %(message)s', \
                    level=logging.INFO,datefmt='%I:%M:%S')
logger = logging.getLogger()

# command line options
_desc = __doc__[1:]
parser = argparse.ArgumentParser(description=_desc)
parser.add_argument('--version', action='version', version=__version__)
_workflow.add_workflow_command_line_group(parser)
args = parser.parse_args()

workflow = _workflow.Workflow(args, 'daily_check')
currDir = os.getcwd()
segDir = os.path.join(currDir,"segments")
dfDirSYR = os.path.join(currDir,"datafindSYR")
dfDirCIT = os.path.join(currDir,"datafindCIT")
dfDirLHO = os.path.join(currDir,"datafindLHO")
dfDirLLO = os.path.join(currDir,"datafindLLO")
dfDirUWM = os.path.join(currDir,"datafindUWM")

print "BEGIN BY GENERATING SCIENCE AND CAT_X VETOES"

def segment_report(sSegs):
    fullLen = 0
    fullNum = 0
    shortLen = 0
    shortNum = 0
    longLen = 0
    longNum = 0
    for ifo in sSegs.keys():
        for seg in sSegs[ifo]:
            fullLen += abs(seg)
            fullNum += 1
            if abs(seg) > 500:
                shortLen+=abs(seg)
                shortNum+=1
            if abs(seg) > 2000:
                longLen+=abs(seg)
                longNum+=1
        print "For ifo %s there is %d seconds of data in %d segments, %d seconds (%d unique segments) in segments longer than 500s and %d seconds (%d unique segments) longer than 2000s." %(ifo, fullLen, fullNum, shortLen, shortNum, longLen, longNum)


scienceSegs, segsList = _workflow.setup_segment_generation(workflow, segDir)

segment_report(scienceSegs)

print "STARTING DF"
print

# Start with SYR comparison
# FIXME: Used to use deecopy here, but now that seems to fail so repeating
#        segment query calls with logging off. This may be slow!
logger.disabled = True
scienceSegsS, _ = _workflow.setup_segment_generation(workflow, segDir)
logger.disabled = False
print "RUNNING DATAFIND FOR SYR"
datafinds, scienceSegsS = _workflow.setup_datafind_workflow(workflow, scienceSegsS,
                     dfDirSYR, segsList, tag="SYR")

segment_report(scienceSegsS)

print 
print
print "RUNNING DATAFIND FOR CIT"
logger.disabled = True
scienceSegsC, _ = _workflow.setup_segment_generation(workflow, segDir)
logger.disabled = False
datafinds, scienceSegsC = _workflow.setup_datafind_workflow(workflow, scienceSegsC,
                       dfDirCIT, segsList, tag="CIT")

segment_report(scienceSegsC)

print "Frames present a SYR and not at CIT:"
for ifo in scienceSegsS.keys():
    print "For ifo", ifo
    if ifo in scienceSegsC.keys():
       print (scienceSegsS[ifo] - scienceSegsC[ifo])
    else:
       print "No science segments for ifo %s at CIT" %(ifo)
    print

print "Frames present at CIT and not at SYR:"

for ifo in scienceSegsC.keys():
    print "For ifo", ifo
    if ifo in scienceSegsS.keys():
       print (scienceSegsC[ifo] - scienceSegsS[ifo])
    else:
       print "No science segments for ifo %s at SYR" %(ifo)
    print

# Next do LHO comparison

print
print "RUNNING DATAFIND FOR LHO"
logger.disabled = True
scienceSegsS, _ = _workflow.setup_segment_generation(workflow, segDir)
logger.disabled = False
datafinds, scienceSegsS = _workflow.setup_datafind_workflow(workflow, scienceSegsS,
                       dfDirLHO, segsList, tag="LHO")

segment_report(scienceSegsS)

print "Frames present at LHO and not at CIT:"
for ifo in scienceSegsS.keys():
    print "For ifo", ifo
    if ifo in scienceSegsC.keys():
       print (scienceSegsS[ifo] - scienceSegsC[ifo])
    else:
       print "No science segments for ifo %s at CIT" %(ifo)
    print

print "Frames present at CIT and not at LHO:"

for ifo in scienceSegsC.keys():
    print "For ifo", ifo
    if ifo in scienceSegsS.keys():
       print (scienceSegsC[ifo] - scienceSegsS[ifo])
    else:
       print "No science segments for ifo %s at LHO" %(ifo)
    print

# Next do LLO comparison

print
print "RUNNING DATAFIND FOR LLO"
logger.disabled = True
scienceSegsS, _ = _workflow.setup_segment_generation(workflow, segDir)
logger.disabled = False
datafinds, scienceSegsS = _workflow.setup_datafind_workflow(workflow, scienceSegsS,
                       dfDirLLO, segsList, tag="LLO")

segment_report(scienceSegsS)

print "Frames present at LLO and not at CIT:"
for ifo in scienceSegsS.keys():
    print "For ifo", ifo
    if ifo in scienceSegsC.keys():
       print (scienceSegsS[ifo] - scienceSegsC[ifo])
    else:
       print "No science segments for ifo %s at CIT" %(ifo)
    print

print "Frames present at CIT and not at LLO:"

for ifo in scienceSegsC.keys():
    print "For ifo", ifo
    if ifo in scienceSegsS.keys():
       print (scienceSegsC[ifo] - scienceSegsS[ifo])
    else:
       print "No science segments for ifo %s at LLO" %(ifo)
    print

# Next do UWM comparison

print
print "RUNNING DATAFIND FOR UWM"
logger.disabled = True
scienceSegsS, _ = _workflow.setup_segment_generation(workflow, segDir)
logger.disabled = False
datafinds, scienceSegsS = _workflow.setup_datafind_workflow(workflow, scienceSegsS,
                       dfDirUWM, segsList, tag="UWM")

segment_report(scienceSegsS)

print "Frames present at UWM and not at CIT:"
for ifo in scienceSegsS.keys():
    print "For ifo", ifo
    if ifo in scienceSegsC.keys():
       print (scienceSegsS[ifo] - scienceSegsC[ifo])
    else:
       print "No science segments for ifo %s at CIT" %(ifo)
    print

print "Frames present at CIT and not at UWM:"

for ifo in scienceSegsC.keys():
    print "For ifo", ifo
    if ifo in scienceSegsS.keys():
       print (scienceSegsC[ifo] - scienceSegsS[ifo])
    else:
       print "No science segments for ifo %s at UWM" %(ifo)
    print

