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

import pycbc.version

import pycbc

__author__ = "Ian Harry <ian.harry@astro.cf.ac.uk>"
__version__ = pycbc.version.git_verbose_msg
__date__ = pycbc.version.date
__program__ = "daily_test"

import argparse
import logging
import os

import pycbc.workflow as _workflow

logging.basicConfig(
    format="%(asctime)s:%(levelname)s : %(message)s",
    level=logging.INFO,
    datefmt="%I:%M:%S",
)
logger = logging.getLogger()

# command line options
_desc = __doc__[1:]
parser = argparse.ArgumentParser(description=_desc)
parser.add_argument("--version", action="version", version=__version__)
_workflow.add_workflow_command_line_group(parser)
args = parser.parse_args()

workflow = _workflow.Workflow(args, "daily_check")
curr_dir = os.getcwd()
seg_dir = os.path.join(curr_dir, "segments")
df_dir_syr = os.path.join(curr_dir, "datafindSYR")
df_dir_cit = os.path.join(curr_dir, "datafindCIT")
df_dir_lho = os.path.join(curr_dir, "datafindLHO")
df_dir_llo = os.path.join(curr_dir, "datafindLLO")
df_dir_uwm = os.path.join(curr_dir, "datafindUWM")

print("BEGIN BY GENERATING SCIENCE AND CAT_X VETOES")


def segment_report(s_segs):
    full_len = 0
    full_num = 0
    short_len = 0
    short_num = 0
    long_len = 0
    long_num = 0
    for ifo in s_segs.keys():
        for seg in s_segs[ifo]:
            full_len += abs(seg)
            full_num += 1
            if abs(seg) > 500:
                short_len += abs(seg)
                short_num += 1
            if abs(seg) > 2000:
                long_len += abs(seg)
                long_num += 1
        print(
            f"For ifo {ifo} there is {int(full_len)} seconds of data in {full_num} "
            f"segments, {int(short_len)} seconds ({short_num} unique segments) in "
            f"segments longer than 500s and {int(long_len)} seconds ({long_num} unique "
            "segments) longer than 2000s."
        )


science_segs, segs_list = _workflow.setup_segment_generation(workflow, seg_dir)

segment_report(science_segs)

print("STARTING DF")
print()

# Start with SYR comparison
# FIXME: Used to use deecopy here, but now that seems to fail so repeating
#        segment query calls with logging off. This may be slow!
logger.disabled = True
science_segs_s, _ = _workflow.setup_segment_generation(workflow, seg_dir)
logger.disabled = False
print("RUNNING DATAFIND FOR SYR")
datafinds, science_segs_s = _workflow.setup_datafind_workflow(
    workflow, science_segs_s, df_dir_syr, segs_list, tag="SYR"
)

segment_report(science_segs_s)

print()
print()
print("RUNNING DATAFIND FOR CIT")
logger.disabled = True
science_segs_c, _ = _workflow.setup_segment_generation(workflow, seg_dir)
logger.disabled = False
datafinds, science_segs_c = _workflow.setup_datafind_workflow(
    workflow, science_segs_c, df_dir_cit, segs_list, tag="CIT"
)

segment_report(science_segs_c)

print("Frames present a SYR and not at CIT:")
for ifo in science_segs_s.keys():
    print("For ifo", ifo)
    if ifo in science_segs_c.keys():
        print(science_segs_s[ifo] - science_segs_c[ifo])
    else:
        print(f"No science segments for ifo {ifo} at CIT")
    print()

print("Frames present at CIT and not at SYR:")

for ifo in science_segs_c.keys():
    print("For ifo", ifo)
    if ifo in science_segs_s.keys():
        print(science_segs_c[ifo] - science_segs_s[ifo])
    else:
        print(f"No science segments for ifo {ifo} at SYR")
    print()

# Next do LHO comparison

print()
print("RUNNING DATAFIND FOR LHO")
logger.disabled = True
science_segs_s, _ = _workflow.setup_segment_generation(workflow, seg_dir)
logger.disabled = False
datafinds, science_segs_s = _workflow.setup_datafind_workflow(
    workflow, science_segs_s, df_dir_lho, segs_list, tag="LHO"
)

segment_report(science_segs_s)

print("Frames present at LHO and not at CIT:")
for ifo in science_segs_s.keys():
    print("For ifo", ifo)
    if ifo in science_segs_c.keys():
        print(science_segs_s[ifo] - science_segs_c[ifo])
    else:
        print(f"No science segments for ifo {ifo} at CIT")
    print()

print("Frames present at CIT and not at LHO:")

for ifo in science_segs_c.keys():
    print("For ifo", ifo)
    if ifo in science_segs_s.keys():
        print(science_segs_c[ifo] - science_segs_s[ifo])
    else:
        print(f"No science segments for ifo {ifo} at LHO")
    print()

# Next do LLO comparison

print()
print("RUNNING DATAFIND FOR LLO")
logger.disabled = True
science_segs_s, _ = _workflow.setup_segment_generation(workflow, seg_dir)
logger.disabled = False
datafinds, science_segs_s = _workflow.setup_datafind_workflow(
    workflow, science_segs_s, df_dir_llo, segs_list, tag="LLO"
)

segment_report(science_segs_s)

print("Frames present at LLO and not at CIT:")
for ifo in science_segs_s.keys():
    print("For ifo", ifo)
    if ifo in science_segs_c.keys():
        print(science_segs_s[ifo] - science_segs_c[ifo])
    else:
        print(f"No science segments for ifo {ifo} at CIT")
    print()

print("Frames present at CIT and not at LLO:")

for ifo in science_segs_c.keys():
    print("For ifo", ifo)
    if ifo in science_segs_s.keys():
        print(science_segs_c[ifo] - science_segs_s[ifo])
    else:
        print(f"No science segments for ifo {ifo} at LLO")
    print()

# Next do UWM comparison

print()
print("RUNNING DATAFIND FOR UWM")
logger.disabled = True
science_segs_s, _ = _workflow.setup_segment_generation(workflow, seg_dir)
logger.disabled = False
datafinds, science_segs_s = _workflow.setup_datafind_workflow(
    workflow, science_segs_s, df_dir_uwm, segs_list, tag="UWM"
)

segment_report(science_segs_s)

print("Frames present at UWM and not at CIT:")
for ifo in science_segs_s.keys():
    print("For ifo", ifo)
    if ifo in science_segs_c.keys():
        print(science_segs_s[ifo] - science_segs_c[ifo])
    else:
        print(f"No science segments for ifo {ifo} at CIT")
    print()

print("Frames present at CIT and not at UWM:")

for ifo in science_segs_c.keys():
    print("For ifo", ifo)
    if ifo in science_segs_s.keys():
        print(science_segs_c[ifo] - science_segs_s[ifo])
    else:
        print(f"No science segments for ifo {ifo} at UWM")
    print()
