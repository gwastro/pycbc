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
Programe for running a single detector ahope analysis up to matched-filtering.
This is designed to mimic the current behaviour of dail_ihope.
"""
import pycbc
import pycbc.version
__author__  = "Ian Harry <ian.harry@astro.cf.ac.uk>"
__version__ = pycbc.version.git_verbose_msg
__date__    = pycbc.version.date
__program__ = "daily_ahope"

import os
import copy
import logging
import urlparse
import optparse
import lal
from glue import pipeline
from glue import segments
import pycbc.ahope as ahope

logging.basicConfig(format='%(asctime)s:%(levelname)s : %(message)s', \
                    level=logging.INFO,datefmt='%I:%M:%S')

# command line options
usage = """usage: %prog [options]"""
_desc = __doc__[1:]
parser = optparse.OptionParser(usage, version=__version__, description=_desc)
parser.add_option("-s", "--start-time", type="int",\
                  help="Time to start analysis from.")
parser.add_option("-i", "--config-file", help="Location of ini file.")
parser.add_option("-d", "--output-dir", help="Path to output directory.")
(opts,args) = parser.parse_args()

# Add check that all of these exist
if not opts.start_time:
    parser.error("Must supply --start-time")
if not opts.config_file:
    parser.error("Must supply --config-file")
if not opts.output_dir:
    parser.error("Must supply --output-dir")

workflow = ahope.Workflow(opts.config_file)
# Get dates and stuff
# This feels hacky!
yestDate = lal.GPSToUTC(opts.start_time)
yestMnight = copy.deepcopy(yestDate)
yestMnight[3] = 0
yestMnight[4] = 0
yestMnight[5] = 0
yestMnightGPS = lal.UTCToGPS(yestMnight)

monthName = '%04d%02d' %(yestMnight[0], yestMnight[1])
dayName = '%04d%02d%02d' %(yestMnight[0], yestMnight[1], yestMnight[2])

workingDir = os.path.join(opts.output_dir, monthName, dayName)

if not os.path.exists(workingDir):
    os.makedirs(workingDir)

os.chdir(workingDir)

if not os.path.exists('log'):
    os.makedirs('log')
if not os.path.exists('logs'):
    os.makedirs('logs')

# THIS CONCLUDES EVERYTHING THAT THE SETUP SCRIPT DID BEFORE RUNNING DAG

# Set start and end times need some padding to assume the whole day
# is analysed
pad_time = 72
start_time = yestMnightGPS - pad_time
end_time = start_time + 60*60*24 + 2*pad_time

# Set the ifos to analyse
ifos = []
for ifo in workflow.cp.options('ahope-ifos'):
    ifos.append(ifo.upper())

# Get segments
scienceSegs, segsDict = ahope.setup_segment_generation(workflow, ifos,
                               start_time, end_time, workingDir, maxVetoCat=4,
                               minSegLength=2000, 
                               generate_coincident_segs=False)

# Get frames, this can be slow, as we ping every frame to check it exists,
# the second option shows how to turn this off
#datafinds, scienceSegs = ahope.setup_datafind_workflow(cp, scienceSegs, dag,\
#                         dfDir)
# This second case will also update the segment list on missing data, not fail
datafinds, scienceSegs = ahope.setup_datafind_workflow(workflow, scienceSegs,
                           workingDir, checkFramesExist='update_times',\
                           checkSegmentGaps='update_times')

# Template bank stuff
banks = ahope.setup_tmpltbank_workflow(workflow, scienceSegs, datafinds, 
                                       workingDir)
# Split bank up
splitBanks = ahope.setup_splittable_workflow(workflow, banks, workingDir)
# Do matched-filtering
insps = ahope.setup_matchedfltr_workflow(workflow, scienceSegs, datafinds,
                                         splitBanks, workingDir)


# Now I start doing things not supported by ahope at present.
# NOTE: llwadd was moved over to ahope functionality

# First we need direct access to the dag and cp objects
dag = workflow.dag
cp = workflow.cp

basename = 'daily_ahope'
# Set up condor jobs for this stuff
llwadd_exe = ahope.LigolwAddExec('llwadd')
llwadd_job = llwadd_exe.create_job(cp, ''.join(scienceSegs.keys()),
                                   out_dir=workingDir)

cp_exec = ahope.Executable('cp')
cp_job = cp_exec.create_job(cp, out_dir=workingDir)

# Hopefully with tags, I would need only one exec and two jobs
si_exec_coarse = ahope.Executable('siclustercoarse')
si_exec_fine = ahope.Executable('siclusterfine')
si_job_coarse = si_exec_coarse.create_job(cp, out_dir=workingDir)
si_job_fine = si_exec_fine.create_job(cp, out_dir=workingDir)

inspstr = 'INSPIRAL'
pageDagParents = []
for inspOutGroup in insps:
    ifo = inspOutGroup.ifo
    analysis_seg = inspOutGroup.segment

    # Create a cache file to hole the input to ligolw_add
    llwadd_node = llwadd_job.create_node(analysis_seg, [inspOutGroup]) 
    # HACK TO OVERWRITE FILE NAMING
    llwadd_node.output_files = ahope.AhopeFileList([])
    output_name = '%s-INSPIRAL_UNCLUSTERED-%d-%d.xml.gz'\
                   %(ifo, analysis_seg[0], abs(analysis_seg))
    output_url = urlparse.urlunparse(['file', 'localhost', 
                                      os.path.join(workingDir,output_name),
                                      None, None, None])
    llwaddFile = ahope.AhopeFile(ifo, 'LLWADD_UNCLUSTERED', analysis_seg,
                              file_url=output_url)
    llwaddFile.node = llwadd_node
    llwadd_node.add_output(llwaddFile, opt='output')
    llwadd_node.add_var_opt('lfn-start-time', analysis_seg[0])
    llwadd_node.add_var_opt('lfn-end-time',analysis_seg[1])
    workflow.add_node(llwadd_node)

    # Finally run 30ms and 16s clustering on the combined files
    clustered_30ms_name = output_name.replace('UNCLUSTERED',\
                                              '30MILLISEC_CLUSTERED')
    clustered_30ms_url = urlparse.urlunparse(['file', 'localhost',
                                 os.path.join(workingDir, clustered_30ms_name),
                                 None, None, None])
    clustered_30ms_file = ahope.AhopeFile(ifo, 'LLWADD_30MS_CLUSTERED',
                            analysis_seg, file_url=clustered_30ms_url)
    clustered_16s_name  = output_name.replace('UNCLUSTERED', '16SEC_CLUSTERED')
    clustered_16s_url = urlparse.urlunparse(['file', 'localhost',
                                 os.path.join(workingDir, clustered_16s_name),
                                 None, None, None])
    clustered_16s_file = ahope.AhopeFile(ifo, 'LLWADD_16S_CLUSTERED',
                            analysis_seg, file_url=clustered_16s_url)
 

    for cfile in [clustered_30ms_file, clustered_16s_file]:
        cpnode = cp_job.create_node()
        cpnode.add_input(llwaddFile, argument=True)
        cpnode.add_output(cfile, argument=True)
        workflow.add_node(cpnode)

        if cfile == clustered_16s_file:
            sinode = si_job_coarse.create_node()
        else:
            sinode = si_job_fine.create_node()

        # FIXME: this node overwrites the input file. Better
        # that this take command line options, remove the cp job and write to
        # a different file
        sinode.add_input(cfile, argument=True)
        workflow.add_node(sinode)
        pageDagParents.append(sinode)

# Now we construct the page_conf.txt for the daily_ihope_page code
pageConts = []
pageConts.append('# Where I live')
pageConts.append('ihope_daily_page = %s' \
                 %(cp.get('executables', 'ihope_daily_page')) )
pageConts.append('# The location and configuration for the glitch page')
pageConts.append('ligolw_cbc_glitch_page = %s' \
                 %(cp.get('executables', 'cbc_glitch_page')) )
pageConts.append('omega_conf = %s' %(cp.get('ahope-omega','omega-conf-file')) )
pageConts.append('omega_frame_dir = %s' \
                 %(cp.get('ahope-omega','omega-frame-dir')) )
pageConts.append('# Location and configuration for the hw injection page')
pageConts.append('ligolw_cbc_hardware_inj_page = %s' \
                 %(cp.get('executables', 'cbc_hardware_inj_page')) )
pageConts.append('hwinj_file = %s' \
                 %(cp.get('ahope-hwinj', 'hwinj-file')) )
pageConts.append('# Location of asset files (.css, .js, etc)')
pageConts.append('asset_dir = %s' \
                 %(cp.get('ahope', 'ahope-asset-dir')) )
pageConts.append('# Veto definer file')
pageConts.append('veto_definer_file = %s' \
                 %(cp.get('ahope-segments','segments-veto-definer-file')) )
pageConts.append('# source dir with triggers')
pageConts.append('trigger_dir = %s' %(workingDir)) 
pageConts.append('# temp directory')
pageConts.append('tmp_dir = %s' %(workingDir))
pageConts.append('# target directory')
htmlBaseDir = cp.get('ahope','ahope-html-basedir')
htmlOutDir = os.path.join(htmlBaseDir, monthName, dayName)
pageConts.append('out_dir = %s' %(htmlOutDir))
pageText = '\n'.join(pageConts)
pageConfFile = os.path.join(workingDir, 'page_conf.txt')
pageConfFP = open(pageConfFile, 'w')
pageConfFP.write(pageText)
pageConfFP.close()

# Run the code to make the daily page dag
ihopePageCmd = []
ihopePageCmd.append(cp.get('executables','ihope_daily_page'))
ihopePageCmd.append('-a')
ihopePageCmd.append('make_dag')
ihopePageCmd.append('-f')
ihopePageCmd.append('%s' %(pageConfFile))
ihopePageCmd.append('-s')
ihopePageCmd.append('%s' %(str(start_time)))
ihopePageCmd.append('-i')
ihopePageCmd.append(','.join(ifos))
ahope.make_external_call(ihopePageCmd, outDir=os.path.join(workingDir,'logs'),\
                         outBaseName='daily_ihope_page_daggen')

# Add this to the workflow and make it a child of all cluster jobs
dailyPageJob = pipeline.CondorDAGManJob('daily_page.dag', workingDir,\
                                        'daily_page.dax') 
dailyPageNode = dailyPageJob.create_node()
for job in pageDagParents:
    dailyPageNode.add_parent(job)
dag.add_node(dailyPageNode)

# One final job to make the output page
# Make sub file for summary page job
summJob = pipeline.CondorDAGJob('vanilla', \
                            cp.get('executables', 'ihope_daily_page'))
summJob.set_stderr_file('logs/summ-page-$(cluster)-$(process).err')
summJob.set_stdout_file('logs/summ-page-$(cluster)-$(process).out')
summJob.set_sub_file('summary_page.sub')
summJob.add_condor_cmd('getenv', 'True')

summNode = summJob.create_node()
summNode.add_var_opt('action', 'make_index_page')
summNode.add_var_opt('config', pageConfFile)
summNode.add_var_opt('gps-start-time', start_time)
summNode.add_var_opt('ifos', ','.join(ifos))
summNode.add_parent(dailyPageNode)
dag.add_node(summNode)

workflow.write_plans()

logging.info("Finished.")
