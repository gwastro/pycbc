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
import optparse
import lal
from glue import pipeline
from glue import segments
import pycbc.ahope as ahope

# Force vanilla universe in legacy code
#ahope.legacy_ihope_job_utils.condorUniverse='vanilla'

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
dag = workflow.dag
cp = workflow.cp
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
for ifo in cp.options('ahope-ifos'):
    ifos.append(ifo.upper())

# Get segments
scienceSegs, segsDict = ahope.setup_segment_generation(workflow, ifos, start_time,
                               end_time, workingDir, maxVetoCat=4,
                               minSegLength=2000)

# Get frames, this can be slow, as we ping every frame to check it exists,
# the second option shows how to turn this off
#datafinds, scienceSegs = ahope.setup_datafind_workflow(cp, scienceSegs, dag,\
#                         dfDir)
# This second case will also update the segment list on missing data, not fail
datafinds, scienceSegs = ahope.setup_datafind_workflow(workflow, scienceSegs,
                           workingDir, updateSegmentTimes=True,
                           checkFramesExist=True, checkSegmentGaps=False)

# Template bank stuff
banks = ahope.setup_tmpltbank_workflow(workflow, scienceSegs, datafinds, 
                                       workingDir)
# Split bank up
splitBanks = ahope.setup_splittable_workflow(workflow, banks, workingDir)
# Do matched-filtering
insps = ahope.setup_matchedfltr_workflow(workflow, scienceSegs, datafinds,
                                         splitBanks, workingDir)


# Now I start doing things not supported by ahope at present.
basename = 'daily_ahope'
# Set up condor jobs for this stuff
cp.add_section('condor')
# Hack to avoid hardcoding in glue
cp.set('condor','ligolw_add', cp.get('executables', 'ligolw_add'))
lwadd_job = pipeline.LigolwAddJob('logs',cp)
lwadd_job.add_condor_cmd('priority', '20')
lwadd_job.add_opt('add-lfn-table','')
lwadd_job.set_sub_file( basename + '_' + ifo + '.ligolwadd.sub' )

cp_job = pipeline.CondorDAGJob('vanilla','/bin/cp')
cp_job.set_stderr_file('logs/cp-$(cluster)-$(process).err')
cp_job.set_stdout_file('logs/cp-$(cluster)-$(process).out')
cp_job.set_sub_file('cp.sub')

si_job_coarse = ahope.Job(cp, 'siclustercoarse', 'vanilla')
si_job_coarse.add_condor_cmd('getenv','True')

si_job_fine = ahope.Job(cp, 'siclusterfine', 'vanilla')
si_job_fine.add_condor_cmd('getenv','True')

inspstr = 'INSPIRAL'
pageDagParents = []
for inspOutGroup in insps:
    ifo = inspOutGroup.ifo
    analysis_seg = inspOutGroup.segment
    fileList = inspOutGroup.paths
    jobList = [f.job for f in fileList]

    # Create a cache file to hole the input to ligolw_add
    out_basename = ifo + '-' + inspstr + '-' + str(analysis_seg[0]) + '-' 
    out_basename += str(abs(analysis_seg))
    insp_cache_file_name = out_basename + '.cache'
    insp_cache_file = open(insp_cache_file_name, 'w')
    fileList.tofile(insp_cache_file)
    insp_cache_file.close()

    # use ligolw_add to join everything back together
    output_name = '%s-INSPIRAL_UNCLUSTERED-%d-%d.xml.gz'\
                   %(ifo, analysis_seg[0], abs(analysis_seg))
    lwadd = pipeline.LigolwAddNode(lwadd_job)
    lwadd.add_var_opt('input-cache',insp_cache_file_name)
    lwadd.set_output(output_name)
    lwadd.add_var_opt('lfn-start-time',analysis_seg[0])
    lwadd.add_var_opt('lfn-end-time',analysis_seg[1])
    for job in jobList:
        lwadd.add_parent(job)
    output_file = lwadd.get_output()
    dag.add_node(lwadd)

    # Finally run 30ms and 16s clustering on the combined files
    clustered_30ms_name = output_name.replace('UNCLUSTERED',\
                                              '30MILLISEC_CLUSTERED')
    clustered_16s_name  = output_name.replace('UNCLUSTERED', '16SEC_CLUSTERED')

    for cname in [clustered_30ms_name, clustered_16s_name]:
        cpnode = pipeline.CondorDAGNode(cp_job)
        cpnode.add_file_arg(output_name)
        cpnode.add_file_arg(cname)
        cpnode.add_parent(lwadd)
        dag.add_node(cpnode)

        if cname == clustered_16s_name:
            sinode = ahope.Node(si_job_coarse)
        else:
            sinode = ahope.Node(si_job_fine)

        sinode.add_file_arg(cname)
        sinode.add_parent(cpnode)
        dag.add_node(sinode)
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
dailyPageNode = pipeline.CondorDAGManNode(dailyPageJob)
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

summNode = pipeline.CondorDAGNode(summJob)
summNode.add_var_opt('action', 'make_index_page')
summNode.add_var_opt('config', pageConfFile)
summNode.add_var_opt('gps-start-time', start_time)
summNode.add_var_opt('ifos', ','.join(ifos))
summNode.add_parent(dailyPageNode)
dag.add_node(summNode)

workflow.write_plan()

logging.info("Finished.")
