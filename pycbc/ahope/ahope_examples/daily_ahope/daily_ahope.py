import os
import logging
from glue import pipeline
from glue import segments
import pycbc.ahope as ahope

logging.basicConfig(format='%(asctime)s:%(levelname)s : %(message)s', \
                    level=logging.INFO,datefmt='%I:%M:%S')

# Parse ini file
cp = ahope.parse_ahope_ini_file('./daily_ahope.ini')

# Make directories for output
currDir = os.getcwd()
segDir = os.path.join(currDir,"segments")
if not os.path.exists(segDir+'/logs'):
    os.makedirs(segDir+'/logs')
dfDir = os.path.join(currDir,"datafind")
if not os.path.exists(dfDir+'/logs'):
    os.makedirs(dfDir+'/logs')
tmpltbankDir = os.path.join(currDir,"tmpltbank")
if not os.path.exists(tmpltbankDir+'/logs'):
    os.makedirs(tmpltbankDir+'/logs')
inspiralDir = os.path.join(currDir,"inspiral")
if not os.path.exists(inspiralDir+'/logs'):
    os.makedirs(inspiralDir+'/logs')

# Set start and end times
start_time = 961545543
end_time = 961594487
ifos = ['H1','L1']
# ALEX: please add V1 to this and commit to repo!

# Initialize the dag
basename = 'ahope_test'
logfile = 'asdtwe.log'
fh = open( logfile, "w" )
fh.close()
dag = pipeline.CondorDAG(logfile,dax=False)
dag.set_dax_file(basename)
dag.set_dag_file(basename)

# Get segments
scienceSegs, segsDict = ahope.setup_segment_generation(cp, ifos, start_time,\
                               end_time, None, segDir)

# Get frames, this can be slow, as we ping every frame to check it exists,
# the second option shows how to turn this off
#datafinds, scienceSegs = ahope.setup_datafind_workflow(cp, scienceSegs, dag,\
#                       dfDir)
# This second case will also update the segment list on missing data, not fail
datafinds, scienceSegs = ahope.setup_datafind_workflow(cp, scienceSegs, dag,\
                       dfDir, checkSegmentGaps=False, checkFramesExist=False,\
                       updateSegmentTimes=True)

# Template bank stuff
banks = ahope.setup_tmpltbank_workflow(cp, scienceSegs, datafinds, dag, \
                                       tmpltbankDir)
# Split bank up
splitBanks = ahope.setup_splittable_workflow(cp, dag, banks, tmpltbankDir)
# Do matched-filtering
insps = ahope.setup_matchedfltr_workflow(cp, scienceSegs, datafinds, dag, \
                                         splitBanks, inspiralDir)

dag.write_sub_files()
dag.write_dag()
#dag.write_abstract_dag()
dag.write_script()
