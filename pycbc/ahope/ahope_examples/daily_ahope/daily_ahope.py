import os
import logging
from glue import pipeline
from glue import segments
import pycbc.ahope as ahope

logging.basicConfig(format='%(asctime)s:%(levelname)s : %(message)s', \
                    level=logging.INFO,datefmt='%I:%M:%S')

workflow = ahope.Workflow('./daily_ahope.ini')

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

# These are Chris' example S6 day of data
start_time = 961545543
end_time = 961594487

# This is S6D chunk 3, an example full-scale ihope analysis time-frame
start_time = 968543943
end_time = 971622087

# Set the ifos to analyse
ifos = ['H1','L1']
# ALEX: please add V1 to this and commit to repo!

# Get segments
scienceSegs, segsDict = ahope.setup_segment_generation(workflow, ifos, start_time,
                               end_time, segDir, maxVetoCat=5,
                               minSegLength=2000)

# Get frames, this can be slow, as we ping every frame to check it exists,
# the second option shows how to turn this off
#datafinds, scienceSegs = ahope.setup_datafind_workflow(cp, scienceSegs, dag,\
#                         dfDir)
# This second case will also update the segment list on missing data, not fail
datafinds, scienceSegs = ahope.setup_datafind_workflow(workflow, scienceSegs,
                       dfDir, checkSegmentGaps=False, checkFramesExist=False,
                       updateSegmentTimes=True)

# Template bank stuff
banks = ahope.setup_tmpltbank_workflow(cp, scienceSegs, datafinds, dag, 
                                       tmpltbankDir)
# Split bank up
splitBanks = ahope.setup_splittable_workflow(cp, dag, banks, tmpltbankDir)
# Do matched-filtering
insps = ahope.setup_matchedfltr_workflow(cp, scienceSegs, datafinds, dag, 
                                         splitBanks, inspiralDir)



logging.info("Finished.")
