import os
import logging
from glue import pipeline
from glue import segments
import pycbc.ahope as ahope

logging.basicConfig(format='%(asctime)s:%(levelname)s : %(message)s', 
                    level=logging.INFO,datefmt='%I:%M:%S')

workflow = ahope.Workflow(['./test_ahope.ini', 'inj.ini'])

# Make directories for output
currDir = os.getcwd()
segDir = os.path.join(currDir, "segments")
if not os.path.exists(segDir+'/logs'):
    os.makedirs(segDir+'/logs')
dfDir = os.path.join(currDir,"datafind")
if not os.path.exists(dfDir+'/logs'):
    os.makedirs(dfDir+'/logs')
if not os.path.exists('coincidence/logs'):
    os.makedirs('coincidence/logs')

# Set start and end times

# These are Chris' example S6 day of data
start_time = 961585543
end_time = 961588487

# This is S6D chunk 3, an example full-scale ihope analysis time-frame
##start_time = 968543943
#end_time = 971622087

# Set the ifos to analyse
ifos = ['H1','L1']

# Get segments and find where the data is
scienceSegs, segsDict = ahope.setup_segment_generation(workflow, ifos, 
                                            start_time, end_time, segDir, 
                                            maxVetoCat=5,
                                            minSegLength=2000)
datafind_files, scienceSegs = ahope.setup_datafind_workflow(workflow, 
                                            scienceSegs, dfDir, 
                                            checkSegmentGaps='update_times',
                                            checkFramesExist='warn')

# Template bank stuff
bank_files = ahope.setup_tmpltbank_workflow(workflow, scienceSegs, 
                                            datafind_files, 'tmpltbank')
splitbank_files = ahope.setup_splittable_workflow(workflow, bank_files, 
                                                            'tmpltbank')

# setup the injection files
inj_files, inj_tags = ahope.setup_injection_workflow(workflow, scienceSegs, 
                                           datafind_files, splitbank_files, 
                                           start_time, end_time, 'inspiral')
tags = ["zero"] + inj_tags
inj_files = [None] + inj_files
all_coincs = []
for inj_file, tag in zip(inj_files, tags):
    insps = ahope.setup_matchedfltr_workflow(workflow, scienceSegs, 
                                           datafind_files, splitbank_files, 
                                           'inspiral', injection_file=inj_file,
                                           tags = [tag])
    coincs = ahope.setup_coincidence_workflow(workflow, scienceSegs, segsDict, 
                                              insps, 'coincidence', tags=[tag],
                                              maxVetoCat = 5)
    all_coincs.append(coincs)
    
workflow.write_plans()
logging.info("Finished.")
