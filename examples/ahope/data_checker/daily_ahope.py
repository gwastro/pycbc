import os
import logging
from glue import pipeline
from glue import segments
import pycbc.ahope as ahope

logging.basicConfig(format='%(asctime)s:%(levelname)s : %(message)s', 
                    level=logging.DEBUG, datefmt='%I:%M:%S')

start_time = 1060881616
end_time = 1061856016

workflow = ahope.Workflow('./daily_er4.ini')
ifos = ['H1','L1','V1']
currDir = os.getcwd()
segDir = os.path.join(currDir,"segments")
if not os.path.exists(segDir+'/logs'):
    os.makedirs(segDir+'/logs')
dfDir = os.path.join(currDir,"datafind")
if not os.path.exists(dfDir+'/logs'):
    os.makedirs(dfDir+'/logs')

scienceSegs, segsDict = ahope.setup_segment_generation(workflow, ifos, 
                                     start_time, end_time, segDir)
datafinds, scienceSegs = ahope.setup_datafind_workflow(workflow, scienceSegs, 
                        dfDir, updateSegmentTimes=True, checkSegmentGaps=False)
