import os
import copy
import logging
from glue import pipeline
from glue import segments
import pycbc.ahope as ahope

logging.basicConfig(format='%(asctime)s:%(levelname)s : %(message)s', \
                    level=logging.INFO,datefmt='%I:%M:%S')

start_time = 1057881568
end_time = 1061856064
end_time = 1058181568

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
scienceSegsS = copy.deepcopy(scienceSegs)
datafinds, scienceSegsS = ahope.setup_datafind_workflow(workflow, scienceSegsS,
                       dfDir, updateSegmentTimes=True, checkSegmentGaps=False,
                       checkFramesExist=False)
os.environ["LIGO_DATAFIND_SERVER"] = """ldr.ligo.caltech.edu:443"""
scienceSegsC = copy.deepcopy(scienceSegs)
datafinds, scienceSegsC = ahope.setup_datafind_workflow(workflow, scienceSegsC,
                       dfDir, updateSegmentTimes=True, checkSegmentGaps=False,
                       checkFramesExist=False)

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

# This does not use the dag anywhere
#dag.write_sub_files()
#dag.write_dag()
#dag.write_abstract_dag()
#dag.write_script()
