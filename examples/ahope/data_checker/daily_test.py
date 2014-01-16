import os
import copy
import logging
from glue import pipeline
from glue import segments
import pycbc.ahope as ahope

logging.basicConfig(format='%(asctime)s:%(levelname)s : %(message)s', \
                    level=logging.INFO,datefmt='%I:%M:%S')

start_time = 1073779216
end_time = 1073838334

workflow = ahope.Workflow('./daily_er5.ini')
ifos = ['H1','L1','V1']
currDir = os.getcwd()
segDir = os.path.join(currDir,"segments")
if not os.path.exists(segDir+'/logs'):
    os.makedirs(segDir+'/logs')
dfDirSYR = os.path.join(currDir,"datafindSYR")
if not os.path.exists(dfDirSYR+'/logs'):
    os.makedirs(dfDirSYR+'/logs')
dfDirCIT = os.path.join(currDir,"datafindCIT")
if not os.path.exists(dfDirCIT+'/logs'):
    os.makedirs(dfDirCIT+'/logs')
dfDirLHO = os.path.join(currDir,"datafindLHO")
if not os.path.exists(dfDirLHO+'/logs'):
    os.makedirs(dfDirLHO+'/logs')
dfDirLLO = os.path.join(currDir,"datafindLLO")
if not os.path.exists(dfDirLLO+'/logs'):
    os.makedirs(dfDirLLO+'/logs')

print "BEGIN BY GENERATING SCIENCE AND CAT_X VETOES"

scienceSegs, segsDict = ahope.setup_segment_generation(workflow, ifos,
                                start_time, end_time, segDir)
print
print

# Start with SYR comparison
scienceSegsS = copy.deepcopy(scienceSegs)
print "RUNNING DATAFIND FOR SYR"
datafinds, scienceSegsS = ahope.setup_datafind_workflow(workflow, scienceSegsS,
                       dfDirSYR,checkSegmentGaps='update_times',\
                       checkFramesExist='no_test')

print 
print
print "RUNNING DATAFIND FOR CIT"
os.environ["LIGO_DATAFIND_SERVER"] = """ldr.ligo.caltech.edu:443"""
scienceSegsC = copy.deepcopy(scienceSegs)
datafinds, scienceSegsC = ahope.setup_datafind_workflow(workflow, scienceSegsC,
                       dfDirCIT, checkSegmentGaps='update_times',\
                       checkFramesExist='no_test')

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
os.environ["LIGO_DATAFIND_SERVER"] = """ldr.ligo-wa.caltech.edu:443"""
scienceSegsS = copy.deepcopy(scienceSegs)
datafinds, scienceSegsS = ahope.setup_datafind_workflow(workflow, scienceSegsS,
                       dfDirLHO,checkSegmentGaps='update_times',\
                       checkFramesExist='no_test')


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
os.environ["LIGO_DATAFIND_SERVER"] = """ldr.ligo-la.caltech.edu:443"""
scienceSegsS = copy.deepcopy(scienceSegs)
datafinds, scienceSegsS = ahope.setup_datafind_workflow(workflow, scienceSegsS,
                       dfDirLLO,checkSegmentGaps='update_times',\
                       checkFramesExist='no_test')

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

