import os
import copy
import logging
from glue import pipeline
from glue import segments
import pycbc.ahope as ahope

logging.basicConfig(format='%(asctime)s:%(levelname)s : %(message)s', \
                    level=logging.INFO,datefmt='%I:%M:%S')

# This is the "official" start according to Stuart
start_time = 1073822416
# End times set as current time - one day (86400s)
end_time= 1076352849

workflow = ahope.Workflow('./daily_er5.ini')
ifos = ['H1','L1','V1']
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


scienceSegs, segsList = ahope.setup_segment_generation(workflow, ifos,
                                start_time, end_time, segDir)

segment_report(scienceSegs)

print
print

# Start with SYR comparison
scienceSegsS = copy.deepcopy(scienceSegs)
print "RUNNING DATAFIND FOR SYR"
datafinds, scienceSegsS = ahope.setup_datafind_workflow(workflow, scienceSegsS,
                     dfDirSYR, segsList, tag="SYR")

segment_report(scienceSegsS)

print 
print
print "RUNNING DATAFIND FOR CIT"
scienceSegsC = copy.deepcopy(scienceSegs)
datafinds, scienceSegsC = ahope.setup_datafind_workflow(workflow, scienceSegsC,
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
scienceSegsS = copy.deepcopy(scienceSegs)
datafinds, scienceSegsS = ahope.setup_datafind_workflow(workflow, scienceSegsS,
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
scienceSegsS = copy.deepcopy(scienceSegs)
datafinds, scienceSegsS = ahope.setup_datafind_workflow(workflow, scienceSegsS,
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
scienceSegsS = copy.deepcopy(scienceSegs)
datafinds, scienceSegsS = ahope.setup_datafind_workflow(workflow, scienceSegsS,
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

