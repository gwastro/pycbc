import os
from glue import pipeline
from glue import segments
import pycbc.ahope as ahope

cp = ahope.parse_ahope_ini_file('./reduced.ini')
seg = segments.segment([853000000,853010000])
scienceSeg = segments.segmentlist([seg])
scienceSegs = {}
scienceSegs['H1'] = scienceSeg
scienceSegs['L1'] = scienceSeg
scienceSegs['H2'] = scienceSeg
basename = 'ahope_test'
logfile = 'asdtwe.log'
fh = open( logfile, "w" )
fh.close()
dag = pipeline.CondorDAG(logfile,dax=False)
dag.set_dax_file(basename)
dag.set_dag_file(basename)
currDir = os.getcwd()
dfDir = os.path.join(currDir,"datafind")
if not os.path.exists(dfDir+'/logs'):
    os.makedirs(dfDir+'/logs')
tmpltbankDir = os.path.join(currDir,"tmpltbank")
if not os.path.exists(tmpltbankDir+'/logs'):
    os.makedirs(tmpltbankDir+'/logs')
inspiralDir = os.path.join(currDir,"inspiral")
if not os.path.exists(inspiralDir+'/logs'):
    os.makedirs(inspiralDir+'/logs')

datafinds, scienceSegs = ahope.setup_datafind_workflow(cp, scienceSegs, dag,\
                       dfDir)
banks = ahope.setup_tmpltbank_workflow(cp, scienceSegs, datafinds, dag, \
                                       tmpltbankDir)
splitBanks = ahope.setup_splittable_workflow(cp, dag, banks, tmpltbankDir)
insps = ahope.setup_matchedfltr_workflow(cp, scienceSegs, datafinds, dag, \
                                         splitBanks, inspiralDir)

for insp in insps:
    for file in insp.get_output():
        print file.path

dag.write_sub_files()
dag.write_dag()
#dag.write_abstract_dag()
dag.write_script()
