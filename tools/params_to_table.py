#!/usr/bin/env python

import time
import os
import sys
 
from optparse import OptionParser

from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils as ligolw_utils
from glue.ligolw.utils import process as ligolw_process
from glue.segmentdb import segmentdb_utils
from glue import pidfile as pidfile
from glue import git_version
from scipy.interpolate import interp1d
import matplotlib
matplotlib.use('Agg')
import pylab
from pylab import arange,pi,sin,cos,sqrt
from numpy import loadtxt

parser = OptionParser(
    version = git_version.verbose_msg,
    usage   = "%prog [OPTIONS]",
    description = "Creates a template bank and writes it to XML." )

parser.add_option("-V", "--verbose", action="store_true", help="print extra debugging information", default=False )
parser.add_option("-p", "--parameter-names",help="list of parameters names")
parser.add_option("-i", "--input-file",  help="file with the list of paramters")
parser.add_option("--type")

options, argv_frame_files = parser.parse_args()

outdoc = ligolw.Document()
outdoc.appendChild(ligolw.LIGO_LW())

proc_id = ligolw_process.register_to_xmldoc(outdoc, 
                "params_to_table.py", options.__dict__, comment="", ifos=[""],
                version=git_version.id, cvs_repository=git_version.branch,
                cvs_entry_time=git_version.date).process_id
                
params = loadtxt(options.input_file)

col_names = []

param_names = str.split(options.parameter_names)

for name in param_names:
    col_names.append(name)

print param_names                
if options.type == "sngl":
     sngl_inspiral_table = lsctables.New(lsctables.SnglInspiralTable,columns= col_names)
elif options.type == "sim":
     sngl_inspiral_table = lsctables.New(lsctables.SimInspiralTable,columns= col_names)

outdoc.childNodes[0].appendChild(sngl_inspiral_table)

for values in params:
    if options.type == "sngl":
        tmplt = lsctables.SnglInspiral()
    elif options.type == "sim":
        tmplt = lsctables.SimInspiral()

    tmplt.process_id = proc_id
    index = 0 
    for value in values:
        if value is 'skip':
            continue
        setattr(tmplt,param_names[index],value)
        index += 1
    sngl_inspiral_table.append(tmplt)

# write the xml doc to disk
proctable = table.get_table(outdoc, lsctables.ProcessTable.tableName)
outname = 'table.xml'
ligolw_utils.write_filename(outdoc, outname)

