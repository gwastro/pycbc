""" A minimal pycbc workflow example """

import argparse
import pycbc
import pycbc.workflow as wf
import os

pycbc.init_logging(True)
parser = argparse.ArgumentParser(description=__doc__[1:])
parser.add_argument("--multilevel", action='store_true', default=False)
wf.add_workflow_command_line_group(parser)
wf.add_workflow_settings_cli(parser)
args = parser.parse_args()

input_file = wf.resolve_url_to_file("test_input.txt")
input_file.add_pfn(os.path.abspath('./test_input.txt'), 'local')

cont = wf.Workflow(args, 'cont')
sub1 = wf.Workflow(args, 'sub1')
sub1_1 = wf.Workflow(args, 'sub1_1')
sub2 = wf.Workflow(args, 'sub2')

exe1 = wf.Executable(cont.cp, 'exe1')

SUBSUB = args.multilevel

# Subworkflow 1: generate file that will be needed later

# PATH1: generate that input in a sub-sub workflow
if SUBSUB:
    # Subworkflow 1: sub-subworkflow
    node = exe1.create_node()
    wf1_out_file = wf.File.from_path(os.path.abspath("test_output.txt.2"))
    node.add_input_arg(input_file)
    node.add_output_arg(wf1_out_file)
    sub1_1 += node
    sub1 += sub1_1
    
#PATH2: generate that input within a sub-workflow
else:
    node = exe1.create_node()
    wf1_out_file = wf.File.from_path(os.path.abspath("test_output.txt.2"))
    node.add_input_arg(input_file)
    node.add_output_arg(wf1_out_file)
    sub1 += node

# Subworkflow 2
# TEST GOAL: Job within subworkflow gets the input file produce in 
# some external workflow
node2 = exe1.create_node()
node2.add_input_arg(wf1_out_file)
node2.add_output_arg(wf.File.from_path(os.path.abspath("test_output.txt.3")))
sub2 += node2

# Regular job in top-level workflow
# Test GOAL: job *directly* in workflow gets file produced by subworkflow in 
# the same workflow
node2 = exe1.create_node()
node2.add_input_arg(wf1_out_file)
node2.add_output_arg(wf.File.from_path(os.path.abspath("test_output.txt.4")))
cont += node2

cont += sub1
cont += sub2
cont.save()
