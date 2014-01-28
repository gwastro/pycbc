import os
import logging
import pycbc.ahope
from pycbc.ahope.jobsetup_utils import *
from pycbc.ahope.matchedfltr_utils import *
from glue import segments

def setup_injection_workflow(workflow, science_segs, datafind_outs, 
                             tmplt_banks, start_time, end_time, output_dir=None):
    '''
    Setup matched filter with injections section of ahope workflow.
    FIXME: ADD MORE DOCUMENTATION
    '''
    logging.info("Entering injection module.")
    inj_exe = LalappsInspinjExec('injection')
    all_sec = workflow.cp.sections()
    sections = [sec for sec in all_sec if sec.startswith('injection-')]
    inj_tags = []
    inj_files = []
    for sec in sections:  
        segment = segments.segment(start_time, end_time)
        inj_tag = sec.split('-')[1]
        inj_job = inj_exe.create_job(workflow.cp, tags=[inj_tag], out_dir=output_dir)
        node = inj_job.create_node(segment)
        workflow.add_node(node)
        injection_file = node.output_files[0]
        inj_files.append(injection_file)
        inj_tags.append(inj_tag)
    logging.info("Leaving injection module.")  
    return inj_files, inj_tags

