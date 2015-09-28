# Copyright (C) 2013  Ian Harry
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

"""This module is responsible for setting up PSD-related jobs in workflows.
"""

from pycbc.workflow.core import FileList, make_analysis_dir, Executable, File
from pycbc.events import segments_to_file
from glue.segments import segmentlist

class CalcPSDExecutable(Executable):
    current_retention_level = Executable.CRITICAL

class MergePSDFiles(Executable):
    current_retention_level = Executable.CRITICAL

def chunks(l, n):
    """ Yield n successive chunks from l.
    """
    newn = int(len(l) / n)
    for i in xrange(0, n-1):
        yield l[i*newn:i*newn+newn]
    yield l[n*newn-newn:]

def merge_psds(workflow, files, ifo, out_dir, tags=None):
    make_analysis_dir(out_dir)
    tags = [] if not tags else tags
    node = MergePSDFiles(workflow.cp, 'merge_psds',
                         ifos=ifo, out_dir=out_dir,
                         tags=tags).create_node()
    node.add_input_list_opt('--psd-files', files)
    node.new_output_file_opt(workflow.analysis_time, '.hdf', '--output-file')
    workflow += node
    return node.output_files[0]        

def setup_psd_calculate(workflow, frame_files, ifo, segments,
                        segment_name, out_dir,
                        gate_files=None, tags=None):
    make_analysis_dir(out_dir)
    tags = [] if not tags else tags
    if workflow.cp.has_option_tags('workflow-psd', 'parallelization-factor', tags=tags):
        num_parts = int(workflow.cp.get_opt_tags('workflow-psd', 
                                                 'parallelization-factor',
                                                 tags=tags))
    else:
        num_parts = 1
        
    segment_lists = list(chunks(segments, num_parts)) 
    
    psd_files = FileList([])
    for i, segs in enumerate(segment_lists):
        seg_file = segments_to_file(segmentlist(segs), 
                               out_dir + '/%s-INSPIRAL_DATA-%s.xml' % (ifo, i), 
                               'INSPIRAL_DATA', ifo=ifo)

        psd_files += [make_psd_file(workflow, frame_files, seg_file,
                                    segment_name, out_dir, 
                                    gate_files=gate_files, 
                                    tags=tags + ['PART%s' % i])]
    
    if num_parts > 1:
        return merge_psds(workflow, psd_files, ifo, out_dir, tags=tags)
    else:
        return psd_files[0]
    
def make_psd_file(workflow, frame_files, segment_file, segment_name, out_dir,
                  gate_files=None, tags=None):
    make_analysis_dir(out_dir)
    tags = [] if not tags else tags
    exe = CalcPSDExecutable(workflow.cp, 'calculate_psd',
                             ifos=segment_file.ifo, out_dir=out_dir,
                             tags=tags)
    node = exe.create_node()
    node.add_input_opt('--analysis-segment-file', segment_file)
    node.add_opt('--segment-name', segment_name)
    
    if gate_files is not None:
        ifo_gate = None
        for gate_file in gate_files:
            if gate_file.ifo == segment_file.ifo:
                ifo_gate = gate_file
        
        if ifo_gate is not None:
            node.add_input_opt('--gating-file', ifo_gate)
    
    if not exe.has_opt('frame-type'):
        node.add_input_list_opt('--frame-files', frame_files)
    node.new_output_file_opt(workflow.analysis_time, '.hdf', '--output-file')
    workflow += node
    return node.output_files[0]

class AvgPSDExecutable(Executable):
    current_retention_level = Executable.FINAL_RESULT

def make_average_psd(workflow, psd_files, out_dir, tags=None,
                     gate_files=None,
                     output_fmt='.txt'):
    make_analysis_dir(out_dir)
    tags = [] if tags is None else tags
    node = AvgPSDExecutable(workflow.cp, 'average_psd', ifos=workflow.ifos,
                            out_dir=out_dir, tags=tags).create_node()
    node.add_input_list_opt('--input-files', psd_files)
    node.new_output_file_opt(workflow.analysis_time, output_fmt,
                             '--detector-avg-file')

    # FIXME should Node have a public method for handling
    # multidetector output options of type --option H1:foo L1:bar?
    node.add_opt('--time-avg-file')
    for ifo in workflow.ifos:
        time_avg_file = File(ifo, node.executable.name, workflow.analysis_time,
                             extension=output_fmt, directory=out_dir,
                             tags=tags)
        multi_ifo_string = ifo + ':' + time_avg_file.name
        node.add_opt(multi_ifo_string)
        node._add_output(time_avg_file)
    
        if gate_files is not None:
            ifo_gate = None
            for gate_file in gate_files:
                if gate_file.ifo == ifo:
                    ifo_gate = gate_file
            
            if ifo_gate is not None:
                node.add_input_opt('--gating-file', ifo_gate)

    workflow += node
    return node.output_files

# keep namespace clean
__all__ = ['make_psd_file', 'make_average_psd', 'setup_psd_calculate', 'merge_psds']
