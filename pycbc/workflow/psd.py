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

from pycbc.workflow.core import FileList, make_analysis_dir, Executable
from pycbc.workflow.core import SegFile
from ligo.segments import segmentlist

class CalcPSDExecutable(Executable):
    current_retention_level = Executable.ALL_TRIGGERS

class MergePSDFiles(Executable):
    current_retention_level = Executable.MERGED_TRIGGERS

def chunks(l, n):
    """ Yield n successive chunks from l.
    """
    newn = int(len(l) / n)
    for i in range(0, n-1):
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
                        segment_name, out_dir, tags=None):
    make_analysis_dir(out_dir)
    tags = [] if not tags else tags
    if workflow.cp.has_option_tags('workflow-psd', 'parallelization-factor', tags=tags):
        num_parts = int(workflow.cp.get_opt_tags('workflow-psd',
                                                 'parallelization-factor',
                                                 tags=tags))
    else:
        num_parts = 1

    # get rid of duplicate segments which happen when splitting the bank
    segments = segmentlist(frozenset(segments))

    segment_lists = list(chunks(segments, num_parts))

    psd_files = FileList([])
    for i, segs in enumerate(segment_lists):
        seg_file = SegFile.from_segment_list('%s_%s' %(segment_name, i),
                         segmentlist(segs), segment_name, ifo,
                         valid_segment=workflow.analysis_time,
                         extension='xml', directory=out_dir)

        psd_files += [make_psd_file(workflow, frame_files, seg_file,
                                    segment_name, out_dir,
                                    tags=tags + ['PART%s' % i])]

    if num_parts > 1:
        return merge_psds(workflow, psd_files, ifo, out_dir, tags=tags)
    else:
        return psd_files[0]

def make_psd_file(workflow, frame_files, segment_file, segment_name, out_dir,
                  tags=None):
    make_analysis_dir(out_dir)
    tags = [] if not tags else tags
    exe = CalcPSDExecutable(workflow.cp, 'calculate_psd',
                             ifos=segment_file.ifo, out_dir=out_dir,
                             tags=tags)
    node = exe.create_node()
    node.add_input_opt('--analysis-segment-file', segment_file)
    node.add_opt('--segment-name', segment_name)

    if frame_files and not exe.has_opt('frame-type'):
        node.add_input_list_opt('--frame-files', frame_files)

    node.new_output_file_opt(workflow.analysis_time, '.hdf', '--output-file')
    workflow += node
    return node.output_files[0]

class AvgPSDExecutable(Executable):
    current_retention_level = Executable.FINAL_RESULT

def make_average_psd(workflow, psd_files, out_dir, tags=None,
                     output_fmt='.txt'):
    make_analysis_dir(out_dir)
    tags = [] if tags is None else tags
    node = AvgPSDExecutable(workflow.cp, 'average_psd', ifos=workflow.ifos,
                            out_dir=out_dir, tags=tags).create_node()
    node.add_input_list_opt('--input-files', psd_files)

    if len(workflow.ifos) > 1:
        node.new_output_file_opt(workflow.analysis_time, output_fmt,
                                 '--detector-avg-file')

    node.new_multiifo_output_list_opt('--time-avg-file', workflow.ifos,
                                 workflow.analysis_time, output_fmt, tags=tags)

    workflow += node
    return node.output_files

# keep namespace clean
__all__ = ['make_psd_file', 'make_average_psd', 'setup_psd_calculate', 'merge_psds']
