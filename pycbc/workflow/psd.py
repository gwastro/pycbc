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


class CalcPSDExecutable(Executable):
    current_retention_level = Executable.CRITICAL

def make_psd_file(workflow, frame_files, segment_file, segment_name, out_dir,
                  tags=None):
    make_analysis_dir(out_dir)
    tags = [] if not tags else tags
    node = CalcPSDExecutable(workflow.cp, 'calculate_psd',
                             ifos=segment_file.ifo, out_dir=out_dir,
                             tags=tags).create_node()
    node.add_input_opt('--analysis-segment-file', segment_file)
    node.add_opt('--segment-name', segment_name)
    node.add_input_list_opt('--frame-files', frame_files)
    node.new_output_file_opt(workflow.analysis_time, '.hdf', '--output-file')
    workflow += node
    return node.output_files[0]

class AvgPSDExecutable(Executable):
    current_retention_level = Executable.FINAL_RESULT

def make_average_psd(workflow, psd_files, out_dir, tags=None):
    make_analysis_dir(out_dir)
    tags = [] if tags is None else tags
    node = AvgPSDExecutable(workflow.cp, 'average_psd', ifos=workflow.ifos,
                            out_dir=out_dir, tags=tags).create_node()
    node.add_input_list_opt('--input-files', psd_files)
    node.new_output_file_opt(workflow.analysis_time, '.txt',
                             '--detector-avg-file')

    # FIXME should Node have a public method for handling
    # multidetector output options of type --option H1:foo L1:bar?
    node.add_opt('--time-avg-file')
    for ifo in workflow.ifos:
        time_avg_file = File(ifo, node.executable.name, workflow.analysis_time,
                             extension='.txt', directory=out_dir, tags=tags)
        multi_ifo_string = ifo + ':' + time_avg_file.name
        node.add_opt(multi_ifo_string)
        node._add_output(time_avg_file)

    workflow += node
    return node.output_files[0]

# keep namespace clean
__all__ = ['make_psd_file', 'make_average_psd']
