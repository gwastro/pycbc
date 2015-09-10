# Copyright (C) 2015 Christopher M. Biwer
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

import logging
from pycbc.workflow.core import Executable, FileList, Node

def setup_minifollowups(workflow, out_dir, frame_files,
                             coinc_file, tmpltbank_file, data_type, tags=None):
    ''' This performs a series of followup jobs on the num_events-th loudest
    events.
    '''

    logging.info('Entering minifollowups module')

    if tags == None: tags = []

    # create a FileList that will contain all output files
    output_filelist = FileList([])

    # check if minifollowups section exists
    # if not then do not do add minifollowup jobs to the workflow
    if not workflow.cp.has_section('workflow-minifollowups'):
      logging.info('There is no [workflow-minifollowups] section in configuration file')
      logging.info('Leaving minifollowups')
      return output_filelist

    # loop over number of loudest events to be followed up
    num_events = int(workflow.cp.get_opt_tags('workflow-minifollowups', 'num-events', ''))
    for num_event in range(num_events):

        # increment by 1 for human readability
        num_event += 1

        # get output directory for this event
        tag_str = '_'.join(tags)
        output_dir = out_dir['result/loudest_event_%d_of_%d_%s'%(num_event, num_events, tag_str)]

        # make a pycbc_mf_table node for this event
        table_exe = MinifollowupsTableExecutable(workflow.cp, 'mf_table',
                        workflow.ifo_string, output_dir, tags=tags)
        table_node = table_exe.create_node(workflow.analysis_time, coinc_file,
                        tmpltbank_file, data_type, num_event)
        workflow.add_node(table_node)
        output_filelist.extend(table_node.output_files)

    logging.info('Leaving minifollowups module')

    return output_filelist

class MinifollowupsTableExecutable(Executable):
    ''' The class responsible for creating jobs for pycbc_mf_table.'''

    current_retention_level = Executable.FINAL_RESULT
    def __init__(self, cp, exe_name,
                 ifo=None, out_dir=None, tags=[], universe=None):
        super(MinifollowupsTableExecutable, self).__init__(cp, exe_name, universe, ifo, out_dir, tags=tags)

    def create_node(self, segment, coinc_file, tmpltbank_file, data_type, loudest_event_number):
        ''' Creates a node for the loudest_event_number-th loudest event in the
        coinc_file. '''

        # make a node
        node = Node(self)

        # add input files
        node.add_input_opt('--coinc-file', coinc_file)
        node.add_input_opt('--tmpltbank-file', tmpltbank_file)

        # add options
        node.add_opt('--data-type', data_type)
        node.add_opt('--loudest-event-number', loudest_event_number)

        # add output file
        node.new_output_file_opt(segment, '.html', '--output-file',
                                 store_file=self.retain_files)

        return node
