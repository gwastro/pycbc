# Copyright (C) 2020 Max Trevor and Derek Davis
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

#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#

import os
import logging
from ligo import segments
from pycbc.workflow.core import (FileList, Executable, Node,
                                 File, SegFile, make_analysis_dir)
from pycbc.workflow.datafind import setup_datafind_workflow

class PyCBCCalculateDQExecutable(Executable):
    current_retention_level = Executable.ALL_TRIGGERS
    def create_node(self, segment, frames):
        start = int(segment[0])
        end = int(segment[1])
        node = Node(self)
        node.add_input_list_opt('--frame-files', frames)
        node.add_opt('--gps-start-time', start)
        node.add_opt('--gps-end-time', end)
        node.new_output_file_opt(segment, '.hdf', '--output-file')
        return node

class PyCBCRerankDQExecutable(Executable):
    current_retention_level = Executable.MERGED_TRIGGERS
    def create_node(self, workflow, ifo, dq_type, dq_files, binned_rate_file):
        node = Node(self)
        node.add_opt('--dq-type', dq_type)
        node.add_opt('--ifo', ifo)
        node.add_input_list_opt('--input-file', dq_files)
        node.add_input_opt('--rate-file', binned_rate_file)
        node.new_output_file_opt(workflow.analysis_time, '.hdf',
                                 '--output-file')
        return node

class PyCBCBinTriggerRatesDQExecutable(Executable):
    current_retention_level = Executable.MERGED_TRIGGERS
    def create_node(self, workflow, ifo, dq_files, trig_file, bank_file):
        node = Node(self)
        node.add_opt('--ifo', ifo)
        node.add_input_opt('--bank-file', bank_file)
        node.add_input_opt('--trig-file', trig_file)
        node.add_input_list_opt('--dq-file', dq_files)
        node.new_output_file_opt(workflow.analysis_time,'.hdf',
                                 '--output-file')
        return node

class PyCBCCalculateDQFlagExecutable(Executable):
    # current_retention_level = Executable.ALL_TRIGGERS
    current_retention_level = Executable.MERGED_TRIGGERS

    def create_node(self, workflow, segment, dq_file, flag):
        node = Node(self)
        # Executable objects are initialized with ifo information
        start = int(segment[0])
        end = int(segment[1])
        node.add_opt('--ifo', self.ifo_string)
        node.add_opt('--flag', flag)
        node.add_opt('--gps-start-time', start)
        node.add_opt('--gps-end-time', end)
        node.add_input_opt('--dq-segments', dq_file)
        node.new_output_file_opt(workflow.analysis_time, '.hdf',
                                 '--output-file')
        return node

def setup_dq_reranking(workflow, dq_label, insps, bank,
                        segs, analyzable_file, dq_file,
                        output_dir=None, tags=None):
    make_analysis_dir(output_dir)
    output = FileList()
    if tags:
        dq_tags = tags + [dq_label]
    else:
        dq_tags = [dq_label]
    dq_type =  workflow.cp.get_opt_tags("workflow-data_quality",
                                         'dq-type', [dq_label])
    if dq_type == 'timeseries':
        if dq_label not in workflow.cp.get_subsections('workflow-datafind'):
            msg = """No workflow-datafind section with dq tag.
              Tags must be used in workflow-datafind sections "
              if more than one source of data is used.
              Strain data source must be tagged
              workflow-datafind-hoft.
              Consult the documentation for more info."""
            raise ValueError(msg)
        dq_ifos = workflow.cp.get_opt_tags("workflow-data_quality",
                                         'ifos', [dq_label])
        dq_ifos = dq_ifos.split(',')
        dq_segs = {}
        dq_segs_for_file = {}
        for ifo in dq_ifos:
            dq_segs[ifo] = segs[ifo]
            dq_segs_for_file[ifo+':'+dq_label] = segs[ifo]
        dq_segs_file = SegFile.from_segment_list_dict(dq_label,
                                          dq_segs_for_file,
                                          extension='.xml',
                                          valid_segment=workflow.analysis_time,
                                          directory=output_dir)
        datafind_files, dq_file, dq_segs, dq_name = \
                                           setup_datafind_workflow(workflow,
                                           dq_segs, "datafind_dq",
                                           seg_file=dq_segs_file,
                                           tags=dq_tags)
        for ifo in dq_ifos:
            ifo_insp = [insp for insp in insps if (insp.ifo == ifo)]
            assert len(ifo_insp)==1
            ifo_insp = ifo_insp[0]

            dq_files = FileList()
            for seg in dq_segs[ifo]:
                seg_frames = datafind_files.find_all_output_in_range(ifo, seg)
                raw_exe  = PyCBCCalculateDQExecutable(workflow.cp,
                                                   'calculate_dq', ifos=ifo,
                                                   out_dir=output_dir,
                                                   tags=dq_tags)
                raw_node = raw_exe.create_node(seg, seg_frames)
                workflow += raw_node
                dq_files += raw_node.output_files

            intermediate_exe = PyCBCBinTriggerRatesDQExecutable(workflow.cp,
                                                   'bin_trigger_rates_dq',
                                                   ifos=ifo,
                                                   out_dir=output_dir,
                                                   tags=dq_tags)
            intermediate_node = intermediate_exe.create_node(workflow, ifo,
                                                             dq_files,
                                                             ifo_insp, bank)
            workflow += intermediate_node
            binned_rate_file = intermediate_node.output_file

            new_exe = PyCBCRerankDQExecutable(workflow.cp,
                                                   'rerank_dq', ifos=ifo,
                                                   out_dir=output_dir,
                                                   tags=dq_tags)
            new_node = new_exe.create_node(workflow, ifo, dq_label,
                                           dq_files, binned_rate_file)
            workflow += new_node
            output += new_node.output_files
    elif dq_type == 'flag':
        flag_str = workflow.cp.get_opt_tags("workflow-data_quality",
                                            'flag-name', dq_tags)
        ifo = flag_str[:2]
        ifo_insp = [insp for insp in insps if (insp.ifo == ifo)]
        assert len(ifo_insp)==1
        ifo_insp = ifo_insp[0]
        flag_name = flag_str
        logging.info("Creating job for flag %s", flag_name)
        for seg in segs[ifo]:
            raw_exe = PyCBCCalculateDQFlagExecutable(workflow.cp,
                                                 'calculate_dqflag', ifos=ifo,
                                                 out_dir=output_dir,
                                                 tags=dq_tags)
            raw_node = raw_exe.create_node(workflow, seg, dq_file,
                                           flag_name)
            workflow += raw_node
            dq_files = raw_node.output_files
        intermediate_exe = PyCBCBinTriggerRatesDQExecutable(workflow.cp,
                                               'bin_trigger_rates_dq',
                                               ifos=ifo,
                                               out_dir=output_dir,
                                               tags=dq_tags)
        intermediate_node = intermediate_exe.create_node(workflow, ifo,
                                                         dq_files,
                                                         ifo_insp, bank)
        workflow += intermediate_node
        binned_rate_file = intermediate_node.output_file

        new_exe = PyCBCRerankDQExecutable(workflow.cp,
                                               'rerank_dq', ifos=ifo,
                                               out_dir=output_dir,
                                               tags=dq_tags)
        new_node = new_exe.create_node(workflow, ifo, dq_label,
                                       dq_files, binned_rate_file)
        workflow += new_node
        output += new_node.output_files
    else:
        msg = """Incorrect DQ type specified.
              Only valid DQ types are 'flag'
              and 'timeseries'.
              Consult the documentation for more info."""
        raise ValueError(msg)

    return output
