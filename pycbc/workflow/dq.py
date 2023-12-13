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

import numpy
from pycbc.workflow.core import (FileList, Executable, Node,
                                 SegFile, make_analysis_dir)
from pycbc.workflow.datafind import setup_datafind_workflow


class PyCBCBinTemplatesDQExecutable(Executable):
    current_retention_level = Executable.MERGED_TRIGGERS

    def create_node(self, workflow, ifo, template_bank_file, trigger_file):
        node = Node(self)
        node.add_opt('--ifo', ifo)
        node.add_input_opt('--template-bank-file', template_bank_file)
        node.add_input_opt('--trigger-file', trigger_file)
        node.new_output_file_opt(workflow.analysis_time, '.hdf', '--output-file')
        return node


class PyCBCDQFlagFromTimeseriesExecutable(Executable):
    current_retention_level = Executable.MERGED_TRIGGERS

    def create_node(self, workflow, dq_label, analyzable_segs, dq_frames):
        node = Node(self)
        node.add_input_opt('--dq-label', dq_label)
        node.add_input_list_opt('--analyzable-segs', analyzable_segs)
        node.add_input_list_opt('--dq-frame-files', dq_frames)
        node.new_output_file_opt(workflow.analysis_time, '.hdf',
                                 '--output-file')
        return node


class PyCBCBinTriggerRatesDQExecutable(Executable):
    current_retention_level = Executable.MERGED_TRIGGERS

    def create_node(self, workflow, ifo, dq_label,
                    flag_file, trig_file, template_bins_file):
        node = Node(self)
        node.add_opt('--ifo', ifo)
        node.add_input_opt('--dq-label', dq_label)
        node.add_input_opt('--template-bins-file', template_bins_file)
        node.add_input_opt('--trig-file', trig_file)
        node.add_input_opt('--flag-file', flag_file)
        node.new_output_file_opt(workflow.analysis_time, '.hdf',
                                 '--output-file')
        return node


def setup_dq_reranking(workflow, insps, bank,
                       analyzable_segs, analyzable_file,
                       dq_seg_file,
                       output_dir=None, tags=None):
    make_analysis_dir(output_dir)
    output_files = FileList()
    output_labels = []

    dq_labels = workflow.cp.get_subsections('workflow-data_quality')
    if tags:
        dq_tags = tags + dq_labels
    else:
        dq_tags = dq_labels

    dq_types = [workflow.cp.get_opt_tags('workflow-data_quality', 'dq-type', [dq_label])
                for dq_label in dq_labels]
    dq_ifos = [workflow.cp.get_opt_tags('workflow-data_quality', 'ifo', [dq_label])
               for dq_label in dq_labels]
    ifos = numpy.unique(dq_ifos)

    for ifo in ifos:
        # get the dq labels and types for this ifo
        ifo_dq_info = [(dq_label, dq_type) for dq_label, dq_type, dq_ifo
                       in zip(dq_labels, dq_types, dq_ifos) if dq_ifo == ifo]
        assert len(ifo_dq_info) == 1, f"Received more than one dq file for {ifo}"
        dq_label, dq_type = ifo_dq_info[0]

        # get triggers for this ifo
        ifo_insp = [insp for insp in insps if (insp.ifo == ifo)]
        assert len(ifo_insp) == 1, f"Received more than one inspiral file for {ifo}"
        ifo_insp = ifo_insp[0]

        # calculate template bins for this ifo
        bin_templates_exe = PyCBCBinTemplatesDQExecutable(
            workflow.cp,
            'bin_templates_dq',
            ifos=ifo,
            out_dir=output_dir,
            tags=dq_tags)
        bin_templates_node = bin_templates_exe.create_node(workflow, ifo, bank, ifo_insp)
        workflow += bin_templates_node
        template_bins_file = bin_templates_node.output_file

        flag_file = None

        # if dq is a timeseries, need to convert it to a flag
        if dq_type == 'timeseries':
            if dq_label not in workflow.cp.get_subsections('workflow-datafind'):
                msg = """No workflow-datafind section with dq tag.
                Tags must be used in workflow-datafind sections "
                if more than one source of data is used.
                Strain data source must be tagged
                workflow-datafind-hoft.
                Consult the documentation for more info."""
                raise ValueError(msg)

            # find timeseries frames
            ifo_seg_dict = {ifo: analyzable_segs[ifo]}
            dq_segs_for_file = {ifo+':'+dq_label: analyzable_segs[ifo]}
            ts_dq_seg_file = SegFile.from_segment_list_dict(
                dq_label,
                dq_segs_for_file,
                extension='.xml',
                valid_segment=workflow.analysis_time,
                directory=output_dir)
            datafind_files, dq_file, dq_segs, dq_name = setup_datafind_workflow(
                workflow,
                ifo_seg_dict,
                "datafind_dq",
                seg_file=ts_dq_seg_file,
                tags=dq_tags)

            # calculate dq flag from timeseries
            calculate_flag_exe = PyCBCDQFlagFromTimeseriesExecutable(
                workflow.cp,
                'calculate_dqflag',
                ifos=ifo,
                out_dir=output_dir,
                tags=dq_tags)
            calculate_flag_node = calculate_flag_exe.create_node(
                workflow,
                analyzable_segs,
                datafind_files)
            workflow += calculate_flag_node
            flag_file = calculate_flag_node.output_file
        else:
            assert dq_type == 'flag', f"Unknown dq type {dq_type}"
            flag_file = dq_seg_file

        # calculate trigger rates during dq flags
        bin_triggers_exe = PyCBCBinTriggerRatesDQExecutable(
            workflow.cp,
            'bin_trigger_rates_dq',
            ifos=ifo,
            out_dir=output_dir,
            tags=dq_tags)
        bin_triggers_node = bin_triggers_exe.create_node(
            workflow,
            ifo,
            flag_file,
            ifo_insp,
            template_bins_file)
        workflow += bin_triggers_node
        output_files += bin_triggers_node.output_files
        output_labels += [dq_label]

    return output_files, output_labels
