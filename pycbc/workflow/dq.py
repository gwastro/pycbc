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

import logging
import numpy
from pycbc.workflow.core import (FileList, Executable, Node, make_analysis_dir)


class PyCBCBinTemplatesDQExecutable(Executable):
    current_retention_level = Executable.MERGED_TRIGGERS

    def create_node(self, workflow, ifo, template_bank_file):
        node = Node(self)
        node.add_opt('--ifo', ifo)
        node.add_input_opt('--bank-file', template_bank_file)
        node.new_output_file_opt(
            workflow.analysis_time, '.hdf', '--output-file')
        return node


class PyCBCBinTriggerRatesDQExecutable(Executable):
    current_retention_level = Executable.MERGED_TRIGGERS

    def create_node(self, workflow, flag_file, flag_name,
                    analysis_segment_file, analysis_segment_name,
                    trig_file, template_bins_file):
        node = Node(self)
        node.add_input_opt('--template-bins-file', template_bins_file)
        node.add_input_opt('--trig-file', trig_file)
        node.add_input_opt('--flag-file', flag_file)
        node.add_opt('--flag-name', flag_name)
        node.add_input_opt('--analysis-segment-file', analysis_segment_file)
        node.add_opt('--analysis-segment-name', analysis_segment_name)
        node.new_output_file_opt(workflow.analysis_time, '.hdf',
                                 '--output-file')
        return node


def setup_dq_reranking(workflow, insps, bank,
                       analyzable_seg_file,
                       analyzable_name,
                       dq_seg_file,
                       output_dir=None, tags=None):
    logging.info("Setting up dq reranking")
    make_analysis_dir(output_dir)
    output_files = FileList()
    output_labels = []
    if tags is None:
        tags = []

    dq_labels = workflow.cp.get_subsections('workflow-data_quality')
    dq_labels = numpy.array(dq_labels)

    dq_types = numpy.array([workflow.cp.get_opt_tags(
            'workflow-data_quality', 'dq-type', [dq_label])
            for dq_label in dq_labels])
    dq_ifos = numpy.array([workflow.cp.get_opt_tags(
            'workflow-data_quality', 'dq-ifo', [dq_label])
            for dq_label in dq_labels])
    dq_names = numpy.array([workflow.cp.get_opt_tags(
            'workflow-data_quality', 'dq-name', [dq_label])
            for dq_label in dq_labels])

    ifos = set(dq_ifos)

    for ifo in ifos:
        # FIXME : make this able to take multiple dq files per ifo

        # get the dq label, type, and name for this ifo
        mask = numpy.array([ifo == dq_ifo for dq_ifo in dq_ifos])
        assert numpy.sum(mask) > 0, f"Received no dq files for {ifo}"
        assert numpy.sum(mask) < 2, f"Received more than one dq file for {ifo}"
        dq_label = str(dq_labels[mask][0])
        dq_type = str(dq_types[mask][0])
        dq_name = str(dq_names[mask][0])

        dq_tags = tags + [dq_label]

        # get triggers for this ifo
        ifo_insp = [insp for insp in insps if (insp.ifo == ifo)]
        assert len(ifo_insp) == 1, \
            f"Received more than one inspiral file for {ifo}"
        ifo_insp = ifo_insp[0]

        # calculate template bins for this ifo
        bin_templates_exe = PyCBCBinTemplatesDQExecutable(
            workflow.cp,
            'bin_templates',
            ifos=ifo,
            out_dir=output_dir,
            tags=tags)
        bin_templates_node = bin_templates_exe.create_node(workflow, ifo, bank)
        workflow += bin_templates_node
        template_bins_file = bin_templates_node.output_file

        if dq_type == 'flag':
            flag_file = dq_seg_file
            flag_name = dq_name
        else:
            raise ValueError(f"{dq_type} dq support not yet implemented")

        # calculate trigger rates during dq flags
        bin_triggers_exe = PyCBCBinTriggerRatesDQExecutable(
            workflow.cp,
            'bin_trigger_rates_dq',
            ifos=ifo,
            out_dir=output_dir,
            tags=dq_tags)
        bin_triggers_node = bin_triggers_exe.create_node(
            workflow,
            flag_file,
            flag_name,
            analyzable_seg_file,
            analyzable_name,
            ifo_insp,
            template_bins_file)
        workflow += bin_triggers_node
        output_files += bin_triggers_node.output_files
        output_labels += [dq_label]

    logging.info("Finished setting up dq reranking")
    return output_files, output_labels
