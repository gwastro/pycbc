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

#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
"""
This module is responsible for setting up the coincidence stage of pycbc
workflows. For details about this module and its capabilities see here:
https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/coincidence.html
"""

import os
import logging
from ligo import segments
from pycbc.workflow.core import FileList, make_analysis_dir, Executable, Node, File

class PyCBCBank2HDFExecutable(Executable):
    """Converts xml tmpltbank to hdf format"""

    current_retention_level = Executable.MERGED_TRIGGERS
    def create_node(self, bank_file):
        node = Node(self)
        node.add_input_opt('--bank-file', bank_file)
        node.new_output_file_opt(bank_file.segment, '.hdf', '--output-file')
        return node


class PyCBCTrig2HDFExecutable(Executable):
    """Converts xml triggers to hdf format, grouped by template hash"""

    current_retention_level = Executable.MERGED_TRIGGERS
    def create_node(self, trig_files, bank_file):
        node = Node(self)
        node.add_input_opt('--bank-file', bank_file)
        node.add_input_list_opt('--trigger-files', trig_files)
        node.new_output_file_opt(trig_files[0].segment, '.hdf',
                                 '--output-file', use_tmp_subdirs=True)
        return node


class PyCBCFitByTemplateExecutable(Executable):
    """Calculates values that describe the background distribution template by template"""

    current_retention_level = Executable.MERGED_TRIGGERS
    def create_node(self, trig_file, bank_file, veto_file, veto_name):
        node = Node(self)
        # Executable objects are initialized with ifo information
        node.add_opt('--ifo', self.ifo_string)
        node.add_input_opt('--trigger-file', trig_file)
        node.add_input_opt('--bank-file', bank_file)
        node.add_input_opt('--veto-file', veto_file)
        node.add_opt('--veto-segment-name', veto_name)
        node.new_output_file_opt(trig_file.segment, '.hdf', '--output')
        return node


class PyCBCFitOverParamExecutable(Executable):
    """Smooths the background distribution parameters over a continuous parameter"""

    current_retention_level = Executable.MERGED_TRIGGERS
    def create_node(self, raw_fit_file, bank_file):
        node = Node(self)
        node.add_input_opt('--template-fit-file', raw_fit_file)
        node.add_input_opt('--bank-file', bank_file)
        node.new_output_file_opt(raw_fit_file.segment, '.hdf', '--output')
        return node


class PyCBCFindCoincExecutable(Executable):
    """Find coinc triggers using a folded interval method"""

    current_retention_level = Executable.ALL_TRIGGERS
    def create_node(self, trig_files, bank_file, stat_files, veto_file,
                    veto_name, template_str, pivot_ifo, fixed_ifo, tags=None):
        if tags is None:
            tags = []
        segs = trig_files.get_times_covered_by_files()
        seg = segments.segment(segs[0][0], segs[-1][1])
        node = Node(self)
        node.add_input_opt('--template-bank', bank_file)
        node.add_input_list_opt('--trigger-files', trig_files)
        if len(stat_files) > 0:
            node.add_input_list_opt('--statistic-files', stat_files)
        if veto_file is not None:
            node.add_input_opt('--veto-files', veto_file)
            node.add_opt('--segment-name', veto_name)
        node.add_opt('--pivot-ifo', pivot_ifo)
        node.add_opt('--fixed-ifo', fixed_ifo)
        node.add_opt('--template-fraction-range', template_str)
        node.new_output_file_opt(seg, '.hdf', '--output-file', tags=tags)
        return node


class PyCBCStatMapExecutable(Executable):
    """Calculate FAP, IFAR, etc"""

    current_retention_level = Executable.MERGED_TRIGGERS
    def create_node(self, coinc_files, ifos, tags=None):
        if tags is None:
            tags = []
        segs = coinc_files.get_times_covered_by_files()
        seg = segments.segment(segs[0][0], segs[-1][1])

        node = Node(self)
        node.add_input_list_opt('--coinc-files', coinc_files)
        node.add_opt('--ifos', ifos)
        node.new_output_file_opt(seg, '.hdf', '--output-file', tags=tags)
        return node


class PyCBCStatMapInjExecutable(Executable):
    """Calculate FAP, IFAR, etc"""

    current_retention_level = Executable.MERGED_TRIGGERS
    def create_node(self, zerolag, full_data,
                    injfull, fullinj, ifos, tags=None):
        if tags is None:
            tags = []
        segs = zerolag.get_times_covered_by_files()
        seg = segments.segment(segs[0][0], segs[-1][1])

        node = Node(self)
        node.add_input_list_opt('--zero-lag-coincs', zerolag)

        if isinstance(full_data, list):
            node.add_input_list_opt('--full-data-background', full_data)
        else:
            node.add_input_opt('--full-data-background', full_data)

        node.add_input_list_opt('--mixed-coincs-inj-full', injfull)
        node.add_input_list_opt('--mixed-coincs-full-inj', fullinj)
        node.add_opt('--ifos', ifos)
        node.new_output_file_opt(seg, '.hdf', '--output-file', tags=tags)
        return node


class PyCBCHDFInjFindExecutable(Executable):
    """Find injections in the hdf files output"""

    current_retention_level = Executable.MERGED_TRIGGERS
    def create_node(self, inj_coinc_file, inj_xml_file, veto_file, veto_name, tags=None):
        if tags is None:
            tags = []
        node = Node(self)
        node.add_input_list_opt('--trigger-file', inj_coinc_file)
        node.add_input_list_opt('--injection-file', inj_xml_file)
        if veto_name is not None:
            node.add_input_opt('--veto-file', veto_file)
        node.add_opt('--segment-name', veto_name)
        node.new_output_file_opt(inj_xml_file[0].segment, '.hdf',
                                 '--output-file', tags=tags)
        return node


class PyCBCDistributeBackgroundBins(Executable):
    """Distribute coinc files among different background bins"""

    current_retention_level = Executable.ALL_TRIGGERS
    def create_node(self, coinc_files, bank_file, background_bins, tags=None):
        if tags is None:
            tags = []
        node = Node(self)
        node.add_input_list_opt('--coinc-files', coinc_files)
        node.add_input_opt('--bank-file', bank_file)
        node.add_opt('--background-bins', ' '.join(background_bins))

        names = [b.split(':')[0] for b in background_bins]

        output_files = [File(coinc_files[0].ifo_list,
                             self.name,
                             coinc_files[0].segment,
                             directory=self.out_dir,
                             tags = tags + ['mbin-%s' % i],
                             extension='.hdf') for i in range(len(background_bins))]
        node.add_output_list_opt('--output-files', output_files)
        node.names = names
        return node


class PyCBCCombineStatmap(Executable):
    """Combine coincs over different bins and apply trials factor"""

    current_retention_level = Executable.MERGED_TRIGGERS
    def create_node(self, statmap_files, tags=None):
        if tags is None:
            tags = []
        node = Node(self)
        node.add_input_list_opt('--statmap-files', statmap_files)
        node.new_output_file_opt(statmap_files[0].segment, '.hdf',
                                 '--output-file', tags=tags)
        return node


class PyCBCAddStatmap(PyCBCCombineStatmap):
    """Combine statmap files and add FARs over different coinc types"""

    current_retention_level = Executable.MERGED_TRIGGERS
    def create_node(self, statmap_files, background_files, tags=None):
        if tags is None:
            tags = []
        node = super(PyCBCAddStatmap, self).create_node(statmap_files,
                                                            tags=tags)
        # Enforce upper case
        ctags = [t.upper() for t in (tags + self.tags)]
        if 'INJECTIONS' in ctags:
            node.add_input_list_opt('--background-files', background_files)

        return node


class PyCBCExcludeZerolag(Executable):
    """ Remove times of zerolag coincidences of all types from exclusive
        background """
    current_retention_level = Executable.MERGED_TRIGGERS
    def create_node(self, statmap_file, other_statmap_files, tags=None):
        if tags is None:
            tags = []
        node = Node(self)
        node.add_input_opt('--statmap-file', statmap_file)
        node.add_input_list_opt('--other-statmap-files',
                                other_statmap_files)
        node.new_output_file_opt(statmap_file.segment, '.hdf',
                                 '--output-file', tags=None)

        return node


class MergeExecutable(Executable):
    current_retention_level = Executable.MERGED_TRIGGERS


class CensorForeground(Executable):
    current_retention_level = Executable.MERGED_TRIGGERS


def make_foreground_censored_veto(workflow, bg_file, veto_file, veto_name,
                                  censored_name, out_dir, tags=None):
    tags = [] if tags is None else tags
    node = CensorForeground(workflow.cp, 'foreground_censor', ifos=workflow.ifos,
                            out_dir=out_dir, tags=tags).create_node()
    node.add_input_opt('--foreground-triggers', bg_file)
    node.add_input_opt('--veto-file', veto_file)
    node.add_opt('--segment-name', veto_name)
    node.add_opt('--output-segment-name', censored_name)
    node.new_output_file_opt(workflow.analysis_time, '.xml', '--output-file')
    workflow += node
    return node.output_files[0]


def merge_single_detector_hdf_files(workflow, bank_file, trigger_files, out_dir, tags=None):
    if tags is None:
        tags = []
    make_analysis_dir(out_dir)
    out = FileList()
    for ifo in workflow.ifos:
        node = MergeExecutable(workflow.cp, 'hdf_trigger_merge',
                        ifos=ifo, out_dir=out_dir, tags=tags).create_node()
        node.add_input_opt('--bank-file', bank_file)
        node.add_input_list_opt('--trigger-files', trigger_files.find_output_with_ifo(ifo))
        node.new_output_file_opt(workflow.analysis_time, '.hdf', '--output-file')
        workflow += node
        out += node.output_files
    return out


def setup_trigger_fitting(workflow, insps, hdfbank, veto_file, veto_name,
                          output_dir=None, tags=None):
    if not workflow.cp.has_option('workflow-coincidence', 'do-trigger-fitting'):
        return FileList()
    else:
        smoothed_fit_files = FileList()
        for i in workflow.ifos:
            ifo_insp = [insp for insp in insps if (insp.ifo == i)]
            assert len(ifo_insp)==1
            ifo_insp = ifo_insp[0]
            raw_exe = PyCBCFitByTemplateExecutable(workflow.cp,
                                                   'fit_by_template', ifos=i,
                                                   out_dir=output_dir,
                                                   tags=tags)
            raw_node = raw_exe.create_node(ifo_insp, hdfbank,
                                           veto_file, veto_name)
            workflow += raw_node
            smooth_exe = PyCBCFitOverParamExecutable(workflow.cp,
                                                     'fit_over_param', ifos=i,
                                                     out_dir=output_dir,
                                                     tags=tags)
            smooth_node = smooth_exe.create_node(raw_node.output_file,
                                                 hdfbank)
            workflow += smooth_node
            smoothed_fit_files += smooth_node.output_files
        return smoothed_fit_files


def find_injections_in_hdf_coinc(workflow, inj_coinc_file, inj_xml_file,
                                 veto_file, veto_name, out_dir, tags=None):
    if tags is None:
        tags = []
    make_analysis_dir(out_dir)
    exe = PyCBCHDFInjFindExecutable(workflow.cp, 'hdfinjfind',
                                    ifos=workflow.ifos,
                                    out_dir=out_dir, tags=tags)
    node = exe.create_node(inj_coinc_file, inj_xml_file, veto_file, veto_name)
    workflow += node
    return node.output_files[0]


def convert_bank_to_hdf(workflow, xmlbank, out_dir, tags=None):
    """Return the template bank in hdf format"""
    if tags is None:
        tags = []
    #FIXME, make me not needed
    if len(xmlbank) > 1:
        raise ValueError('Can only convert a single template bank')

    logging.info('convert template bank to HDF')
    make_analysis_dir(out_dir)
    bank2hdf_exe = PyCBCBank2HDFExecutable(workflow.cp, 'bank2hdf',
                                            ifos=workflow.ifos,
                                            out_dir=out_dir, tags=tags)
    bank2hdf_node = bank2hdf_exe.create_node(xmlbank[0])
    workflow.add_node(bank2hdf_node)
    return bank2hdf_node.output_files


def convert_trig_to_hdf(workflow, hdfbank, xml_trigger_files, out_dir, tags=None):
    """Return the list of hdf5 trigger files outputs"""
    if tags is None:
        tags = []
    #FIXME, make me not needed
    logging.info('convert single inspiral trigger files to hdf5')
    make_analysis_dir(out_dir)

    trig_files = FileList()
    for ifo, insp_group in zip(*xml_trigger_files.categorize_by_attr('ifo')):
        trig2hdf_exe = PyCBCTrig2HDFExecutable(workflow.cp, 'trig2hdf',
                                       ifos=ifo, out_dir=out_dir, tags=tags)
        _, insp_bundles = insp_group.categorize_by_attr('segment')
        for insps in insp_bundles:
            trig2hdf_node =  trig2hdf_exe.create_node(insps, hdfbank[0])
            workflow.add_node(trig2hdf_node)
            trig_files += trig2hdf_node.output_files
    return trig_files


def setup_statmap(workflow, ifos, coinc_files, out_dir, tags=None):
    tags = [] if tags is None else tags

    statmap_exe = PyCBCStatMapExecutable(workflow.cp, 'statmap',
                                         ifos=ifos,
                                         tags=tags, out_dir=out_dir)

    ifolist = ' '.join(ifos)
    stat_node = statmap_exe.create_node(coinc_files, ifolist)
    workflow.add_node(stat_node)
    return stat_node.output_file


def setup_statmap_inj(workflow, ifos, coinc_files, background_file,
                      out_dir, tags=None):
    tags = [] if tags is None else tags

    statmap_exe = PyCBCStatMapInjExecutable(workflow.cp,
                                            'statmap_inj',
                                            ifos=ifos,
                                            tags=tags, out_dir=out_dir)

    ifolist = ' '.join(ifos)
    stat_node = statmap_exe.create_node(FileList(coinc_files['injinj']),
                                        background_file,
                                        FileList(coinc_files['injfull']),
                                        FileList(coinc_files['fullinj']),
                                        ifolist)
    workflow.add_node(stat_node)
    return stat_node.output_files[0]


def setup_interval_coinc_inj(workflow, hdfbank, full_data_trig_files,
                             inj_trig_files, stat_files,
                             background_file, veto_file, veto_name,
                             out_dir, pivot_ifo, fixed_ifo, tags=None):
    """
    This function sets up exact match coincidence for injections
    """
    if tags is None:
        tags = []
    make_analysis_dir(out_dir)
    logging.info('Setting up coincidence for injections')

    # Wall time knob and memory knob
    factor = int(workflow.cp.get_opt_tags('workflow-coincidence',
                                          'parallelization-factor', tags))

    ffiles = {}
    ifiles = {}
    for ifo, ffi in zip(*full_data_trig_files.categorize_by_attr('ifo')):
        ffiles[ifo] = ffi[0]
    for ifo, ifi in zip(*inj_trig_files.categorize_by_attr('ifo')):
        ifiles[ifo] = ifi[0]

    injinj_files = FileList()
    injfull_files = FileList()
    fullinj_files = FileList()
    # For the injfull and fullinj separation we take the pivot_ifo on one side,
    # and the rest that are attached to the fixed_ifo on the other side
    for ifo in ifiles:  # ifiles is keyed on ifo
        if ifo == pivot_ifo:
            injinj_files.append(ifiles[ifo])
            injfull_files.append(ifiles[ifo])
            fullinj_files.append(ffiles[ifo])
        else:
            injinj_files.append(ifiles[ifo])
            injfull_files.append(ffiles[ifo])
            fullinj_files.append(ifiles[ifo])

    combo = [(injinj_files, "injinj"),
             (injfull_files, "injfull"),
             (fullinj_files, "fullinj"),
            ]
    bg_files = {'injinj':[], 'injfull':[], 'fullinj':[]}

    for trig_files, ctag in combo:
        findcoinc_exe = PyCBCFindCoincExecutable(workflow.cp,
                                                 'coinc',
                                                 ifos=ifiles.keys(),
                                                 tags=tags + [ctag],
                                                 out_dir=out_dir)
        for i in range(factor):
            group_str = '%s/%s' % (i, factor)
            coinc_node = findcoinc_exe.create_node(trig_files, hdfbank,
                                                   stat_files,
                                                   veto_file, veto_name,
                                                   group_str,
                                                   pivot_ifo,
                                                   fixed_ifo,
                                                   tags=['JOB'+str(i)])

            bg_files[ctag] += coinc_node.output_files
            workflow.add_node(coinc_node)

    logging.info('...leaving coincidence for injections')

    return setup_statmap_inj(workflow, ifiles.keys(), bg_files,
                             background_file, out_dir,
                             tags=tags + [veto_name])


def setup_interval_coinc(workflow, hdfbank, trig_files, stat_files,
                         veto_file, veto_name, out_dir, pivot_ifo,
                         fixed_ifo, tags=None):
    """
    This function sets up exact match coincidence
    """
    if tags is None:
        tags = []
    make_analysis_dir(out_dir)
    logging.info('Setting up coincidence')

    ifos, _ = trig_files.categorize_by_attr('ifo')
    findcoinc_exe = PyCBCFindCoincExecutable(workflow.cp, 'coinc',
                                             ifos=ifos,
                                             tags=tags, out_dir=out_dir)

    # Wall time knob and memory knob
    factor = int(workflow.cp.get_opt_tags('workflow-coincidence',
                                          'parallelization-factor',
                                          [findcoinc_exe.ifo_string] + tags))

    statmap_files = []
    bg_files = FileList()
    for i in range(factor):
        group_str = '%s/%s' % (i, factor)
        coinc_node = findcoinc_exe.create_node(trig_files, hdfbank,
                                               stat_files,
                                               veto_file, veto_name,
                                               group_str,
                                               pivot_ifo,
                                               fixed_ifo,
                                               tags=['JOB'+str(i)])
        bg_files += coinc_node.output_files
        workflow.add_node(coinc_node)

    statmap_files = setup_statmap(workflow, ifos, bg_files,
                                  out_dir, tags=tags)

    logging.info('...leaving coincidence ')
    return statmap_files


def select_files_by_ifo_combination(ifocomb, insps):
    """
    This function selects single-detector files ('insps') for a given ifo combination
    """
    inspcomb = FileList()
    for ifo, ifile in zip(*insps.categorize_by_attr('ifo')):
        if ifo in ifocomb:
            inspcomb += ifile

    return inspcomb


def get_ordered_ifo_list(ifocomb, ifo_ids):
    """
    This function sorts the combination of ifos (ifocomb) based on the given
    precedence list (ifo_ids dictionary) and returns the first ifo as pivot
    the second ifo as fixed, and the ordered list joined as a string.
    """
    # combination_prec stores precedence info for the detectors in the combination
    combination_prec = {ifo: ifo_ids[ifo] for ifo in ifocomb}
    ordered_ifo_list = sorted(combination_prec, key = combination_prec.get)
    pivot_ifo = ordered_ifo_list[0]
    fixed_ifo = ordered_ifo_list[1]

    return pivot_ifo, fixed_ifo, ''.join(ordered_ifo_list)


def setup_combine_statmap(workflow, final_bg_file_list, bg_file_list,
                          out_dir, tags=None):
    """
    Combine the statmap files into one background file
    """
    if tags is None:
        tags = []
    make_analysis_dir(out_dir)
    logging.info('Setting up combine statmap')

    cstat_exe_name = os.path.basename(workflow.cp.get("executables",
                                                      "combine_statmap"))
    if cstat_exe_name == 'pycbc_combine_statmap':
        cstat_class = PyCBCCombineStatmap
    elif cstat_exe_name == 'pycbc_add_statmap':
        cstat_class = PyCBCAddStatmap
    else:
        raise NotImplementedError('executable should be '
            'pycbc_combine_statmap or pycbc_add_statmap')

    cstat_exe = cstat_class(workflow.cp, 'combine_statmap', ifos=workflow.ifos,
                            tags=tags, out_dir=out_dir)

    if cstat_exe_name == 'pycbc_combine_statmap':
        combine_statmap_node = cstat_exe.create_node(final_bg_file_list)
    elif cstat_exe_name == 'pycbc_add_statmap':
        combine_statmap_node = cstat_exe.create_node(final_bg_file_list,
                                                     bg_file_list)

    workflow.add_node(combine_statmap_node)
    return combine_statmap_node.output_file


def setup_exclude_zerolag(workflow, statmap_file, other_statmap_files,
                          out_dir, ifos, tags=None):
    """
    Exclude single triggers close to zerolag triggers from forming any
    background events
    """
    if tags is None:
        tags = []
    make_analysis_dir(out_dir)
    logging.info('Setting up exclude zerolag')

    exc_zerolag_exe = PyCBCExcludeZerolag(workflow.cp, 'exclude_zerolag',
                                          ifos=ifos, tags=tags,
                                          out_dir=out_dir)
    exc_zerolag_node = exc_zerolag_exe.create_node(statmap_file,
                                                   other_statmap_files,
                                                   tags=None)
    workflow.add_node(exc_zerolag_node)
    return exc_zerolag_node.output_file


def rerank_coinc_followup(workflow, statmap_file, bank_file, out_dir,
                          tags=None,
                          injection_file=None,
                          ranking_file=None):
    if tags is None:
        tags = []

    make_analysis_dir(out_dir)

    if not workflow.cp.has_section("workflow-rerank"):
        logging.info("No reranking done in this workflow")
        return statmap_file
    else:
        logging.info("Setting up reranking of candidates")

    # Generate reduced data files (maybe this could also be used elsewhere?)
    stores = FileList([])
    for ifo in workflow.ifos:
        make_analysis_dir('strain_files')
        node = Executable(workflow.cp, 'strain_data_reduce', ifos=[ifo],
                          out_dir='strain_files', tags=tags).create_node()
        node.add_opt('--gps-start-time', workflow.analysis_time[0])
        node.add_opt('--gps-end-time', workflow.analysis_time[1])
        if injection_file:
            node.add_input_opt('--injection-file', injection_file)

        fil = node.new_output_file_opt(workflow.analysis_time, '.hdf',
                                       '--output-file')
        stores.append(fil)
        workflow += node

    # Generate trigger input file
    node = Executable(workflow.cp, 'rerank_trigger_input', ifos=workflow.ifos,
                      out_dir=out_dir, tags=tags).create_node()
    node.add_input_opt('--statmap-file', statmap_file)
    node.add_input_opt('--bank-file', bank_file)
    trigfil = node.new_output_file_opt(workflow.analysis_time, '.hdf',
                                   '--output-file')
    workflow += node

    # Parallelize coinc trigger followup
    factor = int(workflow.cp.get_opt_tags("workflow-rerank",
                                          "parallelization-factor", tags))
    exe = Executable(workflow.cp, 'coinc_followup', ifos=workflow.ifos,
                     out_dir=out_dir, tags=tags)

    stat_files = FileList([])
    for i in range(factor):
        node = exe.create_node()
        node.new_output_file_opt(workflow.analysis_time, '.hdf',
                                 '--output-file', tags=[str(i)])
        node.add_multiifo_input_list_opt('--hdf-store', stores)
        node.add_input_opt('--input-file', trigfil)
        node.add_opt('--start-index', str(i))
        node.add_opt('--stride', factor)
        workflow += node
        stat_files += node.output_files

    exe = Executable(workflow.cp, 'rerank_coincs', ifos=workflow.ifos,
                     out_dir=out_dir, tags=tags)
    node = exe.create_node()
    node.add_input_list_opt('--stat-files', stat_files)
    node.add_input_opt('--statmap-file', statmap_file)
    node.add_input_opt('--followup-file', trigfil)

    if ranking_file:
        node.add_input_opt('--ranking-file', ranking_file)

    node.new_output_file_opt(workflow.analysis_time, '.hdf',
                             '--output-file')
    workflow += node
    return node.output_file
