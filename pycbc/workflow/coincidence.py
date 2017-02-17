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

import logging
from pycbc.workflow.core import FileList, make_analysis_dir, Executable, Node, File
from glue import segments

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
    file_input_options = ['--statistic-files']
    def create_node(self, trig_files, bank_file, stat_files, veto_file,
                    veto_name, template_str, tags=None):
        if tags is None:
            tags = []
        segs = trig_files.get_times_covered_by_files()
        seg = segments.segment(segs[0][0], segs[-1][1])
        node = Node(self)
        node.set_memory(10000)
        node.add_input_opt('--template-bank', bank_file)
        node.add_input_list_opt('--trigger-files', trig_files)
        if len(stat_files) > 0:
            node.add_input_list_opt('--statistic-files', stat_files)
        if veto_file is not None:
            node.add_input_opt('--veto-files', veto_file)
            node.add_opt('--segment-name', veto_name)
        node.add_opt('--template-fraction-range', template_str)
        node.new_output_file_opt(seg, '.hdf', '--output-file', tags=tags)
        return node

class PyCBCStatMapExecutable(Executable):
    """Calculate FAP, IFAR, etc"""
    current_retention_level = Executable.MERGED_TRIGGERS
    def create_node(self, coinc_files, tags=None):
        if tags is None:
            tags = []
        segs = coinc_files.get_times_covered_by_files()
        seg = segments.segment(segs[0][0], segs[-1][1])

        node = Node(self)
        node.set_memory(5000)
        node.add_input_list_opt('--coinc-files', coinc_files)
        node.new_output_file_opt(seg, '.hdf', '--output-file', tags=tags)
        return node
        
class PyCBCStatMapInjExecutable(Executable):
    """Calculate FAP, IFAR, etc"""
    current_retention_level = Executable.MERGED_TRIGGERS
    def create_node(self, zerolag, full_data, injfull, fullinj, tags=None):
        if tags is None:
            tags = []
        segs = zerolag.get_times_covered_by_files()
        seg = segments.segment(segs[0][0], segs[-1][1])

        node = Node(self)
        node.set_memory(5000)
        node.add_input_list_opt('--zero-lag-coincs', zerolag)
        node.add_input_list_opt('--full-data-background', full_data)
        node.add_input_list_opt('--mixed-coincs-inj-full', injfull)
        node.add_input_list_opt('--mixed-coincs-full-inj', fullinj)
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
        node.new_output_file_opt(inj_xml_file[0].segment, '.hdf', '--output-file', 
                                 tags=tags)
        return node

class PyCBCDistributeBackgroundBins(Executable):
    """Distribute coinc files amoung different background bins"""
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
    current_retention_level = Executable.MERGED_TRIGGERS
    def create_node(self, statmap_files, tags=None):
        if tags is None:
            tags = []
        node = Node(self)
        node.add_input_list_opt('--statmap-files', statmap_files)
        node.new_output_file_opt(statmap_files[0].segment, '.hdf',
                                 '--output-file', tags=tags)
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

def setup_trigger_fitting(workflow, insps, hdfbank, veto_file, veto_name):
    if not workflow.cp.has_option('workflow-coincidence', 'do-trigger-fitting'):
        return FileList()
    else:
        assert len(hdfbank) == 1  # must be a list with exactly 1 bank file
        assert len(veto_file) == 1
        assert len(veto_name) == 1
        smoothed_fit_files = FileList()
        for i in workflow.ifos:
            ifo_insp = [insp for insp in insps if (insp.ifo == i)]
            assert len(ifo_insp)==1
            raw_node = PyCBCFitByTemplateExecutable(workflow.cp,
                'fit_by_template', ifos=i).create_node(ifo_insp[0], hdfbank[0],
                                                    veto_file[0], veto_name[0])
            workflow += raw_node
            smooth_node = PyCBCFitOverParamExecutable(workflow.cp,
                'fit_over_param', ifos=i).create_node(raw_node.output_files[0],
                                                                    hdfbank[0])
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
    node = exe.create_node(inj_coinc_file, inj_xml_file, veto_file, veto_name, tags)
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

    ifos, insp_groups = xml_trigger_files.categorize_by_attr('ifo')
    trig_files = FileList()
    for ifo, insp_group in zip(ifos,  insp_groups):
        trig2hdf_exe = PyCBCTrig2HDFExecutable(workflow.cp, 'trig2hdf',
                                       ifos=ifo, out_dir=out_dir, tags=tags)
        segs, insp_bundles = insp_group.categorize_by_attr('segment')
        for insps in  insp_bundles:
            trig2hdf_node =  trig2hdf_exe.create_node(insps, hdfbank[0])
            workflow.add_node(trig2hdf_node)
            trig_files += trig2hdf_node.output_files
    return trig_files

def setup_statmap(workflow, coinc_files, bank_file, out_dir, tags=None):
    tags = [] if tags is None else tags
    if workflow.cp.has_option_tags('workflow-coincidence', 'background-bins', tags):
        return setup_background_bins(workflow, coinc_files, bank_file, out_dir, tags=tags)
    else:
        return setup_simple_statmap(workflow, coinc_files, out_dir, tags=tags)
        
def setup_simple_statmap(workflow, coinc_files, out_dir, tags=None):
    tags = [] if tags is None else tags

    statmap_exe = PyCBCStatMapExecutable(workflow.cp, 'statmap',
                                              ifos=workflow.ifos,
                                              tags=tags, out_dir=out_dir) 

    stat_node = statmap_exe.create_node(coinc_files, tags=tags)
    workflow.add_node(stat_node)
    return stat_node.output_files[0], stat_node.output_files

def setup_background_bins(workflow, coinc_files, bank_file, out_dir, tags=None):
    tags = [] if tags is None else tags
    
    bins_exe = PyCBCDistributeBackgroundBins(workflow.cp, 'distribute_background_bins', 
                                       ifos=workflow.ifos, tags=tags, out_dir=out_dir)
                                       
    statmap_exe = PyCBCStatMapExecutable(workflow.cp, 'statmap',
                                              ifos=workflow.ifos,
                                              tags=tags, out_dir=out_dir)       

    cstat_exe = PyCBCCombineStatmap(workflow.cp, 'combine_statmap', ifos=workflow.ifos,
                                    tags=tags, out_dir=out_dir)                       
           
    background_bins = workflow.cp.get_opt_tags('workflow-coincidence', 'background-bins', tags).split(' ')             
    background_bins = [x for x in background_bins if x != '']
    bins_node = bins_exe.create_node(coinc_files, bank_file, background_bins)
    workflow += bins_node

    statmap_files = FileList([])
    for i, coinc_file in enumerate(bins_node.output_files):
        statnode = statmap_exe.create_node(FileList([coinc_file]), tags=tags + ['BIN_%s' % i])
        workflow += statnode
        statmap_files.append(statnode.output_files[0])
        statmap_files[i].bin_name = bins_node.names[i]
    
    cstat_node = cstat_exe.create_node(statmap_files, tags=tags)
    workflow += cstat_node
    
    return cstat_node.output_files[0], statmap_files
    
def setup_statmap_inj(workflow, coinc_files, background_file, bank_file, out_dir, tags=None):
    tags = [] if tags is None else tags
    if workflow.cp.has_option_tags('workflow-coincidence', 'background-bins', tags):
        return setup_background_bins_inj(workflow, coinc_files, background_file, bank_file, out_dir, tags=tags)
    else:
        return setup_simple_statmap_inj(workflow, coinc_files, background_file, out_dir, tags=tags)
        
def setup_simple_statmap_inj(workflow, coinc_files, background_file, out_dir, tags=None):
    tags = [] if tags is None else tags

    statmap_exe = PyCBCStatMapInjExecutable(workflow.cp, 'statmap_inj',
                                              ifos=workflow.ifos,
                                              tags=tags, out_dir=out_dir) 

    stat_node = statmap_exe.create_node(FileList(coinc_files['injinj']), background_file, 
                                     FileList(coinc_files['injfull']), FileList(coinc_files['fullinj']))
    workflow.add_node(stat_node)
    return stat_node.output_files[0]

def setup_background_bins_inj(workflow, coinc_files, background_file, bank_file, out_dir, tags=None):
    tags = [] if tags is None else tags
    
    bins_exe = PyCBCDistributeBackgroundBins(workflow.cp, 'distribute_background_bins', 
                                       ifos=workflow.ifos, tags=tags, out_dir=out_dir)
                                       
    statmap_exe = PyCBCStatMapInjExecutable(workflow.cp, 'statmap_inj',
                                              ifos=workflow.ifos,
                                              tags=tags, out_dir=out_dir)   
    
    cstat_exe = PyCBCCombineStatmap(workflow.cp, 'combine_statmap', ifos=workflow.ifos,
                                    tags=tags, out_dir=out_dir)                       
           
    background_bins = workflow.cp.get_opt_tags('workflow-coincidence', 'background-bins', tags).split(' ')   
    background_bins = [x for x in background_bins if x != '']
    
    for inj_type in ['injinj', 'injfull', 'fullinj']:          
        bins_node = bins_exe.create_node(FileList(coinc_files[inj_type]), bank_file, background_bins, tags=tags + [inj_type])
        workflow += bins_node
        coinc_files[inj_type] = bins_node.output_files
    
    statmap_files = FileList([])
    for i in range(len(background_bins)):
        statnode = statmap_exe.create_node(FileList([coinc_files['injinj'][i]]), FileList([background_file[i]]), 
                                     FileList([coinc_files['injfull'][i]]), FileList([coinc_files['fullinj'][i]]), 
                                     tags=tags + ['BIN_%s' % i])
        workflow += statnode
        statmap_files.append(statnode.output_files[0])
        
    cstat_node = cstat_exe.create_node(statmap_files, tags=tags)
    workflow += cstat_node
    
    return cstat_node.output_files[0]

def setup_interval_coinc_inj(workflow, hdfbank, full_data_trig_files, inj_trig_files,
              stat_files, background_file, veto_file, veto_name, out_dir, tags=None):
    """
    This function sets up exact match coincidence and background estimation 

    using a folded interval technique.
    """
    if tags is None:
        tags = []
    make_analysis_dir(out_dir)
    logging.info('Setting up coincidence for injection')

    if len(hdfbank) > 1:
        raise ValueError('This coincidence method only supports a '
                         'pregenerated template bank')
    hdfbank = hdfbank[0]

    if len(workflow.ifos) > 2:
        raise ValueError('This coincidence method only supports two ifo searches')

    # Wall time knob and memory knob
    factor = int(workflow.cp.get_opt_tags('workflow-coincidence', 'parallelization-factor', tags))
    
    ffiles = {}
    ifiles = {}
    ifos, files = full_data_trig_files.categorize_by_attr('ifo')
    for ifo, file in zip(ifos, files):
        ffiles[ifo] = file[0]
    ifos, files = inj_trig_files.categorize_by_attr('ifo')
    for ifo, file in zip(ifos, files):
        ifiles[ifo] = file[0]
    ifo0, ifo1 = ifos[0], ifos[1]
    combo = [(FileList([ifiles[ifo0], ifiles[ifo1]]), "injinj"),
             (FileList([ifiles[ifo0], ffiles[ifo1]]), "injfull"),
             (FileList([ifiles[ifo1], ffiles[ifo0]]), "fullinj"),
            ]
    bg_files = {'injinj':[],'injfull':[],'fullinj':[]}

    for trig_files, ctag in combo:
        findcoinc_exe = PyCBCFindCoincExecutable(workflow.cp, 'coinc',
                                              ifos=workflow.ifos,
                                              tags=tags + [ctag], out_dir=out_dir)
        for i in range(factor):
            group_str = '%s/%s' % (i, factor)
            coinc_node = findcoinc_exe.create_node(trig_files, hdfbank,
                                                   stat_files,
                                                   veto_file, veto_name,
                                                   group_str,
                                                   tags=[str(i)])
            bg_files[ctag] += coinc_node.output_files
            workflow.add_node(coinc_node)

    return setup_statmap_inj(workflow, bg_files, background_file, hdfbank, out_dir, tags=tags)
    
def setup_interval_coinc(workflow, hdfbank, trig_files, stat_files,
                         veto_files, veto_names, out_dir, tags=None):
    """
    This function sets up exact match coincidence and background estimation

    using a folded interval technique.
    """
    if tags is None:
        tags = []
    make_analysis_dir(out_dir)
    logging.info('Setting up coincidence')

    if len(hdfbank) != 1:
        raise ValueError('Must use exactly 1 bank file for this coincidence '
                         'method, I got %i !' % len(hdfbank))
    hdfbank = hdfbank[0]

    if len(workflow.ifos) > 2:
        raise ValueError('This coincidence method only supports two ifo searches')

    findcoinc_exe = PyCBCFindCoincExecutable(workflow.cp, 'coinc',
                                             ifos=workflow.ifos,
                                             tags=tags, out_dir=out_dir)
                                         
    # Wall time knob and memory knob
    factor = int(workflow.cp.get_opt_tags('workflow-coincidence', 'parallelization-factor', tags))

    statmap_files = []
    for veto_file, veto_name in zip(veto_files, veto_names):
        bg_files = FileList()
        for i in range(factor):
            group_str = '%s/%s' % (i, factor)
            coinc_node = findcoinc_exe.create_node(trig_files, hdfbank,
                                                   stat_files,
                                                   veto_file, veto_name,
                                                   group_str,
                                                   tags=[veto_name, str(i)])
            bg_files += coinc_node.output_files
            workflow.add_node(coinc_node)
             
        statmap_files += [setup_statmap(workflow, bg_files, hdfbank, out_dir, tags=tags + [veto_name])]

    logging.info('...leaving coincidence ')
    return statmap_files
