# Copyright (C) 2015  Alex Nitz
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
This module is responsible for setting up plotting jobs.
https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/NOTYETCREATED.html
"""
from pycbc.workflow.core import FileList, makedir, Executable, Node
 
class PlotExecutable(Executable):
    """ This converts xml tmpltbank to an hdf format
    """
    current_retention_level = Executable.FINAL_RESULT
 
def make_segments_plot(workflow, seg_files, out_dir, tags=[]):
    makedir(out_dir)
    node = PlotExecutable(workflow.cp, 'plot_segments', ifos=workflow.ifos,
                         out_dir=out_dir, tags=tags).create_node()
    node.add_input_list_opt('--segment-files', seg_files)
    node.new_output_file_opt(workflow.analysis_time, '.html', '--output-file')
    workflow += node
        
def make_foreground_table(workflow, trig_file, bank_file, ftag, out_dir, tags=[]):
    makedir(out_dir)
    node = PlotExecutable(workflow.cp, 'page_foreground', ifos=workflow.ifos,
                    out_dir=out_dir, tags=tags).create_node()
    node.add_input_opt('--bank-file', bank_file)
    node.add_opt('--foreground-tag', ftag)
    node.add_input_opt('--trigger-file', trig_file)
    node.new_output_file_opt(bank_file.segment, '.html', '--output-file')
    workflow += node

def make_sensitivity_plot(workflow, inj_file, out_dir, tags=[]):
    makedir(out_dir)   
    for tag in workflow.cp.get_subsections('plot_sensitivity'):
        node = PlotExecutable(workflow.cp, 'plot_sensitivity', ifos=workflow.ifos,
                    out_dir=out_dir, tags=[tag] + tags).create_node()
        node.add_input_opt('--injection-file', inj_file)
        node.new_output_file_opt(inj_file.segment, '.png', '--output-file')
        workflow += node

def make_coinc_snrchi_plot(workflow, inj_file, inj_trig, stat_file, trig_file, out_dir, tags=[]):
    makedir(out_dir)    
    for tag in workflow.cp.get_subsections('plot_coinc_snrchi'):
        node = PlotExecutable(workflow.cp, 'plot_coinc_snrchi', ifos=inj_trig.ifo,
                    out_dir=out_dir, tags=[tag] + tags).create_node()
        node.add_input_opt('--found-injection-file', inj_file)
        node.add_input_opt('--single-injection-file', inj_trig)
        node.add_input_opt('--coinc-statistic-file', stat_file)
        node.add_input_opt('--single-trigger-file', trig_file)
        node.new_output_file_opt(inj_file.segment, '.png', '--output-file')
        workflow += node

def make_inj_table(workflow, inj_file, out_dir, tags=[]):
    makedir(out_dir)
    node = PlotExecutable(workflow.cp, 'page_injections', ifos=workflow.ifos,
                    out_dir=out_dir, tags=tags).create_node()
    node.add_input_opt('--injection-file', inj_file)
    node.new_output_file_opt(inj_file.segment, '.html', '--output-file')
    workflow += node   

def make_snrchi_plot(workflow, trig_files, veto_file, out_dir, tags=[]):
    makedir(out_dir)    
    for tag in workflow.cp.get_subsections('plot_snrchi'):
        for trig_file in trig_files:
            node = PlotExecutable(workflow.cp, 'plot_snrchi',
                        ifos=trig_file.ifo, 
                        out_dir=out_dir, 
                        tags=[tag] + tags).create_node()

            node.set_memory(15000)
            node.add_input_opt('--trigger-file', trig_file)
            node.add_input_opt('--veto-file', veto_file)
            node.new_output_file_opt(trig_file.segment, '.png', '--output-file')
            workflow += node  

def make_foundmissed_plot(workflow, inj_file, out_dir, tags=[]):
    makedir(out_dir)   
    for tag in workflow.cp.get_subsections('plot_foundmissed'):
        exe = PlotExecutable(workflow.cp, 'plot_foundmissed', ifos=workflow.ifos,
                    out_dir=out_dir, tags=[tag] + tags)
        node = exe.create_node()        
        ext = '.html' if exe.has_opt('dynamic') else '.png'
        node.add_input_opt('--injection-file', inj_file)
        node.new_output_file_opt(inj_file.segment, ext, '--output-file')
        workflow += node   
    
def make_snrifar_plot(workflow, bg_file, out_dir, tags=[]):
    makedir(out_dir)
    node = PlotExecutable(workflow.cp, 'plot_snrifar', ifos=workflow.ifos,
                out_dir=out_dir, tags=tags).create_node()
    node.add_input_opt('--trigger-file', bg_file)
    node.new_output_file_opt(bg_file.segment, '.png', '--output-file')
    workflow += node
    
def make_results_web_page(workflow, results_dir):
    import pycbc.results
    template_path =  pycbc.results.__path__[0] + '/templates/concrete.html'

    out_dir = workflow.cp.get('results_page', 'output-path')
    makedir(out_dir)
    node = PlotExecutable(workflow.cp, 'results_page', ifos=workflow.ifos,
                out_dir=out_dir).create_node()
    node.add_opt('--plots-dir', results_dir)
    node.add_opt('--template-file', template_path)
    workflow += node

def make_singles_plot(workflow, trig_files, bank_file, veto_file, out_dir, tags=[]):
    makedir(out_dir)    
    for tag in workflow.cp.get_subsections('plot_singles'):
        for trig_file in trig_files:
            node = PlotExecutable(workflow.cp, 'plot_singles',
                        ifos=trig_file.ifo, 
                        out_dir=out_dir, 
                        tags=[tag] + tags).create_node()

            node.set_memory(15000)
            node.add_input_opt('--bank-file', bank_file)
            node.add_input_opt('--veto-file', veto_file)
            node.add_opt('--detector', trig_file.ifo)
            node.add_input_opt('--single-trig-file', trig_file)
            node.new_output_file_opt(trig_file.segment, '.png', '--output-file')
            workflow += node  
