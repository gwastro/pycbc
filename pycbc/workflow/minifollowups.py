# Copyright (C) 2015 Christopher M. Biwer, Alexander Harvey Nitz
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
from pycbc.workflow.core import Executable, FileList, Node, makedir
from pycbc.workflow.plotting import PlotExecutable, requirestr, excludestr
from itertools import izip_longest

def grouper(iterable, n, fillvalue=None):
    """ Create a list of n length tuples
    """
    args = [iter(iterable)] * n
    return izip_longest(*args, fillvalue=fillvalue)

def setup_minifollowups(workflow, coinc_file, single_triggers, tmpltbank_file, 
                       insp_segs, insp_seg_name, out_dir, tags=None):
    """ Create plots that followup the Nth loudest coincident injection
    from a statmap produced HDF file.
    
    Parameters
    ----------
    workflow: pycbc.workflow.Workflow
        The core workflow instance we are populating
    coinc_file: 
    single_triggers: list of pycbc.workflow.File
        A list cointaining the file objects associated with the merged
        single detector trigger files for each ifo.
    tmpltbank_file: pycbc.workflow.File
        The file object pointing to the HDF format template bank
    insp_segs: dict
        A dictionary, keyed by ifo name, of the data read by each inspiral job.
    insp_segs_name: str 
        The name of the segmentlist to read from the inspiral segment file
    out_dir: path
        The directory to store minifollowups result plots and files
    tags: {None, optional}
        Tags to add to the minifollowups executables
    
    Returns
    -------
    layout: list
        A list of tuples which specify the displayed file layout for the 
        minifollops plots.
    """
    logging.info('Entering minifollowups module')

    # check if minifollowups section exists
    # if not then do not do add minifollowup jobs to the workflow
    if not workflow.cp.has_section('workflow-minifollowups'):
        logging.info('There is no [workflow-minifollowups] section in configuration file')
        logging.info('Leaving minifollowups')
        return []
        
    tags = [] if tags is None else tags
    makedir(out_dir)

    # create a FileList that will contain all output files
    layout = []

    # loop over number of loudest events to be followed up
    num_events = int(workflow.cp.get_opt_tags('workflow-minifollowups', 'num-events', ''))
    for num_event in range(num_events):
        files = FileList([])
        layout += (make_coinc_info(workflow, single_triggers, tmpltbank_file,
                                  coinc_file, num_event, 
                                  out_dir, tags=tags + [str(num_event)]),)        
        files += make_trigger_timeseries(workflow, single_triggers,
                                  coinc_file, num_event, 
                                  out_dir, tags=tags + [str(num_event)])
        files += make_single_template_plots(workflow, insp_segs,
                                  insp_seg_name, coinc_file, tmpltbank_file,
                                  num_event, out_dir, tags=tags + [str(num_event)])
        layout += list(grouper(files, 2))
        num_event += 1
    logging.info('Leaving minifollowups module')

    return layout

def make_single_template_plots(workflow, segs, seg_name, coinc, bank, num, out_dir, 
                               exclude=None, require=None, tags=None):
    tags = [] if tags is None else tags
    makedir(out_dir)
    name = 'single_template_plot'
    secs = requirestr(workflow.cp.get_subsections(name), require)  
    secs = excludestr(secs, exclude)
    files = FileList([])
    for tag in secs:
        for ifo in workflow.ifos:
            # Reanalyze the time around the trigger in each detector
            node = PlotExecutable(workflow.cp, 'single_template', ifos=[ifo],
                                 out_dir=out_dir, tags=[tag] + tags).create_node()
            node.add_input_opt('--statmap-file', coinc)
            node.add_opt('--n-loudest', str(num))
            node.add_input_opt('--inspiral-segments', segs[ifo])
            node.add_opt('--segment-name', seg_name)
            node.add_input_opt('--bank-file', bank)
            node.new_output_file_opt(workflow.analysis_time, '.hdf', '--output-file')
            data = node.output_files[0]
            workflow += node
            
            # Make the plot for this trigger and detector
            node = PlotExecutable(workflow.cp, name, ifos=[ifo],
                                  out_dir=out_dir, tags=[tag] + tags).create_node()
            node.add_input_opt('--single-template-file', data)
            node.new_output_file_opt(workflow.analysis_time, '.png', '--output-file')
            workflow += node
            files += node.output_files
    return files      

def make_coinc_info(workflow, singles, bank, coinc, num, out_dir,
                    exclude=None, require=None, tags=None):
    tags = [] if tags is None else tags
    makedir(out_dir)
    name = 'page_coincinfo'
    secs = requirestr(workflow.cp.get_subsections(name), require)  
    secs = excludestr(secs, exclude)
    files = FileList([])
    node = PlotExecutable(workflow.cp, name, ifos=workflow.ifos,
                              out_dir=out_dir, tags=tags).create_node()
    node.add_input_list_opt('--single-trigger-files', singles)
    node.add_input_opt('--statmap-file', coinc)
    node.add_input_opt('--bank-file', bank)
    node.add_opt('--n-loudest', str(num))
    node.new_output_file_opt(workflow.analysis_time, '.html', '--output-file')
    workflow += node
    files += node.output_files
    return files
    
def make_trigger_timeseries(workflow, singles, coinc, num, out_dir,
                            exclude=None, require=None, tags=None):
    tags = [] if tags is None else tags
    makedir(out_dir)
    name = 'plot_trigger_timeseries'
    secs = requirestr(workflow.cp.get_subsections(name), require)  
    secs = excludestr(secs, exclude)
    files = FileList([])
    for tag in secs:
        node = PlotExecutable(workflow.cp, name, ifos=workflow.ifos,
                              out_dir=out_dir, tags=[tag] + tags).create_node()
        node.add_input_list_opt('--single-trigger-files', singles)
        node.add_input_opt('--statmap-file', coinc)
        node.add_opt('--n-loudest', str(num))
        node.new_output_file_opt(workflow.analysis_time, '.png', '--output-file')
        workflow += node
        files += node.output_files
    return files



