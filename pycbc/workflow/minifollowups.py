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

import logging, os.path
from pycbc.workflow.core import Executable, FileList, Node, makedir, File, Workflow
from pycbc.workflow.plotting import PlotExecutable, requirestr, excludestr
from itertools import izip_longest
from Pegasus import DAX3 as dax
from pycbc.workflow import pegasus_workflow as wdax

def grouper(iterable, n, fillvalue=None):
    """ Create a list of n length tuples
    """
    args = [iter(iterable)] * n
    return izip_longest(*args, fillvalue=fillvalue)

def setup_foreground_minifollowups(workflow, coinc_file, single_triggers, tmpltbank_file, 
                       insp_segs, insp_seg_name, dax_output, out_dir, tags=None):
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
    
    if not workflow.cp.has_section('workflow-minifollowups'):
        logging.info('There is no [workflow-minifollowups] section in configuration file')
        logging.info('Leaving minifollowups')
    
    tags = [] if tags is None else tags
    makedir(dax_output)
    
    # turn the config file into a File class
    config_path = dax_output + '/' + '_'.join(tags) + 'foreground_minifollowup.ini'
    workflow.cp.write(open(config_path, 'w'))
    
    config_file = wdax.File(os.path.basename(config_path))
    config_file.PFN(config_path, 'local')
    
    exe = Executable(workflow.cp, 'foreground_minifollowup', ifos=workflow.ifos, out_dir=dax_output)
    
    node = exe.create_node()
    node.add_input_opt('--config-files', config_file)
    node.add_input_opt('--bank-file', tmpltbank_file)
    node.add_input_opt('--statmap-file', coinc_file)
    node.add_multiifo_input_list_opt('--single-detector-triggers', single_triggers)
    node.add_multiifo_input_list_opt('--inspiral-segments', insp_segs.values())
    node.add_opt('--inspiral-segment-name', insp_seg_name)
    node.new_output_file_opt(workflow.analysis_time, '.dax', '--output-file', tags=tags)
    
    name = node.output_files[0].name
    map_loc = name + '.map'
    node.add_opt('--workflow-name', name)
    node.add_opt('--output-dir', out_dir)
    node.add_opt('--output-map', map_loc)
    
    workflow += node
    
    # execute this is a sub-workflow
    fil = node.output_files[0]
    
    ## FIXME not clear why I have to set the id here, pegasus should do this!
    job = dax.DAX(fil, id=id(node))
    Workflow.set_job_properties(job, map_loc)
    workflow._adag.addJob(job)    
    dep = dax.Dependency(parent=node._dax_node, child=job)
    
    logging.info('Leaving minifollowups module')

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



