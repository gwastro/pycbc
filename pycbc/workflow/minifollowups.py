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
        return
    
    tags = [] if tags is None else tags
    makedir(dax_output)
    
    # turn the config file into a File class
    config_path = os.path.abspath(dax_output + '/' + '_'.join(tags) + 'foreground_minifollowup.ini')
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
    
    job = dax.DAX(fil)
    job.addArguments('--basename %s' % os.path.splitext(os.path.basename(name))[0])
    Workflow.set_job_properties(job, map_loc)
    workflow._adag.addJob(job)
    dep = dax.Dependency(parent=node._dax_node, child=job)
    workflow._adag.addDependency(dep)
    logging.info('Leaving minifollowups module')

def setup_single_det_minifollowups(workflow, single_trig_file, tmpltbank_file,
                                  insp_segs, insp_seg_name, dax_output,
                                  out_dir, veto_file=None,
                                  veto_segment_name=None, tags=None):
    """ Create plots that followup the Nth loudest clustered single detector
    triggers from a merged single detector trigger HDF file.
    
    Parameters
    ----------
    workflow: pycbc.workflow.Workflow
        The core workflow instance we are populating
    single_trig_file: pycbc.workflow.File
        The File class holding the single detector triggers.
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
        msg = 'There is no [workflow-minifollowups] section in '
        msg += 'configuration file'
        logging.info(msg)
        logging.info('Leaving minifollowups')
        return

    tags = [] if tags is None else tags
    makedir(dax_output)

    # turn the config file into a File class
    curr_ifo = single_trig_file.ifo
    config_path = os.path.abspath(dax_output + '/' + curr_ifo + \
                                   '_'.join(tags) + 'singles_minifollowup.ini')
    workflow.cp.write(open(config_path, 'w'))

    config_file = wdax.File(os.path.basename(config_path))
    config_file.PFN(config_path, 'local')

    exe = Executable(workflow.cp, 'singles_minifollowup',
                     ifos=curr_ifo, out_dir=dax_output)

    node = exe.create_node()
    node.add_input_opt('--config-files', config_file)
    node.add_input_opt('--bank-file', tmpltbank_file)
    node.add_input_opt('--single-detector-file', single_trig_file)
    node.add_input_opt('--inspiral-segments', insp_segs[curr_ifo])
    node.add_opt('--inspiral-segment-name', insp_seg_name)
    node.add_opt('--instrument', curr_ifo)
    if veto_file is not None:
        assert(veto_segment_name is not None)
        node.add_input_opt('--veto-file', veto_file)
        node.add_opt('--veto-segment-name', veto_segment_name)
    node.new_output_file_opt(workflow.analysis_time, '.dax', '--output-file',
                             tags=tags)

    name = node.output_files[0].name
    map_loc = name + '.map'
    node.add_opt('--workflow-name', name)
    node.add_opt('--output-dir', out_dir)
    node.add_opt('--output-map', map_loc)

    workflow += node

    # execute this is a sub-workflow
    fil = node.output_files[0]

    job = dax.DAX(fil)
    job.addArguments('--basename %s' \
                     % os.path.splitext(os.path.basename(name))[0])
    Workflow.set_job_properties(job, map_loc)
    workflow._adag.addJob(job)
    dep = dax.Dependency(parent=node._dax_node, child=job)
    workflow._adag.addDependency(dep)
    logging.info('Leaving minifollowups module')


def setup_injection_minifollowups(workflow, injection_file, inj_xml_file,
                                  single_triggers, 
                                  tmpltbank_file, insp_segs, insp_seg_name,
                                  dax_output, out_dir, tags=None):
    """ Create plots that followup the closest missed injections
    
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
    logging.info('Entering injection minifollowups module')
    
    if not workflow.cp.has_section('workflow-injection_minifollowups'):
        logging.info('There is no [workflow-injection_minifollowups] section in configuration file')
        logging.info('Leaving minifollowups')
        return
    
    tags = [] if tags is None else tags
    makedir(dax_output)
    
    # turn the config file into a File class
    config_path = os.path.abspath(dax_output + '/' + '_'.join(tags) + 'injection_minifollowup.ini')
    workflow.cp.write(open(config_path, 'w'))
    
    config_file = wdax.File(os.path.basename(config_path))
    config_file.PFN(config_path, 'local')
    
    exe = Executable(workflow.cp, 'injection_minifollowup', ifos=workflow.ifos, out_dir=dax_output)
    
    node = exe.create_node()
    node.add_input_opt('--config-files', config_file)
    node.add_input_opt('--bank-file', tmpltbank_file)
    node.add_input_opt('--injection-file', injection_file)
    node.add_input_opt('--injection-xml-file', inj_xml_file)
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
    
    job = dax.DAX(fil)
    job.addArguments('--basename %s' % os.path.splitext(os.path.basename(name))[0])
    Workflow.set_job_properties(job, map_loc)
    workflow._adag.addJob(job)
    dep = dax.Dependency(parent=node._dax_node, child=job)
    workflow._adag.addDependency(dep)
    logging.info('Leaving injection minifollowups module')

def make_single_template_plots(workflow, segs, seg_name, params,
                                   out_dir, inj_file=None, exclude=None,
                                   require=None, tags=None, params_str=None):
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
            node.add_opt('--mass1', "%.6f" % params['mass1'])
            node.add_opt('--mass2', "%.6f" % params['mass2'])
            node.add_opt('--spin1z',"%.6f" % params['spin1z'])
            node.add_opt('--spin2z',"%.6f" % params['spin2z'])
            # str(numpy.float64) restricts to 2d.p. BE CAREFUL WITH THIS!!!
            str_trig_time = '%.6f' %(params[ifo + '_end_time'])
            node.add_opt('--trigger-time', str_trig_time)
            node.add_input_opt('--inspiral-segments', segs[ifo])
            if inj_file is not None:
                node.add_input_opt('--injection-file', inj_file)
            node.add_opt('--segment-name', seg_name)
            node.new_output_file_opt(workflow.analysis_time, '.hdf',
                                     '--output-file', store_file=False)
            data = node.output_files[0]
            workflow += node
            # Make the plot for this trigger and detector
            node = PlotExecutable(workflow.cp, name, ifos=[ifo],
                              out_dir=out_dir, tags=[tag] + tags).create_node()
            node.add_input_opt('--single-template-file', data)
            node.new_output_file_opt(workflow.analysis_time, '.png',
                                     '--output-file')
            title="'%s SNR and chi^2 timeseries" %(ifo) 
            if params_str is not None:
                title+= " using %s" %(params_str)
            title+="'"
            node.add_opt('--plot-title', title)
            caption = "'The SNR and chi^2 timeseries around the injection"
            if params_str is not None:
                caption += " using %s" %(params_str)
            caption += ". The template used has the following parameters: "
            caption += "mass1=%s, mass2=%s, spin1z=%s, spin2z=%s'"\
                       %(params['mass1'], params['mass2'], params['spin1z'],
                         params['spin2z'])
            node.add_opt('--plot-caption', caption)
            workflow += node
            files += node.output_files
    return files

def make_inj_info(workflow, injection_file, injection_index, num, out_dir,
                  tags=None):
    tags = [] if tags is None else tags
    makedir(out_dir)
    name = 'page_injinfo'
    files = FileList([])
    node = PlotExecutable(workflow.cp, name, ifos=workflow.ifos,
                              out_dir=out_dir, tags=tags).create_node()
    node.add_input_opt('--injection-file', injection_file)
    node.add_opt('--injection-index', str(injection_index))
    node.add_opt('--n-nearest', str(num))
    node.new_output_file_opt(workflow.analysis_time, '.html', '--output-file')
    workflow += node
    files += node.output_files
    return files

def make_coinc_info(workflow, singles, bank, coinc, num, out_dir, tags=None):
    tags = [] if tags is None else tags
    makedir(out_dir)
    name = 'page_coincinfo'
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

def make_sngl_ifo(workflow, sngl_file, bank_file, num, out_dir, ifo,
                  veto_file=None, veto_segment_name=None, tags=None):
    """Setup a job to create sngl detector sngl ifo html summary snippet.
    """
    tags = [] if tags is None else tags
    makedir(out_dir)
    name = 'page_snglinfo'
    files = FileList([])
    node = PlotExecutable(workflow.cp, name, ifos=[ifo],
                              out_dir=out_dir, tags=tags).create_node()
    node.add_input_opt('--single-trigger-file', sngl_file)
    node.add_input_opt('--bank-file', bank_file)
    if veto_file is not None:
        assert(veto_segment_name is not None)
        node.add_input_opt('--veto-file', veto_file)
        node.add_opt('--veto-segment-name', veto_segment_name)
    node.add_opt('--n-loudest', str(num))
    node.add_opt('--instrument', ifo)
    node.new_output_file_opt(workflow.analysis_time, '.html', '--output-file')
    workflow += node
    files += node.output_files
    return files


def make_trigger_timeseries(workflow, singles, ifo_times, out_dir, special_tids=None,
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
        node.add_multiifo_input_list_opt('--single-trigger-files', singles)
        node.add_opt('--times', ifo_times)
        node.new_output_file_opt(workflow.analysis_time, '.png', '--output-file')
        
        if special_tids is not None:
            node.add_opt('--special-trigger-ids', special_tids)
        
        workflow += node
        files += node.output_files
    return files

    
def make_singles_timefreq(workflow, single, bank_file, start, end, out_dir,
                          veto_file=None, tags=None):
    tags = [] if tags is None else tags
    makedir(out_dir)
    name = 'plot_singles_timefreq'

    node = PlotExecutable(workflow.cp, name, ifos=workflow.ifos,
                          out_dir=out_dir, tags=tags).create_node()
    node.add_input_opt('--trig-file', single)
    node.add_input_opt('--bank-file', bank_file)
    node.add_opt('--gps-start-time', int(start))
    node.add_opt('--gps-end-time', int(end))
    
    if veto_file:
        node.add_input_opt('--veto-file', veto_file)
        
    node.add_opt('--detector', single.ifo)
    node.new_output_file_opt(workflow.analysis_time, '.png', '--output-file')
    workflow += node
    return node.output_files


