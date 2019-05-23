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
from six.moves.urllib.request import pathname2url
from six.moves.urllib.parse import urljoin
import distutils.spawn
from pycbc.workflow.core import Executable, FileList, Node, makedir, File, Workflow
from pycbc.workflow.plotting import PlotExecutable, requirestr, excludestr
try:
    # Python 3
    from itertools import zip_longest
except ImportError:
    # Python 2
    from itertools import izip_longest as zip_longest
from Pegasus import DAX3 as dax
from pycbc.workflow import pegasus_workflow as wdax

def grouper(iterable, n, fillvalue=None):
    """ Create a list of n length tuples
    """
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)

def setup_foreground_minifollowups(workflow, coinc_file, single_triggers,
                       tmpltbank_file, insp_segs, insp_data_name,
                       insp_anal_name, dax_output, out_dir, tags=None):
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
    insp_segs: SegFile
       The segment file containing the data read and analyzed by each inspiral
       job.
    insp_data_name: str
        The name of the segmentlist storing data read.
    insp_anal_name: str
        The name of the segmentlist storing data analyzed.
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
    config_file.PFN(urljoin('file:', pathname2url(config_path)), site='local')

    exe = Executable(workflow.cp, 'foreground_minifollowup', ifos=workflow.ifos, out_dir=dax_output)

    node = exe.create_node()
    node.add_input_opt('--config-files', config_file)
    node.add_input_opt('--bank-file', tmpltbank_file)
    node.add_input_opt('--statmap-file', coinc_file)
    node.add_multiifo_input_list_opt('--single-detector-triggers', single_triggers)
    node.add_input_opt('--inspiral-segments', insp_segs)
    node.add_opt('--inspiral-data-read-name', insp_data_name)
    node.add_opt('--inspiral-data-analyzed-name', insp_anal_name)
    node.new_output_file_opt(workflow.analysis_time, '.dax', '--output-file', tags=tags)
    node.new_output_file_opt(workflow.analysis_time, '.dax.map', '--output-map', tags=tags)
    node.new_output_file_opt(workflow.analysis_time, '.tc.txt', '--transformation-catalog', tags=tags)

    name = node.output_files[0].name
    map_file = node.output_files[1]
    tc_file = node.output_files[2]

    node.add_opt('--workflow-name', name)
    node.add_opt('--output-dir', out_dir)

    workflow += node

    # execute this in a sub-workflow
    fil = node.output_files[0]

    # determine if a staging site has been specified
    try:
        staging_site = workflow.cp.get('workflow-foreground_minifollowups',
                                       'staging-site')
    except:
        staging_site = None

    job = dax.DAX(fil)
    job.addArguments('--basename %s' % os.path.splitext(os.path.basename(name))[0])
    Workflow.set_job_properties(job, map_file, tc_file, staging_site=staging_site)
    workflow._adag.addJob(job)
    dep = dax.Dependency(parent=node._dax_node, child=job)
    workflow._adag.addDependency(dep)
    logging.info('Leaving minifollowups module')

def setup_single_det_minifollowups(workflow, single_trig_file, tmpltbank_file,
                                  insp_segs, insp_data_name, insp_anal_name,
                                  dax_output, out_dir, veto_file=None,
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
    insp_segs: SegFile
       The segment file containing the data read by each inspiral job.
    insp_data_name: str
        The name of the segmentlist storing data read.
    insp_anal_name: str
        The name of the segmentlist storing data analyzed.
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

    if not workflow.cp.has_section('workflow-sngl_minifollowups'):
        msg = 'There is no [workflow-sngl_minifollowups] section in '
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
    config_file.PFN(urljoin('file:', pathname2url(config_path)), site='local')

    exe = Executable(workflow.cp, 'singles_minifollowup',
                     ifos=curr_ifo, out_dir=dax_output, tags=tags)

    wikifile = curr_ifo + '_'.join(tags) + 'loudest_table.txt'

    node = exe.create_node()
    node.add_input_opt('--config-files', config_file)
    node.add_input_opt('--bank-file', tmpltbank_file)
    node.add_input_opt('--single-detector-file', single_trig_file)
    node.add_input_opt('--inspiral-segments', insp_segs)
    node.add_opt('--inspiral-data-read-name', insp_data_name)
    node.add_opt('--inspiral-data-analyzed-name', insp_anal_name)
    node.add_opt('--instrument', curr_ifo)
    node.add_opt('--wiki-file', wikifile)
    if veto_file is not None:
        assert(veto_segment_name is not None)
        node.add_input_opt('--veto-file', veto_file)
        node.add_opt('--veto-segment-name', veto_segment_name)
    node.new_output_file_opt(workflow.analysis_time, '.dax', '--output-file', tags=tags)
    node.new_output_file_opt(workflow.analysis_time, '.dax.map', '--output-map', tags=tags)
    node.new_output_file_opt(workflow.analysis_time, '.tc.txt', '--transformation-catalog', tags=tags)

    name = node.output_files[0].name
    map_file = node.output_files[1]
    tc_file = node.output_files[2]

    node.add_opt('--workflow-name', name)
    node.add_opt('--output-dir', out_dir)

    workflow += node

    # execute this in a sub-workflow
    fil = node.output_files[0]

    # determine if a staging site has been specified
    try:
        staging_site = workflow.cp.get('workflow-sngl_minifollowups',
                                       'staging-site')
    except:
        staging_site = None

    job = dax.DAX(fil)
    job.addArguments('--basename %s' \
                     % os.path.splitext(os.path.basename(name))[0])
    Workflow.set_job_properties(job, map_file, tc_file, staging_site=staging_site)
    workflow._adag.addJob(job)
    dep = dax.Dependency(parent=node._dax_node, child=job)
    workflow._adag.addDependency(dep)
    logging.info('Leaving minifollowups module')


def setup_injection_minifollowups(workflow, injection_file, inj_xml_file,
                                  single_triggers, tmpltbank_file,
                                  insp_segs, insp_data_name, insp_anal_name,
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
    insp_segs: SegFile
       The segment file containing the data read by each inspiral job.
    insp_data_name: str
        The name of the segmentlist storing data read.
    insp_anal_name: str
        The name of the segmentlist storing data analyzed.
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
    config_file.PFN(urljoin('file:', pathname2url(config_path)), site='local')

    exe = Executable(workflow.cp, 'injection_minifollowup', ifos=workflow.ifos, out_dir=dax_output)

    node = exe.create_node()
    node.add_input_opt('--config-files', config_file)
    node.add_input_opt('--bank-file', tmpltbank_file)
    node.add_input_opt('--injection-file', injection_file)
    node.add_input_opt('--injection-xml-file', inj_xml_file)
    node.add_multiifo_input_list_opt('--single-detector-triggers', single_triggers)
    node.add_input_opt('--inspiral-segments', insp_segs)
    node.add_opt('--inspiral-data-read-name', insp_data_name)
    node.add_opt('--inspiral-data-analyzed-name', insp_anal_name)
    node.new_output_file_opt(workflow.analysis_time, '.dax', '--output-file', tags=tags)
    node.new_output_file_opt(workflow.analysis_time, '.dax.map', '--output-map', tags=tags)
    node.new_output_file_opt(workflow.analysis_time, '.tc.txt', '--transformation-catalog', tags=tags)

    name = node.output_files[0].name
    map_file = node.output_files[1]
    tc_file = node.output_files[2]

    node.add_opt('--workflow-name', name)
    node.add_opt('--output-dir', out_dir)

    workflow += node

    # execute this in a sub-workflow
    fil = node.output_files[0]

    # determine if a staging site has been specified
    try:
        staging_site = workflow.cp.get('workflow-injection_minifollowups',
                                       'staging-site')
    except:
        staging_site = None

    job = dax.DAX(fil)
    job.addArguments('--basename %s' % os.path.splitext(os.path.basename(name))[0])
    Workflow.set_job_properties(job, map_file, tc_file, staging_site=staging_site)
    workflow._adag.addJob(job)
    dep = dax.Dependency(parent=node._dax_node, child=job)
    workflow._adag.addDependency(dep)
    logging.info('Leaving injection minifollowups module')


class SingleTemplateExecutable(PlotExecutable):
    """Class to be used for to create workflow.Executable instances for the
    pycbc_single_template executable. Basically inherits directly from
    PlotExecutable but adds the file_input_options.
    """
    file_input_options = ['--gating-file']


class SingleTimeFreqExecutable(PlotExecutable):
    """Class to be used for to create workflow.Executable instances for the
    pycbc_plot_singles_timefreq executable. Basically inherits directly from
    PlotExecutable but adds the file_input_options.
    """
    file_input_options = ['--gating-file']

class PlotQScanExecutable(PlotExecutable):
    """Class to be used for to create workflow.Executable instances for the
    pycbc_plot_qscan executable. Basically inherits directly from
    PlotExecutable but adds the file_input_options.
    """
    file_input_options = ['--gating-file']


def make_single_template_plots(workflow, segs, data_read_name, analyzed_name,
                                  params, out_dir, inj_file=None, exclude=None,
                                  require=None, tags=None, params_str=None,
                                  use_exact_inj_params=False):
    """Function for creating jobs to run the pycbc_single_template code and
    to run the associated plotting code pycbc_single_template_plots and add
    these jobs to the workflow.

    Parameters
    -----------
    workflow : workflow.Workflow instance
        The pycbc.workflow.Workflow instance to add these jobs to.
    segs : workflow.File instance
        The pycbc.workflow.File instance that points to the XML file containing
        the segment lists of data read in and data analyzed.
    data_read_name : str
        The name of the segmentlist containing the data read in by each
        inspiral job in the segs file.
    analyzed_name : str
        The name of the segmentlist containing the data analyzed by each
        inspiral job in the segs file.
    params : dictionary
        A dictionary containing the parameters of the template to be used.
        params[ifo+'end_time'] is required for all ifos in workflow.ifos.
        If use_exact_inj_params is False then also need to supply values for
        [mass1, mass2, spin1z, spin2x]. For precessing templates one also
        needs to supply [spin1y, spin1x, spin2x, spin2y, inclination]
        additionally for precession one must supply u_vals or
        u_vals_+ifo for all ifos. u_vals is the ratio between h_+ and h_x to
        use when constructing h(t). h(t) = (h_+ * u_vals) + h_x.
    out_dir : str
        Directory in which to store the output files.
    inj_file : workflow.File (optional, default=None)
        If given send this injection file to the job so that injections are
        made into the data.
    exclude : list (optional, default=None)
        If given, then when considering which subsections in the ini file to
        parse for options to add to single_template_plot, only use subsections
        that *do not* match strings in this list.
    require : list (optional, default=None)
        If given, then when considering which subsections in the ini file to
        parse for options to add to single_template_plot, only use subsections
        matching strings in this list.
    tags : list (optional, default=None)
        Add this list of tags to all jobs.
    params_str : str (optional, default=None)
        If given add this string to plot title and caption to describe the
        template that was used.
    use_exact_inj_params : boolean (optional, default=False)
        If True do not use masses and spins listed in the params dictionary
        but instead use the injection closest to the filter time as a template.

    Returns
    --------
    output_files : workflow.FileList
        The list of workflow.Files created in this function.
    """
    tags = [] if tags is None else tags
    makedir(out_dir)
    name = 'single_template_plot'
    secs = requirestr(workflow.cp.get_subsections(name), require)
    secs = excludestr(secs, exclude)
    files = FileList([])
    for tag in secs:
        for ifo in workflow.ifos:
            if params['%s_end_time' % ifo] == -1.0:
                continue
            # Reanalyze the time around the trigger in each detector
            node = SingleTemplateExecutable(workflow.cp, 'single_template',
                                            ifos=[ifo], out_dir=out_dir,
                                            tags=[tag] + tags).create_node()
            if use_exact_inj_params:
                node.add_opt('--use-params-of-closest-injection')
            else:
                node.add_opt('--mass1', "%.6f" % params['mass1'])
                node.add_opt('--mass2', "%.6f" % params['mass2'])
                node.add_opt('--spin1z',"%.6f" % params['spin1z'])
                node.add_opt('--spin2z',"%.6f" % params['spin2z'])
                node.add_opt('--template-start-frequency',
                             "%.6f" % params['f_lower'])
                # Is this precessing?
                if 'u_vals' in params or 'u_vals_%s' % ifo in params:
                    node.add_opt('--spin1x',"%.6f" % params['spin1x'])
                    node.add_opt('--spin1y',"%.6f" % params['spin1y'])
                    node.add_opt('--spin2x',"%.6f" % params['spin2x'])
                    node.add_opt('--spin2y',"%.6f" % params['spin2y'])
                    node.add_opt('--inclination',"%.6f" % params['inclination'])
                    try:
                        node.add_opt('--u-val',"%.6f" % params['u_vals'])
                    except:
                        node.add_opt('--u-val',
                                     "%.6f" % params['u_vals_%s' % ifo])

            # str(numpy.float64) restricts to 2d.p. BE CAREFUL WITH THIS!!!
            str_trig_time = '%.6f' %(params[ifo + '_end_time'])
            node.add_opt('--trigger-time', str_trig_time)
            node.add_input_opt('--inspiral-segments', segs)
            if inj_file is not None:
                node.add_input_opt('--injection-file', inj_file)
            node.add_opt('--data-read-name', data_read_name)
            node.add_opt('--data-analyzed-name', analyzed_name)
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
            if use_exact_inj_params:
                caption += ". The injection itself was used as the template.'"
            else:
                caption += ". The template used has the following parameters: "
                caption += "mass1=%s, mass2=%s, spin1z=%s, spin2z=%s'"\
                       %(params['mass1'], params['mass2'], params['spin1z'],
                         params['spin2z'])
            node.add_opt('--plot-caption', caption)
            workflow += node
            files += node.output_files
    return files

def make_plot_waveform_plot(workflow, params, out_dir, ifos, exclude=None,
                            require=None, tags=None):
    """ Add plot_waveform jobs to the workflow.
    """
    tags = [] if tags is None else tags
    makedir(out_dir)
    name = 'single_template_plot'
    secs = requirestr(workflow.cp.get_subsections(name), require)
    secs = excludestr(secs, exclude)
    files = FileList([])
    for tag in secs:
        node = PlotExecutable(workflow.cp, 'plot_waveform', ifos=ifos,
                              out_dir=out_dir, tags=[tag] + tags).create_node()
        node.add_opt('--mass1', "%.6f" % params['mass1'])
        node.add_opt('--mass2', "%.6f" % params['mass2'])
        node.add_opt('--spin1z',"%.6f" % params['spin1z'])
        node.add_opt('--spin2z',"%.6f" % params['spin2z'])
        if 'u_vals' in params:
            # Precessing options
            node.add_opt('--spin1x',"%.6f" % params['spin1x'])
            node.add_opt('--spin2x',"%.6f" % params['spin2x'])
            node.add_opt('--spin1y',"%.6f" % params['spin1y'])
            node.add_opt('--spin2y',"%.6f" % params['spin2y'])
            node.add_opt('--inclination',"%.6f" % params['inclination'])
            node.add_opt('--u-val', "%.6f" % params['u_vals'])
        node.new_output_file_opt(workflow.analysis_time, '.png',
                                     '--output-file')
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

def make_coinc_info(workflow, singles, bank, coinc, out_dir,
                    n_loudest=None, trig_id=None, file_substring=None,
                    tags=None):
    tags = [] if tags is None else tags
    makedir(out_dir)
    name = 'page_coincinfo'
    files = FileList([])
    node = PlotExecutable(workflow.cp, name, ifos=workflow.ifos,
                              out_dir=out_dir, tags=tags).create_node()
    node.add_input_list_opt('--single-trigger-files', singles)
    node.add_input_opt('--statmap-file', coinc)
    node.add_input_opt('--bank-file', bank)
    if n_loudest is not None:
        node.add_opt('--n-loudest', str(n_loudest))
    if trig_id is not None:
        node.add_opt('--trigger-id', str(trig_id))
    if file_substring is not None:
        node.add_opt('--statmap-file-subspace-name', file_substring)
    node.new_output_file_opt(workflow.analysis_time, '.html', '--output-file')
    workflow += node
    files += node.output_files
    return files

def make_sngl_ifo(workflow, sngl_file, bank_file, trigger_id, out_dir, ifo,
                  tags=None, rank=None):
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
    node.add_opt('--trigger-id', str(trigger_id))
    if rank is not None:
        node.add_opt('--n-loudest', str(rank))
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

def make_qscan_plot(workflow, ifo, trig_time, out_dir, injection_file=None,
                    data_segments=None, time_window=100, tags=None):
    """ Generate a make_qscan node and add it to workflow.

    This function generates a single node of the singles_timefreq executable
    and adds it to the current workflow. Parent/child relationships are set by
    the input/output files automatically.

    Parameters
    -----------
    workflow: pycbc.workflow.core.Workflow
        The workflow class that stores the jobs that will be run.
    ifo: str
        Which interferometer are we using?
    trig_time: int
        The time of the trigger being followed up.
    out_dir: str
        Location of directory to output to
    injection_file: pycbc.workflow.File (optional, default=None)
        If given, add the injections in the file to strain before making the
        plot.
    data_segments: ligo.segments.segmentlist (optional, default=None)
        The list of segments for which data exists and can be read in. If given
        the start/end times given to singles_timefreq will be adjusted if
        [trig_time - time_window, trig_time + time_window] does not completely
        lie within a valid data segment. A ValueError will be raised if the
        trig_time is not within a valid segment, or if it is not possible to
        find 2*time_window (plus the padding) of continuous data around the
        trigger. This **must** be coalesced.
    time_window: int (optional, default=None)
        The amount of data (not including padding) that will be read in by the
        singles_timefreq job. The default value of 100s should be fine for most
        cases.
    tags: list (optional, default=None)
        List of tags to add to the created nodes, which determine file naming.
    """
    tags = [] if tags is None else tags
    makedir(out_dir)
    name = 'plot_qscan'

    curr_exe = PlotQScanExecutable(workflow.cp, name, ifos=[ifo],
                          out_dir=out_dir, tags=tags)
    node = curr_exe.create_node()

    # Determine start/end times, using data segments if needed.
    # Begin by choosing "optimal" times
    start = trig_time - time_window
    end = trig_time + time_window
    # Then if data_segments is available, check against that, and move if
    # needed
    if data_segments is not None:
        # Assumes coalesced, so trig_time can only be within one segment
        for seg in data_segments:
            if trig_time in seg:
                data_seg = seg
                break
            elif trig_time == -1.0:
                node.add_opt('--gps-start-time', int(trig_time))
                node.add_opt('--gps-end-time', int(trig_time))
                node.add_opt('--center-time', trig_time)
                caption_string = "'No trigger in %s'" % ifo
                node.add_opt('--plot-caption', caption_string)
                node.new_output_file_opt(workflow.analysis_time, '.png', '--output-file')
                workflow += node
                return node.output_files
        else:
            err_msg = "Trig time {} ".format(trig_time)
            err_msg += "does not seem to lie within any data segments. "
            err_msg += "This shouldn't be possible, please ask for help!"
            raise ValueError(err_msg)
        # Check for pad-data
        if curr_exe.has_opt('pad-data'):
            pad_data = int(curr_exe.get_opt('pad-data'))
        else:
            pad_data = 0
        # We only read data that's available. The code must handle the case
        # of not much data being available.
        if end > (data_seg[1] - pad_data):
            end = data_seg[1] - pad_data
        if start < (data_seg[0] + pad_data):
            start = data_seg[0] + pad_data

    node.add_opt('--gps-start-time', int(start))
    node.add_opt('--gps-end-time', int(end))
    node.add_opt('--center-time', trig_time)

    if injection_file is not None:
        node.add_input_opt('--injection-file', injection_file)

    node.new_output_file_opt(workflow.analysis_time, '.png', '--output-file')
    workflow += node
    return node.output_files

def make_singles_timefreq(workflow, single, bank_file, trig_time, out_dir,
                          veto_file=None, time_window=10, data_segments=None,
                          tags=None):
    """ Generate a singles_timefreq node and add it to workflow.

    This function generates a single node of the singles_timefreq executable
    and adds it to the current workflow. Parent/child relationships are set by
    the input/output files automatically.

    Parameters
    -----------
    workflow: pycbc.workflow.core.Workflow
        The workflow class that stores the jobs that will be run.
    single: pycbc.workflow.core.File instance
        The File object storing the single-detector triggers to followup.
    bank_file: pycbc.workflow.core.File instance
        The File object storing the template bank.
    trig_time: int
        The time of the trigger being followed up.
    out_dir: str
        Location of directory to output to
    veto_file: pycbc.workflow.core.File (optional, default=None)
        If given use this file to veto triggers to determine the loudest event.
        FIXME: Veto files *should* be provided a definer argument and not just
        assume that all segments should be read.
    time_window: int (optional, default=None)
        The amount of data (not including padding) that will be read in by the
        singles_timefreq job. The default value of 10s should be fine for most
        cases.
    data_segments: ligo.segments.segmentlist (optional, default=None)
        The list of segments for which data exists and can be read in. If given
        the start/end times given to singles_timefreq will be adjusted if
        [trig_time - time_window, trig_time + time_window] does not completely
        lie within a valid data segment. A ValueError will be raised if the
        trig_time is not within a valid segment, or if it is not possible to
        find 2*time_window (plus the padding) of continuous data around the
        trigger. This **must** be coalesced.
    tags: list (optional, default=None)
        List of tags to add to the created nodes, which determine file naming.
    """
    tags = [] if tags is None else tags
    makedir(out_dir)
    name = 'plot_singles_timefreq'

    curr_exe = SingleTimeFreqExecutable(workflow.cp, name, ifos=[single.ifo],
                          out_dir=out_dir, tags=tags)
    node = curr_exe.create_node()
    node.add_input_opt('--trig-file', single)
    node.add_input_opt('--bank-file', bank_file)

    # Determine start/end times, using data segments if needed.
    # Begin by choosing "optimal" times
    start = trig_time - time_window
    end = trig_time + time_window
    # Then if data_segments is available, check against that, and move if
    # needed
    if data_segments is not None:
        # Assumes coalesced, so trig_time can only be within one segment
        for seg in data_segments:
            if trig_time in seg:
                data_seg = seg
                break
            elif trig_time == -1.0:
                node.add_opt('--gps-start-time', int(trig_time))
                node.add_opt('--gps-end-time', int(trig_time))
                node.add_opt('--center-time', trig_time)

                if veto_file:
                    node.add_input_opt('--veto-file', veto_file)

                node.add_opt('--detector', single.ifo)
                node.new_output_file_opt(workflow.analysis_time, '.png', '--output-file')
                workflow += node
                return node.output_files
        else:
            err_msg = "Trig time {} ".format(trig_time)
            err_msg += "does not seem to lie within any data segments. "
            err_msg += "This shouldn't be possible, please ask for help!"
            raise ValueError(err_msg)
        # Check for pad-data
        if curr_exe.has_opt('pad-data'):
            pad_data = int(curr_exe.get_opt('pad-data'))
        else:
            pad_data = 0
        if abs(data_seg) < (2 * time_window + 2 * pad_data):
            tl = 2 * time_window + 2 * pad_data
            err_msg = "I was asked to use {} seconds of data ".format(tl)
            err_msg += "to run a plot_singles_timefreq job. However, I have "
            err_msg += "only {} seconds available.".format(abs(data_seg))
            raise ValueError(err_msg)
        if data_seg[0] > (start - pad_data):
            start = data_seg[0] + pad_data
            end = start + 2 * time_window
        if data_seg[1] < (end + pad_data):
            end = data_seg[1] - pad_data
            start = end - 2 * time_window
        # Sanity check, shouldn't get here!
        if data_seg[0] > (start - pad_data):
            err_msg = "I shouldn't be here! Go ask Ian what he broke."
            raise ValueError(err_msg)

    node.add_opt('--gps-start-time', int(start))
    node.add_opt('--gps-end-time', int(end))
    node.add_opt('--center-time', trig_time)

    if veto_file:
        node.add_input_opt('--veto-file', veto_file)

    node.add_opt('--detector', single.ifo)
    node.new_output_file_opt(workflow.analysis_time, '.png', '--output-file')
    workflow += node
    return node.output_files

def create_noop_node():
    """
    Creates a noop node that can be added to a DAX doing nothing. The reason
    for using this is if a minifollowups dax contains no triggers currently
    the dax will contain no jobs and be invalid. By adding a noop node we
    ensure that such daxes will actually run if one adds one such noop node.
    Adding such a noop node into a workflow *more than once* will cause a
    failure.
    """
    exe = wdax.Executable('NOOP')
    pfn = distutils.spawn.find_executable('true')
    exe.add_pfn(pfn)
    node = wdax.Node(exe)
    return node
