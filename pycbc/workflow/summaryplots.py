# Copyright (C) 2014 Christopher M. Biwer
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
import itertools
import logging
import os
import urllib

from pycbc.workflow.core import Executable, Node, File, FileList, make_external_call

import Pegasus.DAX3 as dax

def setup_plotinspiral(workflow, input_files, cache_filename,
                       inspiral_cachepattern, output_dir, tags=[], **kwargs):
    """
    This function sets up the nodes that will generate summary plots given a
    list of inspiral files.

    Parameters
    -----------
    Workflow : ahope.Workflow
        The ahope workflow instance that the coincidence jobs will be added to.
    input_files : ahope.FileList
        An FileList of files that are used as input at this stage.
    cache_filename : str
        Filename of the ihope cache.
    inspiral_cachepattern : str
        The pattern that will be used to find inspiral filenames in the cache.
    output_dir : path
        The directory in which output files will be stored.
    tags : list of strings (optional, default = [])
        A list of the tagging strings that will be used for all jobs created
        by this call to the workflow. An example might be ['full_data'].
        This will be used in output names and directories.

    Returns
    --------
    plot_files : ahope.FileList
        A list of the output files from this stage.
    """

    plot_files = FileList([])

    # create executable
    plotinspiral_job = Executable(workflow.cp, 'plotinspiral', 'vanilla',
                                  workflow.ifos, output_dir, tags)

    for tag in tags:
        for ifo in plotinspiral_job.ifo_list:
            # create node
            node = Node(plotinspiral_job)
            node.add_opt('--gps-start-time', workflow.analysis_time[0])
            node.add_opt('--gps-end-time', workflow.analysis_time[1])
            node.add_opt('--cache-file', cache_filename)
            node.add_opt('--ifo-times', ifo)
            node.add_opt('--ifo-tag', 'FIRST_'+ifo)
            node.add_opt('--user-tag', tag.upper()+'_SUMMARY_PLOTS')
            node.add_opt('--output-path', output_dir)
            node.add_opt('--trig-pattern', inspiral_cachepattern)
            node.add_opt('--enable-output')

            # add node to workflow
            workflow.add_node(node)

            # make all input_files parents
            #for f in input_files:
            #    dep = dax.Dependency(parent=f.node._dax_node, child=node._dax_node)
            #    workflow._adag.addDependency(dep)

    return plot_files


def setup_plotnumtemplates(workflow, input_files, cache_filename,
                           tmpltbank_cachepattern, output_dir, tags=[],
                           **kwargs):
    """
    This function sets up the nodes that will generate a plot of the number
    of templates against time.

    Parameters
    -----------
    Workflow : ahope.Workflow
        The ahope workflow instance that the coincidence jobs will be added to.
    input_files : ahope.FileList
        An FileList of files that are used as input at this stage.
    cache_filename : str
        Filename of the ihope cache.
    tmpltbank_cachepattern : str
        The pattern that will be used to find template_bank filenames in the cache.
    output_dir : path
        The directory in which output files will be stored.
    tags : list of strings (optional, default = [])
        A list of the tagging strings that will be used for all jobs created
        by this call to the workflow. An example might be ['full_data'].
        This will be used in output names and directories.

    Returns
    --------
    plot_files : ahope.FileList
        A list of the output files from this stage.
    """

    plot_files = FileList([])

    # create executable
    plotnumtemplates_job = Executable(workflow.cp, 'plotnumtemplates',
                                  'vanilla', workflow.ifos, output_dir, tags)

    for tag in tags:
        # create node
        node = Node(plotnumtemplates_job)
        node.add_opt('--gps-start-time', workflow.analysis_time[0])
        node.add_opt('--gps-end-time', workflow.analysis_time[1])
        node.add_opt('--cache-file', cache_filename)
        node.add_opt('--ifo-times', node.executable.ifo_string)
        node.add_opt('--user-tag', tag.upper()+'_SUMMARY_PLOTS')
        node.add_opt('--output-path', output_dir)
        node.add_opt('--bank-pattern', tmpltbank_cachepattern)
        node.add_opt('--enable-output')

        # add node to workflow
        workflow.add_node(node)

        # make all input_files parents
        #for f in input_files:
        #    dep = dax.Dependency(parent=f.node._dax_node, child=node._dax_node)
        #    workflow._adag.addDependency(dep)

    return plot_files


def setup_plotthinca(workflow, input_files, cache_filename, coinc_cachepattern,
                     slide_cachepattern, output_dir, tags=[], **kwargs):
    """
    This function sets up the nodes that will generate summary from a list of
    thinca files.

    Parameters
    -----------
    Workflow : ahope.Workflow
        The ahope workflow instance that the coincidence jobs will be added to.
    input_files : ahope.FileList
        An FileList of files that are used as input at this stage.
    cache_filename : str
        Filename of the ihope cache.
    coinc_cachepattern : str
        The pattern that will be used to find zero-lag coincidence filenames in the cache.
    slide_cachepattern : str
        The pattern that will be used to find time slide filenames in the cache.
    output_dir : path
        The directory in which output files will be stored.
    tags : list of strings (optional, default = [])
        A list of the tagging strings that will be used for all jobs created
        by this call to the workflow. An example might be ['full_data'].
        This will be used in output names and directories.

    Returns
    --------
    plot_files : ahope.FileList
        A list of the output files from this stage.
    """

    plot_files = FileList([])

    # create executable
    plotthinca_job = Executable(workflow.cp, 'plotthinca', 'vanilla',
                                workflow.ifos, output_dir, tags)

    # get all ifo combinations of at least 2 coincident ifos
    ifo_combos = []
    for n in xrange(len(plotthinca_job.ifo_list)+1):
        for ifo_list in itertools.combinations(plotthinca_job.ifo_list, n+2):
            ifo_combos.append(ifo_list)

    for tag in tags:
        for ifo_list in ifo_combos:
            ifo_string = ''.join(ifo_list)

            # create node
            node = Node(plotthinca_job)
            node.add_opt('--gps-start-time', workflow.analysis_time[0])
            node.add_opt('--gps-end-time', workflow.analysis_time[1])
            node.add_opt('--cache-file', cache_filename)
            node.add_opt('--ifo-times', ifo_string)
            node.add_opt('--ifo-tag', 'SECOND_'+ifo_string)
            for ifo in ifo_list:
                node.add_opt('--%s-triggers'%ifo.lower(), '')
            node.add_opt('--user-tag', tag.upper()+'_SUMMARY_PLOTS')
            node.add_opt('--output-path', output_dir)
            node.add_opt('--coinc-pattern', coinc_cachepattern)
            node.add_opt('--slide-pattern', slide_cachepattern)
            node.add_opt('--enable-output')

            # add node to workflow
            workflow.add_node(node)

            # make all input_files parents
            #for f in input_files:
            #    dep = dax.Dependency(parent=f.node._dax_node, child=node._dax_node)
            #    workflow._adag.addDependency(dep)

    return plot_files

def setup_plotinspiralrange(workflow, input_files, cache_filename,
                            tmpltbank_cachepattern, inspiral_cachepattern,
                            output_dir, tags=[], **kwargs):
    """
    This function sets up the nodes that will generate the inspiral horizon distance
    plots.

    Parameters
    -----------
    Workflow : ahope.Workflow
        The ahope workflow instance that the coincidence jobs will be added to.
    input_files : ahope.FileList
        An FileList of files that are used as input at this stage.
    cache_filename : str
        Filename of the ihope cache.
    tmpltbank_cachepattern : str
        The pattern that will be used to find template_bank filenames in the cache.
    inspiral_cachepattern : str
        The pattern that will be used to find inspiral filenames in the cache.
    output_dir : path
        The directory in which output files will be stored.
    tags : list of strings (optional, default = [])
        A list of the tagging strings that will be used for all jobs created
        by this call to the workflow. An example might be ['full_data'].
        This will be used in output names and directories.

    Returns
    --------
    plot_files : ahope.FileList
        A list of the output files from this stage.
    """

    plot_files = FileList([])

    # create executable
    plotinspiralrange_job = Executable(workflow.cp, 'plotinspiralrange',
                                   'vanilla', workflow.ifos, output_dir, tags)

    for tag in tags:
        # create node
        node = Node(plotinspiralrange_job)
        node.add_opt('--gps-start-time', workflow.analysis_time[0])
        node.add_opt('--gps-end-time', workflow.analysis_time[1])
        node.add_opt('--cache-file', cache_filename)
        node.add_opt('--ifo-times', node.executable.ifo_string)
        node.add_opt('--user-tag', tag.upper()+'_SUMMARY_PLOTS')
        node.add_opt('--output-path', output_dir)
        node.add_opt('--bank-pattern', tmpltbank_cachepattern)
        node.add_opt('--trig-pattern', inspiral_cachepattern)
        node.add_opt('--enable-output')

        # add node to workflow
        workflow.add_node(node)

        # make all input_files parents
        #for f in input_files:
        #    dep = dax.Dependency(parent=f.node._dax_node, child=node._dax_node)
        #    workflow._adag.addDependency(dep)

    return plot_files

def setup_summary_plots(workflow, input_files, cache_filename,
                        tmpltbank_cachepattern, inspiral_cachepattern,
                        #coinc_cachepattern, slide_cachepattern,
                        output_dir, tags=[], **kwargs):
    """
    This function sets up the summary plots jobs.

    Parameters
    -----------
    Workflow : ahope.Workflow
        The ahope workflow instance that the coincidence jobs will be added to.
    input_files : ahope.FileList
        An FileList of files that are used as input at this stage.
    cache_filename : str
        Filename of the ihope cache.
    tmpltbank_cachepattern : str
        The pattern that will be used to find template_bank filenames in the cache.
    inspiral_cachepattern : str
        The pattern that will be used to find inspiral filenames in the cache.
    coinc_cachepattern : str (currently not implemented)
        The pattern that will be used to find zero-lag coincidence filenames in the cache.
    slide_cachepattern : str (currently not implemented)
        The pattern that will be used to find time slide filenames in the cache.
    output_dir : path
        The directory in which output files will be stored.
    tags : list of strings (optional, default = [])
        A list of the tagging strings that will be used for all jobs created
        by this call to the workflow. An example might be ['full_data'].
        This will be used in output names and directories.

    Returns
    --------
    plot_files : ahope.FileList
        A list of the output files from this stage.
    """

    plot_files = FileList([])

    # make summary plot dir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # setup summary plot jobs
    plotinspiral_files = setup_plotinspiral(workflow, input_files,
                       cache_filename, inspiral_cachepattern, output_dir, tags)
    plotinspiralrange_files = setup_plotinspiralrange(workflow, input_files,
                                       cache_filename, tmpltbank_cachepattern,
                                       inspiral_cachepattern, output_dir, tags)
    plotnumtemplates_files = setup_plotnumtemplates(workflow, input_files,
                      cache_filename, tmpltbank_cachepattern, output_dir, tags)
    #plotthinca_files = setup_plotthinca(workflow, input_files, cache_filename,
    #                                        coinc_cachepattern, slide_cachepattern,
    #                                        output_dir, tags)

    # concatenate plot files
    plot_files += plotinspiral_files
    plot_files += plotinspiralrange_files
    plot_files += plotnumtemplates_files
    #plot_files += plotthinca_files

    return plot_files

def setup_hardware_injection_page(workflow, input_files, cache_filename,
                         inspiral_cachepattern, output_dir, tags=[], **kwargs):
    """
    This function sets up the nodes that will create the hardware injection page.

    Parameters
    -----------
    Workflow : ahope.Workflow
        The ahope workflow instance that the coincidence jobs will be added to.
    input_files : ahope.FileList
        An FileList of files that are used as input at this stage.
    cache_filename : str
        Filename of the ihope cache.
    inspiral_cachepattern : str
        The pattern that will be used to find inspiral filenames in the cache.
    output_dir : path
        The directory in which output files will be stored.
    tags : list of strings (optional, default = [])
        A list of the tagging strings that will be used for all jobs created
        by this call to the workflow. An example might be ['full_data'].
        This will be used to search the cache.

    Returns
    --------
    plot_files : ahope.FileList
        A list of the output files from this stage.
    """

    logging.info("Entering hardware injection page setup.")

    out_files = FileList([])

    # check if hardware injection section exists
    # if not then do not do add hardware injection job to the workflow
    if not workflow.cp.has_section('workflow-hardware-injections'):
      msg  = "There is no workflow-hardware-injections section. "
      msg += "The hardware injection page will not be added to the workflow."
      logging.info(msg)
      logging.info("Leaving hardware injection page setup.")
      return out_files

    # make the output dir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # create executable
    hwinjpage_job = Executable(workflow.cp, 'hardware_injection_page',
                                    'vanilla', workflow.ifos, output_dir, tags)

    # retrieve hardware injection file
    hwinjDefUrl = workflow.cp.get_opt_tags('workflow-hardware-injections',
                                                     'hwinj-definer-url', tags)
    hwinjDefBaseName = os.path.basename(hwinjDefUrl)
    hwinjDefNewPath = os.path.join(output_dir, hwinjDefBaseName)
    urllib.urlretrieve (hwinjDefUrl, hwinjDefNewPath)

    # update hwinj definer file location
    workflow.cp.set("workflow-hardware-injections", "hwinj-definer-file",
                                                               hwinjDefNewPath)

    # query for the hardware injection segments
    get_hardware_injection_segment_files(workflow, output_dir, hwinjDefNewPath)

    # create node
    node = Node(hwinjpage_job)
    node.add_opt('--gps-start-time', workflow.analysis_time[0])
    node.add_opt('--gps-end-time', workflow.analysis_time[1])
    node.add_opt('--source-xml', hwinjDefNewPath)
    node.add_opt('--segment-dir', output_dir)
    node.add_opt('--cache-file', cache_filename)
    node.add_opt('--cache-pattern', inspiral_cachepattern)
    node.add_opt('--analyze-injections', '')
    for ifo in workflow.ifos:
        node.add_opt('--%s-injections'%ifo.lower(), '')
    outfile = File(node.executable.ifo_string, 'HWINJ_SUMMARY',
                workflow.analysis_time, extension='html', directory=output_dir)
    node.add_opt('--outfile', outfile.storage_path)

    # add node to workflow
    workflow.add_node(node)

    # make all input_files parents
    #for f in input_files:
    #    dep = dax.Dependency(parent=f.node._dax_node, child=node._dax_node)
    #    workflow._adag.addDependency(dep)

    out_files += node.output_files

    logging.info("Leaving hardware injection page setup.")

    return out_files

def get_hardware_injection_segment_files(workflow, output_dir, hwinjDefPath,
                                         tag=None):
    """
    This function queries the segment database for the hardware
    injection segments and saves them to the output_dir.

    Parameters
    -----------
    workflow : ahope.Workflow
        The ahope workflow instance that the coincidence jobs will be added to.
    output_dir : path
        The directory in which output files will be stored.
    hwinjDefPath : path
        The path to the hardware injection definer file.
    tag : string, optional (default=None)
        Use this to specify a tag. This can be used if this module is being
        called more than once to give call specific configuration (by setting
        options in [workflow-datafind-${TAG}] rather than [workflow-datafind]).
        This is also used to tag the Files returned by the class to uniqueify
        the Files and uniqueify the actual filename.
    """

    # create log dir
    log_dir = os.path.join(output_dir, 'logs')
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    # get hardware injection segments
    # ligolw_cbc_hardware_inj_page expects seperate XML files for each IFO
    for ifo in workflow.ifos:
        hwinjSegName = workflow.cp.get_opt_tags('workflow-hardware-injections',
                                 'segments-%s-hwinj-name' % (ifo.lower()), [tag])

        output_filename = '-'.join(map(str, [ifo, 'HWINJ_SEGMENTS',
                                       workflow.analysis_time[0],
                                       abs(workflow.analysis_time)])) + '.xml'
        output_path = os.path.join(output_dir, output_filename)

        segFindCall = [ workflow.cp.get("executables","segment_query"),
            "--query-segments",
            "--gps-start-time", str(workflow.analysis_time[0]),
            "--gps-end-time", str(workflow.analysis_time[1]),
            "--include-segments", hwinjSegName,
            "--output-file", output_path,
            "--segment-url", workflow.cp.get("workflow-hardware-injections",
                                                 "segments-database-url")]

        make_external_call(segFindCall, out_dir=log_dir,
                                out_basename='%s-hwinj-call' %(ifo.lower()) )
