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
This module is responsible for setting up the segment generation stage of
workflows. For details about this module and its capabilities see here:
https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/ahope/segments.html
"""

import os, sys, shutil, stat, copy, itertools
import logging
from six.moves.urllib.request import pathname2url
from six.moves.urllib.parse import urljoin, urlunparse
import lal
from ligo import segments
from ligo.segments import utils as segmentsUtils
from glue.ligolw import table, lsctables, ligolw
from pycbc.workflow.core import Executable, FileList, Node, SegFile, make_analysis_dir, make_external_call, File
from pycbc.workflow.core import resolve_url
from pycbc.workflow.jobsetup import LigolwAddExecutable, LigoLWCombineSegsExecutable

def get_science_segments(workflow, out_dir, tags=None):
    """
    Get the analyzable segments after applying ini specified vetoes.

    Parameters
    -----------
    workflow : Workflow object
        Instance of the workflow object
    out_dir : path
        Location to store output files
    tags : list of strings
        Used to retrieve subsections of the ini file for
        configuration options.

    Returns
    --------
    sci_seg_file : workflow.core.SegFile instance
        The segment file combined from all ifos containing the science segments.
    sci_segs : Ifo keyed dict of ligo.segments.segmentlist instances
        The science segs for each ifo, keyed by ifo
    sci_seg_name : str
        The name with which science segs are stored in the output XML file.
    """
    if tags is None:
        tags = []
    logging.info('Starting generation of science segments')

    make_analysis_dir(out_dir)
    start_time = workflow.analysis_time[0]
    end_time = workflow.analysis_time[1]

    # NOTE: Should this be overrideable in the config file?
    sci_seg_name = "SCIENCE"
    sci_segs = {}
    sci_seg_dict = segments.segmentlistdict()
    sci_seg_summ_dict = segments.segmentlistdict()

    for ifo in workflow.ifos:
        curr_sci_segs, curr_sci_xml, curr_seg_name = get_sci_segs_for_ifo(ifo,
                              workflow.cp, start_time, end_time, out_dir, tags)
        sci_seg_dict[ifo + ':' + sci_seg_name] = curr_sci_segs
        sci_segs[ifo] = curr_sci_segs
        sci_seg_summ_dict[ifo + ':' + sci_seg_name] = \
                          curr_sci_xml.seg_summ_dict[ifo + ':' + curr_seg_name]
    sci_seg_file = SegFile.from_segment_list_dict(sci_seg_name,
                                          sci_seg_dict, extension='xml',
                                          valid_segment=workflow.analysis_time,
                                          seg_summ_dict=sci_seg_summ_dict,
                                          directory=out_dir, tags=tags)
    logging.info('Done generating science segments')
    return sci_seg_file, sci_segs, sci_seg_name

def get_files_for_vetoes(workflow, out_dir,
                         runtime_names=None, in_workflow_names=None, tags=None):
    """
    Get the various sets of veto segments that will be used in this analysis.

    Parameters
    -----------
    workflow : Workflow object
        Instance of the workflow object
    out_dir : path
        Location to store output files
    runtime_names : list
        Veto category groups with these names in the [workflow-segment] section
        of the ini file will be generated now.
    in_workflow_names : list
        Veto category groups with these names in the [workflow-segment] section
        of the ini file will be generated in the workflow. If a veto category
        appears here and in runtime_names, it will be generated now.
    tags : list of strings
        Used to retrieve subsections of the ini file for
        configuration options.

    Returns
    --------
    veto_seg_files : FileList
        List of veto segment files generated
    """
    if tags is None:
        tags = []
    if runtime_names is None:
        runtime_names = []
    if in_workflow_names is None:
        in_workflow_names = []
    logging.info('Starting generating veto files for analysis')
    make_analysis_dir(out_dir)
    start_time = workflow.analysis_time[0]
    end_time = workflow.analysis_time[1]
    save_veto_definer(workflow.cp, out_dir, tags)

    now_cat_sets = []
    for name in runtime_names:
        cat_sets = parse_cat_ini_opt(workflow.cp.get_opt_tags(
                                              'workflow-segments', name, tags))
        now_cat_sets.extend(cat_sets)

    now_cats = set()
    for cset in now_cat_sets:
        now_cats = now_cats.union(cset)

    later_cat_sets = []
    for name in in_workflow_names:
        cat_sets = parse_cat_ini_opt(workflow.cp.get_opt_tags(
                                              'workflow-segments', name, tags))
        later_cat_sets.extend(cat_sets)

    later_cats = set()
    for cset in later_cat_sets:
        later_cats = later_cats.union(cset)
        # Avoid duplication
        later_cats = later_cats - now_cats

    veto_gen_job = create_segs_from_cats_job(workflow.cp, out_dir,
                                             workflow.ifo_string, tags=tags)

    cat_files = FileList()
    for ifo in workflow.ifos:
        for category in now_cats:
            cat_files.append(get_veto_segs(workflow, ifo,
                                        cat_to_veto_def_cat(category),
                                        start_time, end_time, out_dir,
                                        veto_gen_job, execute_now=True,
                                        tags=tags))

        for category in later_cats:
            cat_files.append(get_veto_segs(workflow, ifo,
                                        cat_to_veto_def_cat(category),
                                        start_time, end_time, out_dir,
                                        veto_gen_job, tags=tags,
                                        execute_now=False))

    logging.info('Done generating veto segments')
    return cat_files

def get_analyzable_segments(workflow, sci_segs, cat_files, out_dir, tags=None):
    """
    Get the analyzable segments after applying ini specified vetoes and any
    other restrictions on the science segs, e.g. a minimum segment length, or
    demanding that only coincident segments are analysed.

    Parameters
    -----------
    workflow : Workflow object
        Instance of the workflow object
    sci_segs : Ifo-keyed dictionary of glue.segmentlists
        The science segments for each ifo to which the vetoes, or any other
        restriction, will be applied.
    cat_files : FileList of SegFiles
        The category veto files generated by get_veto_segs
    out_dir : path
        Location to store output files
    tags : list of strings
        Used to retrieve subsections of the ini file for
        configuration options.

    Returns
    --------
    sci_ok_seg_file : workflow.core.SegFile instance
        The segment file combined from all ifos containing the analyzable
        science segments.
    sci_ok_segs : Ifo keyed dict of ligo.segments.segmentlist instances
        The analyzable science segs for each ifo, keyed by ifo
    sci_ok_seg_name : str
        The name with which analyzable science segs are stored in the output
        XML file.
    """
    if tags is None:
        tags = []
    logging.info('Starting reducing to analysable science segments')

    make_analysis_dir(out_dir)
    # NOTE: Should this be overrideable in the config file?
    sci_ok_seg_name = "SCIENCE_OK"
    sci_ok_seg_dict = segments.segmentlistdict()
    sci_ok_segs = {}

    cat_sets = parse_cat_ini_opt(workflow.cp.get_opt_tags('workflow-segments',
                                                'segments-science-veto', tags))
    if len(cat_sets) > 1:
        raise ValueError('Provide only 1 category group to determine'
                         ' analyzable segments')
    cat_set = cat_sets[0]

    for ifo in workflow.ifos:
        curr_segs = copy.copy(sci_segs[ifo])
        files = cat_files.find_output_with_ifo(ifo)
        for category in cat_set:
            veto_def_cat = cat_to_veto_def_cat(category)
            file_list = files.find_output_with_tag('VETO_CAT%d' %(veto_def_cat))
            if len(file_list) > 1:
                err_msg = "Found more than one veto file for %s " %(ifo,)
                err_msg += "and category %s." %(category,)
                raise ValueError(err_msg)
            if len(file_list) == 0:
                err_msg = "Found no veto files for %s " %(ifo,)
                err_msg += "and category %s." %(category,)
                raise ValueError(err_msg)
            curr_veto_file = file_list[0]
            cat_segs = curr_veto_file.return_union_seglist()
            curr_segs -= cat_segs
            curr_segs.coalesce()
        sci_ok_seg_dict[ifo + ':' + sci_ok_seg_name] = curr_segs

    sci_ok_seg_file = SegFile.from_segment_list_dict(sci_ok_seg_name,
                                          sci_ok_seg_dict, extension='xml',
                                          valid_segment=workflow.analysis_time,
                                          directory=out_dir, tags=tags)


    if workflow.cp.has_option_tags("workflow-segments",
                          "segments-minimum-segment-length", tags):
        min_seg_length = int( workflow.cp.get_opt_tags("workflow-segments",
                              "segments-minimum-segment-length", tags) )
        sci_ok_seg_file.remove_short_sci_segs(min_seg_length)

    # FIXME: Another test we can do is limit to coinc time +/- some window
    #        this should *not* be set through segments-method, but currently
    #        is not implemented
    #segments_method = workflow.cp.get_opt_tags("workflow-segments",
    #                                  "segments-method", tags)
    #if segments_method == 'ALL_SINGLE_IFO_TIME':
    #    pass
    #elif segments_method == 'COINC_TIME':
    #    cum_segs = None
    #    for ifo in sci_segs:
    #        if cum_segs is not None:
    #            cum_segs = (cum_segs & sci_segs[ifo]).coalesce()
    #        else:
    #            cum_segs = sci_segs[ifo]
    #
    #    for ifo in sci_segs:
    #        sci_segs[ifo] = cum_segs
    #else:
    #    raise ValueError("Invalid segments-method, %s. Options are "
    #                     "ALL_SINGLE_IFO_TIME and COINC_TIME" % segments_method)

    for ifo in workflow.ifos:
        sci_ok_segs[ifo] = \
                      sci_ok_seg_file.segment_dict[ifo + ':' + sci_ok_seg_name]

    logging.info('Done generating analyzable science segments')
    return sci_ok_seg_file, sci_ok_segs, sci_ok_seg_name


def get_cumulative_veto_group_files(workflow, option, cat_files,
                                    out_dir, execute_now=True, tags=None):
    """
    Get the cumulative veto files that define the different backgrounds
    we want to analyze, defined by groups of vetos.

    Parameters
    -----------
    workflow : Workflow object
        Instance of the workflow object
    option : str
        ini file option to use to get the veto groups
    cat_files : FileList of SegFiles
        The category veto files generated by get_veto_segs
    out_dir : path
        Location to store output files
    execute_now : Boolean
        If true outputs are generated at runtime. Else jobs go into the workflow
        and are generated then.
    tags : list of strings
        Used to retrieve subsections of the ini file for
        configuration options.

    Returns
    --------
    seg_files : workflow.core.FileList instance
        The cumulative segment files for each veto group.
    names : list of strings
        The segment names for the corresponding seg_file
    cat_files : workflow.core.FileList instance
        The list of individual category veto files
    """
    if tags is None:
        tags = []
    logging.info("Starting generating vetoes for groups in %s" %(option))
    make_analysis_dir(out_dir)

    cat_sets = parse_cat_ini_opt(workflow.cp.get_opt_tags('workflow-segments',
                                            option, tags))

    cum_seg_files = FileList()
    names = []
    for cat_set in cat_sets:
        segment_name = "CUMULATIVE_CAT_%s" % (''.join(sorted(cat_set)))
        logging.info('getting information for %s' % segment_name)
        categories = [cat_to_veto_def_cat(c) for c in cat_set]

        cum_seg_files += [get_cumulative_segs(workflow, categories, cat_files,
                          out_dir, execute_now=execute_now,
                          segment_name=segment_name, tags=tags)]
        names.append(segment_name)

    logging.info("Done generating vetoes for groups in %s" %(option))

    return cum_seg_files, names, cat_files

def setup_segment_generation(workflow, out_dir, tag=None):
    """
    This function is the gateway for setting up the segment generation steps in a
    workflow. It is designed to be able to support multiple ways of obtaining
    these segments and to combine/edit such files as necessary for analysis.
    The current modules have the capability to generate files at runtime or to
    generate files that are not needed for workflow generation within the workflow.

    Parameters
    -----------
    workflow : pycbc.workflow.core.Workflow
        The workflow instance that the coincidence jobs will be added to.
        This instance also contains the ifos for which to attempt to obtain
        segments for this analysis and the start and end times to search for
        segments over.
    out_dir : path
        The directory in which output will be stored.
    tag : string, optional (default=None)
        Use this to specify a tag. This can be used if this module is being
        called more than once to give call specific configuration (by setting
        options in [workflow-datafind-${TAG}] rather than [workflow-datafind]). This
        is also used to tag the Files returned by the class to uniqueify
        the Files and uniqueify the actual filename.
        FIXME: Filenames may not be unique with current codes!

    Returns
    -------
    segsToAnalyse : dictionay of ifo-keyed glue.segment.segmentlist instances
        This will contain the times that your code should analyse. By default this
        is science time - CAT_1 vetoes. (This default could be changed if desired)
    segFilesList : pycbc.workflow.core.FileList of SegFile instances
        These are representations of the various segment files that were constructed
        at this stage of the workflow and may be needed at later stages of the
        analysis (e.g. for performing DQ vetoes). If the file was generated at
        run-time the segment lists contained within these files will be an attribute
        of the instance. (If it will be generated in the workflow it will not be
        because I am not psychic).
    """
    logging.info("Entering segment generation module")
    make_analysis_dir(out_dir)

    cp = workflow.cp

    # Parse for options in ini file
    segmentsMethod = cp.get_opt_tags("workflow-segments",
                                     "segments-method", [tag])
    # These only needed if calling setup_segment_gen_mixed
    if segmentsMethod in ['AT_RUNTIME','CAT2_PLUS_DAG','CAT3_PLUS_DAG',
                          'CAT4_PLUS_DAG']:
        veto_cats = cp.get_opt_tags("workflow-segments",
                                    "segments-veto-categories", [tag])
        max_veto_cat = max([int(c) for c in veto_cats.split(',')])
        veto_categories = range(1, max_veto_cat + 1)
        if cp.has_option_tags("workflow-segments",
                              "segments-generate-coincident-segments", [tag]):
            generate_coincident_segs = True
        else:
            generate_coincident_segs = False
        # Need to curl the veto-definer file
        vetoDefUrl = cp.get_opt_tags("workflow-segments",
                                     "segments-veto-definer-url", [tag])
        vetoDefBaseName = os.path.basename(vetoDefUrl)
        vetoDefNewPath = os.path.join(out_dir, vetoDefBaseName)
        resolve_url(vetoDefUrl,out_dir)
        # and update location
        cp.set("workflow-segments", "segments-veto-definer-file",
                vetoDefNewPath)


    if cp.has_option_tags("workflow-segments",
                          "segments-minimum-segment-length", [tag]):
        minSegLength = int( cp.get_opt_tags("workflow-segments",
                           "segments-minimum-segment-length", [tag]) )
    else:
        minSegLength = 0

    if segmentsMethod == "AT_RUNTIME":
        max_veto = 1000
    elif segmentsMethod == "CAT2_PLUS_DAG":
        max_veto = 1
    elif segmentsMethod == "CAT3_PLUS_DAG":
        max_veto = 2
    elif segmentsMethod == "CAT4_PLUS_DAG":
        max_veto = 3
    else:
        msg = "Entry segments-method in [workflow-segments] does not have "
        msg += "expected value. Valid values are AT_RUNTIME, CAT4_PLUS_DAG, "
        msg += "CAT2_PLUS_DAG or CAT3_PLUS_DAG."
        raise ValueError(msg)

    logging.info("Generating segments with setup_segment_gen_mixed")
    segFilesList = setup_segment_gen_mixed(workflow, veto_categories,
                             out_dir, max_veto, tag=tag,
                             generate_coincident_segs=generate_coincident_segs)
    logging.info("Segments obtained")

    # This creates the segsToAnalyse from the segFilesList. Currently it uses
    # the 'SCIENCE_OK' segFilesList, which is science - CAT_1 in
    # setup_segment_gen_mixed.
    # This also applies the minimum science length
    segsToAnalyse = {}
    for ifo in workflow.ifos:
        analSegs = segFilesList.find_output_with_ifo(ifo)
        analSegs = analSegs.find_output_with_tag('SCIENCE_OK')
        assert len(analSegs) == 1
        analSegs = analSegs[0]
        if analSegs.segment_list:
            if minSegLength:
                analSegs.remove_short_sci_segs(minSegLength)
                analSegs.to_segment_xml(override_file_if_exists=True)
            segsToAnalyse[ifo] = analSegs.segment_list
        else:
            msg = "No science segments found for ifo %s. " %(ifo)
            msg += "If this is unexpected check the files that were dumped "
            msg += "in the %s directory. Also the " %(out_dir)
            msg += "commands that can be used to reproduce some of these "
            msg += "in %s/*.sh" %(os.path.join(out_dir,'logs'))
            logging.warn(msg)

    logging.info("Leaving segment generation module")
    return segsToAnalyse, segFilesList

def setup_segment_gen_mixed(workflow, veto_categories, out_dir,
                            maxVetoAtRunTime, tag=None,
                            generate_coincident_segs=True):
    """
    This function will generate veto files for each ifo and for each veto
    category.
    It can generate these vetoes at run-time or in the workflow (or do some at
    run-time and some in the workflow). However, the CAT_1 vetoes and science
    time must be generated at run time as they are needed to plan the workflow.
    CATs 2 and higher *may* be needed for other workflow construction.
    It can also combine these files to create a set of cumulative,
    multi-detector veto files, which can be used in ligolw_thinca and in
    pipedown. Again these can be created at run time or within the workflow.

    Parameters
    -----------
    workflow : pycbc.workflow.core.Workflow
        The Workflow instance that the coincidence jobs will be added to.
        This instance also contains the ifos for which to attempt to obtain
        segments for this analysis and the start and end times to search for
        segments over.
    veto_categories : list of ints
        List of veto categories to generate segments for. If this stops being
        integers, this can be changed here.
    out_dir : path
        The directory in which output will be stored.
    maxVetoAtRunTime : int
        Generate veto files at run time up to this category. Veto categories
        beyond this in veto_categories will be generated in the workflow.
        If we move to a model where veto
        categories are not explicitly cumulative, this will be rethought.
    tag : string, optional (default=None)
        Use this to specify a tag. This can be used if this module is being
        called more than once to give call specific configuration (by setting
        options in [workflow-datafind-${TAG}] rather than [workflow-datafind]). This
        is also used to tag the Files returned by the class to uniqueify
        the Files and uniqueify the actual filename.
        FIXME: Filenames may not be unique with current codes!
    generate_coincident_segs : boolean, optional (default = True)
        If given this module will generate a set of coincident, cumulative veto
        files that can be used with ligolw_thinca and pipedown.

    Returns
    -------
    segFilesList : dictionary of pycbc.workflow.core.SegFile instances
        These are representations of the various segment files that were
        constructed
        at this stage of the workflow and may be needed at later stages of the
        analysis (e.g. for performing DQ vetoes). If the file was generated at
        run-time the segment lists contained within these files will be an
        attribute
        of the instance. (If it will be generated in the workflow it will
        not be because I am not psychic).
    """
    cp = workflow.cp
    segFilesList = FileList([])
    start_time = workflow.analysis_time[0]
    end_time = workflow.analysis_time[1]
    segValidSeg = workflow.analysis_time
    # Will I need to add some jobs to the workflow?
    vetoGenJob = create_segs_from_cats_job(cp, out_dir, workflow.ifo_string)

    for ifo in workflow.ifos:
        logging.info("Generating science segments for ifo %s" %(ifo))
        currSciSegs, currSciXmlFile, _ = get_sci_segs_for_ifo(ifo, cp,
                                        start_time, end_time, out_dir, tags=tag)
        segFilesList.append(currSciXmlFile)

        for category in veto_categories:
            if category > maxVetoAtRunTime:
                msg = "Adding creation of CAT_%d segments " %(category)
                msg += "for ifo %s to workflow." %(ifo)
                logging.info(msg)
                execute_status = False

            if category <= maxVetoAtRunTime:
                logging.info("Generating CAT_%d segments for ifo %s." \
                             %(category,ifo))
                execute_status = True

            currVetoXmlFile = get_veto_segs(workflow, ifo, category,
                                                start_time, end_time, out_dir,
                                                vetoGenJob,
                                                execute_now=execute_status)

            segFilesList.append(currVetoXmlFile)
            # Store the CAT_1 veto segs for use below
            if category == 1:
                cat1Segs = currVetoXmlFile.return_union_seglist()

        analysedSegs = currSciSegs - cat1Segs
        analysedSegs.coalesce()
        analysedSegDict = segments.segmentlistdict()
        analysedSegDict[ifo + ':SCIENCE_OK'] = analysedSegs
        analysedXmlFile = os.path.join(out_dir,
                             "%s-SCIENCE_OK_SEGMENTS.xml" %(ifo.upper()) )
        currUrl = urlunparse(['file', 'localhost', analysedXmlFile,
                              None, None, None])
        if tag:
            currTags = [tag, 'SCIENCE_OK']
        else:
            currTags = ['SCIENCE_OK']
        currFile = SegFile(ifo, 'SEGMENTS', analysedSegs,
                           segment_dict=analysedSegDict, file_url=currUrl,
                           tags=currTags)
        segFilesList.append(currFile)
        currFile.to_segment_xml()


    if generate_coincident_segs:
        # Need to make some combined category veto files to use when vetoing
        # segments and triggers.
        ifo_string = workflow.ifo_string
        categories = []
        cum_cat_files = []
        for category in veto_categories:
            categories.append(category)
            # Set file name in workflow standard
            if tag:
                currTags = [tag]
            else:
                currTags = []

            # And actually make the file (or queue it in the workflow)
            logging.info("Generating combined, cumulative CAT_%d segments."\
                             %(category))
            if category <= maxVetoAtRunTime:
                execute_status = True
            else:
                execute_status = False
            currSegFile = get_cumulative_segs(workflow, categories,
                                segFilesList, out_dir,
                                execute_now=execute_status, tags=currTags)

            segFilesList.append(currSegFile)
            cum_cat_files.append(currSegFile)
        # Create a combined file
        # Set file tag in workflow standard
        if tag:
            currTags = [tag, 'COMBINED_CUMULATIVE_SEGMENTS']
        else:
            currTags = ['COMBINED_CUMULATIVE_SEGMENTS']

        combined_veto_file = os.path.join(out_dir,
                               '%s-CUMULATIVE_ALL_CATS_SEGMENTS.xml' \
                               %(ifo_string) )
        curr_url = urlunparse(['file', 'localhost',
                               combined_veto_file, None, None, None])
        curr_file = SegFile(ifo_string, 'SEGMENTS', segValidSeg,
                            file_url=curr_url, tags=currTags)

        for category in veto_categories:
            if category <= maxVetoAtRunTime:
                execute_status = True
                break
        else:
            execute_status = False
        add_cumulative_files(workflow, curr_file, cum_cat_files, out_dir,
                             execute_now=execute_status)
        segFilesList.append(curr_file)

    return segFilesList

def get_sci_segs_for_ifo(ifo, cp, start_time, end_time, out_dir, tags=None):
    """
    Obtain science segments for the selected ifo

    Parameters
    -----------
    ifo : string
        The string describing the ifo to obtain science times for.
    start_time : gps time (either int/LIGOTimeGPS)
        The time at which to begin searching for segments.
    end_time : gps time (either int/LIGOTimeGPS)
        The time at which to stop searching for segments.
    out_dir : path
        The directory in which output will be stored.
    tag : string, optional (default=None)
        Use this to specify a tag. This can be used if this module is being
        called more than once to give call specific configuration (by setting
        options in [workflow-datafind-${TAG}] rather than [workflow-datafind]).
        This is also used to tag the Files returned by the class to uniqueify
        the Files and uniqueify the actual filename.

    Returns
    --------
    sci_segs : ligo.segments.segmentlist
        The segmentlist generated by this call
    sci_xml_file : pycbc.workflow.core.SegFile
        The workflow File object corresponding to this science segments file.
    out_sci_seg_name : string
        The name of the output segment list in the output XML file.

    """
    if tags is None:
        tags = []
    seg_valid_seg = segments.segment([start_time,end_time])
    sci_seg_name = cp.get_opt_tags(
        "workflow-segments", "segments-%s-science-name" %(ifo.lower()), tags)
    sci_seg_url = cp.get_opt_tags(
        "workflow-segments", "segments-database-url", tags)
    # NOTE: ligolw_segment_query returns slightly strange output. The output
    #       segment list is put in with name "RESULT". So this is hardcoded here
    out_sci_seg_name = "RESULT"
    if tags:
        sci_xml_file_path = os.path.join(
            out_dir, "%s-SCIENCE_SEGMENTS_%s.xml" \
                     %(ifo.upper(), '_'.join(tags)))
        tag_list=tags + ['SCIENCE']
    else:
        sci_xml_file_path = os.path.join(
            out_dir, "%s-SCIENCE_SEGMENTS.xml" %(ifo.upper()) )
        tag_list = ['SCIENCE']

    if file_needs_generating(sci_xml_file_path, cp, tags=tags):
        seg_find_call = [ resolve_url(cp.get("executables","segment_query"),
                permissions=stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR),
            "--query-segments",
            "--segment-url", sci_seg_url,
            "--gps-start-time", str(start_time),
            "--gps-end-time", str(end_time),
            "--include-segments", sci_seg_name,
            "--output-file", sci_xml_file_path ]

        make_external_call(seg_find_call, out_dir=os.path.join(out_dir,'logs'),
                                out_basename='%s-science-call' %(ifo.lower()) )

    # Yes its yucky to generate a file and then read it back in.
    sci_xml_file_path = os.path.abspath(sci_xml_file_path)
    sci_xml_file = SegFile.from_segment_xml(sci_xml_file_path, tags=tag_list,
                                        valid_segment=seg_valid_seg)
    # NOTE: ligolw_segment_query returns slightly strange output. The output
    #       segment_summary output does not use RESULT. Therefore move the
    #       segment_summary across.
    sci_xml_file.seg_summ_dict[ifo.upper() + ":" + out_sci_seg_name] = \
             sci_xml_file.seg_summ_dict[':'.join(sci_seg_name.split(':')[0:2])]

    sci_segs = sci_xml_file.return_union_seglist()
    return sci_segs, sci_xml_file, out_sci_seg_name

def get_veto_segs(workflow, ifo, category, start_time, end_time, out_dir,
                  veto_gen_job, tags=None, execute_now=False):
    """
    Obtain veto segments for the selected ifo and veto category and add the job
    to generate this to the workflow.

    Parameters
    -----------
    workflow: pycbc.workflow.core.Workflow
        An instance of the Workflow class that manages the workflow.
    ifo : string
        The string describing the ifo to generate vetoes for.
    category : int
        The veto category to generate vetoes for.
    start_time : gps time (either int/LIGOTimeGPS)
        The time at which to begin searching for segments.
    end_time : gps time (either int/LIGOTimeGPS)
        The time at which to stop searching for segments.
    out_dir : path
        The directory in which output will be stored.
    vetoGenJob : Job
        The veto generation Job class that will be used to create the Node.
    tag : string, optional (default=None)
        Use this to specify a tag. This can be used if this module is being
        called more than once to give call specific configuration (by setting
        options in [workflow-datafind-${TAG}] rather than [workflow-datafind]).
        This is also used to tag the Files returned by the class to uniqueify
        the Files and uniqueify the actual filename.
        FIXME: Filenames may not be unique with current codes!
    execute_now : boolean, optional
        If true, jobs are executed immediately. If false, they are added to the
        workflow to be run later.

    Returns
    --------
    veto_def_file : pycbc.workflow.core.SegFile
        The workflow File object corresponding to this DQ veto file.
    """
    if tags is None:
        tags = []
    seg_valid_seg = segments.segment([start_time,end_time])
    # FIXME: This job needs an internet connection and X509_USER_PROXY
    #        For internet connection, it may need a headnode (ie universe local)
    #        For X509_USER_PROXY, I don't know what pegasus is doing
    node = Node(veto_gen_job)
    node.add_opt('--veto-categories', str(category))
    node.add_opt('--ifo-list', ifo)
    node.add_opt('--gps-start-time', str(start_time))
    node.add_opt('--gps-end-time', str(end_time))
    if tags:
        veto_xml_file_name = "%s-VETOTIME_CAT%d_%s-%d-%d.xml" \
                               %(ifo, category, '_'.join(tags), start_time,
                                 end_time-start_time)
    else:
        veto_xml_file_name = "%s-VETOTIME_CAT%d-%d-%d.xml" \
                         %(ifo, category, start_time, end_time-start_time)
    veto_xml_file_path = os.path.abspath(os.path.join(out_dir,
                                         veto_xml_file_name))
    curr_url = urlunparse(['file', 'localhost',
                           veto_xml_file_path, None, None, None])
    if tags:
        curr_tags = tags + ['VETO_CAT%d' %(category)]
    else:
        curr_tags = ['VETO_CAT%d' %(category)]

    if file_needs_generating(veto_xml_file_path, workflow.cp, tags=tags):
        if execute_now:
            workflow.execute_node(node, verbatim_exe = True)
            veto_xml_file = SegFile.from_segment_xml(veto_xml_file_path,
                                                 tags=curr_tags,
                                                 valid_segment=seg_valid_seg)
        else:
            veto_xml_file = SegFile(ifo, 'SEGMENTS', seg_valid_seg,
                                    file_url=curr_url, tags=curr_tags)
            node._add_output(veto_xml_file)
            workflow.add_node(node)
    else:
        node.executed = True
        for fil in node._outputs:
            fil.node = None
        veto_xml_file = SegFile.from_segment_xml(veto_xml_file_path,
                                                 tags=curr_tags,
                                                 valid_segment=seg_valid_seg)
    return veto_xml_file

def create_segs_from_cats_job(cp, out_dir, ifo_string, tags=None):
    """
    This function creates the CondorDAGJob that will be used to run
    ligolw_segments_from_cats as part of the workflow

    Parameters
    -----------
    cp : pycbc.workflow.configuration.WorkflowConfigParser
        The in-memory representation of the configuration (.ini) files
    out_dir : path
        Directory in which to put output files
    ifo_string : string
        String containing all active ifos, ie. "H1L1V1"
    tag : list of strings, optional (default=None)
        Use this to specify a tag(s). This can be used if this module is being
        called more than once to give call specific configuration (by setting
        options in [workflow-datafind-${TAG}] rather than [workflow-datafind]).
        This is also used to tag the Files returned by the class to uniqueify
        the Files and uniqueify the actual filename.
        FIXME: Filenames may not be unique with current codes!

    Returns
    --------
    job : Job instance
        The Job instance that will run segments_from_cats jobs
    """
    if tags is None:
        tags = []

    seg_server_url = cp.get_opt_tags("workflow-segments",
                                   "segments-database-url", tags)
    veto_def_file = cp.get_opt_tags("workflow-segments",
                                  "segments-veto-definer-file", tags)

    job = Executable(cp, 'segments_from_cats', universe='local',
                               ifos=ifo_string, out_dir=out_dir, tags=tags)
    job.add_opt('--separate-categories')
    job.add_opt('--segment-url', seg_server_url)

    job.add_opt('--veto-file', veto_def_file)
    # FIXME: Would like the proxy in the Workflow instance
    # FIXME: Explore using the x509 condor commands
    # If the user has a proxy set in the environment, add it to the job
    return job

def get_cumulative_segs(workflow, categories, seg_files_list, out_dir,
                        tags=None, execute_now=False, segment_name=None):
    """
    Function to generate one of the cumulative, multi-detector segment files
    as part of the workflow.

    Parameters
    -----------
    workflow: pycbc.workflow.core.Workflow
        An instance of the Workflow class that manages the workflow.
    categories : int
        The veto categories to include in this cumulative veto.
    seg_files_list : Listionary of SegFiles
        The list of segment files to be used as input for combining.
    out_dir : path
        The directory to write output to.
    tags : list of strings, optional
        A list of strings that is used to identify this job
    execute_now : boolean, optional
        If true, jobs are executed immediately. If false, they are added to the
        workflow to be run later.
    segment_name : str
        The name of the combined, cumulative segments in the output file.
    """
    if tags is None:
        tags = []
    add_inputs = FileList([])
    valid_segment = workflow.analysis_time
    if segment_name is None:
        segment_name = 'VETO_CAT%d_CUMULATIVE' % (categories[-1])
    cp = workflow.cp
    # calculate the cumulative veto files for a given ifo
    for ifo in workflow.ifos:
        cum_job = LigoLWCombineSegsExecutable(cp, 'ligolw_combine_segments',
                       out_dir=out_dir, tags=[segment_name]+tags, ifos=ifo)
        inputs = []
        files = seg_files_list.find_output_with_ifo(ifo)
        for category in categories:
            file_list = files.find_output_with_tag('VETO_CAT%d' %(category))
            inputs+=file_list

        cum_node  = cum_job.create_node(valid_segment, inputs, segment_name)
        if file_needs_generating(cum_node.output_files[0].cache_entry.path,
                                 workflow.cp, tags=tags):
            if execute_now:
                workflow.execute_node(cum_node)
            else:
                workflow.add_node(cum_node)
        else:
            cum_node.executed = True
            for fil in cum_node._outputs:
                fil.node = None
                fil.PFN(urljoin('file:', pathname2url(fil.storage_path)),
                        site='local')
        add_inputs += cum_node.output_files

    # add cumulative files for each ifo together
    name = '%s_VETO_SEGMENTS' %(segment_name)
    outfile = File(workflow.ifos, name, workflow.analysis_time,
                                            directory=out_dir, extension='xml',
                                            tags=[segment_name] + tags)
    add_job = LigolwAddExecutable(cp, 'llwadd', ifos=ifo, out_dir=out_dir,
                                  tags=tags)
    add_node = add_job.create_node(valid_segment, add_inputs, output=outfile)
    if file_needs_generating(add_node.output_files[0].cache_entry.path,
                             workflow.cp, tags=tags):
        if execute_now:
            workflow.execute_node(add_node)
        else:
            workflow.add_node(add_node)
    else:
        add_node.executed = True
        for fil in add_node._outputs:
            fil.node = None
            fil.PFN(urljoin('file:', pathname2url(fil.storage_path)),
                    site='local')
    return outfile

def add_cumulative_files(workflow, output_file, input_files, out_dir,
                         execute_now=False, tags=None):
    """
    Function to combine a set of segment files into a single one. This function
    will not merge the segment lists but keep each separate.

    Parameters
    -----------
    workflow: pycbc.workflow.core.Workflow
        An instance of the Workflow class that manages the workflow.
    output_file: pycbc.workflow.core.File
        The output file object
    input_files: pycbc.workflow.core.FileList
        This list of input segment files
    out_dir : path
        The directory to write output to.
    execute_now : boolean, optional
        If true, jobs are executed immediately. If false, they are added to the
        workflow to be run later.
    tags : list of strings, optional
        A list of strings that is used to identify this job
    """
    if tags is None:
        tags = []
    llwadd_job = LigolwAddExecutable(workflow.cp, 'llwadd',
                       ifo=output_file.ifo_list, out_dir=out_dir, tags=tags)
    add_node = llwadd_job.create_node(output_file.segment, input_files,
                                   output=output_file)
    if file_needs_generating(add_node.output_files[0].cache_entry.path,
                             workflow.cp, tags=tags):
        if execute_now:
            workflow.execute_node(add_node)
        else:
            workflow.add_node(add_node)
    else:
        add_node.executed = True
        for fil in add_node._outputs:
            fil.node = None
            fil.PFN(urljoin('file:', pathname2url(fil.storage_path)),
                    site='local')
    return add_node.output_files[0]

def find_playground_segments(segs):
    '''Finds playground time in a list of segments.

      Playground segments include the first 600s of every 6370s stride starting
      at GPS time 729273613.

      Parameters
      ----------
      segs : segmentfilelist
          A segmentfilelist to find playground segments.

      Returns
      -------
      outlist : segmentfilelist
          A segmentfilelist with all playground segments during the input
          segmentfilelist (ie. segs).
    '''

    # initializations
    start_s2 = 729273613
    playground_stride = 6370
    playground_length = 600
    outlist = segments.segmentlist()

    # loop over segments
    for seg in segs:
        start = seg[0]
        end = seg[1]

        # the first playground segment whose end is after the start of seg
        playground_start = start_s2 + playground_stride * ( 1 + \
                     int(start-start_s2-playground_length) / playground_stride)

        while playground_start < end:
            # find start of playground segment
            if playground_start > start:
                ostart = playground_start
            else:
                ostart = start

            playground_end = playground_start + playground_length

            # find end of playground segment
            if playground_end < end:
                oend = playground_end
            else:
                oend = end

            # append segment
            x = segments.segment(ostart, oend)
            outlist.append(x)

            # increment
            playground_start = playground_start + playground_stride

    return outlist

def get_triggered_coherent_segment(workflow, sciencesegs):
    """
    Construct the coherent network on and off source segments. Can switch to
    construction of segments for a single IFO search when coherent segments
    are insufficient for a search.

    Parameters
    -----------
    workflow : pycbc.workflow.core.Workflow
        The workflow instance that the calculated segments belong to.
    sciencesegs : dict
        Dictionary of all science segments within analysis time.

    Returns
    --------
    onsource : ligo.segments.segmentlistdict
        A dictionary containing the on source segments for network IFOs

    offsource : ligo.segments.segmentlistdict
        A dictionary containing the off source segments for network IFOs
    """

    # Load parsed workflow config options
    cp = workflow.cp
    triggertime = int(os.path.basename(cp.get('workflow', 'trigger-time')))
    minduration = int(os.path.basename(cp.get('workflow-exttrig_segments',
                                              'min-duration')))
    maxduration = int(os.path.basename(cp.get('workflow-exttrig_segments',
                                              'max-duration')))
    onbefore = int(os.path.basename(cp.get('workflow-exttrig_segments',
                                           'on-before')))
    onafter = int(os.path.basename(cp.get('workflow-exttrig_segments',
                                          'on-after')))
    padding = int(os.path.basename(cp.get('workflow-exttrig_segments',
                                          'pad-data')))
    if cp.has_option("workflow-condition_strain", "do-gating"):
        padding += int(os.path.basename(cp.get("condition_strain",
                                               "pad-data")))
    quanta = int(os.path.basename(cp.get('workflow-exttrig_segments',
                                         'quanta')))

    # Check available data segments meet criteria specified in arguments
    commonsegs = sciencesegs.extract_common(sciencesegs.keys())
    offsrclist = commonsegs[commonsegs.keys()[0]]
    if len(offsrclist) > 1:
        logging.info("Removing network segments that do not contain trigger "
                     "time")
        for seg in offsrclist:
            if triggertime in seg:
                offsrc = seg
    else:
        offsrc = offsrclist[0]

    if abs(offsrc) < minduration + 2 * padding:
        fail = segments.segment([triggertime - minduration / 2. - padding,
                                 triggertime + minduration / 2. + padding])
        logging.warning("Available network segment shorter than minimum "
                        "allowed duration.")
        return None, fail

    # Will segment duration be the maximum desired length or not?
    if abs(offsrc) >= maxduration + 2 * padding:
        logging.info("Available network science segment duration (%ds) is "
                     "greater than the maximum allowed segment length (%ds). "
                     "Truncating..." % (abs(offsrc), maxduration))
    else:
        logging.info("Available network science segment duration (%ds) is "
                     "less than the maximum allowed segment length (%ds)."
                     % (abs(offsrc), maxduration))

    logging.info("%ds of padding applied at beginning and end of segment."
                 % padding)


    # Construct on-source
    onstart = triggertime - onbefore
    onend = triggertime + onafter
    oncentre = onstart + ((onbefore + onafter) / 2)
    onsrc = segments.segment(onstart, onend)
    logging.info("Constructed ON-SOURCE: duration %ds (%ds before to %ds after"
                 " trigger)."
                 % (abs(onsrc), triggertime - onsrc[0],
                    onsrc[1] - triggertime))
    onsrc = segments.segmentlist([onsrc])

    # Maximal, centred coherent network segment
    idealsegment = segments.segment(int(oncentre - padding -
                                    0.5 * maxduration),
                                    int(oncentre + padding +
                                    0.5 * maxduration))

    # Construct off-source
    if (idealsegment in offsrc):
        offsrc = idealsegment

    elif idealsegment[1] not in offsrc:
        offsrc &= segments.segment(offsrc[1] - maxduration - 2 * padding,
                                   offsrc[1])

    elif idealsegment[0] not in offsrc:
        offsrc &= segments.segment(offsrc[0],
                                   offsrc[0] + maxduration + 2 * padding)

    # Trimming off-source
    excess = (abs(offsrc) - 2 * padding) % quanta
    if excess != 0:
        logging.info("Trimming %ds excess time to make OFF-SOURCE duration a "
                     "multiple of %ds" % (excess, quanta))
        offset = (offsrc[0] + abs(offsrc) / 2.) - oncentre
        if 2 * abs(offset) > excess:
            if offset < 0:
                offsrc &= segments.segment(offsrc[0] + excess,
                                           offsrc[1])
            elif offset > 0:
                offsrc &= segments.segment(offsrc[0],
                                           offsrc[1] - excess)
            assert abs(offsrc) % quanta == 2 * padding
        else:
            logging.info("This will make OFF-SOURCE symmetrical about trigger "
                         "time.")
            start = int(offsrc[0] - offset + excess / 2)
            end = int(offsrc[1] - offset - round(float(excess) / 2))
            offsrc = segments.segment(start, end)
            assert abs(offsrc) % quanta == 2 * padding

    logging.info("Constructed OFF-SOURCE: duration %ds (%ds before to %ds "
                 "after trigger)."
                 % (abs(offsrc) - 2 * padding,
                    triggertime - offsrc[0] - padding,
                    offsrc[1] - triggertime - padding))
    offsrc = segments.segmentlist([offsrc])

    # Put segments into segmentlistdicts
    onsource = segments.segmentlistdict()
    offsource = segments.segmentlistdict()
    ifos = ''
    for iifo in sciencesegs.keys():
        ifos += str(iifo)
        onsource[iifo] = onsrc
        offsource[iifo] = offsrc

    return onsource, offsource

def generate_triggered_segment(workflow, out_dir, sciencesegs):
    cp = workflow.cp

    if cp.has_option("workflow", "allow-single-ifo-search"):
        min_ifos = 1
    else:
        min_ifos = 2

    triggertime = int(os.path.basename(cp.get('workflow', 'trigger-time')))
    minbefore = int(os.path.basename(cp.get('workflow-exttrig_segments',
                                            'min-before')))
    minafter = int(os.path.basename(cp.get('workflow-exttrig_segments',
                                           'min-after')))
    minduration = int(os.path.basename(cp.get('workflow-exttrig_segments',
                                              'min-duration')))
    onbefore = int(os.path.basename(cp.get('workflow-exttrig_segments',
                                           'on-before')))
    onafter = int(os.path.basename(cp.get('workflow-exttrig_segments',
                                          'on-after')))
    padding = int(os.path.basename(cp.get('workflow-exttrig_segments',
                                          'pad-data')))
    if cp.has_option("workflow-condition_strain", "do-gating"):
        padding += int(os.path.basename(cp.get("condition_strain",
                                               "pad-data")))

    # How many IFOs meet minimum data requirements?
    min_seg = segments.segment(triggertime - onbefore - minbefore - padding,
                               triggertime + onafter + minafter + padding)
    scisegs = segments.segmentlistdict({ifo: sciencesegs[ifo]
            for ifo in sciencesegs.keys() if min_seg in sciencesegs[ifo]
            and abs(sciencesegs[ifo]) >= minduration})

    # Find highest number of IFOs that give an acceptable coherent segment
    num_ifos = len(scisegs.keys())
    while num_ifos >= min_ifos:
        # Consider all combinations for a given number of IFOs
        ifo_combos = itertools.combinations(scisegs.keys(), num_ifos)
        onsource = {}
        offsource = {}
        for ifo_combo in ifo_combos:
            ifos = "".join(ifo_combo)
            logging.info("Calculating optimal segment for %s.", ifos)
            segs = segments.segmentlistdict({ifo: scisegs[ifo]
                                             for ifo in ifo_combo})
            onsource[ifos], offsource[ifos] = get_triggered_coherent_segment(\
                    workflow, segs)

        # Which combination gives the longest coherent segment?
        valid_combs = [iifos for iifos in onsource.keys()
                       if onsource[iifos] is not None]

        if len(valid_combs) == 0:
            # If none, offsource dict will contain segments showing criteria
            # that have not been met, for use in plotting
            if len(offsource.keys()) > 1:
                seg_lens = {ifos: abs(next(offsource[ifos].itervalues())[0])
                            for ifos in offsource.keys()}
                best_comb = max(seg_lens.iterkeys(),
                                key=(lambda key: seg_lens[key]))
            else:
                best_comb = offsource.keys()[0]
            logging.info("No combination of %d IFOs with suitable science "
                         "segment.", num_ifos)
        else:
            # Identify best analysis segment
            if len(valid_combs) > 1:
                seg_lens = {ifos: abs(next(offsource[ifos].itervalues())[0])
                            for ifos in valid_combs}
                best_comb = max(seg_lens.iterkeys(),
                                key=(lambda key: seg_lens[key]))
            else:
                best_comb = valid_combs[0]
            logging.info("Calculated science segments.")

            offsourceSegfile = os.path.join(out_dir, "offSourceSeg.txt")
            segmentsUtils.tosegwizard(open(offsourceSegfile, "w"),
                                      next(offsource[best_comb].itervalues()))

            onsourceSegfile = os.path.join(out_dir, "onSourceSeg.txt")
            segmentsUtils.tosegwizard(file(onsourceSegfile, "w"),
                                      next(onsource[best_comb].itervalues()))

            bufferleft = int(cp.get('workflow-exttrig_segments',
                                    'num-buffer-before'))
            bufferright = int(cp.get('workflow-exttrig_segments',
                                     'num-buffer-after'))
            onlen = onbefore + onafter
            bufferSegment = segments.segment(\
                    triggertime - onbefore - bufferleft * onlen,
                    triggertime + onafter + bufferright * onlen)
            bufferSegfile = os.path.join(out_dir, "bufferSeg.txt")
            segmentsUtils.tosegwizard(file(bufferSegfile, "w"),
                                      segments.segmentlist([bufferSegment]))

            return onsource[best_comb], offsource[best_comb]

        num_ifos -= 1

    logging.warning("No suitable science segments available.")
    try:
        return None, offsource[best_comb]
    except UnboundLocalError:
        return None, min_seg

def save_veto_definer(cp, out_dir, tags=None):
    """ Retrieve the veto definer file and save it locally

    Parameters
    -----------
    cp : ConfigParser instance
    out_dir : path
    tags : list of strings
        Used to retrieve subsections of the ini file for
        configuration options.
    """
    if tags is None:
        tags = []
    make_analysis_dir(out_dir)
    veto_def_url = cp.get_opt_tags("workflow-segments",
                                 "segments-veto-definer-url", tags)
    veto_def_base_name = os.path.basename(veto_def_url)
    veto_def_new_path = os.path.abspath(os.path.join(out_dir,
                                        veto_def_base_name))
    # Don't need to do this if already done
    resolve_url(veto_def_url,out_dir)

    # and update location
    cp.set("workflow-segments", "segments-veto-definer-file", veto_def_new_path)
    return veto_def_new_path

def parse_cat_ini_opt(cat_str):
    """ Parse a cat str from the ini file into a list of sets """
    if cat_str == "":
        return []

    cat_groups = cat_str.split(',')
    cat_sets = []
    for group in cat_groups:
        group = group.strip()
        cat_sets += [set(c for c in group)]
    return cat_sets

def cat_to_veto_def_cat(val):
    """ Convert a category character to the corresponding value in the veto
    definer file.

    Parameters
    ----------
    str : single character string
        The input category character

    Returns
    -------
    pipedown_str : str
        The pipedown equivelant notation that can be passed to programs
    that expect this definition.
    """
    if val == '1':
        return 1
    if val == '2':
        return 2
    if val == '3':
        return 4
    if val == 'H':
        return 3
    else:
        raise ValueError('Invalid Category Choice')

def file_needs_generating(file_path, cp, tags=None):
    """
    This job tests the file location and determines if the file should be
    generated now or if an error should be raised. This uses the
    generate_segment_files variable, global to this module, which is described
    above and in the documentation.

    Parameters
    -----------
    file_path : path
        Location of file to check
    cp : ConfigParser
        The associated ConfigParser from which the
        segments-generate-segment-files variable is returned.
        It is recommended for most applications to use the default option by
        leaving segments-generate-segment-files blank, which will regenerate
        all segment files at runtime. Only use this facility if you need it.
        Choices are
        * 'always' : DEFAULT: All files will be generated even if they already exist.
        * 'if_not_present': Files will be generated if they do not already exist. Pre-existing files will be read in and used.
        * 'error_on_duplicate': Files will be generated if they do not already exist. Pre-existing files will raise a failure.
        * 'never': Pre-existing files will be read in and used. If no file exists the code will fail.

    Returns
    --------
    int
        1 = Generate the file. 0 = File already exists, use it. Other cases
        will raise an error.
    """
    if tags is None:
        tags = []
    if cp.has_option_tags("workflow-segments",
                          "segments-generate-segment-files", tags):
        value = cp.get_opt_tags("workflow-segments",
                                   "segments-generate-segment-files", tags)
        generate_segment_files = value
    else:
        generate_segment_files = 'always'

    # Does the file exist
    if os.path.isfile(file_path):
        if generate_segment_files in ['if_not_present', 'never']:
            return 0
        elif generate_segment_files == 'always':
            err_msg = "File %s already exists. " %(file_path,)
            err_msg += "Regenerating and overwriting."
            logging.warn(err_msg)
            return 1
        elif generate_segment_files == 'error_on_duplicate':
            err_msg = "File %s already exists. " %(file_path,)
            err_msg += "Refusing to overwrite file and exiting."
            raise ValueError(err_msg)
        else:
            err_msg = 'Global variable generate_segment_files must be one of '
            err_msg += '"always", "if_not_present", "error_on_duplicate", '
            err_msg += '"never". Got %s.' %(generate_segment_files,)
            raise ValueError(err_msg)
    else:
        if generate_segment_files in ['always', 'if_not_present',
                                      'error_on_duplicate']:
            return 1
        elif generate_segment_files == 'never':
            err_msg = 'File %s does not exist. ' %(file_path,)
            raise ValueError(err_msg)
        else:
            err_msg = 'Global variable generate_segment_files must be one of '
            err_msg += '"always", "if_not_present", "error_on_duplicate", '
            err_msg += '"never". Got %s.' %(generate_segment_files,)
            raise ValueError(err_msg)

def get_segments_file(workflow, name, option_name, out_dir):
    """Get cumulative segments from option name syntax for each ifo.

    Use syntax of configparser string to define the resulting segment_file
    e.x. option_name = +up_flag1,+up_flag2,+up_flag3,-down_flag1,-down_flag2
    Each ifo may have a different string and is stored separately in the file.
    Flags which add time must precede flags which subtract time.

    Parameters
    ----------
    workflow: pycbc.workflow.Workflow
    name: string
        Name of the segment list being created
    option_name: str
        Name of option in the associated config parser to get the flag list

    returns
    --------
    seg_file: pycbc.workflow.SegFile
        SegFile intance that points to the segment xml file on disk.
    """
    from pycbc.dq import query_str
    make_analysis_dir(out_dir)
    cp = workflow.cp
    start = workflow.analysis_time[0]
    end = workflow.analysis_time[1]

    # Check for veto definer file
    veto_definer = None
    if cp.has_option("workflow-segments", "segments-veto-definer-url"):
        veto_definer = save_veto_definer(workflow.cp, out_dir, [])

    # Check for provided server
    server = "segments.ligo.org"
    if cp.has_option("workflow-segments", "segments-database-url"):
        server = cp.get("workflow-segments",
                                 "segments-database-url")

    segs = {}
    for ifo in workflow.ifos:
        flag_str = cp.get_opt_tags("workflow-segments", option_name, [ifo])
        key = ifo + ':' + name
        segs[key] = query_str(ifo, flag_str, start, end,
                              server=server,
                              veto_definer=veto_definer)
        logging.info("%s: got %s flags", ifo, option_name)

    return SegFile.from_segment_list_dict(name, segs,
                                          extension='.xml',
                                          valid_segment=workflow.analysis_time,
                                          directory=out_dir)
