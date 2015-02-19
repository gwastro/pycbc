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
https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/NOTYETCREATED.html
"""


from __future__ import division

import re
import os
import os.path
import logging
from glue import segments
from glue.ligolw import lsctables,ligolw
from glue.ligolw import utils as ligolw_utils
from pycbc.workflow.core import FileList, make_analysis_dir, Executable, Node
from pycbc.workflow.core import get_random_label
from pycbc.workflow.jobsetup import LigolwAddExecutable, LigolwSSthincaExecutable, SQLInOutExecutable
from pylal import ligolw_cafe

class ContentHandler(ligolw.LIGOLWContentHandler):
        pass

lsctables.use_in(ContentHandler)

def setup_coincidence_workflow(workflow, segsList, timeSlideFiles,
                               inspiral_outs, output_dir, maxVetoCat=5,
                               tags=[], timeSlideTags=None):
    '''
    This function aims to be the gateway for setting up a set of coincidence
    jobs in a workflow. The goal is that this function can support a
    number of different ways/codes that could be used for doing this.
    For now it only supports ligolw_sstinca.

    Parameters
    -----------
    workflow : pycbc.workflow.core.Workflow
        The Workflow instance that the coincidence jobs will be added to.
    segsList : pycbc.workflow.core.FileList
        The list of files returned by workflow's segment module that contains
        pointers to all the segment files generated in the workflow. If the
        coincidence code will be applying the data quality vetoes, then this
        will be used to ensure that the codes get the necessary input to do
        this.
    timeSlideFiles : pycbc.workflow.core.FileList
        An FileList of the timeSlide input files that are needed to
        determine what time sliding needs to be done if the coincidence code
        will be running time slides to facilitate background computations later
        in the workflow.
    inspiral_outs : pycbc.workflow.core.FileList
        An FileList of the matched-filter module output that is used as
        input to the coincidence codes running at this stage.
    output_dir : path
        The directory in which coincidence output will be stored.
    maxVetoCat : int (optional, default=5)
        The maximum veto category that will be applied. If this takes the
        default value the code will run data quality at cumulative categories
        1, 2, 3, 4 and 5. Note that if we change the flag definitions to be
        non-cumulative then this option will need to be revisited,
    tags : list of strings (optional, default = [])
        A list of the tagging strings that will be used for all jobs created
        by this call to the workflow. An example might be ['BNSINJECTIONS'] or
        ['NOINJECTIONANALYSIS']. This will be used in output names.
    timeSlideTags : list of strings (optional, default = [])
        A list of the tags corresponding to the timeSlideFiles that are to be
        used in this call to the module. This can be used to ensure that the
        injection runs do no time sliding, but the no-injection runs do perform
        time slides (or vice-versa if you prefer!)
    Returns
    --------
    coinc_outs : pycbc.workflow.core.FileList
        A list of the *final* outputs of the coincident stage. This *does not*
        include any intermediate products produced within the workflow. If you
        require access to intermediate products call the various sub-functions
        in this module directly.
    '''
    logging.info('Entering coincidence setup module.')
    make_analysis_dir(output_dir)

    # Parse for options in .ini file
    coincidenceMethod = workflow.cp.get_opt_tags("workflow-coincidence",
                                        "coincidence-method", tags)

    # Scope here for adding different options/methods here. For now we only
    # have the single_stage ihope method which consists of using ligolw_add
    # to create a large job for coincidence and then running ligolw_thinca
    # on that output.
    if coincidenceMethod == "WORKFLOW_DISCRETE_SLIDES":
        # If I am doing exact match I can parallelize these jobs and reduce
        # memory footprint. This will require all input inspiral jobs to have
        # a JOB%d tag to distinguish between them.
        if workflow.cp.has_option_tags("workflow-coincidence",
                             "coincidence-exact-match-parallelize", tags):
            parallelize_split_input = True
        else:
            parallelize_split_input = False

        # If you want the ligolw_add outputs, call this function directly
        coinc_outs, other_outs = setup_coincidence_workflow_ligolw_thinca(
                     workflow,
                     segsList, timeSlideFiles, inspiral_outs,
                     output_dir, maxVetoCat=maxVetoCat, tags=tags,
                     timeSlideTags=timeSlideTags,
                     parallelize_split_input=parallelize_split_input)
    else:
        errMsg = "Coincidence method not recognized. Must be one of "
        errMsg += "WORKFLOW_DISCRETE_SLIDES (currently only one option)."
        raise ValueError(errMsg)

    logging.info('Leaving coincidence setup module.')

    return coinc_outs, other_outs

def setup_coincidence_workflow_ligolw_thinca(workflow, segsList,
                                             timeSlideFiles, inspiral_outs,
                                             output_dir, maxVetoCat=5, tags=[],
                                             timeSlideTags=None,
                                             parallelize_split_input=False):
    """
    This function is used to setup a single-stage ihope style coincidence stage
    of the workflow using ligolw_sstinca (or compatible code!).

    Parameters
    -----------
    workflow : pycbc.workflow.core.Workflow
        The workflow instance that the coincidence jobs will be added to.
    segsList : pycbc.workflow.core.FileList
        The list of files returned by workflow's segment module that contains
        pointers to all the segment files generated in the workflow. If the
        coincidence code will be applying the data quality vetoes, then this
        will be used to ensure that the codes get the necessary input to do
        this.
    timeSlideFiles : pycbc.workflow.core.FileList
        An FileList of the timeSlide input files that are needed to
        determine what time sliding needs to be done. One of the timeSlideFiles
        will normally be "zero-lag only", the others containing time slides
        used to facilitate background computations later in the workflow.
    inspiral_outs : pycbc.workflow.core.FileList
        An FileList of the matched-filter module output that is used as
        input to the coincidence codes running at this stage.
    output_dir : path
        The directory in which coincidence output will be stored.
    maxVetoCat : int (optional, default=5)
        The maximum veto category that will be applied. If this takes the
        default value the code will run data quality at cumulative categories
        1, 2, 3, 4 and 5. Note that if we change the flag definitions to be
        non-cumulative then this option will need to be revisited,
    tags : list of strings (optional, default = [])
        A list of the tagging strings that will be used for all jobs created
        by this call to the workflow. An example might be ['BNSINJECTIONS'] or
        ['NOINJECTIONANALYSIS']. This will be used in output names.
    timeSlideTags : list of strings (optional, default = [])
        A list of the tags corresponding to the timeSlideFiles that are to be
        used in this call to the module. This can be used to ensure that the
        injection runs do no time sliding, but the no-injection runs do perform
        time slides (or vice-versa if you prefer!)
    Returns
    --------
    ligolwThincaOuts : pycbc.workflow.core.FileList
        A list of the output files generated from ligolw_sstinca.
    ligolwAddOuts : pycbc.workflow.core.FileList
        A list of the output files generated from ligolw_add.
    """
    logging.debug("Entering coincidence module.")
    veto_categories = range(1,maxVetoCat+1)

    # setup code for each veto_category

    ligolwThincaOuts = FileList([])
    other_outs = {}

    if not timeSlideTags:
        # Get all sections by looking in ini file, use all time slide files.
        timeSlideTags = [(sec.split('-')[-1]).upper()
                  for sec in workflow.cp.sections() if sec.startswith('tisi-')]

    if parallelize_split_input:
        # Want to split all input jobs according to their JOB%d tag.
        # This matches any string that is the letters JOB followed by some
        # numbers and nothing else.
        inspiral_outs_dict = {}
        regex_match = re.compile('JOB([0-9]+)\Z')
        for file in inspiral_outs:
            matches = [regex_match.match(tag) for tag in file.tags]
            # Remove non matching entries
            matches = [i for i in matches if i is not None]
            # Must have one entry
            if len(matches) == 0:
                warn_msg = "I was asked to parallelize over split inspiral "
                warn_msg += "files at the coincidence stage, but at least one "
                warn_msg += "input file does not have a JOB\%d tag indicating "
                warn_msg += "that it was split. Assuming that I do not have "
                warn_msg += "split input files and turning "
                warn_msg += "parallelize_split_input off."
                logging.warn(warn_msg)
                parallelize_split_input = False
                break
            if len(matches) > 1:
                err_msg = "One of my input files has two tags fitting JOB\%d "
                err_msg += "this means I cannot tell which split job this "
                err_msg += "file is from."
                raise ValueError(err_msg)
            # Extract the job ID
            id = int(matches[0].string[3:])
            if not inspiral_outs_dict.has_key(id):
                inspiral_outs_dict[id] = FileList([])
            inspiral_outs_dict[id].append(file)
        else:
            # If I got through all the files I want to sort the dictionaries so
            # that file with key a and index 3 is the same file as key b and
            # index 3 other than the tag is JOBA -> JOBB ... ie. it has used
            # a different part of the template bank.
            sort_lambda = lambda x: (x.ifo_string, x.segment,
                                     x.tagged_description)
            for key in inspiral_outs_dict.keys():
                inspiral_outs_dict[id].sort(key = sort_lambda)
            # These should be in ascending order, so I can assume the existence
            # of a JOB0 tag
            inspiral_outs = inspiral_outs_dict[0]
            for index, file in enumerate(inspiral_outs):
                # Store the index in the file for quicker mapping later
                file.thinca_index = index
    else:
        inspiral_outs_dict = None

    for timeSlideTag in timeSlideTags:
        # Get the time slide file from the inputs
        tisiOutFile = timeSlideFiles.find_output_with_tag(timeSlideTag)
        if not len(tisiOutFile) == 1:
            errMsg = "If you are seeing this, something batshit is going on!"
            if len(tisiOutFile) == 0:
                errMsg = "No time slide files found matching %s." \
                                                                %(timeSlideTag)
            if len(tisiOutFile) > 1:
                errMsg = "More than one time slide files match %s." \
                                                                %(timeSlideTag)
            raise ValueError(errMsg)
        tisiOutFile = tisiOutFile[0]

        # Next we run ligolw_cafe. This is responsible for
        # identifying what times will be used for the ligolw_thinca jobs and
        # what files are needed for each. If doing time sliding there
        # will be some triggers read into multiple jobs
        cacheInspOuts = inspiral_outs.convert_to_lal_cache()
        if workflow.cp.has_option_tags("workflow-coincidence", 
                                       "maximum-extent", tags):
            max_extent = float( workflow.cp.get_opt_tags(
                              "workflow-coincidence", "maximum-extent", tags) )
        else:
            # hard-coded default value for extent of time in a single job
            max_extent = 3600
        logging.debug("Calling into cafe.")
        time_slide_table = lsctables.TimeSlideTable.get_table(\
                ligolw_utils.load_filename(tisiOutFile.storage_path,
                                 gz=tisiOutFile.storage_path.endswith(".gz"),
                                 contenthandler=ContentHandler,
                                 verbose=False))
        time_slide_table.sync_next_id()
        time_slide_dict = time_slide_table.as_dict()
        cafe_seglists, cafe_caches = ligolw_cafe.ligolw_cafe(cacheInspOuts,
            time_slide_dict.values(), extentlimit=max_extent, verbose=False)
        logging.debug("Done with cafe.")


        # Now loop over vetoes
        for category in veto_categories:
            logging.debug("Preparing %s %s" %(timeSlideTag,category))
            # FIXME: There are currently 3 names to say cumulative cat_3
            dqSegFile = segsList.find_output_with_tag(
                                               'CUMULATIVE_CAT_%d' %(category))
            if not len(dqSegFile) == 1:
                errMsg = "Did not find exactly 1 data quality file."
                raise ValueError(errMsg)
            dqVetoName = 'VETO_CAT%d_CUMULATIVE' %(category)
            # FIXME: Here we set the dqVetoName to be compatible with pipedown
            pipedownDQVetoName = 'CAT_%d_VETO' %(category)
            # FIXME: For pipedown must put the slide identifier first and
            # dqVetoName last.
            curr_thinca_job_tags = [timeSlideTag] + tags + [pipedownDQVetoName]

            logging.debug("Starting workflow")
            currLigolwThincaOuts, currOtherOuts = \
                  setup_snglveto_workflow_ligolw_thinca(workflow,
                               dqSegFile, tisiOutFile, dqVetoName,
                               cafe_seglists, cafe_caches, output_dir,
                               tags=curr_thinca_job_tags,
                               parallelize_split_input=parallelize_split_input,
                               insp_files_dict=inspiral_outs_dict)
            logging.debug("Done")
            ligolwThincaOuts.extend(currLigolwThincaOuts)
            for key, file_list in currOtherOuts.items():
                if other_outs.has_key(key):
                    other_outs[key].extend(currOtherOuts[key])
                else:
                    other_outs[key] = currOtherOuts[key]
    return ligolwThincaOuts, other_outs

def setup_snglveto_workflow_ligolw_thinca(workflow, dqSegFile, tisiOutFile,
                                          dqVetoName, cafe_seglists,
                                          cafe_caches, output_dir, tags=[],
                                          parallelize_split_input=False,
                                          insp_files_dict=None):
    '''
    This function is used to setup a single-stage ihope style coincidence stage
    of the workflow for a given dq category and for a given timeslide file
    using ligolw_sstinca (or compatible code!).

    Parameters
    -----------
    workflow : pycbc.workflow.core.Workflow
        The workflow instance that the coincidence jobs will be added to.
    dqSegFile : pycbc.workflow.core.SegFile
        The file that contains the data-quality segments to be applied to jobs
        setup by this call to this function.
    tisiOutFile : pycbc.workflow.core.File
        The file that contains the time-slides that will be performed in jobs
        setup by this call to this function. A file containing only one,
        zero-lag, time slide is still a valid "time slide" file.
    cafe_seglists : List of glue.segments.segmentlistdicts
        The times that each sstinca job will be valid for. Each entry represents
        a unique job
    cafe_caches : List of glue.lal.Cache objects
        The files that are inputs to each of the sstinca jobs.
    output_dir : path
        The directory in which coincidence output will be stored.
    tags : list of strings (optional, default = [])
        A list of the tagging strings that will be used for all jobs created
        by this call to the workflow. An example might be ['BNSINJECTIONS'] or
        ['NOINJECTIONANALYSIS']. This will be used in output names. At this
        stage a tag describing the time slide file, and a tag describing the
        data quality veto should be included.
    Returns
    --------
    ligolwThincaOuts : pycbc.workflow.core.FileList
        A list of the output files generated from ligolw_sstinca.
    ligolwAddOuts : pycbc.workflow.core.FileList
        A list of the output files generated from ligolw_add.

    '''
    cp = workflow.cp
    ifoString = workflow.ifo_string

    # This flag will add a clustering job after ligolw_thinca
    if workflow.cp.has_option_tags("workflow-coincidence",
                                             "coincidence-post-cluster", tags):
        coinc_post_cluster = True
    else:
        coinc_post_cluster = False


    # Set up jobs for ligolw_add and ligolw_thinca
    ligolwadd_job = LigolwAddExecutable(cp, 'llwadd', ifo=ifoString,
                                     out_dir=output_dir, tags=tags)
    ligolwthinca_job = LigolwSSthincaExecutable(cp, 'thinca', ifo=ifoString,
                                     out_dir=output_dir,
                                     dqVetoName=dqVetoName, tags=tags)
    if coinc_post_cluster:
        cluster_job = SQLInOutExecutable(cp, 'pycbccluster',
                                     ifo=ifoString, out_dir=output_dir,
                                     tags=tags)

    # Set up the nodes to do the coincidence analysis
    ligolwAddOuts = FileList([])
    ligolwThincaOuts = FileList([])
    ligolwThincaLikelihoodOuts = FileList([])
    ligolwClusterOuts = FileList([])
    for idx, cafe_cache in enumerate(cafe_caches):
        if not len(cafe_cache.objects):
            raise ValueError("One of the cache objects contains no files!")

        # Determine segments to accept coincidences.
        # If cache is not the first or last in the timeseries, check if the
        # two closes caches in the timeseries and see if their extent
        # match. If they match, they're adjacent and use the time where
        # they meet as a bound for accepting coincidences. If they're not
        # adjacent, then there is no bound for accepting coincidences.
        coincStart, coincEnd = None, None
        if idx and (cafe_cache.extent[0] == cafe_caches[idx-1].extent[1]):
            coincStart = cafe_cache.extent[0]
        if idx + 1 - len(cafe_caches) and \
                        (cafe_cache.extent[1] == cafe_caches[idx+1].extent[0]):
            coincEnd = cafe_cache.extent[1]
        coincSegment = (coincStart, coincEnd)

        # Need to create a list of the File(s) contained in the cache.
        # Assume that if we have partitioned input then if *one* job in the
        # partitioned input is an input then *all* jobs will be.
        if not parallelize_split_input:
            job_label = get_random_label()
            inputTrigFiles = FileList([])
            for object in cafe_cache.objects:
                inputTrigFiles.append(object.workflow_file)

            llw_files = inputTrigFiles + dqSegFile + [tisiOutFile]

            # Now we can create the nodes
            node = ligolwadd_job.create_node(cafe_cache.extent, llw_files)
            node.add_profile('pegasus', 'label', job_label)
            ligolwAddFile = node.output_files[0]
            ligolwAddOuts.append(ligolwAddFile)
            workflow.add_node(node)
            node = ligolwthinca_job.create_node(cafe_cache.extent,
                                                   coincSegment, ligolwAddFile)
            node.add_profile('pegasus', 'label', job_label)
            ligolwThincaOuts += \
                        node.output_files.find_output_without_tag('DIST_STATS')
            ligolwThincaLikelihoodOuts += \
                           node.output_files.find_output_with_tag('DIST_STATS')
            workflow.add_node(node)
            if coinc_post_cluster:
                node = cluster_job.create_node(cafe_cache.extent,
                                               ligolwThincaOuts[-1])
                node.add_profile('pegasus', 'label', job_label)
                ligolwClusterOuts += node.output_files
                workflow.add_node(node)
        else:
            for key in insp_files_dict.keys():
                job_label = get_random_label()
                curr_tags = ["JOB%d" %(key)]
                curr_list = insp_files_dict[key]
                inputTrigFiles = FileList([])
                for object in cafe_cache.objects:
                    inputTrigFiles.append(
                                     curr_list[object.workflow_file.thinca_index])

                llw_files = inputTrigFiles + dqSegFile + [tisiOutFile]

                # Now we can create the nodes
                node = ligolwadd_job.create_node(cafe_cache.extent, llw_files,
                                                 tags=curr_tags)
                node.add_profile('pegasus', 'label', job_label)
                ligolwAddFile = node.output_files[0]
                ligolwAddOuts.append(ligolwAddFile)
                workflow.add_node(node)
                if workflow.cp.has_option_tags("workflow-coincidence",
                                         "coincidence-write-likelihood", tags):
                    write_likelihood=True
                else:
                    write_likelihood=False
                node = ligolwthinca_job.create_node(cafe_cache.extent,
                                   coincSegment, ligolwAddFile, tags=curr_tags,
                                   write_likelihood=write_likelihood)
                node.add_profile('pegasus', 'label', job_label)
                ligolwThincaOuts += \
                       node.output_files.find_output_without_tag('DIST_STATS')
                ligolwThincaLikelihoodOuts += \
                          node.output_files.find_output_with_tag('DIST_STATS')
                workflow.add_node(node)
                if coinc_post_cluster:
                    node = cluster_job.create_node(cafe_cache.extent,
                                                   ligolwThincaOuts[-1])
                    node.add_profile('pegasus', 'label', job_label)
                    ligolwClusterOuts += node.output_files
                    workflow.add_node(node)

    other_returns = {}
    other_returns['LIGOLW_ADD'] = ligolwAddOuts
    other_returns['DIST_STATS'] = ligolwThincaLikelihoodOuts

    if coinc_post_cluster:
        main_return = ligolwClusterOuts
        other_returns['THINCA'] = ligolwThincaOuts
    else:
        main_return = ligolwThincaOuts

    return main_return, other_returns

class PyCBCBank2HDFExecutable(Executable):
    """ This converts xml tmpltbank to an hdf format
    """
    current_retention_level = Executable.CRITICAL
    def create_node(self, bank_file):
        node = Node(self)
        node.add_input_opt('--bank-file', bank_file)
        node.new_output_file_opt(bank_file.segment, '.hdf', '--output-file')
        return node

class PyCBCTrig2HDFExecutable(Executable):
    """ This converts xml triggers to an hdf format, grouped by template hash
    """
    current_retention_level = Executable.CRITICAL
    def create_node(self, trig_files, bank_file):
        node = Node(self)
        node.add_input_opt('--bank-file', bank_file)
        node.add_input_list_opt('--trigger-files', trig_files)
        node.new_output_file_opt(trig_files[0].segment, '.hdf', '--output-file')
        return node

class PyCBCFindCoincExecutable(Executable):
    """ Find coinc triggers using a folded interval method
    """
    current_retention_level = Executable.CRITICAL
    def create_node(self, trig_files, bank_file, veto_files, template_str, tags=[]):
        segs = trig_files.get_times_covered_by_files()
        seg = segments.segment(segs[0][0], segs[-1][1])
        node = Node(self)
        node.set_memory(10000)
        node.add_input_opt('--template-bank', bank_file)
        node.add_input_list_opt('--trigger-files', trig_files)
        if len(veto_files) != 0:
            node.add_input_list_opt('--veto-files', veto_files)
        node.add_opt('--template-fraction-range', template_str)
        node.new_output_file_opt(seg, '.hdf', '--output-file', tags=tags)
        return node

class PyCBCStatMapExecutable(Executable):
    """ Calculate FAP, IFAR, etc
    """
    current_retention_level = Executable.CRITICAL
    def create_node(self, coinc_files, external_background=None, tags=[]):
        segs = coinc_files.get_times_covered_by_files()
        seg = segments.segment(segs[0][0], segs[-1][1])

        node = Node(self)
        node.set_memory(5000)
        node.add_input_list_opt('--coinc-files', coinc_files)
        node.new_output_file_opt(seg, '.hdf', '--output-file', tags=tags)
        if external_background:
            node.add_input_opt('--external-background', external_background)
        return node
        
def get_subsections(cp, section_name):
    sections = cp.sections()   
    subsections = [sec.split('-')[1] for sec in sections if sec.startswith(section_name + '-')]   
    if len(subsections) > 0:
        return subsections
    else:
        return ['']
 
def make_segments_plot(workflow, seg_files, out_dir, tags=[]):
    make_analysis_dir(out_dir)
    node = Node(Executable(workflow.cp, 'plot_segments', ifos=workflow.ifos,
                    out_dir=out_dir, tags=tags))
    node.add_input_list_opt('--segment-files', seg_files)
    node.new_output_file_opt(workflow.segment, '.html', '--output-file')
    workflow += node
        
def merge_single_detector_hdf_files(workflow, trigger_files, out_dir, tags=[]):
    make_analysis_dir(out_dir)
    out = FileList()
    for ifo in workflow.ifos:
        node = Node(Executable(workflow.cp, 'hdf_trigger_merge', 
                        ifos=ifo, out_dir=out_dir, tags=tags))  
        node.add_input_list_opt('--trigger-files', trigger_files.find_output_with_ifo(ifo))
        node.new_output_file_opt(workflow.analysis_time, '.hdf', '--output-file')
        workflow += node
        out += node.output_files
    return out
        
def make_foreground_table(workflow, trig_file, bank_file, ftag, out_dir, tags=[]):
    make_analysis_dir(out_dir)
    node = Node(Executable(workflow.cp, 'page_foreground', ifos=workflow.ifos,
                    out_dir=out_dir, tags=tags))
    node.add_input_opt('--bank-file', bank_file)
    node.add_opt('--foreground-tag', ftag)
    node.add_input_opt('--trigger-file', trig_file)
    node.new_output_file_opt(bank_file.segment, '.html', '--output-file')
    workflow += node

def make_sensitivity_plot(workflow, inj_file, out_dir, tags=[]):
    make_analysis_dir(out_dir)
    
    for tag in get_subsections(workflow.cp, 'plot_sensitivity'):
        node = Node(Executable(workflow.cp, 'plot_sensitivity', ifos=workflow.ifos,
                    out_dir=out_dir, tags=[tag] + tags))
        node.add_input_opt('--injection-file', inj_file)
        node.new_output_file_opt(inj_file.segment, '.png', '--output-file')
        workflow += node

def make_coinc_snrchi_plot(workflow, inj_file, inj_trig, stat_file, trig_file, out_dir, tags=[]):
    make_analysis_dir(out_dir)
    
    for tag in get_subsections(workflow.cp, 'plot_coinc_snrchi'):
        node = Node(Executable(workflow.cp, 'plot_coinc_snrchi', ifos=inj_trig.ifo,
                    out_dir=out_dir, tags=[tag] + tags))
        node.add_input_opt('--found-injection-file', inj_file)
        node.add_input_opt('--single-injection-file', inj_trig)
        node.add_input_opt('--coinc-statistic-file', stat_file)
        node.add_input_opt('--single-trigger-file', trig_file)
        node.new_output_file_opt(inj_file.segment, '.png', '--output-file')
        workflow += node

def make_inj_table(workflow, inj_file, out_dir, tags=[]):
    make_analysis_dir(out_dir)
    node = Node(Executable(workflow.cp, 'page_injections', ifos=workflow.ifos,
                    out_dir=out_dir, tags=tags))
    node.add_input_opt('--injection-file', inj_file)
    node.new_output_file_opt(bank_file.segment, '.html', '--output-file')
    workflow += node   

def make_snrchi_plot(workflow, trig_files, veto_file, out_dir, tags=[]):
    make_analysis_dir(out_dir)
    
    for tag in get_subsections(workflow.cp, 'plot_snrchi'):
        for trig_file in trig_files:
            node = Node(Executable(workflow.cp, 'plot_snrchi',
                        ifos=trig_file.ifo, 
                        out_dir=out_dir, 
                        tags=[tag] + tags))

            node.set_memory(15000)
            node.add_input_opt('--trigger-file', trig_file)
            node.add_input_opt('--veto-file', veto_file)
            node.new_output_file_opt(trig_file.segment, '.png', '--output-file')
            workflow += node  

def make_foundmissed_plot(workflow, inj_file, inj_tag, out_dir, tags=[]):
    make_analysis_dir(out_dir)
    
    for tag in get_subsections(workflow.cp, 'plot_foundmissed'):
        node = Node(Executable(workflow.cp, 'plot_foundmissed', ifos=workflow.ifos,
                    out_dir=out_dir, tags=[tag] + tags))
        node.add_opt('--injection-tag', inj_tag)
        node.add_input_opt('--injection-file', inj_file)
        node.new_output_file_opt(inj_file.segment, '.html', '--output-file')
        workflow += node   
    
def make_snrifar_plot(workflow, bg_file, out_dir, tags=[]):
    make_analysis_dir(out_dir)
    node = Node(Executable(workflow.cp, 'plot_snrifar', ifos=workflow.ifos,
                out_dir=out_dir, tags=tags))
    node.add_input_opt('--trigger-file', bg_file)
    node.new_output_file_opt(bg_file.segment, '.png', '--output-file')
    workflow += node
 
class PyCBCHDFInjFindExecutable(Executable):
    """ Find injections in the hdf files output
    """
    current_retention_level = Executable.CRITICAL
    def create_node(self, inj_coinc_file, inj_xml_file, veto_file, tags=[]):
        node = Node(self)
        node.add_input_list_opt('--trigger-file', inj_coinc_file)
        node.add_input_list_opt('--injection-file', inj_xml_file)
        node.add_input_opt('--veto-file', veto_file)
        node.new_output_file_opt(inj_xml_file[0].segment, '.hdf', '--output-file', 
                                 tags=tags)
        return node

def find_injections_in_hdf_coinc(workflow, inj_coinc_file, inj_xml_file, 
                                 veto_file, out_dir, tags=[]):
    make_analysis_dir(out_dir)
    exe = PyCBCHDFInjFindExecutable(workflow.cp, 'hdfinjfind', 
                                    ifos=workflow.ifos, 
                                    out_dir=out_dir, tags=tags)
    node = exe.create_node(inj_coinc_file, inj_xml_file, veto_file, tags)
    workflow += node
    return node.output_files[0]     

def convert_bank_to_hdf(workflow, xmlbank, out_dir, tags=[]):
    """Return the template bank in hdf format
    """
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

def convert_trig_to_hdf(workflow, hdfbank, xml_trigger_files, out_dir, tags=[]):
    """Return the list of hdf5 trigger files outpus
    """
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

def setup_interval_coinc_inj(workflow, hdfbank, trig_files,
                           background_file, out_dir, tags=[]):
    """
    This function sets up exact match coincidence and background estimation
    using a folded interval technique.
    """
    make_analysis_dir(out_dir)
    logging.info('Setting up coincidence for injection')

    if len(hdfbank) > 1:
        raise ValueError('This coincidence method only supports a '
                         'pregenerated template bank')
    hdfbank = hdfbank[0]

    if len(workflow.ifos) > 2:
        raise ValueError('This coincidence method only supports two ifo searches')

    findcoinc_exe = PyCBCFindCoincExecutable(workflow.cp, 'coinc',
                                              ifos=workflow.ifos,
                                              tags=tags, out_dir=out_dir)

    combinecoinc_exe = PyCBCStatMapExecutable(workflow.cp, 'statmap',
                                              ifos=workflow.ifos,
                                              tags=tags, out_dir=out_dir)

    # Wall time knob and memory knob
    factor = int(workflow.cp.get_opt_tags('workflow-coincidence', 'parallelization-factor', tags))

    bg_files = FileList()
    for i in range(factor):
        group_str = '%s/%s' % (i, factor)
        coinc_node = findcoinc_exe.create_node(trig_files, hdfbank, [],
                                           group_str, tags=[str(i)])
        bg_files += coinc_node.output_files
        workflow.add_node(coinc_node)

    combine_node = combinecoinc_exe.create_node(bg_files,
                                           external_background=background_file[0])
    workflow.add_node(combine_node)

    logging.info('...leaving coincidence ')
    return combine_node.output_files[0]


def setup_interval_coinc(workflow, hdfbank, trig_files,
                         veto_files, out_dir, tags=[]):
    """
    This function sets up exact match coincidence and background estimation
    using a folded interval technique.
    """
    make_analysis_dir(out_dir)
    logging.info('Setting up coincidence')

    if len(hdfbank) > 1:
        raise ValueError('This coincidence method only supports a '
                         'pregenerated template bank')
    hdfbank = hdfbank[0]

    if len(workflow.ifos) > 2:
        raise ValueError('This coincidence method only supports two ifo searches')

    findcoinc_exe = PyCBCFindCoincExecutable(workflow.cp, 'coinc',
                                              ifos=workflow.ifos,
                                              tags=tags, out_dir=out_dir)

    combinecoinc_exe = PyCBCStatMapExecutable(workflow.cp, 'statmap',
                                              ifos=workflow.ifos,
                                              tags=tags, out_dir=out_dir)
                                              
    # Wall time knob and memory knob
    factor = int(workflow.cp.get_opt_tags('workflow-coincidence', 'parallelization-factor', tags))

    tags, veto_file_groups = veto_files.categorize_by_attr('tags')
    stat_files = FileList()
    for tag, veto_files in zip(tags, veto_file_groups):
        bg_files = FileList()
        for i in range(factor):
            group_str = '%s/%s' % (i, factor)
            coinc_node = findcoinc_exe.create_node(trig_files, hdfbank, veto_files,
                                                   group_str,
                                                   tags= tag + [str(i)])
            bg_files += coinc_node.output_files
            workflow.add_node(coinc_node)
             
        combine_node = combinecoinc_exe.create_node(bg_files, tags=tag)
        workflow.add_node(combine_node)
        stat_files += combine_node.output_files
        
        make_snrifar_plot(workflow, combine_node.output_files[0],
                          'plots/background', tags=tag)

    return stat_files
    logging.info('...leaving coincidence ')

