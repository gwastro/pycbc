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


from __future__ import division

import re
import logging
from glue import segments
from glue.ligolw import lsctables,ligolw
from glue.ligolw import utils as ligolw_utils
from pycbc.workflow.core import FileList, make_analysis_dir, Executable, Node, File
from pycbc.workflow.jobsetup import LigolwAddExecutable, LigolwSSthincaExecutable, SQLInOutExecutable

class ContentHandler(ligolw.LIGOLWContentHandler):
    pass

lsctables.use_in(ContentHandler)

def setup_coincidence_workflow(workflow, segsList, timeSlideFiles,
                               inspiral_outs, output_dir, veto_cats=[2,3,4],
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
    veto_cats : list of ints (optional, default = [2,3,4])
        Veto categories that will be applied in the coincidence jobs. If this
        takes the default value the code will run data quality at cumulative 
        categories 2, 3 and 4. Note that if we change the flag definitions to
        be non-cumulative then this option will need to be revisited.
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
                     output_dir, veto_cats=veto_cats, tags=tags,
                     timeSlideTags=timeSlideTags,
                     parallelize_split_input=parallelize_split_input)
    else:
        errMsg = "Coincidence method not recognized. Must be one of "
        errMsg += "WORKFLOW_DISCRETE_SLIDES (currently only one option)."
        raise ValueError(errMsg)

    logging.info('Leaving coincidence setup module.')

    return coinc_outs, other_outs

def setup_coincidence_workflow_ligolw_thinca(
        workflow, segsList, timeSlideFiles, inspiral_outs, output_dir,
        veto_cats=[2,3,4], tags=[], timeSlideTags=None,
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
    veto_cats : list of ints (optional, default = [2,3,4])
        Veto categories that will be applied in the coincidence jobs. If this
        takes the default value the code will run data quality at cumulative 
        categories 2, 3 and 4. Note that if we change the flag definitions to
        be non-cumulative then this option will need to be revisited.
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
    from pylal import ligolw_cafe

    logging.debug("Entering coincidence module.")
    cp = workflow.cp
    ifoString = workflow.ifo_string

    # setup code for each veto_category

    coinc_outs = FileList([])
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

        # Take the combined seglist file
        dqSegFile=segsList.find_output_with_tag('COMBINED_CUMULATIVE_SEGMENTS')
        if not len(dqSegFile) == 1:
            errMsg = "Did not find exactly 1 data quality file."
            print len(dqSegFile), dqSegFile
            raise ValueError(errMsg)
        dqSegFile=dqSegFile[0]

        # Set up llwadd job
        llwadd_tags = [timeSlideTag] + tags 
        ligolwadd_job = LigolwAddExecutable(cp, 'llwadd', ifo=ifoString,
                                          out_dir=output_dir, tags=llwadd_tags)
        ligolwAddOuts = FileList([])

        # Go global setup at each category
        # This flag will add a clustering job after ligolw_thinca
        if workflow.cp.has_option_tags("workflow-coincidence",
                                      "coincidence-post-cluster", llwadd_tags):
            coinc_post_cluster = True
        else:
            coinc_post_cluster = False

        # Go global setup at each category
        ligolwthinca_job = {}
        cluster_job = {}
        thinca_tags = {}
        for category in veto_cats:
            logging.debug("Preparing %s %s" %(timeSlideTag,category))
            dqVetoName = 'VETO_CAT%d_CUMULATIVE' %(category)
            # FIXME: Should we resolve this now?
            # FIXME: Here we set the dqVetoName to be compatible with pipedown
            #        For pipedown must put the slide identifier first and
            #        dqVetoName last.
            pipedownDQVetoName = 'CAT_%d_VETO' %(category)
            curr_thinca_job_tags = [timeSlideTag] + tags + [pipedownDQVetoName]
            thinca_tags[category]=curr_thinca_job_tags
            # Set up jobs for ligolw_thinca
            ligolwthinca_job[category] = LigolwSSthincaExecutable(cp, 'thinca',
                                             ifo=ifoString, out_dir=output_dir,
                                             dqVetoName=dqVetoName,
                                             tags=curr_thinca_job_tags)
            if coinc_post_cluster:
                cluster_job[category] = SQLInOutExecutable(cp, 'pycbccluster',
                                             ifo=ifoString, out_dir=output_dir,
                                             tags=curr_thinca_job_tags)
        
        for idx, cafe_cache in enumerate(cafe_caches):
            ligolwAddOuts = FileList([])
            ligolwThincaOuts = FileList([])
            ligolwThincaLikelihoodOuts = FileList([])
            ligolwClusterOuts = FileList([])

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
                inputTrigFiles = FileList([])
                for object in cafe_cache.objects:
                    inputTrigFiles.append(object.workflow_file)
        
                llw_files = inputTrigFiles + [dqSegFile] + [tisiOutFile]
        
                # Now we can create the nodes
                node = ligolwadd_job.create_node(cafe_cache.extent, llw_files)
                ligolwAddFile = node.output_files[0]
                ligolwAddOuts.append(ligolwAddFile)
                workflow.add_node(node)
                for category in veto_cats:
                    node = ligolwthinca_job[category].create_node(\
                                cafe_cache.extent, coincSegment, ligolwAddFile)
                    ligolwThincaOuts += \
                        node.output_files.find_output_without_tag('DIST_STATS')
                    ligolwThincaLikelihoodOuts += \
                           node.output_files.find_output_with_tag('DIST_STATS')
                    workflow.add_node(node)
                    if coinc_post_cluster:
                        node = cluster_job[category].create_node(\
                                       cafe_cache.extent, ligolwThincaOuts[-1])
                        ligolwClusterOuts += node.output_files
                        workflow.add_node(node)
            else:
                for key in inspiral_outs_dict.keys():
                    curr_tags = ["JOB%d" %(key)]
                    curr_list = inspiral_outs_dict[key]
                    inputTrigFiles = FileList([])
                    for object in cafe_cache.objects:
                        inputTrigFiles.append(
                                  curr_list[object.workflow_file.thinca_index])
        
                    llw_files = inputTrigFiles + [dqSegFile] + [tisiOutFile]

                    # Now we can create the nodes
                    node = ligolwadd_job.create_node(cafe_cache.extent,
                                                     llw_files, tags=curr_tags)
                    ligolwAddFile = node.output_files[0]
                    ligolwAddOuts.append(ligolwAddFile)
                    workflow.add_node(node)
                    if workflow.cp.has_option_tags("workflow-coincidence",
                          "coincidence-write-likelihood",curr_thinca_job_tags):
                        write_likelihood=True
                    else:
                        write_likelihood=False
                    for category in veto_cats:
                        node = ligolwthinca_job[category].create_node(\
                             cafe_cache.extent, coincSegment, ligolwAddFile,
                             tags=curr_tags, write_likelihood=write_likelihood)
                        ligolwThincaOuts += \
                               node.output_files.find_output_without_tag(\
                                                                  'DIST_STATS')
                        ligolwThincaLikelihoodOuts += \
                              node.output_files.find_output_with_tag(\
                                                                  'DIST_STATS')
                        workflow.add_node(node)
                        if coinc_post_cluster:
                            node = cluster_job[category].create_node(\
                                       cafe_cache.extent, ligolwThincaOuts[-1])
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
        
            logging.debug("Done")
            coinc_outs.extend(main_return)
            for key, file_list in other_returns.items():
                if other_outs.has_key(key):
                    other_outs[key].extend(other_returns[key])
                else:
                    other_outs[key] = other_returns[key]
    return coinc_outs, other_outs


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
        node.new_output_file_opt(trig_files[0].segment, '.hdf',
                                 '--output-file', use_tmp_subdirs=True)
        return node

class PyCBCFindCoincExecutable(Executable):
    """ Find coinc triggers using a folded interval method
    """
    current_retention_level = Executable.CRITICAL
    def create_node(self, trig_files, bank_file, veto_file, veto_name, template_str, tags=[]):
        segs = trig_files.get_times_covered_by_files()
        seg = segments.segment(segs[0][0], segs[-1][1])
        node = Node(self)
        node.set_memory(10000)
        node.add_input_opt('--template-bank', bank_file)
        node.add_input_list_opt('--trigger-files', trig_files)
        if veto_file is not None:
            node.add_input_opt('--veto-files', veto_file)
            node.add_opt('--segment-name', veto_name)
        node.add_opt('--template-fraction-range', template_str)
        node.new_output_file_opt(seg, '.hdf', '--output-file', tags=tags)
        return node

class PyCBCStatMapExecutable(Executable):
    """ Calculate FAP, IFAR, etc
    """
    current_retention_level = Executable.CRITICAL
    def create_node(self, coinc_files, tags=[]):
        segs = coinc_files.get_times_covered_by_files()
        seg = segments.segment(segs[0][0], segs[-1][1])

        node = Node(self)
        node.set_memory(5000)
        node.add_input_list_opt('--coinc-files', coinc_files)
        node.new_output_file_opt(seg, '.hdf', '--output-file', tags=tags)
        return node
        
class PyCBCStatMapInjExecutable(Executable):
    """ Calculate FAP, IFAR, etc
    """
    current_retention_level = Executable.CRITICAL
    def create_node(self, zerolag, full_data, injfull, fullinj, tags=[]):
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
    """ Find injections in the hdf files output
    """
    current_retention_level = Executable.CRITICAL
    def create_node(self, inj_coinc_file, inj_xml_file, veto_file, veto_name, tags=[]):
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
    """ Distribute coinc files amoung different background bins """
    current_retention_level = Executable.CRITICAL
    def create_node(self, coinc_files, bank_file,  background_bins, tags=[]):
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
    current_retention_level = Executable.CRITICAL
    def create_node(self, stat_files, tags=[]):
        node = Node(self)
        node.add_input_list_opt('--statmap-files', stat_files)
        node.new_output_file_opt(stat_files[0].segment, '.hdf', '--output-file', tags=tags)
        return node
 
class MergeExecutable(Executable):
    current_retention_level = Executable.CRITICAL

class CensorForeground(Executable):
    current_retention_level = Executable.CRITICAL

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

def merge_single_detector_hdf_files(workflow, bank_file, trigger_files, out_dir, tags=[]):
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

def find_injections_in_hdf_coinc(workflow, inj_coinc_file, inj_xml_file, 
                                 veto_file, veto_name, out_dir, tags=[]):
    make_analysis_dir(out_dir)
    exe = PyCBCHDFInjFindExecutable(workflow.cp, 'hdfinjfind', 
                                    ifos=workflow.ifos, 
                                    out_dir=out_dir, tags=tags)
    node = exe.create_node(inj_coinc_file, inj_xml_file, veto_file, veto_name, tags)
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
    bins_node = bins_exe.create_node(coinc_files, bank_file, background_bins)
    workflow += bins_node
    
    stat_files = FileList([])
    for i, coinc_file in enumerate(bins_node.output_files):
        statnode = statmap_exe.create_node(FileList([coinc_file]), tags=tags + ['BIN_%s' % i])
        workflow += statnode
        stat_files.append(statnode.output_files[0])
        stat_files[i].bin_name = bins_node.names[i]
    
    cstat_node = cstat_exe.create_node(stat_files, tags=tags)
    workflow += cstat_node
    
    return cstat_node.output_files[0], stat_files
    
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
    
    
    for inj_type in ['injinj', 'injfull', 'fullinj']:          
        bins_node = bins_exe.create_node(FileList(coinc_files[inj_type]), bank_file, background_bins, tags=tags + [inj_type])
        workflow += bins_node
        coinc_files[inj_type] = bins_node.output_files
    
    stat_files = FileList([])
    for i in range(len(background_bins)):
        statnode = statmap_exe.create_node(FileList([coinc_files['injinj'][i]]), FileList([background_file[i]]), 
                                     FileList([coinc_files['injfull'][i]]), FileList([coinc_files['fullinj'][i]]), 
                                     tags=tags + ['BIN_%s' % i])
        workflow += statnode
        stat_files.append(statnode.output_files[0])
        
    cstat_node = cstat_exe.create_node(stat_files, tags=tags)
    workflow += cstat_node
    
    return cstat_node.output_files[0]

def setup_interval_coinc_inj(workflow, hdfbank, full_data_trig_files, inj_trig_files,
                           background_file, veto_file, veto_name, out_dir, tags=[]):
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
                                           veto_file, veto_name,
                                           group_str, tags=([str(i)]))
            bg_files[ctag] += coinc_node.output_files
            workflow.add_node(coinc_node)

    return setup_statmap_inj(workflow, bg_files, background_file, hdfbank, out_dir, tags=tags)
    
def setup_interval_coinc(workflow, hdfbank, trig_files,
                         veto_files, veto_names, out_dir, tags=[]):
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
                                         
    # Wall time knob and memory knob
    factor = int(workflow.cp.get_opt_tags('workflow-coincidence', 'parallelization-factor', tags))

    stat_files = []
    for veto_file, veto_name in zip(veto_files, veto_names):
        bg_files = FileList()
        for i in range(factor):
            group_str = '%s/%s' % (i, factor)
            coinc_node = findcoinc_exe.create_node(trig_files, hdfbank, 
                                                   veto_file, veto_name,
                                                   group_str,
                                                   tags= [veto_name, str(i)])
            bg_files += coinc_node.output_files
            workflow.add_node(coinc_node)
             
        stat_files += [setup_statmap(workflow, bg_files, hdfbank, out_dir, tags=tags + [veto_name])]

    return stat_files
    logging.info('...leaving coincidence ')
