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
This module is responsible for setting up the coincidence stage of ahope
workflows. For details about this module and its capabilities see here:
https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/NOTYETCREATED.html
"""


from __future__ import division

import os
import os.path
from glue import segments
from pycbc.ahope.ahope_utils import *
from pycbc.ahope.jobsetup_utils import *
from pylal import ligolw_cafe, ligolw_tisi

def setup_coincidence_workflow(workflow, segsList, timeSlideFiles,
                               inspiral_outs, output_dir, maxVetoCat=5,
                               tags=[], timeSlideTags=None):
    '''
    This function aims to be the gateway for setting up a set of coincidence
    jobs in an ahope workflow. The goal is that this function can support a
    number of different ways/codes that could be used for doing this.
    For now it only supports ligolw_sstinca.

    Parameters
    -----------
    Workflow : ahope.Workflow
        The ahope workflow instance that the coincidence jobs will be added to.
    segsList : ahope.AhopeFileList
        The list of files returned by ahope's segment module that contains
        pointers to all the segment files generated in the workflow. If the
        coincidence code will be applying the data quality vetoes, then this
        will be used to ensure that the codes get the necessary input to do
        this.
    timeSlideFiles : ahope.AhopeFileList
        An AhopeFileList of the timeSlide input files that are needed to
        determine what time sliding needs to be done if the coincidence code
        will be running time slides to facilitate background computations later
        in the workflow.
    inspiral_outs : ahope.AhopeFileList
        An AhopeFileList of the matched-filter module output that is used as
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
    coinc_outs : ahope.AhopeFileList
        A list of the *final* outputs of the coincident stage. This *does not*
        include any intermediate products produced within the workflow. If you
        require access to intermediate products call the various sub-functions
        in this module directly.
    '''
    logging.info('Entering coincidence setup module.')
    make_analysis_dir(output_dir)

    # Parse for options in .ini file
    coincidenceMethod = workflow.cp.get_opt_tags("ahope-coincidence",
                                        "coincidence-method", tags)
    
    # Scope here for adding different options/methods here. For now we only
    # have the single_stage ihope method which consists of using ligolw_add
    # to create a large job for coincidence and then running ligolw_thinca
    # on that output.
    if coincidenceMethod == "WORKFLOW_DISCRETE_SLIDES":
        # If you want the ligolw_add outputs, call this function directly
        coinc_outs, _ = setup_coincidence_workflow_ligolw_thinca(workflow,
                     segsList, timeSlideFiles, inspiral_outs,
                     output_dir, maxVetoCat=maxVetoCat, tags=tags,
                     timeSlideTags=timeSlideTags)
    else:
        errMsg = "Coincidence method not recognized. Must be one of "
        errMsg += "WORKFLOW_DISCRETE_SLIDES (currently only one option)."
        raise ValueError(errMsg)

    logging.info('Leaving coincidence setup module.')

    return coinc_outs

def setup_coincidence_workflow_ligolw_thinca(workflow, segsList,
                                             timeSlideFiles, inspiral_outs, 
                                             output_dir, maxVetoCat=5, tags=[],
                                             timeSlideTags=None):
    """
    This function is used to setup a single-stage ihope style coincidence stage
    of the workflow using ligolw_sstinca (or compatible code!).

    Parameters
    -----------
    Workflow : ahope.Workflow
        The ahope workflow instance that the coincidence jobs will be added to.
    segsList : ahope.AhopeFileList
        The list of files returned by ahope's segment module that contains
        pointers to all the segment files generated in the workflow. If the
        coincidence code will be applying the data quality vetoes, then this
        will be used to ensure that the codes get the necessary input to do
        this.
    timeSlideFiles : ahope.AhopeFileList
        An AhopeFileList of the timeSlide input files that are needed to
        determine what time sliding needs to be done. One of the timeSlideFiles
        will normally be "zero-lag only", the others containing time slides
        used to facilitate background computations later
        in the workflow.
    inspiral_outs : ahope.AhopeFileList
        An AhopeFileList of the matched-filter module output that is used as
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
    ligolwThincaOuts : ahope.AhopeFileList
        A list of the output files generated from ligolw_sstinca.
    ligolwAddOuts : ahope.AhopeFileList
        A list of the output files generated from ligolw_add.
    """
    veto_categories = range(1,maxVetoCat+1)

    # setup code for each veto_category

    ligolwThincaOuts = AhopeFileList([])
    ligolwAddOuts = AhopeFileList([])

    if not timeSlideTags:
        # Get all sections by looking in ini file, use all time slide files.
        timeSlideTags = [(sec.split('-')[-1]).upper() \
                  for sec in workflow.cp.sections() if sec.startswith('tisi-')]

    for timeSlideTag in timeSlideTags:
        # Get the time slide file from the inputs
        tisiOutFile = timeSlideFiles.find_output_with_tag(timeSlideTag)
        if not len(tisiOutFile) == 1:
            errMsg = "If you are seeing this, something batshit is going on!"
            if len(tisiOutFile) == 0:
                errMsg = "No time slide files found matching %s." %(tag)
            if len(tisiOutFile) > 1:
                errMsg = "More than one time slide files match %s." %(tag)
            raise ValueError(errMsg)
        tisiOutFile = tisiOutFile[0]

        # Now loop over vetoes
        for category in veto_categories:
            # FIXME: There are currently 3 names to say cumulative cat_3
            dqSegFile = segsList.find_output_with_tag(\
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

            currLigolwThincaOuts, currLigolwAddOuts = \
                  setup_snglveto_workflow_ligolw_thinca(workflow, 
                                        dqSegFile, tisiOutFile,
                                        dqVetoName, inspiral_outs, output_dir,
                                        tags=curr_thinca_job_tags)
            ligolwAddOuts.extend(currLigolwAddOuts)
            ligolwThincaOuts.extend(currLigolwThincaOuts)
    return ligolwThincaOuts, ligolwAddOuts

def setup_snglveto_workflow_ligolw_thinca(workflow, dqSegFile,
                                          tisiOutFile, dqVetoName,
                                          inspiral_outs, output_dir, tags=[]):
    '''
    This function is used to setup a single-stage ihope style coincidence stage
    of the workflow for a given dq category and for a given timeslide file
    using ligolw_sstinca (or compatible code!).

    Parameters
    -----------
    Workflow : ahope.Workflow
        The ahope workflow instance that the coincidence jobs will be added to.
    dqSegFile : ahope.AhopeSegFile
        The file that contains the data-quality segments to be applied to jobs
        setup by this call to this function.
    tisiOutFile : ahope.AhopeFile
        The file that contains the time-slides that will be performed in jobs
        setup by this call to this function. A file containing only one,
        zero-lag, time slide is still a valid "time slide" file.
    inspiral_outs : ahope.AhopeFileList
        An AhopeFileList of the matched-filter module output that is used as
        input to the coincidence codes running at this stage.
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
    ligolwThincaOuts : ahope.AhopeFileList
        A list of the output files generated from ligolw_sstinca.
    ligolwAddOuts : ahope.AhopeFileList
        A list of the output files generated from ligolw_add.

    '''
    cp = workflow.cp
    ifoString = workflow.ifoString

    # Next we run ligolw_cafe. This is responsible for
    # identifying what times will be used for the ligolw_thinca jobs and what
    # files are needed for each. If doing time sliding there
    # will be some triggers read into multiple jobs
    cacheInspOuts = inspiral_outs.convert_to_lal_cache()

    cafe_seglists, cafe_caches = ligolw_cafe.ligolw_cafe(cacheInspOuts,
        ligolw_tisi.load_time_slides(tisiOutFile.path,
            gz = tisiOutFile.path.endswith(".gz")).values(),
        extentlimit = 10000, verbose=False)

    # FIXME: This is currently hardcoded to ligolw_thinca and ligolw_add. This
    # is deliberate as different coincidence methods will not just be swapping
    # codes in and out. However, if this needs to be overwritten it can be!

    # Set up jobs for ligolw_add and ligolw_thinca
    ligolwadd_exe = LigolwAddExec('llwadd')
    ligolwthinca_exe = LigolwSSthincaExec('thinca')
    ligolwadd_job = ligolwadd_exe.create_job(cp, ifoString, 
                                     out_dir=output_dir, tags=tags)
    ligolwthinca_job = ligolwthinca_exe.create_job(cp, ifoString, 
                                     out_dir=output_dir, 
                                     dqVetoName=dqVetoName, tags=tags)

    # Set up the nodes to do the coincidence analysis
    ligolwAddOuts = AhopeFileList([])
    ligolwThincaOuts = AhopeFileList([])
    for idx, cafe_cache in enumerate(cafe_caches):
        if not len(cafe_cache.objects):
            raise ValueError("One of the cache objects contains no files!")

        # Need to create a list of the AhopeFile contained in the cache.
        # Assume that if we have partitioned input then if *one* job in the
        # partitioned input is an input then *all* jobs will be.
        inputTrigFiles = AhopeFileList([])
        for object in cafe_cache.objects:
            inputTrigFiles.append(object.ahope_file)
 
        llw_files = inputTrigFiles + dqSegFile + [tisiOutFile]
        # Now we can create the nodes
        node = ligolwadd_job.create_node(cafe_cache.extent, llw_files)
        ligolwAddFile = node.output_files[0]
        ligolwAddOuts.append(ligolwAddFile)
        workflow.add_node(node)
        node = ligolwthinca_job.create_node(cafe_cache.extent, ligolwAddFile)
        ligolwThincaOuts += node.output_files
        workflow.add_node(node)

    return ligolwThincaOuts, ligolwAddOuts
