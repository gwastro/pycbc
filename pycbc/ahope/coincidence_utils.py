from __future__ import division

import os
import os.path
from glue import segments
from pycbc.ahope.ahope_utils import *
from pycbc.ahope.jobsetup_utils import *
from pylal import ligolw_cafe, ligolw_tisi

def setup_coincidence_workflow(workflow, science_segs, segsDict, timeSlideFiles,
                               inspiral_outs, output_dir, maxVetoCat=5,
                               tags=[], timeSlideTags=None):
    '''
    Setup coincidence stage of ahope workflow.
    FIXME: ADD MORE DOCUMENTATION
    '''
    logging.info('Entering coincidence setup module.')
    
    if not os.path.isabs(output_dir):
        output_dir = os.path.join(os.getcwd(), output_dir) 
    
    # Scope here for adding different options/methods here. For now we only
    # have the single_stage ihope method which consists of using ligolw_add
    # to create a large job for coincidence and then running ligolw_thinca
    # on that output.

    # If you want the ligolw_add outputs, call this function directly
    coinc_outs, _ = setup_coincidence_workflow_ligolw_thinca(workflow,
                     science_segs, segsDict, timeSlideFiles, inspiral_outs,
                     output_dir, maxVetoCat=maxVetoCat, tags=tags,
                     timeSlideTags=timeSlideTags)

    logging.info('Leaving coincidence setup module.')

    return coinc_outs

def setup_coincidence_workflow_ligolw_thinca(workflow, science_segs, segsDict,
                                             timeSlideFiles, inspiral_outs, 
                                             output_dir, tags=[], maxVetoCat=5,
                                             timeSlideTags=None):
    """
    ADD DOCUMENTATION
    """
    ifoList = science_segs.keys()
    ifoString = ''.join(ifoList)
    veto_categories = range(1,maxVetoCat+1)

    # setup code for each veto_category

    ligolwThincaOuts = AhopeFileList([])
    ligolwAddOuts = AhopeFileList([])

    if not timeSlideTags:
        # Get all sections by looking in ini file, use all time slide files.
        timeSlideTags = [sec.split('-')[-1] for sec in workflow.cp.sections() \
                                  if sec.startswith('tisi-')]

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
            dqSegFile = segsDict[ifoString]['CUMULATIVE_CAT_%d' %(category)]
            dqVetoName = 'VETO_CAT%d_CUMULATIVE' %(category)
            # FIXME: Here we set the dqVetoName to be compatible with pipedown
            # FIXME: Yes pipedown has a trailing "_", yes I know its stoopid
            pipedownDQVetoName = 'CAT_%d_VETO_' %(category)
            # FIXME: For pipedown must put the slide identifier first and
            # dqVetoName last.
            curr_thinca_job_tags = [timeSlideTag] + tags + [pipedownDQVetoName]

            currLigolwThincaOuts, currLigolwAddOuts = \
                  setup_snglveto_workflow_ligolw_thinca(workflow, 
                                        science_segs, dqSegFile, tisiOutFile,
                                        dqVetoName, inspiral_outs, output_dir,
                                        tags=curr_thinca_job_tags)
            ligolwAddOuts.extend(currLigolwAddOuts)
            ligolwThincaOuts.extend(currLigolwThincaOuts)
    return ligolwThincaOuts, ligolwAddOuts

def setup_snglveto_workflow_ligolw_thinca(workflow, science_segs, dqSegFile,
                                          tisiOutFile, dqVetoName, inspiral_outs,
                                          output_dir, tags=[]):
    '''
    ADD DOCUMENTATION
    '''
    cp = workflow.cp
    ifoString = ''.join(science_segs.keys())

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
        for file in inspiral_outs:
            currCache = file.cache_entries[0]
            if currCache in cafe_cache.objects:
                inputTrigFiles.append(file)
 
        # Now we can create the nodes
        node = ligolwadd_job.create_node(cafe_cache.extent, inputTrigFiles,
                                  tisiOutFile, dqSegFile=dqSegFile)
        ligolwAddFile = node.output_files
        ligolwAddOuts += ligolwAddFile
        workflow.add_node(node)
        node = ligolwthinca_job.create_node(cafe_cache.extent, ligolwAddFile)
        ligolwThincaOuts += node.output_files
        workflow.add_node(node)

    return ligolwThincaOuts, ligolwAddOuts
