from __future__ import division

import os
from glue import segments
from pycbc.ahope.ahope_utils import *
from pycbc.ahope.jobsetup_utils import *
from pylal import ligolw_cafe, ligolw_tisi

def setup_coincidence_workflow(workflow, science_segs, segsDict, inspiral_outs,
                               output_dir, maxVetoCat=5):
    '''
    Setup coincidence stage of ahope workflow.
    FIXME: ADD MORE DOCUMENTATION
    '''
    logging.info('Entering coincidence setup module.')
    
    # Scope here for adding different options/methods here. For now we only
    # have the single_stage ihope method which consists of using ligolw_add
    # to create a large job for coincidence and then running ligolw_thinca
    # on that output.

    # If you want the ligolw_add outputs, call this function directly
    coinc_outs, _ = setup_coincidence_workflow_ligolw_thinca(workflow,
                     science_segs, segsDict, inspiral_outs, output_dir,
                     maxVetoCat=maxVetoCat)

    logging.info('Leaving coincidence setup module.')

    return coinc_outs

def setup_coincidence_workflow_ligolw_thinca(workflow, science_segs, segsDict,
                                     inspiral_outs, output_dir, maxVetoCat=5):
    """
    ADD DOCUMENTATION
    """
    ifoString = ''.join(science_segs.keys())
    veto_categories = range(1,maxVetoCat+1)

    # setup code for each veto_category
    # FIXME: Needs tag feature to avoid overwriting output!

    ligolwThincaOuts = AhopeFileList([])
    ligolwAddOuts = AhopeFileList([])

    for category in veto_categories:
        dqSegFile = segsDict[ifoString]['CUMULATIVE_CAT_%d' %(category)]
        dqVetoName = 'VETO_CAT%d_CUMULATIVE' %(category)
        currLigolwThincaOuts, currLigolwAddOuts = \
              setup_snglveto_workflow_ligolw_thinca(workflow, 
                                        science_segs, dqSegFile,
                                        dqVetoName, inspiral_outs, output_dir)
        ligolwAddOuts.extend(currLigolwAddOuts)
        ligolwThincaOuts.extend(currLigolwThincaOuts)
    return ligolwThincaOuts, ligolwAddOuts

def setup_snglveto_workflow_ligolw_thinca(workflow, science_segs, dqSegFile,
                                dqVetoName, inspiral_outs, output_dir):
    '''
    ADD DOCUMENTATION
    '''
    cp = workflow.cp
    ifoString = ''.join(science_segs.keys())

    # Get full analysis segment
    extents = [science_segs[ifo].extent() for ifo in science_segs.keys()]
    min, max = extents[0]
    for lo, hi in extents:
        if min > lo:
            min = lo
        if max < hi:
            max = hi
    fullSegment = segments.segment(min, max)

    # First we need to run ligolw_tisi to make the necessary time slide input
    # xml files
    tisiOutPath = "%s/TISI_ZEROLAG.xml.gz" %(output_dir)
    tisiOutUrl = urlparse.urlunparse(['file', 'localhost', tisiOutPath,
                                       None, None, None])
    tisiOutFile = AhopeFile(ifoString, 'TIMESLIDES_ZEROLAG', fullSegment,
                            file_url=tisiOutUrl)
    ligolw_tisi_call = ['ligolw_tisi',
        "-v",
        "-i", "H1=0:0:0",
        "-i", "L1=0:100:5",
        "-i", "V1=0:100:10",
        tisiOutFile.path]
    make_external_call(ligolw_tisi_call, outDir=os.path.join(output_dir,'logs'),
                           outBaseName='%s-ligolw_tisi-call' %('ZEROLAG') )

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
    ligolwadd_job = ligolwadd_exe.create_job(cp, ifoString, out_dir=output_dir)
    ligolwthinca_job = ligolwthinca_exe.create_job(cp, ifoString, 
                                     out_dir=output_dir, dqVetoName=dqVetoName)

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
