import os
import logging
from glue import segments
from pycbc.ahope.ahope_utils import *
from pycbc.ahope.configparserutils import get_opt_ifo

def setup_timeslides_workflow(workflow, science_segs, output_dir=None, tags=[],
                              timeSlideSectionName='ligolw_tisi'):
    '''
    Setup generation of time_slide input files in the ahope workflow. Currently used
    only with ligolw_tisi to generate files containing the list of slides to be
    performed in each time slide job.
    '''
    logging.info("Entering time slides setup module.")
    make_analysis_dir(output_dir)
    ifoList = science_segs.keys()
    ifoString = ''.join(ifoList)

    # Get full analysis segment for output file naming
    extents = [science_segs[ifo].extent() for ifo in science_segs.keys()]
    min, max = extents[0]
    for lo, hi in extents:
        if min > lo:
            min = lo
        if max < hi:
            max = hi
    fullSegment = segments.segment(min, max)

    timeSlideOuts = AhopeFileList([])

    # FIXME: Add ability to specify different exes

    # FIXME: Here I think I would prefer to setup a node, like normal, and then either
    # FIXME: add it to the workflow, *or* generate it at runtime.

    # Get all sections by looking in ini file
    timeSlideTags = [sec.split('-')[-1] for sec in workflow.cp.sections() \
                              if sec.startswith('tisi-')]

    # Make the timeSlideFiles
    for timeSlideTag in timeSlideTags:
        # First we need to run ligolw_tisi to make the necessary time slide
        # input xml files
        tisiOutFile = AhopeFile(ifoString, 'TIMESLIDES', fullSegment,
                                directory=output_dir, extension=".xml.gz",
                                tags=[timeSlideTag])
        ligolw_tisi_call = [workflow.cp.get('executables', 'tisi'), "-v"]
        # FIXME: I *really* want a new front end here so I don't need all this!
        subString = 'tisi-%s' %(timeSlideTag.lower())
        if workflow.cp.has_option('tisi', 'inspiral-num-slides'):
            ligolw_tisi_call.append("--inspiral-num-slides")
            ligolw_tisi_call.append(\
                    workflow.cp.get('tisi', 'inspiral-num-slides'))
        elif workflow.cp.has_option(subString, 'inspiral-num-slides'):
            ligolw_tisi_call.append("--inspiral-num-slides")
            ligolw_tisi_call.append(\
                    workflow.cp.get(subString, 'inspiral-num-slides'))
        else:
            for ifo in ifoList:
                ifoSlideStart = get_opt_ifo(workflow.cp, 'tisi',
                                        '%s-slide-start' %(ifo), timeSlideTag)
                ifoSlideEnd = get_opt_ifo(workflow.cp, 'tisi',
                                        '%s-slide-end' %(ifo), timeSlideTag)
                ifoSlideStep = get_opt_ifo(workflow.cp, 'tisi',
                                        '%s-slide-step' %(ifo), timeSlideTag)
                ligolw_tisi_call.append("-i")
                optionString = ':'.join([ifoSlideStart,ifoSlideEnd,ifoSlideStep])
                optionString = '%s=%s' %(ifo.upper(), optionString)
                ligolw_tisi_call.append(optionString)
        if workflow.cp.has_option('tisi', 'remove-zero-lag') or\
                   workflow.cp.has_option(subString, 'remove-zero-lag'):
            ligolw_tisi_call.append("--remove-zero-lag")
        ligolw_tisi_call.append(tisiOutFile.path)
        make_external_call(ligolw_tisi_call,
                            outDir=os.path.join(output_dir,'logs'),
                            outBaseName='%s-ligolw_tisi-call' %(timeSlideTag) )
        timeSlideOuts.append(tisiOutFile)

    return timeSlideOuts



