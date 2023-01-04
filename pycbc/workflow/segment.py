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

import os, shutil, itertools
import logging
from ligo import segments
from ligo.segments import utils as segmentsUtils
from pycbc.workflow.core import SegFile, make_analysis_dir
from pycbc.workflow.core import resolve_url

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
    server = "https://segments.ligo.org"
    if cp.has_option("workflow-segments", "segments-database-url"):
        server = cp.get("workflow-segments",
                                 "segments-database-url")

    if cp.has_option("workflow-segments", "segments-source"):
        source = cp.get("workflow-segments", "segments-source")
    else:
        source = "any"

    if source == "file":
        local_file_path = \
            resolve_url(cp.get("workflow-segments", option_name+"-file"))
        pfn = os.path.join(out_dir, os.path.basename(local_file_path))
        shutil.move(local_file_path, pfn)
        return SegFile.from_segment_xml(pfn)

    segs = {}
    for ifo in workflow.ifos:
        flag_str = cp.get_opt_tags("workflow-segments", option_name, [ifo])
        key = ifo + ':' + name
        
        if flag_str.upper() == "OFF":
            segs[key] = segments.segmentlist([])
        elif flag_str.upper() == "ON":
            all_seg = segments.segment([start, end])
            segs[key] = segments.segmentlist([all_seg])
        else:
            segs[key] = query_str(ifo, flag_str, start, end,
                                  source=source, server=server,
                                  veto_definer=veto_definer)
        logging.info("%s: got %s flags", ifo, option_name)

    return SegFile.from_segment_list_dict(name, segs,
                                          extension='.xml',
                                          valid_segment=workflow.analysis_time,
                                          directory=out_dir)


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
    offsrclist = commonsegs[tuple(commonsegs.keys())[0]]
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
    num_ifos = len(scisegs)
    while num_ifos >= min_ifos:
        # Consider all combinations for a given number of IFOs
        ifo_combos = itertools.combinations(scisegs, num_ifos)
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
                seg_lens = {ifos: abs(next(offsource[ifos].values())[0])
                            for ifos in offsource.keys()}
                best_comb = max(seg_lens.iterkeys(),
                                key=(lambda key: seg_lens[key]))
            else:
                best_comb = tuple(offsource.keys())[0]
            logging.info("No combination of %d IFOs with suitable science "
                         "segment.", num_ifos)
        else:
            # Identify best analysis segment
            if len(valid_combs) > 1:
                seg_lens = {ifos: abs(next(offsource[ifos].values())[0])
                            for ifos in valid_combs}
                best_comb = max(seg_lens.iterkeys(),
                                key=(lambda key: seg_lens[key]))
            else:
                best_comb = valid_combs[0]
            logging.info("Calculated science segments.")

            offsourceSegfile = os.path.join(out_dir, "offSourceSeg.txt")
            segmentsUtils.tosegwizard(open(offsourceSegfile, "w"),
                                      list(offsource[best_comb].values())[0])

            onsourceSegfile = os.path.join(out_dir, "onSourceSeg.txt")
            segmentsUtils.tosegwizard(open(onsourceSegfile, "w"),
                                      list(onsource[best_comb].values())[0])

            bufferleft = int(cp.get('workflow-exttrig_segments',
                                    'num-buffer-before'))
            bufferright = int(cp.get('workflow-exttrig_segments',
                                     'num-buffer-after'))
            onlen = onbefore + onafter
            bufferSegment = segments.segment(\
                    triggertime - onbefore - bufferleft * onlen,
                    triggertime + onafter + bufferright * onlen)
            bufferSegfile = os.path.join(out_dir, "bufferSeg.txt")
            segmentsUtils.tosegwizard(open(bufferSegfile, "w"),
                                      segments.segmentlist([bufferSegment]))

            return onsource[best_comb], offsource[best_comb]

        num_ifos -= 1

    logging.warning("No suitable science segments available.")
    try:
        return None, offsource[best_comb]
    except UnboundLocalError:
        return None, min_seg

def get_flag_segments_file(workflow, name, option_name, out_dir):
     """Get segments from option name syntax for each ifo for indivudal flags.

     Use syntax of configparser string to define the resulting segment_file
     e.x. option_name = +up_flag1,+up_flag2,+up_flag3,-down_flag1,-down_flag2
     Each ifo may have a different string and is stored separately in the file.
     Each flag is stored separately in the file.
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
     server = "https://segments.ligo.org"
     if cp.has_option("workflow-segments", "segments-database-url"):
         server = cp.get("workflow-segments", "segments-database-url")

     source = "any"
     if cp.has_option("workflow-segments", "segments-source"):
         source = cp.get("workflow-segments", "segments-source")
     if source == "file":
         local_file_path = \
             resolve_url(cp.get("workflow-segments", option_name+"-file"))
         pfn = os.path.join(out_dir, os.path.basename(local_file_path))
         shutil.move(local_file_path, pfn)
         return SegFile.from_segment_xml(pfn)

     segs = {}
     for ifo in workflow.ifos:
         if cp.has_option_tags("workflow-segments", option_name, [ifo]):
              flag_str = cp.get_opt_tags("workflow-segments", option_name, [ifo])
              flag_list = flag_str.split(',')
              for flag in flag_list:
                  flag_name = flag[1:]
                  key = flag_name
                  if len(key.split(':')) > 2:
                      key = ':'.join(key.split(':')[:2])
                  segs[key] = query_str(ifo, flag, start, end,
                                        source=source, server=server,
                                        veto_definer=veto_definer)
                  logging.info("%s: got %s segments", ifo, flag_name)
         else:
             logging.info("%s: no segments requested", ifo)

     return SegFile.from_segment_list_dict(name, segs,
                                           extension='.xml',
                                           valid_segment=workflow.analysis_time,
                                           directory=out_dir)
