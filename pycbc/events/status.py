# Copyright (C) 2018 Alex Nitz

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
""" Utitlities to query the status of gravitational-wave instruments
from public sources as well as dqsegdb.
"""

import numpy
import json
import urllib
from glue.segments import segmentlist, segment

def parse_veto_definer(veto_def_filename):
    """ Parse a veto definer file from the filename and return a dictionary 
    indexed by ifo and veto definer category level.
    
    Parameters
    ----------
    veto_def_filename: str
        The path to the veto definer file
    
    Returns:
        parsed_definition: dict   
            Returns a dictionary first indexed by ifo, then category level, and
            finally a list of veto definitions.
    """
    from glue.ligolw import ligolw, table, lsctables, utils as ligolw_utils
    from glue.ligolw.ligolw import LIGOLWContentHandler as h
    lsctables.use_in(h)
    
    indoc = ligolw_utils.load_filename(veto_def_filename, False, contenthandler=h)
    veto_table = table.get_table(indoc, 'veto_definer')
    
    ifo = veto_table.getColumnByName('ifo')
    name = veto_table.getColumnByName('name')
    version = numpy.array(veto_table.getColumnByName('version'))
    category = numpy.array(veto_table.getColumnByName('category'))
    start = numpy.array(veto_table.getColumnByName('start_time'))
    end = numpy.array(veto_table.getColumnByName('end_time'))
    start_pad = numpy.array(veto_table.getColumnByName('start_pad'))
    end_pad = numpy.array(veto_table.getColumnByName('end_pad'))
    
    data = {}
    for i in range(len(veto_table)):
        if ifo[i] not in data:
            data[ifo[i]] = {}
        
        # The veto-definer categories are weird! Hardware injections are stored
        # in "3" and numbers above that are bumped up by one (although not
        # often used any more). So we remap 3 to H and anything above 3 to
        # N-1. 2 and 1 correspond to 2 and 1 (YAY!)
        if category[i] > 3:
            curr_cat = "CAT_{}".format(category[i]-1)
        elif category[i] == 3:
            curr_cat = "CAT_H"
        else:
            curr_cat = "CAT_{}".format(category[i])

        if curr_cat not in data[ifo[i]]:
            data[ifo[i]][curr_cat] = []
            
        veto_info = {'name': name[i],
                     'version': version[i],
                     'start': start[i],
                     'end': end[i],
                     'start_pad': start_pad[i],
                     'end_pad': end_pad[i],
                     }
        data[ifo[i]][curr_cat].append(veto_info)
    return data

def query_flag(ifo, segment_name, start_time, end_time,
               source='any', server="segments.ligo.org",
               veto_definer=None):
    """Return the times where the flag is active 

    Parameters
    ----------
    ifo: string
        The interferometer to query (H1, L1).
    segment_name: string
        The status flag to query from LOSC.
    start_time: int 
        The starting gps time to begin querying from LOSC
    end_time: int 
        The end gps time of the query
    source: str, Optional
        Choice between "GWOSC" or "dqsegdb". If dqsegdb, the server option may
        also be given. The default is to try GWOSC first then try dqsegdb.
    server: str, Optional
        The server path. Only used with dqsegdb atm.
    veto_definer: str, Optional
        The path to a veto definer to define groups of flags which
        themselves define a set of segments.

    Returns
    ---------
    segments: glue.segments.segmentlist
        List of segments
    """
    if source == 'GWOSC' or source == 'any':
        ### Special cases as the LOSC convention is backwards from normal
        ### LIGO / Virgo operation!!!!
        if (('_HW_INJ' in segment_name and not 'NO' in segment_name) or 
           'VETO' in segment_name):
            data = query_flag(ifo, 'DATA', start_time, end_time)
             
            if '_HW_INJ' in segment_name:
                name = 'NO_' + segment_name    
            else:
                name = segment_name.replace('_VETO', '')

            negate = query_flag(ifo, name, start_time, end_time)
            return (data - negate).coalesce()

        duration = end_time - start_time
        url = 'https://www.gw-openscience.org/timeline/segments/json/O1/{}_{}/{}/{}/'
        url = url.format(ifo, segment_name, start_time, duration)
        
        try:
            flag_segments = json.loads(urllib.urlopen(url).read())['segments']
        except:
            msg = "Unable to find segments in GWOSC, check flag name or times"
            if source == 'any':
                return query_flag(ifo, segment_name, start_time, end_time,
                                  source='dqsegdb', server=server,
                                  veto_definer=veto_definer)
            else:
                raise ValueError(msg)
    elif source == 'dqsegdb':
        # Let's not hard require dqsegdb to be installed if we never get here.
        try:
            from dqsegdb.apicalls import dqsegdbQueryTimes as query
        except ImportError:
            raise ValueError("Could not query flag. DQSegdb is not installed" 
                             "to additionally try there")

        # The veto definer will allow the use of MACRO names
        # These directly correspond the name defined in the veto definer file.
        if veto_definer is not None:
            veto_def = parse_veto_definer(veto_definer)

        # We treat the veto definer name as if it were its own flag and
        # a process the flags in the veto definer
        flag_segments = segmentlist([])
        if veto_definer is not None and segment_name in veto_def[ifo]:
            for flag in veto_def[ifo][segment_name]:
                segs = query("https", server, ifo, flag['name'],
                             flag['version'], 'active',
                             start_time, end_time)[0]['active']
                for rseg in segs:
                    s, e = rseg[0] + flag['start_pad'], rseg[1] + flag['end_pad']
                    flag_segments.append(segment(s, e))
        else: # Standard case just query directly.
            try:
                name, version = segment_name.split(':')
                segs = query("https", server, ifo, name, version, 
                      'active', start_time, end_time)[0]['active']
                for rseg in segs:
                    s, e = rseg[0], rseg[1]
                    flag_segments.append(segment(s, e))
            except:
                raise ValueError("Could not query flag, check name (%s) or times",
                                 segment_name)

    else:
        err_msg = "source must be dgseqdb or GWOSC. "
        err_msg += "Got {}".format(source)
        raise ValueError(err_msg)

    return segmentlist(flag_segments).coalesce()

def query_cumulative_flags(ifo, segment_names, start_time, end_time, 
                           source='any', server="segments.ligo.org",
                           veto_definer=None):
    """Return the times where any flag is active 

    Parameters
    ----------
    ifo: string
        The interferometer to query (H1, L1).
    segment_name: list of strings
        The status flag to query from LOSC.
    start_time: int 
        The starting gps time to begin querying from LOSC
    end_time: int 
        The end gps time of the query
    source: str, Optional
        Choice between "GWOSC" or "dqsegdb". If dqsegdb, the server option may
        also be given. The default is to try GWOSC first then try dqsegdb.
    server: str, Optional
        The server path. Only used with dqsegdb atm.
    veto_definer: str, Optional
        The path to a veto definer to define groups of flags which
        themselves define a set of segments.

    Returns
    ---------
    segments: glue.segments.segmentlist
        List of segments
    """
    total_segs = segmentlist([])
    for flag_name in segment_names:
        segs = query_flag(ifo, flag_name, start_time, end_time,                                       
                          source=source, server=server,
                          veto_definer=veto_definer)
        total_segs = (total_segs + segs).coalesce()
    return total_segs

def query_combined_flags(ifo, up_flags, start_time, end_time, down_flags=None,               
                        source='any', server="segments.ligo.org",
                        veto_definer=None):
    """Return the times where any "up" flag is active minus "down" flag times

    Parameters
    ----------
    ifo: string
        The interferometer to query (H1, L1).
    up flags: list of strings
        The status flag to query from LOSC.
    start_time: int 
        The starting gps time to begin querying from LOSC
    end_time: int 
        The end gps time of the query
    down_flags: list of strings
        Flags which indicate times to subtract from the combined segments
    source: str, Optional
        Choice between "GWOSC" or "dqsegdb". If dqsegdb, the server option may
        also be given. The default is to try GWOSC first then try dqsegdb.
    server: str, Optional
        The server path. Only used with dqsegdb atm.
    veto_definer: str, Optional
        The path to a veto definer to define groups of flags which
        themselves define a set of segments.

    Returns
    ---------
    segments: glue.segments.segmentlist
        List of segments
    """
    segments = query_cumulative_flags(ifo, up_flags, start_time, end_time,
                                      source=source, server=server,
                                      veto_definer=veto_definer)
    down_flags = [] if down_flags is None else down_flags
    for flag_name in down_flags:
        mseg = query_flag(ifo, flag_name, start_time, end_time,
                          source=source, server=server,
                          veto_definer=veto_definer)
        segments = (segments - mseg).coalesce()
    return segments
