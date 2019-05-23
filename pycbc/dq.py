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
""" Utilities to query archival instrument status information of
gravitational-wave instruments from public sources as well as dqsegdb.
"""

import json
import numpy
from astropy.utils.data import download_file
from ligo.segments import segmentlist, segment
from pycbc.frame.losc import get_run

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
    from glue.ligolw import table, lsctables, utils as ligolw_utils
    from glue.ligolw.ligolw import LIGOLWContentHandler as h
    lsctables.use_in(h)

    indoc = ligolw_utils.load_filename(veto_def_filename, False,
                                       contenthandler=h)
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


GWOSC_URL = 'https://www.gw-openscience.org/timeline/segments/json/{}/{}_{}/{}/{}/'


def query_flag(ifo, name, start_time, end_time,
               source='any', server="segments.ligo.org",
               veto_definer=None, cache=False):
    """Return the times where the flag is active

    Parameters
    ----------
    ifo: string
        The interferometer to query (H1, L1).
    name: string
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
    cache: bool
        If true cache the query. Default is not to cache

    Returns
    ---------
    segments: glue.segments.segmentlist
        List of segments
    """
    info = name.split(':')
    if len(info) == 2:
        segment_name, version = info
    elif len(info) == 1:
        segment_name = info[0]
        version = 1

    flag_segments = segmentlist([])

    if source in ['GWOSC', 'any']:
        # Special cases as the LOSC convention is backwards from normal
        # LIGO / Virgo operation!!!!
        if (('_HW_INJ' in segment_name and 'NO' not in segment_name) or
           'VETO' in segment_name):
            data = query_flag(ifo, 'DATA', start_time, end_time)

            if '_HW_INJ' in segment_name:
                name = 'NO_' + segment_name
            else:
                name = segment_name.replace('_VETO', '')

            negate = query_flag(ifo, name, start_time, end_time, cache=cache)
            return (data - negate).coalesce()

        duration = end_time - start_time
        try:
            url = GWOSC_URL.format(get_run(start_time + duration/2),
                                   ifo, segment_name,
                                   int(start_time), int(duration))

            fname = download_file(url, cache=cache)
            data = json.load(open(fname, 'r'))
            if 'segments' in data:
                flag_segments = data['segments']

        except Exception as e:
            msg = "Unable to find segments in GWOSC, check flag name or times"
            print(e)
            if source != 'any':
                raise ValueError(msg)
            else:
                print("Tried and failed GWOSC {}, trying dqsegdb", name)


            return query_flag(ifo, segment_name, start_time, end_time,
                              source='dqsegdb', server=server,
                              veto_definer=veto_definer)

    elif source == 'dqsegdb':
        # Let's not hard require dqsegdb to be installed if we never get here.
        try:
            from dqsegdb.apicalls import dqsegdbQueryTimes as query
        except ImportError:
            raise ValueError("Could not query flag. Install dqsegdb"
                             ":'pip install dqsegdb'")

        # The veto definer will allow the use of MACRO names
        # These directly correspond the name defined in the veto definer file.
        if veto_definer is not None:
            veto_def = parse_veto_definer(veto_definer)

        # We treat the veto definer name as if it were its own flag and
        # a process the flags in the veto definer
        if veto_definer is not None and segment_name in veto_def[ifo]:
            for flag in veto_def[ifo][segment_name]:
                partial = segmentlist([])
                segs = query("https", server, ifo, flag['name'],
                             flag['version'], 'active',
                             int(start_time), int(end_time))[0]['active']

                # Apply padding to each segment
                for rseg in segs:
                    seg_start = rseg[0] + flag['start_pad']
                    seg_end = rseg[1] + flag['end_pad']
                    partial.append(segment(seg_start, seg_end))

                # Limit to the veto definer stated valid region of this flag
                send = segmentlist([segment([flag['start'],
                                             flag['end']])])
                flag_segments += (partial.coalesce() & send)

        else:  # Standard case just query directly.
            try:
                segs = query("https", server, ifo, name, version,
                             'active', int(start_time),
                             int(end_time))[0]['active']
                for rseg in segs:
                    flag_segments.append(segment(rseg[0], rseg[1]))
            except Exception as e:
                print("Could not query flag, check name "
                      " (%s) or times" % segment_name)
                raise e

    else:
        raise ValueError("Source must be dqsegdb or GWOSC."
                         " Got {}".format(source))

    return segmentlist(flag_segments).coalesce()


def query_cumulative_flags(ifo, segment_names, start_time, end_time,
                           source='any', server="segments.ligo.org",
                           veto_definer=None,
                           bounds=None,
                           padding=None,
                           override_ifos=None,
                           cache=False):
    """Return the times where any flag is active

    Parameters
    ----------
    ifo: string or dict
        The interferometer to query (H1, L1). If a dict, an element for each
        flag name must be provided.
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
    bounds: dict, Optional
        Dict containing start end tuples keyed by the flag name which
    indicated places which should have a distinct time period to be active.
    padding: dict, Optional
        Dict keyed by the flag name. Each element is a tuple
    (start_pad, end_pad) which indicates how to change the segment boundaries.
    override_ifos: dict, Optional
        A dict keyed by flag_name to override the ifo option on a per flag
    basis.

    Returns
    ---------
    segments: glue.segments.segmentlist
        List of segments
    """
    total_segs = segmentlist([])
    for flag_name in segment_names:
        ifo_name = ifo
        if override_ifos is not None and flag_name in override_ifos:
            ifo_name = override_ifos[flag_name]

        segs = query_flag(ifo_name, flag_name, start_time, end_time,
                          source=source, server=server,
                          veto_definer=veto_definer,
                          cache=cache)

        if padding and flag_name in padding:
            s, e = padding[flag_name]
            segs2 = segmentlist([])
            for seg in segs:
                segs2.append(segment(seg[0] + s, seg[1] + e))
            segs = segs2

        if bounds is not None and flag_name in bounds:
            s, e = bounds[flag_name]
            valid = segmentlist([segment([s, e])])
            segs = (segs & valid).coalesce()


        total_segs = (total_segs + segs).coalesce()
    return total_segs

def parse_flag_str(flag_str):
    """ Parse a dq flag query string

    Parameters
    ----------
    flag_str: str
        String needing to be parsed

    Returns
    -------
    flags: list of strings
        List of reduced name strings which can be passed to lower level
    query commands
    signs: dict
        Dict of bools indicated if this should add positively to the segmentlist
    ifos: dict
        Ifo specified for the given flag
    bounds: dict
        The boundary of a given flag
    padding: dict
        Any padding that should be applied to the segments in for a given
        flag.
    """
    flags = flag_str.replace(' ', '').strip().split(',')

    signs = {}
    ifos = {}
    bounds = {}
    padding = {}
    bflags = []

    for flag in flags:
        # Check if the flag should add or subtract time
        sign = flag[0] == '+'
        flag = flag[1:]

        ifo = pad = bound = None

        # Check for non-default IFO
        if len(flag.split(':')[0]) == 2:
            ifo = flag.split(':')[0]
            flag = flag[3:]

        # Check for padding options
        if '<' in flag:
            popt = flag.split('<')[1].split('>')[0]
            spad, epad = popt.split(':')
            pad = (float(spad), float(epad))
            flag = flag.replace(popt, '').replace('<>', '')

        # Check if there are bounds on the flag
        if '[' in flag:
            bopt = flag.split('[')[1].split(']')[0]
            start, end = bopt.split(':')
            bound = (int(start), int(end))
            flag = flag.replace(bopt, '').replace('[]', '')

        if ifo:
            ifos[flag] = ifo
        if pad:
            padding[flag] = pad
        if bound:
            bounds[flag] = bound
        bflags.append(flag)
        signs[flag] = sign

    return bflags, signs, ifos, bounds, padding

def query_str(ifo, flag_str, start_time, end_time, server="segments.ligo.org",
              veto_definer=None):
    """ Query for flags based on a special str syntax

    Parameters
    ----------
    ifo: str
        The ifo to be mainly quering for. (may be overriden in syntax)
    flag_str: str
        Specification of how to do the query. Ex. +H1:DATA:1<-8,8>[0,100000000]
        would return H1 time for the DATA available flag with version 1. It
        would then apply an 8 second padding and only return times within
        the chosen range 0,1000000000.
    start_time: int
        The start gps time. May be overriden for individual flags with the
        flag str bounds syntax
    end_time: int
        The end gps time. May be overriden for individual flags with the
        flag str bounds syntax

    Returns
    -------
    segs: segmentlist
        A list of segments corresponding to the flag query string
    """
    flags, sign, ifos, bounds, padding = parse_flag_str(flag_str)
    up = [f for f in flags if sign[f]]
    down = [f for f in flags if not sign[f]]

    if len(up) + len(down) != len(flags):
        raise ValueError('Not all flags could be parsed, check +/- prefix')
    segs = query_cumulative_flags(ifo, up, start_time, end_time,
                                  server=server,
                                  veto_definer=veto_definer,
                                  bounds=bounds,
                                  padding=padding,
                                  override_ifos=ifos)

    mseg = query_cumulative_flags(ifo, down, start_time, end_time,
                                  server=server,
                                  veto_definer=veto_definer,
                                  bounds=bounds,
                                  padding=padding,
                                  override_ifos=ifos)

    segs = (segs - mseg).coalesce()
    return segs
