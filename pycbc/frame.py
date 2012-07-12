# Copyright (C) 2012 Andrew Miller
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
This modules uses PyLAL to read in a frame or cache, and returns a pycbc TimeSeries.
"""

import pycbc
import pycbc.types
from glue import lal
from pylal import Fr, frutils
import swiglal as _swiglal
import numpy as _numpy

def ReadFrame(filename, channels, start=None, end=None):
'''
This function will read a Time-Series  or list of Time-Series in from a LIGO frame file. 
It accepts a filename (str), a channel (str) or list of channels, and a start and end time 
(these nust be integers, they can be LIGOTimeGPS instances, or simply  numbers)
'''
    if start is not None and end is not None:
        start = float(start)
        end = float(end)
        span = end - start
        if span < 0:
            raise ValueError('beginning must be before end')
        if start%1 != 0 or end%1 != 0:
            raise ValueError('start and end times must be integers')
    else:
        start = -1
        span = -1

    if type(channels) is list:
        ts = []
        for channel in channels:
            frdata = Fr.frgetvect1d(filename, channel, start, span)
            ts.append(pycbc.types.TimeSeries(frdata[0],frdata[3], _swiglal.LIGOTimeGPS(frdata[1]), copy=False))
    else:
        frdata = Fr.frgetvect1d(filename, channels, start, span)
        ts = pycbc.types.TimeSeries(frdata[0],frdata[3], _swiglal.LIGOTimeGPS(frdata[1]), copy=False)
    return ts

def ReadCache(filename, channels, start, end):
'''
This function will read a Time-Series  or list of Time-Series in from a LIGO cache file. 
It accepts a filename (str), a channel (str) or list of channels, and a start and end time 
(these nust be integers, they can be LIGOTimeGPS instances, or simply  numbers)
'''
    start = float(start)
    end = float(end)
    if start%1 != 0 or end%1 != 0:
        raise ValueError('start and end times must be integers')
    with open(filename,'r') as f:
        lal_cache = lal.Cache.fromfile(f)
    frdata = frutils.FrameCache(lal_cache)
    if type(channels) is list:
        ts = []
        for channel in channels:
            data = frdata.fetch(channel,start,end)
            dt = data.metadata.dt
            ts.append(pycbc.types.TimeSeries(data.A,dt, _swiglal.LIGOTimeGPS(data.metadata.segments[0][0]),copy=False))
    else:
        data = frdata.fetch(channels,start,end)
        dt = data.metadata.dt
        ts = pycbc.types.TimeSeries(data.A,dt, _swiglal.LIGOTimeGPS(data.metadata.segments[0][0]),copy=False)
    return ts
        
        
    
