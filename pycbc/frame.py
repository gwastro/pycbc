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
import lal as _lal
import numpy as _numpy

def read_frame(filename, channels, start=None, end=None):
    if start is not None and end is not None:
        if type(start) is _lal.LIGOTimeGPS:
            if start.gpsNanoSeconds != 0:
                raise ValueError('start and end times must be integer valued')
            else:
                start = start.gpsSeconds
        else:
            if int(start) != start:
                raise ValueError('start and end times must be integer valued')
            else:
                start = int(start)                
        if type(end) is _lal.LIGOTimeGPS:
            if end.gpsNanoSeconds != 0:
                raise ValueError('start and end times must be integer valued')
            else:
                end = end.gpsSeconds
        else:
            if int(end) != end:
                raise ValueError('start and end times must be integer valued')
            else:
                end = int(end)             
        span = end - start
        if span <= 0:
            raise ValueError('beginning must be before end')
    else:
        start = -1
        span = -1

    if type(channels) is list:
        ts = []
        for channel in channels:
            frdata = Fr.frgetvect1d(filename, channel, start, span)
            ts.append(pycbc.types.TimeSeries(initial_array=frdata[0],
                                            delta_t=frdata[3],
                                            epoch=_lal.LIGOTimeGPS(frdata[1]),
                                            copy=False))
    else:
        frdata = Fr.frgetvect1d(filename, channels, start, span)
        ts = pycbc.types.TimeSeries(initial_array=frdata[0],
                                    delta_t=frdata[3],
                                    epoch=_lal.LIGOTimeGPS(frdata[1]),
                                    copy=False)
    return ts

def read_cache(filename, channels, start, end):
    if type(start) is _lal.LIGOTimeGPS:
        if start.gpsNanoSeconds != 0:
            raise ValueError('start and end times must be integer valued')
        else:
            start = start.gpsSeconds
    else:
        if int(start) != start:
            raise ValueError('start and end times must be integer valued')
        else:
            start = int(start)                
    if type(end) is _lal.LIGOTimeGPS:
        if end.gpsNanoSeconds != 0:
            raise ValueError('start and end times must be integer valued')
        else:
            end = end.gpsSeconds
    else:
        if int(end) != end:
            raise ValueError('start and end times must be integer valued')
        else:
            end = int(end)
    span = end - start
    if span <= 0:
        raise ValueError('beginning must be before end')
    with open(filename,'r') as f:
        lal_cache = lal.Cache.fromfile(f)
    frdata = frutils.FrameCache(lal_cache)
    if type(channels) is list:
        ts = []
        for channel in channels:
            data = frdata.fetch(channel,start,end)
            dt = data.metadata.dt
            # Now we actually create the new TimeSeries, and append it to the list
            ts.append(pycbc.types.TimeSeries(initial_array=data.A,
                                                delta_t=dt,
                                                epoch=_lal.LIGOTimeGPS(data.metadata.segments[0][0]),
                                                copy=False))
    else:
        data = frdata.fetch(channels,start,end)
        dt = data.metadata.dt
        ts = pycbc.types.TimeSeries(initial_array=data.A,
                                    delta_t=dt,
                                    epoch=_lal.LIGOTimeGPS(data.metadata.segments[0][0]),
                                    copy=False)
    return ts
        
        
    
