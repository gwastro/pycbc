# Copyright (C) 2012  Alex Nitz
#
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
"""This module provides utilities for injection signals into data
"""

import lal
import lalsimulation as sim
from pycbc.waveform import get_td_waveform
from glue.ligolw import utils as ligolw_utils
from glue.ligolw import table, lsctables
from pycbc.types import float64, float32, TimeSeries
from pycbc.detector import Detector

def inject_td_waveform_into_data(data, hp, hc, detector, ra, dec, pol):
    """ Inject waveform into the data for a given detector.
    """

    signal = sim.XLALSimDetectorStrainREAL8TimeSeries(hp.lal(), hc.lal(), ra, 
                                             dec, spi, detector.FrDetector)
                                             
    laldata = data.lal()                                             
    XLALSimAddInjectionREAL8TimeSeries(laldata, signal, NULL);
    data.data[:]=laldata.data.data[:]    
    
class InjectionSet(object):
    def __init__(self, sim_file, **kwds):
        self.indoc = ligolw_utils.load_filename(sim_file, False)
        self.table = table.get_table(self.indoc, lsctables.SimInspiralTable.tableName)
        self.extra_args = kwds
        
    def get_injections_within_time_window(self, start_time, end_time):
        return []
    
    def inject_into_data(self, strain, detector_name):
        lalstrain = strain.lal()    
        injections = self.get_injections_within_time_window(strain.start_time, strain.end_time)  
        detector = Detector(detector_name)
        
        for injection in injections:
            end_time = lal.LIGOTimeGPS(injection.geocent_end_time, injection.geocent_end_time_ns)
            hp, hc = get_td_waveform(injection, **self.extra_args)
            hp._epoch = end_time
            hc._epoch = end_time
            signal =  sim.XLALSimDetectorStrainREAL8TimeSeries(hp.astype(float64).lal(), hc.astype(float64).lal(), injection.longitude, 
                                             injection.latitude, injection.polarizations, detector.FrDetector)  
                                             
            if strain.dtype is float64:                                  
                XLALSimAddInjectionREAL8TimeSeries(lalstrain, signal.astype(float64), NULL)
            elif strain.dtype is float32:                                  
                XLALSimAddInjectionREAL4TimeSeries(lalstrain, signal.astype(float32), NULL)
            else:
                raise TypeError("Invalid strain dtype: Must be float32, or float64")              
                    
        return TimeSeries(lalstrain, epoch=strain.epoch, delta_t=strain.delta_t, copy=False)

