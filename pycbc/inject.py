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

import lalsimulation as sim

def inject_td_waveform_into_data(data, hp, hc, detector, ra, dec, pol):
    """ Inject waveform into the data for a given detector.
    """

    signal = sim.XLALSimDetectorStrainREAL8TimeSeries(hp.lal(), hc.lal(), ra, 
                                             dec, spi, detector.FrDetector)
                                             
    laldata = data.lal()                                             
    XLALSimAddInjectionREAL8TimeSeries(laldata, signal, NULL);
    data.data[:]=laldata.data.data[:]    
    
class InjectionSet(object):
    def __init__(self, sim_file):
        pass
    
    def inject_into_data(strain):
        strain_with_injections = strain.
        return strain_with_injections

