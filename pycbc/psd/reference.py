#!/usr/bin/python
# Copyright (C) 2012 Alex Nitz
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
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
""" Various reference psds from LALSimulation and utilites to read psds from 
ascii files.
"""

import numpy
from pycbc.types import FrequencySeries,zeros
from scipy.interpolate import interp1d
from pycbc.filter import get_cutoff_indices
import lalsimulation as sim

def psd_from_asd_file(filename,length,delta_f,lower_frequency_cutoff):
    """Returns the psd from an ascii file containing an asd.
    """
    fpsd = numpy.loadtxt(filename)          
    freq_data=fpsd[:,0]
    psd_data=fpsd[:,1]
    psd_data = psd_data **2
    psd_interp= interp1d(freq_data,psd_data) 
    N = (length-1) * 2  

    psd = numpy.zeros(length)
    for k in range(0,kmax,1):
        if (k<kmin):
            psd[k]=0
        else:
            psd[k]=float(psd_interp( k* delta_f ) )                   
   
    return FrequencySeries(psd,delta_f=delta_f,copy=False)


_lalsim_psd_functions = { 'aLIGOBHBH20Deg':sim.SimNoisePSDaLIGOBHBH20Deg,                
    'aLIGOQuantumZeroDetLowPower':sim.SimNoisePSDaLIGOQuantumZeroDetLowPower,
    'AdvVirgo':sim.SimNoisePSDAdvVirgo,                 
    'aLIGOHighFrequency':sim.SimNoisePSDaLIGOHighFrequency,           
    'aLIGOThermal':sim.SimNoisePSDaLIGOThermal,
    'PSDGEO':sim.SimNoisePSDGEO,                  
    'aLIGONSNSOpt':sim.SimNoisePSDaLIGONSNSOpt,                 
    'aLIGOZeroDetHighPower':sim.SimNoisePSDaLIGOZeroDetHighPower,
    'KAGRA':sim.SimNoisePSDKAGRA,                         
    'NoSRMHighPower':sim.SimNoisePSDaLIGONoSRMHighPower,          
    'aLIGOZeroDetLowPower':sim.SimNoisePSDaLIGOZeroDetLowPower,
    'MirrorTherm':sim.SimNoisePSDMirrorTherm,                   
    'aLIGONoSRMLowPower':sim.SimNoisePSDaLIGONoSRMLowPower,            
    'eLIGOModel':sim.SimNoisePSDeLIGOModel,
    'Quantum':sim.SimNoisePSDQuantum,                       
    'QuantumBHBH20Deg':sim.SimNoisePSDaLIGOQuantumBHBH20Deg,         
    'eLIGOShot':sim.SimNoisePSDeLIGOShot,
    'Seismic':sim.SimNoisePSDSeismic,                       
    'aLIGOQuantumHighFrequency':sim.SimNoisePSDaLIGOQuantumHighFrequency,     
    'iLIGOModel':sim.SimNoisePSDiLIGOModel,
    'Shot':sim.SimNoisePSDShot,                          
    'QuantumNSNSOpt':sim.SimNoisePSDaLIGOQuantumNSNSOpt,           
    'iLIGOSRD':sim.SimNoisePSDiLIGOSRD,
    'SuspTherm':sim.SimNoisePSDSuspTherm,                    
    'aLIGOQuantumNoSRMHighPower':sim.SimNoisePSDaLIGOQuantumNoSRMHighPower,    
    'iLIGOSeismic':sim.SimNoisePSDiLIGOSeismic,
    'TAMA':sim.SimNoisePSDTAMA,                          
    'QuantumNoSRMLowPower':sim.SimNoisePSDaLIGOQuantumNoSRMLowPower,     
    'iLIGOShot':sim.SimNoisePSDiLIGOShot,
    'Virgo':sim.SimNoisePSDVirgo,                        
    'aLIGOQuantumZeroDetHighPower':sim.SimNoisePSDaLIGOQuantumZeroDetHighPower,  
    'iLIGOThermal':sim.SimNoisePSDiLIGOThermal
   }

def list_psds():
    """ Return the available reference psds.
    """
    for psd_name in _lalsim_psd_functions.keys():
        print psd_name


def reference_psd(psd_name,length,delta_f,lower_frequency_cutoff):
    """ Return a FrequencySeries containing a reference psd
    """
    psd = FrequencySeries(zeros(length), delta_f=delta_f)
    kmin = lower_frequency_cutoff / delta_f
    for k in range(kmin,length):
        psd[k] = _lalsim_psd_functions[psd_name]( k * delta_f)
    # Set the last element to zero (assuming it is Nyquist)
    psd[length] = 0 
    return psd













