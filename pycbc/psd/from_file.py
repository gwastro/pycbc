"""
Helper functions to load psds and asds from ascii files.
"""

from numpy import loadtxt,zeros
from pycbc import DYN_RANGE_FAC
from pycbc.types import FrequencySeries

def psd_from_file(filename,delta_f,length):
    """Returns the psd from an ascii file containing an asd
    """
    fpsd = loadtxt(options.asd_file)          
    freq_data=fpsd[:,0]
    psd_data=fpsd[:,1]*DYN_RANGE_FAC
    psd_interp= interp1d(freq_data,psd_data) 

    psd = zeros(length)
    for k in range(0,length,1):
        if (k<self.k_min):
            psd[k]=0
        else:
            psd[k]=float(psd_interp( k* delta_f )            
   
    return FrequencySeries(psd,delta_f=delta_f)
  
def asd_from_file(filename,delta_f,length)
    """Returns the interpolated asd from an ascii file containing a asd series
    """
    asd = psd_from_file(filename,delta_f,length)
    return asd ** 2
    
