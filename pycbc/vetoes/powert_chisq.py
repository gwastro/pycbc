# Copyright (C) 2013 Alex Nitz
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
from pycbc.fft import ifft
from pycbc.fft.lalfft import set_measure_level, get_measure_level
from pycbc.types import zeros, float32, complex64, TimeSeries
from math import sqrt
import numpy

def time_matched_filter(template, psd, data, snrs, indices, num_bins, fmin):
    n = len(data)
    N = (n-1)*2
    dt = 1.0 / (N * data.delta_f)
    
    time_data = TimeSeries(zeros(N), delta_t=dt, dtype=float32)
    time_template = TimeSeries(zeros(N), delta_t=dt, dtype=complex64)
    
    asd = psd ** 0.5  
    kmin = fmin / template.delta_f  
    # Separately whiten the template and the data
    tmp = template / asd
    tmp[0:kmin] = 0 
    
    data = data / asd
    data[0:kmin] = 0
    
    tmp.resize(N)
    
    # Convert them to the time domain so we can do the filtering there
    ifft(data, time_data)
    ifft(tmp, time_template)
    
    # Calculate where the bins should be
    power = time_template.squared_norm().cumsum() * 4 * dt
    sigmasq = power[-1]
    edge_vec = numpy.arange(0, num_bins) * sigmasq / num_bins
    bins = numpy.searchsorted(power.numpy(), edge_vec, side='right')
    bins = numpy.append(bins, len(time_template))
    
    norm = (4 * dt / sqrt(sigmasq)) ** 2
    vals = []
    # Calculate the chisq for only the given points in the snr time series
    
    for i, snr in zip(indices, snrs):
        dat = time_data * 1
        dat.roll(-i)
        chisq = 0
        for j in range(len(bins)-1): 
            b = int(bins[j])
            e = int(bins[j+1])           
            time_template._epoch = dat._epoch       
            m = time_template[b:e].inner(dat[b:e]) 
            chisq += (m.conj() * m).real
        chisq *= num_bins * norm     
        chisq -=  (snr.conj() * snr).real
        
        vals.append(chisq) 
    return vals

class PowerTChisq(object):
    def __init__ (self, num_bins, segs, psd, fmin):
        self.nbins = num_bins
        self.kmin = int(fmin / psd.delta_f)
        self.asd = psd ** 0.5  
        self.segs = []
        
        self.N = (len(psd)-1)*2
        self.dt = 1.0 / (self.N * psd.delta_f)
         
        # Pregenerate the whitened h(t) segments 
        for stilde in segs:
            ht = TimeSeries(zeros(self.N), delta_t=self.dt, dtype=float32)    
            stilde = stilde / self.asd
            stilde[0:self.kmin] = 0     
            ifft(stilde, ht)
            self.segs.append(ht)
           
    def template(self, template):
        # Generate the whitened template
        tmp = template / self.asd
        tmp[0:self.kmin] = 0 
        tmp.resize(self.N)
        tmplt = TimeSeries(zeros(self.N), delta_t=self.dt, dtype=complex64)
        ifft(tmp, tmplt)
        
        # Calculate the location of the bins (equal power in each bin)
        power =  tmplt.squared_norm().cumsum() * 4.0 * self.dt
        sigmasq = power[-1]
        edge_vec = numpy.arange(0, self.nbins) * sigmasq / self.nbins
        bins = numpy.searchsorted(power.numpy(), edge_vec, side='right')
        bins = numpy.append(bins, len(tmplt))
    
        norm = (4 * self.dt / sqrt(sigmasq)) ** 2
        return tmplt, bins, norm     
        
    def chisq(self, tmplt, bins, norm, seg_idx, snrs, indices): 
        # For each snr value calculate the chisq value   
        vals = []
        for i, snr in zip(indices, snrs):
            # Move the whitened h(t) to the right time
            dat = self.segs[seg_idx] * 1
            dat.roll(-i)
            
            #Calculate the chisq for that point in time
            chisq = 0
            for j in range(len(bins)-1): 
                b = int(bins[j])
                e = int(bins[j+1])  
                         
                tmplt._epoch = dat._epoch       
                #m = tmplt[b:e].inner(dat[b:e]) 
                m = numpy.vdot(tmplt[b:e].data, dat[b:e].data)
                chisq += (m.conj() * m).real
                
            chisq *= self.nbins * norm     
            chisq -=  (snr.conj() * snr).real      
            vals.append(chisq) 
        return vals

