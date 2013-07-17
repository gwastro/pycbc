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
import numpy
from pycbc.fft import ifft
from pycbc.types import zeros, TimeSeries, real_same_precision_as

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
            ht = TimeSeries(zeros(self.N), delta_t=self.dt, dtype=self.asd.dtype)    
            stilde = stilde / self.asd
            stilde[0:self.kmin].clear()     
            ifft(stilde, ht)
            self.segs.append(ht)
           
    def template(self, template):
        # Generate the whitened template
        tmp = template / self.asd
        tmp[0:self.kmin] = 0 
        tmp.resize(self.N)
        tmplt = TimeSeries(zeros(self.N), delta_t=self.dt, dtype=template.dtype)
        ifft(tmp, tmplt)
        
        # Calculate the location of the bins (equal power in each bin)
        power =  tmplt.squared_norm().cumsum() * 4.0 * self.dt
        sigmasq = power[-1]
        edge_vec = numpy.arange(0, self.nbins) * sigmasq / self.nbins
        bins = numpy.searchsorted(power.numpy(), edge_vec, side='right')
        bins = numpy.append(bins, len(tmplt))
    
        norm = (4 * self.dt / numpy.sqrt(sigmasq)) ** 2
        return tmplt, bins, norm     
        
    def chisq(self, tmplt, bins, norm, seg_idx, snrs, indices): 
        # For each snr value calculate the chisq value   
        vals = []
        for i, snr in zip(indices, snrs):
            # Move the whitened h(t) to the right time
            dat = self.segs[seg_idx]

            #Calculate the chisq for that point in time
            chisq = 0
            for j in range(len(bins)-1): 
                # Get the bin (b)eginning and (e)nd
                b = int(bins[j])
                e = int(bins[j+1])  
                
                # Get the time shifted bin to index the whitened h(t)
                bh = b + i
                eh = e + i
                
                # If the end is beyond a boundary shift it back
                if bh >= len(dat):
                    bh -= len(dat)
                    
                if eh > len(dat):
                    eh -= len(dat)
                        
                # If bin is within the segment go ahead and calculate it
                # else, calculate it as two separate pieces.
                if eh > bh:      
                    m = tmplt[b:e].vdot(dat[bh:eh])
                else:
                    c = b + (len(dat) - bh)
                    m = tmplt[b:c].vdot(dat[bh:len(dat)])
                    m += tmplt[c:e].vdot(dat[0:eh])
                    
                chisq += (m.conj() * m).real
                
            chisq *= self.nbins * norm     
            chisq -=  (snr.conj() * snr).real      
            vals.append(chisq) 
        return vals
