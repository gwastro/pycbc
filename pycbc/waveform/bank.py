# Copyright (C) 2012  Alex Nitz, Josh Willis, Andrew Miller
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
This module provides classes that describe banks of waveforms
"""
from pycbc.types import zeros, complex64
from glue.ligolw import utils as ligolw_utils
from glue.ligolw import table, lsctables
import pycbc.waveform
from pycbc.types import FrequencySeries
from pycbc import DYN_RANGE_FAC        

class TemplateBank(object):
    def __init__(self, filename, approximant, filter_length, delta_f, f_lower,  dtype, psd=None, out=None, **kwds):
        if out:
            self.out = out
        else:
            self.out = zeros(filter_length, dtype=dtype)
        self.dtype = dtype
        self.f_lower = f_lower
        self.approximant = approximant
        self.filename = filename
        self.delta_f = delta_f
        self.filter_length = filter_length
        self.kmin = int(f_lower / delta_f)
        
        self.indoc = ligolw_utils.load_filename(filename, False)     
        self.table = table.get_table(self.indoc, lsctables.SnglInspiralTable.tableName) 
        self.extra_args = kwds
        self.psd = psd
        
        #If we can for this template pregenerate the sigmasq vector 
        self.sigmasq_vec = None
        if (psd is not None) and pycbc.waveform.waveform_norm_exists(approximant):
                self.sigmasq_vec = pycbc.waveform.get_waveform_filter_norm(
                                     approximant, self.psd, filter_length, 
                                     self.delta_f, self.f_lower)             
        
    def __len__(self):
        return len(self.table)
            
    def __iter__(self):
        self.index=-1
        return self
        
    def current_amplitude_norm(self):
        amp_norm = pycbc.waveform.get_template_amplitude_norm(self.table[self.index], 
                              approximant=self.approximant, **self.extra_args)
        if amp_norm:
            return amp_norm
        else:
            return 1

    def current_f_end(self):
        f_end = pycbc.waveform.get_waveform_end_frequency(self.table[self.index], 
                              approximant=self.approximant, **self.extra_args) 
        
        if f_end is None:
            return self.filter_length * self.delta_f
        
        if f_end > self.filter_length * self.delta_f:
            return (self.filter_length-1) * self.delta_f
        else:
            return f_end

    def next(self):
        kmax = int(self.current_f_end() / self.delta_f)
        self.index +=1
        if self.index == len(self.table):
            raise StopIteration
        else:
            poke  = self.out.data
            self.out[self.kmin:kmax].clear()
            htilde = pycbc.waveform.get_waveform_filter(self.out[0:self.filter_length], self.table[self.index], 
                                    approximant=self.approximant, f_lower=self.f_lower, delta_f=self.delta_f, **self.extra_args)
                                    
            htilde = htilde.astype(self.dtype)
            
            htilde.end_frequency = self.current_f_end()
            htilde.end_idx = int(htilde.end_frequency / htilde.delta_f) 
            htilde.params = self.table[self.index]
            htilde.amplitude_norm = self.current_amplitude_norm()
            
            if self.psd is not None:
                if self.sigmasq_vec is not None:
                    
                    htilde.sigmasq = self.sigmasq_vec[htilde.end_idx]
                else:
                    htilde.sigmasq = sigmasq(htilde, self.psd, low_frequency_cutoff=self.f_low) 
 
            return htilde

        
    





