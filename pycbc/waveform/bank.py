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
from pycbc.filter import sigmasq
from pycbc import DYN_RANGE_FAC        

class TemplateBank(object):
    def __init__(self, filename, approximant, filter_length, delta_f, f_lower,  dtype, psd=None, out=None, **kwds):
        self.out = out
        self.dtype = dtype
        self.f_lower = f_lower
        self.approximant = approximant
        self.filename = filename
        self.delta_f = delta_f
        self.N = (filter_length - 1 ) * 2
        self.delta_t = 1.0 / (self.N * self.delta_f)
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
    
    def __getitem__(self, index):
        # Make new memory for templates if we aren't given output memory
        if self.out is None:
            self.out = zeros(self.filter_length, dtype=self.dtype)
    
        # Get the end of the waveform if applicable (only for SPAtmplt atm)
        f_end = pycbc.waveform.get_waveform_end_frequency(self.table[index], 
                              approximant=self.approximant, **self.extra_args) 
        
        if f_end is None or f_end >= (self.filter_length * self.delta_f):
            f_end = (self.filter_length-1) * self.delta_f
        
        poke  = self.out.data
        # Clear the storage memory
        self.out.clear()        
        
        # Get the waveform filter
        distance = 1.0 / DYN_RANGE_FAC
        htilde = pycbc.waveform.get_waveform_filter(self.out[0:self.filter_length], 
                            self.table[index], approximant=self.approximant, 
                            f_lower=self.f_lower, delta_f=self.delta_f, delta_t=self.delta_t, 
                            distance=distance, **self.extra_args)
            
        length_in_time = None      
        if hasattr(htilde, 'length_in_time'):
            length_in_time = htilde.length_in_time      
            
        # Make sure it is the desired type
        htilde = htilde.astype(self.dtype)                   

        htilde.end_frequency = f_end
        htilde.end_idx = int(htilde.end_frequency / htilde.delta_f) 
        htilde.params = self.table[index]
        
        # If we were given a psd, calculate sigmasq so we have it for later
        if self.psd is not None:
            if self.sigmasq_vec is not None: 
                
                # Get an amplitude normalization (mass dependant constant norm)
                amp_norm = pycbc.waveform.get_template_amplitude_norm(self.table[index], 
                                  approximant=self.approximant, **self.extra_args)    
                if amp_norm is None:
                    amp_norm = 1
                scale = DYN_RANGE_FAC * amp_norm
                
                htilde.sigmasq = self.sigmasq_vec[htilde.end_idx] * (scale) **2
            else:
                htilde.sigmasq = sigmasq(htilde, self.psd, low_frequency_cutoff=self.f_lower) 

        if length_in_time is not None:
            htilde.length_in_time = length_in_time  
            self.table[index].template_duration = length_in_time
       
        return htilde
