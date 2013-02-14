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
from pycbc import DYN_RANGE_FAC

class TemplateBank(object):
    def __init__(self, filename, approximant, filter_length, delta_f, f_lower,  dtype, out=None, **kwds):
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
        
        self.indoc = ligolw_utils.load_filename(filename, False)
        
        try :
            self.table = table.get_table(self.indoc, lsctables.SnglInspiralTable.tableName) 
        except ValueError:
            self.table = table.get_table(self.indoc, lsctables.SimInspiralTable.tableName)

        self.extra_args = kwds
        
    def __len__(self):
        return len(self.table)
            
    def __iter__(self):
        self.index=-1
        return self
        
    def current_f_end(self):
        f_end = pycbc.waveform.get_waveform_end_frequency(self.table[self.index], approximant=self.approximant, **self.extra_args) 
        
        if f_end > self.filter_length * self.delta_f:
            return (self.filter_length-1) * self.delta_f
        else:
            return f_end

    def next(self):
        self.index +=1
        if self.index == len(self.table):
            raise StopIteration
        else:
            poke  = self.out.data
            self.out.clear()
            htilde = pycbc.waveform.get_waveform_filter(self.out[0:self.filter_length], self.table[self.index], 
                                    approximant=self.approximant, f_lower=self.f_lower, delta_f=self.delta_f, **self.extra_args)
            return htilde.astype(self.dtype) 

        
    





