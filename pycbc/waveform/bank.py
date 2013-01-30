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
from glue.ligolw import utils as ligolw_utils
from glue.ligolw import table, lsctables
from pycbc.waveform import get_waveform_filter, get_waveform_filter_precondition
from pycbc import DYN_RANGE_FAC

class TemplateBank(object):
    def __init__(self, filename, filter_length, approximant, **kwds):
        self.approximant = approximant
        self.filename = filename
        self.filter_length = filter_length
        self.indoc = ligolw_utils.load_filename(filename, False)
        try :
            self.table = table.get_table(self.indoc, lsctables.SnglInspiralTable.tableName) 
        except ValueError:
            self.table = table.get_table(self.indoc, lsctables.SimInspiralTable.tableName)

        self.index = -1
        self.extra_args = kwds
    
    def precondition(self, stilde):
        """ Return the fourier domain strain preconditioned for the approximant
            used in the template bank.
        """ 
        precondition_factor = get_waveform_filter_precondition(self.approximant, self.filter_length)
        if precondition_factor is not None:
            return precondition_factor * stilde
        else:
            return strain
        

    def __len__(self):
        return len(self.table)
            
    def __iter__(self):
        return self

    def next(self):
        self.index +=1
        if self.index == len(self.table):
            raise StopIteration
        else:
            return get_waveform_filter(self.length, self.table[self.index], **self.extra_args)

        
    





