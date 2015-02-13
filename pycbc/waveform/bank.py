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
import types
import pycbc.waveform
from pycbc.types import zeros
from glue.ligolw import ligolw, table, lsctables, utils as ligolw_utils
from pycbc.filter import sigmasq
from pycbc import DYN_RANGE_FAC

def sigma_cached(self, psd):
    """ Cache sigma calculate for use in tandem with the FilterBank class
    """
    key = id(psd)
    if key not in self._sigmasq:
        # If possible, we precalculate the sigmasq vector for all possible waveforms
        if pycbc.waveform.waveform_norm_exists(self.approximant):
            if not hasattr(psd, 'sigmasq_vec'):
                psd.sigmasq_vec = pycbc.waveform.get_waveform_filter_norm(
                     self.approximant, psd, len(psd), psd.delta_f, self.f_lower)
                
            # Get an amplitude normalization (mass dependant constant norm)
            amp_norm = pycbc.waveform.get_template_amplitude_norm(
                                 self.params, approximant=self.approximant)
            amp_norm = 1 if amp_norm is None else amp_norm
            scale = DYN_RANGE_FAC * amp_norm
            self._sigmasq[key] = psd.sigmasq_vec[self.end_idx] * (scale) **2

        else:
            self._sigmasq[key] = sigmasq(self, psd, low_frequency_cutoff=self.f_lower)                    
    return self._sigmasq[key]
    
# dummy class needed for loading LIGOLW files
class LIGOLWContentHandler(ligolw.LIGOLWContentHandler):
    pass
lsctables.use_in(LIGOLWContentHandler)

class FilterBank(object):
    def __init__(self, filename, approximant, filter_length, delta_f, f_lower,
                 dtype, out=None, **kwds):
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

        self.indoc = ligolw_utils.load_filename(
            filename, False, contenthandler=LIGOLWContentHandler)
        self.table = table.get_table(
            self.indoc, lsctables.SnglInspiralTable.tableName)
        self.extra_args = kwds  

    def __len__(self):
        return len(self.table)

    def __getitem__(self, index):
        # Make new memory for templates if we aren't given output memory
        if self.out is None:
            tempout = zeros(self.filter_length, dtype=self.dtype)
        else:
            tempout = self.out

        # Get the end of the waveform if applicable (only for SPAtmplt atm)
        f_end = pycbc.waveform.get_waveform_end_frequency(self.table[index],
                              approximant=self.approximant, **self.extra_args)

        if f_end is None or f_end >= (self.filter_length * self.delta_f):
            f_end = (self.filter_length-1) * self.delta_f

        poke  = tempout.data
        # Clear the storage memory
        tempout.clear()

        # Get the waveform filter
        distance = 1.0 / DYN_RANGE_FAC
        htilde = pycbc.waveform.get_waveform_filter(
            tempout[0:self.filter_length], self.table[index],
            approximant=self.approximant, f_lower=self.f_lower, f_final=f_end,
            delta_f=self.delta_f, delta_t=self.delta_t, distance=distance,
            **self.extra_args)

        # For time domain templates, record the total duration (which may
        # include ringdown) and the duration up to merger since they will be 
        # erased by the type conversion below
        length_in_time = None
        chirp_length = None
        if hasattr(htilde, 'length_in_time'):
            length_in_time = htilde.length_in_time
        if hasattr(htilde, 'chirp_length'):
            chirp_length = htilde.chirp_length

        htilde = htilde.astype(self.dtype)
        htilde.f_lower = self.f_lower
        htilde.end_frequency = f_end
        htilde.end_idx = int(htilde.end_frequency / htilde.delta_f)
        htilde.params = self.table[index]
        htilde.approximant = self.approximant
        
        # Add sigmasq as a method of this instance
        htilde.sigmasq = types.MethodType(sigma_cached, htilde)
        htilde._sigmasq = {}

        # For time domain templates, assign 'ttotal' to the total template
        # duration (which may include ringdown) and 'template_duration' to
        # the duration up to merger as in previous IMR searches
        if length_in_time is not None:
            htilde.length_in_time = length_in_time
            self.table[index].ttotal = length_in_time
        if chirp_length is not None:
            htilde.chirp_length = chirp_length
            self.table[index].template_duration = chirp_length
        return htilde
