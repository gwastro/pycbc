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
import types, numpy, logging, os.path, math, h5py
import pycbc.waveform
from pycbc.types import zeros
from glue.ligolw import ligolw, table, lsctables, utils as ligolw_utils
from pycbc.filter import sigmasq
from pycbc import DYN_RANGE_FAC
from pycbc.pnutils import nearest_larger_binary_number
from pycbc.io import FieldArray

def sigma_cached(self, psd):
    """ Cache sigma calculate for use in tandem with the FilterBank class
    """
    key = id(psd)
    if key not in self._sigmasq:
        # If possible, we precalculate the sigmasq vector for all possible waveforms
        if pycbc.waveform.waveform_norm_exists(self.approximant):
            if not hasattr(psd, 'sigmasq_vec'):
                psd.sigmasq_vec = {}
            
            if self.approximant not in psd.sigmasq_vec:
                psd.sigmasq_vec[self.approximant] = pycbc.waveform.get_waveform_filter_norm(
                     self.approximant, psd, len(psd), psd.delta_f, self.f_lower)
                
            # Get an amplitude normalization (mass dependant constant norm)
            amp_norm = pycbc.waveform.get_template_amplitude_norm(
                                 self.params, approximant=self.approximant)
            amp_norm = 1 if amp_norm is None else amp_norm
            scale = DYN_RANGE_FAC * amp_norm
            self._sigmasq[key] = psd.sigmasq_vec[self.approximant][self.end_idx] * (scale) **2

        else:
            self._sigmasq[key] = sigmasq(self, psd, low_frequency_cutoff=self.f_lower)                    
    return self._sigmasq[key]
    
# dummy class needed for loading LIGOLW files
class LIGOLWContentHandler(ligolw.LIGOLWContentHandler):
    pass
lsctables.use_in(LIGOLWContentHandler)

class TemplateBank(object):
    """ Class to provide some basic helper functions and information
    about elements of an xml template bank.
    """
    def __init__(self, filename, approximant=None, **kwds):
        ext = os.path.basename(filename)
        if 'xml' in ext:
            self.indoc = ligolw_utils.load_filename(
                filename, False, contenthandler=LIGOLWContentHandler)
            self.table = table.get_table(
                self.indoc, lsctables.SnglInspiralTable.tableName)
            self.table = FieldArray.from_ligolw_table(self.table)

            # inclination stored in xml alpha3 column
            names = list(self.table.dtype.names)
            names = tuple([n if n != 'alpha3' else 'inclination' for n in names]) 
            self.table.dtype.names = names    

        elif 'hdf' in ext:
            f = h5py.File(filename, 'r')
            dtype = []
            data = {}
            for key in f.keys():
                try:
                    data[str(key)] = f[key][:]
                    dtype.append((str(key), data[key].dtype))
                except:
                    pass

            num = len(data[data.keys()[0]])
            self.table = FieldArray(num, dtype=dtype)
            for key in data:
                self.table[key] = data[key]
        else:
            raise ValueError("Unsupported template bank file extension %s" % ext)

        if not hasattr(self.table, 'template_duration'):
            self.table = self.table.add_fields(numpy.zeros(len(self.table),
                                     dtype=numpy.float32), 'template_duration') 
        self.extra_args = kwds  
        self.approximant_str = approximant

    @staticmethod
    def parse_option(row, arg):
        safe_dict = {}
        safe_dict.update(row.__dict__)
        safe_dict.update(math.__dict__)
        return eval(arg, {"__builtins__":None}, safe_dict)

    def end_frequency(self, index):
        """ Return the end frequency of the waveform at the given index value
        """
        from pycbc.waveform.waveform import props

        return pycbc.waveform.get_waveform_end_frequency(self.table[index],
                              approximant=self.approximant(index),
                              **self.extra_args)      

    def approximant(self, index):
        """ Return the name of the approximant ot use at the given index
        """
        if self.approximant_str is not None:
            if 'params' in self.approximant_str:
                t = type('t', (object,), {'params' : self.table[index]})
                approximant = str(self.parse_option(t, self.approximant_str)) 
            else:
                approximant = self.approximant_str
        else:
            raise ValueError("Reading approximant from template bank not yet supported")

        return approximant

    def __len__(self):
        return len(self.table)

class LiveFilterBank(TemplateBank):
    def __init__(self, filename, f_lower, sample_rate, minimum_buffer,
                       approximant=None,
                       **kwds):

        self.f_lower = f_lower
        self.filename = filename
        self.sample_rate = sample_rate
        self.minimum_buffer = minimum_buffer

        super(LiveFilterBank, self).__init__(filename, approximant=approximant, **kwds)

        self.table = sorted(self.table, key=lambda t: t.mchirp)        

    def __getitem__(self, index):
        approximant = self.approximant(index)
        f_end = self.end_frequency(index)

        # Determine the length of time of the filter, rounded up to
        # nearest power of two
        min_buffer = 1.0 + self.minimum_buffer
    
        from pycbc.waveform.waveform import props
        buff_size = pycbc.waveform.get_waveform_filter_length_in_time(approximant, f_lower=self.f_lower, 
                                                                      **props(self.table[index]))
        tlen = nearest_larger_binary_number((buff_size + min_buffer) * self.sample_rate)
        flen = tlen / 2 + 1

        delta_f = self.sample_rate / float(tlen)

        if f_end is None or f_end >= (flen * delta_f):
            f_end = (flen-1) * delta_f

        logging.info("Generating %s, %ss, %i" % (approximant, 1.0/delta_f, index))

        # Get the waveform filter
        distance = 1.0 / DYN_RANGE_FAC
        htilde = pycbc.waveform.get_waveform_filter(
            zeros(flen, dtype=numpy.complex64), self.table[index],
            approximant=approximant, f_lower=self.f_lower, f_final=f_end,
            delta_f=delta_f, delta_t=1.0/self.sample_rate, distance=distance,
            **self.extra_args)

        # If available, record the total duration (which may
        # include ringdown) and the duration up to merger since they will be 
        # erased by the type conversion below.
        ttotal = template_duration = -1
        if hasattr(htilde, 'length_in_time'):
                ttotal = htilde.length_in_time
        if hasattr(htilde, 'chirp_length'):
                template_duration = htilde.chirp_length

        self.table[index].template_duration = template_duration        

        htilde = htilde.astype(numpy.complex64)
        htilde.f_lower = self.f_lower
        htilde.end_frequency = f_end
        htilde.end_idx = int(htilde.end_frequency / htilde.delta_f)
        htilde.params = self.table[index]
        htilde.approximant = approximant
        htilde.chirp_length = template_duration
        htilde.length_in_time = ttotal
        
        # Add sigmasq as a method of this instance
        htilde.sigmasq = types.MethodType(sigma_cached, htilde)
        htilde._sigmasq = {}

        return htilde

class FilterBank(TemplateBank):
    def __init__(self, filename, filter_length, delta_f, f_lower, dtype,
                 out=None, max_template_length=None,
                 approximant=None,
                 **kwds):
        self.out = out
        self.dtype = dtype
        self.f_lower = f_lower
        self.filename = filename
        self.delta_f = delta_f
        self.N = (filter_length - 1 ) * 2
        self.delta_t = 1.0 / (self.N * self.delta_f)
        self.filter_length = filter_length
        self.kmin = int(f_lower / delta_f)
        self.max_template_length = max_template_length

        super(FilterBank, self).__init__(filename, approximant=approximant, **kwds)

    def __getitem__(self, index):
        # Make new memory for templates if we aren't given output memory
        if self.out is None:
            tempout = zeros(self.filter_length, dtype=self.dtype)
        else:
            tempout = self.out

        approximant = self.approximant(index)
        f_end = self.end_frequency(index)
        if f_end is None or f_end >= (self.filter_length * self.delta_f):
            f_end = (self.filter_length-1) * self.delta_f

        # Find the start frequency, if variable
        if self.max_template_length is not None:
            f_low = find_variable_start_frequency(approximant,
                                                  self.table[index],
                                                  self.f_lower,
                                                  self.max_template_length)
        else:
            f_low = self.f_lower

        logging.info('%s: generating %s from %s Hz' % (index, approximant, f_low))

        # Clear the storage memory
        poke  = tempout.data
        tempout.clear()

        # Get the waveform filter
        distance = 1.0 / DYN_RANGE_FAC
        htilde = pycbc.waveform.get_waveform_filter(
            tempout[0:self.filter_length], self.table[index],
            approximant=approximant, f_lower=f_low, f_final=f_end,
            delta_f=self.delta_f, delta_t=self.delta_t, distance=distance,
            **self.extra_args)

        # If available, record the total duration (which may
        # include ringdown) and the duration up to merger since they will be 
        # erased by the type conversion below.
        ttotal = template_duration = None
        if hasattr(htilde, 'length_in_time'):
                ttotal = htilde.length_in_time
        if hasattr(htilde, 'chirp_length'):
                template_duration = htilde.chirp_length

        self.table[index].template_duration = template_duration        

        htilde = htilde.astype(self.dtype)
        htilde.f_lower = self.f_lower
        htilde.end_frequency = f_end
        htilde.end_idx = int(htilde.end_frequency / htilde.delta_f)
        htilde.params = self.table[index]
        htilde.approximant = approximant
        htilde.chirp_length = template_duration
        htilde.length_in_time = ttotal
        
        # Add sigmasq as a method of this instance
        htilde.sigmasq = types.MethodType(sigma_cached, htilde)
        htilde._sigmasq = {}

        return htilde

def find_variable_start_frequency(approximant, parameters, f_start, max_length, 
                                  delta_f = 1):
    """ Find a frequency value above the starting frequency that results in a 
    waveform shorter than max_length.
    """
    l = max_length + 1
    f = f_start - delta_f
    while l > max_length:
        f += delta_f
        l = pycbc.waveform.get_waveform_filter_length_in_time(approximant, parameters, f_lower=f)
    return f


class FilterBankSkyMax(TemplateBank):
    def __init__(self, filename, filter_length, delta_f, f_lower,
                 dtype, out_plus=None, out_cross=None,
                 max_template_length=None, **kwds):
        self.out_plus = out_plus
        self.out_cross = out_cross
        self.dtype = dtype
        self.f_lower = f_lower
        self.filename = filename
        self.delta_f = delta_f
        self.N = (filter_length - 1 ) * 2
        self.delta_t = 1.0 / (self.N * self.delta_f)
        self.filter_length = filter_length
        self.kmin = int(f_lower / delta_f)
        self.max_template_length = max_template_length

        super(FilterBankSkyMax, self).__init__(filename, **kwds)

    def __getitem__(self, index):
        # Make new memory for templates if we aren't given output memory
        if self.out_plus is None:
            tempoutplus = zeros(self.filter_length, dtype=self.dtype)
        else:
            tempoutplus = self.out_plus
        if self.out_cross is None:
            tempoutcross = zeros(self.filter_length, dtype=self.dtype)
        else:
            tempoutcross = self.out_cross

        approximant = self.approximant(index)

        # Get the end of the waveform if applicable (only for SPAtmplt atm)
        f_end = self.end_frequency(index)
        if f_end is None or f_end >= (self.filter_length * self.delta_f):
            f_end = (self.filter_length-1) * self.delta_f

        # Find the start frequency, if variable
        if self.max_template_length is not None:
            f_low = find_variable_start_frequency(approximant,
                                                  self.table[index],
                                                  self.f_lower,
                                                  self.max_template_length)
        else:
            f_low = self.f_lower

        logging.info('%s: generating %s from %s Hz' % (index, approximant, f_low
))

        # What does this do???
        poke1 = tempoutplus.data
        poke2 = tempoutcross.data

        # Clear the storage memory
        tempoutplus.clear()
        tempoutcross.clear()

        # Get the waveform filter
        distance = 1.0 / DYN_RANGE_FAC
        hplus, hcross = pycbc.waveform.get_two_pol_waveform_filter(
            tempoutplus[0:self.filter_length],
            tempoutcross[0:self.filter_length], self.table[index],
            approximant=approximant, f_lower=f_low,
            f_final=f_end, delta_f=self.delta_f, delta_t=self.delta_t,
            distance=distance, **self.extra_args)

        if hasattr(hplus, 'chirp_length') and hplus.chirp_length is not None:
            self.table[index].template_duration = hplus.chirp_length

        hplus = hplus.astype(self.dtype)
        hcross = hcross.astype(self.dtype)
        hplus.f_lower = f_low
        hcross.f_lower = f_low
        hplus.end_frequency = f_end
        hcross.end_frequency = f_end
        hplus.end_idx = int(hplus.end_frequency / hplus.delta_f)
        hcross.end_idx = int(hplus.end_frequency / hplus.delta_f)
        hplus.params = self.table[index]
        hcross.params = self.table[index]
        hplus.approximant = approximant
        hcross.approximant = approximant

        # Add sigmasq as a method of this instance
        hplus.sigmasq = types.MethodType(sigma_cached, hplus)
        hplus._sigmasq = {}
        hcross.sigmasq = types.MethodType(sigma_cached, hcross)
        hcross._sigmasq = {}

        return hplus, hcross
