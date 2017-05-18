# Copyright (C) 2012  Josh Willis, Andrew Miller
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
This package provides a front-end to various fast Fourier transform
implementations within PyCBC.
"""

from pycbc.types import Array as _Array
from pycbc.types import TimeSeries as _TimeSeries
from pycbc.types import FrequencySeries as _FrequencySeries

# The following helper function is in this top-level module because it
# is used by the scheme-dependent files to write their version of the
# _available_backends() function. It cannot go in backend_support as
# that woulc cause circular imports

def _list_available(possible_list, possible_dict):
    # It possibly is strange that we have both a list and a dict.
    # The reason for this is that the name the user specfies for a
    # backend, e.g. 'numpy', may be something that the fft submodule
    # cannot be called, so we need a dict mapping those names to the
    # actual names of the modules.  However when we iterate, it must
    # be in the defined order, because that represents a preference
    # for which backends are most likely to be preferrable.  As a
    # dict is unordered, we cannot simply use its keys for this purpose.
    available_list = []
    available_dict = {}
    for backend in possible_list:
        try:
            mod = __import__('pycbc.fft.' + possible_dict[backend], fromlist = ['pycbc.fft'])
            available_dict.update({backend:mod})
            available_list.append(backend)
        except (ImportError, OSError):
            pass
    return available_list, available_dict

# The main purpose of the top-level module is to present a
# uniform interface for a forward and reverse FFT, independent of
# the underlying backend.  We perform sanity checking here, at the
# top-level, and then don't worry about it in submodules. To
# facilitate this checking, we define dicts mapping the numpy dtype
# to the corresponding precisions and types.

def _check_fft_args(invec, outvec):
    if not isinstance(invec,_Array):
        raise TypeError("Input is not a PyCBC Array")
    if not isinstance(outvec,_Array):
        raise TypeError("Output is not a PyCBC Array")

    if isinstance(invec,_TimeSeries) and not isinstance(
        outvec,_FrequencySeries):
        raise TypeError(
            "When input is TimeSeries output must be FrequencySeries")
    if isinstance(outvec,_TimeSeries) and not isinstance(
        invec,_FrequencySeries):
        raise TypeError(
            "When output is TimeSeries input must be FrequencySeries")
    if isinstance(invec,_FrequencySeries) and not isinstance(
        outvec,_TimeSeries):
        raise TypeError(
            "When input is FrequencySeries output must be TimeSeries")
    if isinstance(outvec,_FrequencySeries) and not isinstance(
        invec,_TimeSeries):
        raise TypeError(
            "When output is FrequencySeries input must be TimeSeries")

    iprec = invec.precision
    oprec = outvec.precision
    if iprec != oprec:
        raise ValueError("Input and output precisions must agree")

    itype = invec.kind
    otype = outvec.kind
    return [iprec,itype,otype]

def _check_fwd_args(invec, itype, outvec, otype, nbatch, size):
    ilen = len(invec)
    olen = len(outvec)
    if nbatch < 1:
        raise ValueError("nbatch must be >= 1")
    if (nbatch > 1) and size is not None:
        raise ValueError("When nbatch > 1, size cannot be 'None'")
    if size is None:
        size = ilen
    inplace = (invec.ptr == outvec.ptr)
    if (ilen % nbatch) != 0:
        raise ValueError("Input length must be divisible by nbatch")
    if (olen % nbatch) != 0:
        raise ValueError("Output length must be divisible by nbatch")
    if itype == 'complex' and otype == 'complex':
        if (ilen/nbatch) != size:
            raise ValueError("For C2C FFT, len(invec) must be nbatch*size")
        if (olen/nbatch) != size:
            raise ValueError("For C2C FFT, len(outvec) must be nbatch*size")
    elif itype == 'real' and otype == 'complex':
        if (olen/nbatch) != (size/2 + 1):
            raise ValueError("For R2C FFT, len(outvec) must be nbatch*(size/2 + 1)")
        if inplace:
            if (ilen/nbatch) != 2*(size/2 + 1):
                raise ValueError("For R2C in-place FFT, len(invec) must be nbatch*2*(size/2+1)")
        else:
            if (ilen/nbatch) != size:
                raise ValueError("For R2C out-of-place FFT, len(invec) must be nbatch*size")
    else:
        raise ValueError("Inconsistent dtypes for forward FFT")

def _check_inv_args(invec, itype, outvec, otype, nbatch, size):
    ilen = len(invec)
    olen = len(outvec)
    if nbatch < 1:
        raise ValueError("nbatch must be >= 1")
    if (nbatch > 1) and size is None:
        raise ValueError("When nbatch > 1, size cannot be 'None'")
    if size is None:
        size = olen
    inplace = (invec.ptr == outvec.ptr)
    if (ilen % nbatch) != 0:
        raise ValueError("Input length must be divisible by nbatch")
    if (olen % nbatch) != 0:
        raise ValueError("Output length must be divisible by nbatch")
    if itype == 'complex' and otype == 'complex':
        if (ilen/nbatch) != size:
            raise ValueError("For C2C IFFT, len(invec) must be nbatch*size")
        if (olen/nbatch) != size:
            raise ValueError("For C2C IFFT, len(outvec) must be nbatch*size")
    elif itype == 'complex' and otype == 'real':
        if (ilen/nbatch) != (size/2 + 1):
            raise ValueError("For C2R IFFT, len(invec) must be nbatch*(size/2 + 1)")
        if inplace:
            if (olen/nbatch) != 2*(size/2 + 1):
                raise ValueError("For C2R in-place IFFT, len(outvec) must be nbatch*2*(size/2+1)")
        else:
            if (olen/nbatch) != size:
                raise ValueError("For C2R out-of-place IFFT, len(outvec) must be nbatch*size")

# The class-based approach requires the following:


# The classes below should serve as the parent for all schemed classes.
# In part, these classes should serve as the location for
# all documentation of the class and its methods, though that is not
# yet implemented.  Perhaps something along the lines of:
#
#    http://stackoverflow.com/questions/2025562/inherit-docstrings-in-python-class-inheritance
#
# will work? Is there a better way?
#
# Unlike some other places within PyCBC, however, the __init__ method of these classes do
# nontrivial work and hence should be called inside the __init__ method of all child classes,
# before anything else.

class _BaseFFT(object):
    def __init__(self, invec, outvec, nbatch, size):
        _, itype, otype = _check_fft_args(invec, outvec)
        _check_fwd_args(invec, itype, outvec, otype, nbatch, size)
        self.forward = True
        self.invec = invec
        self.outvec = outvec
        self.inplace = (self.invec.ptr == self.outvec.ptr)
        self.nbatch = nbatch
        if nbatch > 1:
            self.size = size
        else:
            self.size = len(invec)
        # Whether we are complex-to-complex or real-to-complex is determined
        # by itype:
        if itype == 'complex':
            # Complex-to-complex case:
            self.idist = self.size
            self.odist = self.size
        else:
            # Real-to-complex case:
            self.odist = (self.size/2 + 1)
            if self.inplace:
                self.idist = 2*(self.size/2 + 1)
            else:
                self.idist = self.size

        # For a forward FFT, the length of the *input* vector is the length
        # we should divide by, whether C2C or R2HC transform
        if isinstance(self.invec, _TimeSeries):
            self.outvec._epoch = self.invec._epoch
            self.outvec._delta_f = 1.0/(self.invec._delta_t * len(self.invec))
            self.scale = self.invec._delta_t
        elif isinstance(self.invec, _FrequencySeries):
            self.outvec._epoch = self.invec._epoch
            self.outvec._delta_t = 1.0/(self.invec._delta_f * len(self.invec))
            self.scale = self.invec._delta_f

    def execute(self):
        """
        Compute the (forward) FFT of the input vector specified at object
        instantiation, putting the output into the output vector specified
        at objet instantiation. The intention is that this method should
        be called many times, with the contents of the input vector
        changing between invocations, but not the locations in memory or
        length of either input or output vector.

        *Unlike* the function based API, the class based API does NOT rescale
        its output by the input vector's delta_t (when input is a TimeSeries)
        or delta_f (when input is a FrequencySeries).
        """
        pass

class _BaseIFFT(object):
    def __init__(self, invec, outvec, nbatch, size):
        _, itype, otype = _check_fft_args(invec, outvec)
        _check_inv_args(invec, itype, outvec, otype, nbatch, size)
        self.forward = False
        self.invec = invec
        self.outvec = outvec
        self.inplace = (self.invec.ptr == self.outvec.ptr)
        self.nbatch = nbatch
        if nbatch > 1:
            self.size = size
        else:
            self.size = len(outvec)
        # Whether we are complex-to-complex or complex-to-real is determined
        # by otype:
        if otype == 'complex':
            # Complex-to-complex case:
            self.idist = self.size
            self.odist = self.size
        else:
            # Complex-to-real case:
            self.idist = (self.size/2 + 1)
            if self.inplace:
                self.odist = 2*(self.size/2 + 1)
            else:
                self.odist = self.size

        # For an inverse FFT, the length of the *output* vector is the length
        # we should divide by, whether C2C or HC2R transform
        if isinstance(self.invec, _TimeSeries):
            self.outvec._epoch = self.invec._epoch
            self.outvec._delta_f = 1.0/(self.invec._delta_t * len(self.outvec))
            self.scale = self.invec._delta_t
        elif isinstance(self.invec, _FrequencySeries):
            self.outvec._epoch = self.invec._epoch
            self.outvec._delta_t = 1.0/(self.invec._delta_f * len(self.outvec))
            self.scale = self.invec._delta_f

    def execute(self):
        """
        Compute the (backward) FFT of the input vector specified at object
        instantiation, putting the output into the output vector specified
        at objet instantiation. The intention is that this method should
        be called many times, with the contents of the input vector
        changing between invocations, but not the locations in memory or
        length of either input or output vector.

        *Unlike* the function based API, the class based API does NOT rescale
        its output by the input vector's delta_t (when input is a TimeSeries)
        or delta_f (when input is a FrequencySeries).
        """
        pass

