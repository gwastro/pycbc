# Copyright (C) 2012  Josh Willis
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
This module provides the SwigLAL backend of the fast Fourier transform
for the PyCBC package.
"""

import lal
from lal import CreateForwardCOMPLEX16FFTPlan as _CreateForwardCOMPLEX16FFTPlan
from lal import CreateForwardCOMPLEX8FFTPlan as _CreateForwardCOMPLEX8FFTPlan
from lal import CreateForwardREAL8FFTPlan as _CreateForwardREAL8FFTPlan
from lal import CreateForwardREAL4FFTPlan as _CreateForwardREAL4FFTPlan
from lal import CreateReverseCOMPLEX16FFTPlan as _CreateReverseCOMPLEX16FFTPlan
from lal import CreateReverseCOMPLEX8FFTPlan as _CreateReverseCOMPLEX8FFTPlan
from lal import CreateReverseREAL8FFTPlan as _CreateReverseREAL8FFTPlan
from lal import CreateReverseREAL4FFTPlan as _CreateReverseREAL4FFTPlan

from lal import COMPLEX16VectorFFT as _COMPLEX16VectorFFT
from lal import COMPLEX8VectorFFT as _COMPLEX8VectorFFT
from lal import REAL8ForwardFFT as _REAL8ForwardFFT
from lal import REAL4ForwardFFT as _REAL4ForwardFFT
from lal import REAL8ReverseFFT as _REAL8ReverseFFT
from lal import REAL4ReverseFFT as _REAL4ReverseFFT


# Measure level.  By default 1, which does some but not much planning,
# but we provide functions to read and set it

_default_measurelvl = 1
def get_measure_level():
    return _default_measurelvl

def set_measure_level(mlvl):
    global _default_measurelvl
    if mlvl not in (0,1,2,3):
        raise ValueError("Measure level can only be one of 0, 1, 2, or 3")
    _default_measurelvl = mlvl

_forward_plans = {}
_reverse_plans = {}
_forward_plan_fn_dict = {('double','complex','complex') : _CreateForwardCOMPLEX16FFTPlan,
                    ('single','complex','complex') : _CreateForwardCOMPLEX8FFTPlan,
                    ('double','real','complex') : _CreateForwardREAL8FFTPlan,
                    ('single','real','complex') : _CreateForwardREAL4FFTPlan}
_reverse_plan_fn_dict = {('double','complex','complex') : _CreateReverseCOMPLEX16FFTPlan,
                    ('single','complex','complex') : _CreateReverseCOMPLEX8FFTPlan,
                    ('double','complex','real') : _CreateReverseREAL8FFTPlan,
                    ('single','complex','real') : _CreateReverseREAL4FFTPlan}
_forward_fft_fn_dict = {('double','complex','complex') : _COMPLEX16VectorFFT,
                    ('single','complex','complex') : _COMPLEX8VectorFFT,
                    ('double','real','complex') : _REAL8ForwardFFT,
                    ('single','real','complex') : _REAL4ForwardFFT}
_reverse_fft_fn_dict = {('double','complex','complex') : _COMPLEX16VectorFFT,
                    ('single','complex','complex') : _COMPLEX8VectorFFT,
                    ('double','complex','real') : _REAL8ReverseFFT,
                    ('single','complex','real') : _REAL4ReverseFFT}

_newlalvec = {('single','real') : lal.CreateREAL4Vector,
              ('double','real') : lal.CreateREAL8Vector,
              ('single','complex') : lal.CreateCOMPLEX8Vector,
              ('double','complex') : lal.CreateCOMPLEX16Vector}

def _get_fwd_plan(prec,itype,otype,inlen):
    try:
        theplan, mlvl = _forward_plans[(prec,itype,otype,inlen)]
        # If the measure level at which we created the plan is less than the current
        # global measure level, we must regenerate the plan:
        if (mlvl < _default_measurelvl):
            theplan = _forward_plan_fn_dict[(prec,itype,otype)](inlen,_default_measurelvl)
            _forward_plans.update({(prec,itype,otype,inlen) : (theplan,_default_measurelvl)})
    except KeyError:
        theplan = _forward_plan_fn_dict[(prec,itype,otype)](inlen,_default_measurelvl)
        _forward_plans.update({(prec,itype,otype,inlen) : (theplan,_default_measurelvl)})
    return theplan

def _get_inv_plan(prec,itype,otype,outlen):
    try:
        theplan, mlvl = _reverse_plans[(prec,itype,otype,outlen)]
        # If the measure level at which we created the plan is less than the current
        # global measure level, we must regenerate the plan:
        if (mlvl < _default_measurelvl):
            theplan = _reverse_plan_fn_dict[(prec,itype,otype)](outlen,_default_measurelvl)
            _reverse_plans.update({(prec,itype,otype,outlen) : (theplan,_default_measurelvl)})
    except KeyError:
        theplan = _reverse_plan_fn_dict[(prec,itype,otype)](outlen,_default_measurelvl)
        _reverse_plans.update({(prec,itype,otype,outlen) : (theplan,_default_measurelvl)})
    return theplan


def fft(invec,outvec,prec,itype,otype):
    if invec.ptr == outvec.ptr:
        raise NotImplementedError("lal backend of pycbc.fft does not support in-place transforms.")
    theplan = _get_fwd_plan(prec,itype,otype,len(invec))
    inlal = _newlalvec[(prec,itype)](len(invec))
    outlal = _newlalvec[(prec,otype)](len(outvec))
    inlal.data[:] = invec.numpy()
    _forward_fft_fn_dict[(prec,itype,otype)](outlal,inlal,theplan)
    outvec._data[:] = outlal.data
    del inlal
    del outlal

def ifft(invec,outvec,prec,itype,otype):
    if invec.ptr == outvec.ptr:
        raise NotImplementedError("lal backend of pycbc.fft does not support in-place transforms.")
    theplan = _get_inv_plan(prec,itype,otype,len(outvec))
    inlal = _newlalvec[(prec,itype)](len(invec))
    outlal = _newlalvec[(prec,otype)](len(outvec))
    inlal.data[:] = invec.numpy()
    _reverse_fft_fn_dict[(prec,itype,otype)](outlal,inlal,theplan)
    outvec._data[:] = outlal.data
    del inlal
    del outlal

def insert_fft_options(optgroup):
    """
    Inserts the options that affect the behavior of this backend

    Parameters
    ----------
    optgroup: fft_option
       OptionParser argument group whose options are extended
    """
    optgroup.add_argument("--lalfft-measure-level",
                      help="Determines the measure level used in planning "
                           "LAL-wrapped FFTs; allowed values are: " + str([0,1,2,3]),
                      type=int, default=_default_measurelvl)

def verify_fft_options(opt,parser):
    """Parses the FFT options and verifies that they are
       reasonable.


    Parameters
    ----------
    opt : object
        Result of parsing the CLI with OptionParser, or any object with the
        required attributes.
    parser : object
        OptionParser instance.
    """
    if opt.lalfft_measure_level not in [0,1,2,3]:
        parser.error("{0} is not a valid lalfft measure level.".format(opt.lalfft_measure_level))

def from_cli(opt):
    set_measure_level(opt.lalfft_measure_level)
