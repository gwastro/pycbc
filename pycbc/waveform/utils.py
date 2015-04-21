# Copyright (C) 2013  Alex Nitz
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
"""This module contains convenience utilities for manipulating waveforms
"""
from pycbc.types import TimeSeries, float32, float64
import lal
import lalsimulation as sim
from math import frexp
import numpy
import copy

def ceilpow2(n):
    """convenience function to determine a power-of-2 upper frequency limit"""
    signif,exponent = frexp(n)
    if (signif < 0):
        return 1;
    if (signif == 0.5):
        exponent -= 1;
    return (1) << exponent;

def unwrap_phase(vec, discont, offset):
    """Return a new vector that is unwrapped.

    Return a new vector that is free from discontinuities caused by 
    periodic boundary conditions.
    Where the vector jumps by a value greater than or equal to `discont`,
    it is incremented by the offset value. This function is intended to
    remove boundaries placed by using a cyclic function.

    Parameters
    ----------
    vec : array_like
        Vectors that can be converted to an array as defined by numpy.
        PyCBC types, numpy types, lists, etc can be provided.
    discont : float
        Minimum size of discontinuity that triggers adding an offset to the
        data.  Due to precision on many functions consider a value smaller
        than the maximum (e.g. smaller than 2*pi for a function of an angular
        variable).
    offset : float
        Quantity to be added to/subtracted from the vector at each
        discontinuity.

    Returns
    -------
    nvec : array_like
        New unwrapped vector, same type as vec.
    """
    nvec = copy.deepcopy(vec)
    total_offset = 0
    pval = None
    index = 0
    for val in vec:
        if pval is None:
            pass
        elif val-pval>=discont:
            total_offset -= offset
        elif pval-val>=discont:
            total_offset += offset

        nvec[index] += total_offset
        index += 1
        pval = val

    return nvec

def phase_from_polarizations(h_plus, h_cross, remove_start_phase=True):
    """Return gravitational wave phase

    Return the gravitation-wave phase from the h_plus and h_cross 
    polarizations of the waveform. The returned phase is always
    positive and increasing with an initial phase of 0.

    Parameters
    ----------
    h_plus : TimeSeries
        An PyCBC TmeSeries vector that contains the plus polarization of the
        gravitational waveform.
    h_cross : TimeSeries
        A PyCBC TmeSeries vector that contains the cross polarization of the
        gravitational waveform.

    Returns
    -------
    GWPhase : TimeSeries
        A TimeSeries containing the gravitational wave phase.

    Examples
    --------s
    >>> from pycbc.waveform import get_td_waveform, phase_from_polarizations
    >>> hp, hc = get_td_waveform(approximant="TaylorT4", mass1=10, mass2=10, 
                         f_lower=30, delta_t=1.0/4096)
    >>> phase = phase_from_polarizations(hp, hc)

    """
    p_wrapped = numpy.arctan2(h_plus, h_cross)
    p = unwrap_phase(p_wrapped, 2*lal.PI*0.7, lal.PI*2)
    if remove_start_phase:
        p += -p[0]    
    return TimeSeries(abs(p), delta_t=h_plus.delta_t, epoch=h_plus.start_time)

def amplitude_from_polarizations(h_plus, h_cross):
    """Return gravitational wave amplitude

    Return the gravitation-wave amplitude from the h_plus and h_cross 
    polarizations of the waveform. 

    Parameters
    ----------
    h_plus : TimeSeries
        An PyCBC TmeSeries vector that contains the plus polarization of the
        gravitational waveform.
    h_cross : TimeSeries
        A PyCBC TmeSeries vector that contains the cross polarization of the
        gravitational waveform.

    Returns
    -------
    GWAmplitude : TimeSeries
        A TimeSeries containing the gravitational wave amplitude.

    Examples
    --------
    >>> from pycbc.waveform import get_td_waveform, phase_from_polarizations
    >>> hp, hc = get_td_waveform(approximant="TaylorT4", mass1=10, mass2=10, 
                         f_lower=30, delta_t=1.0/4096)
    >>> amp = amplitude_from_polarizations(hp, hc)

    """
    amp = (h_plus.squared_norm() + h_cross.squared_norm()) ** (0.5)
    return TimeSeries(amp, delta_t=h_plus.delta_t, epoch=h_plus.start_time)

def frequency_from_polarizations(h_plus, h_cross):
    """Return gravitational wave frequency

    Return the gravitation-wave frequency as a function of time
    from the h_plus and h_cross polarizations of the waveform. 
    It is 1 bin shorter than the input vectors and the sample times
    are advanced half a bin.

    Parameters
    ----------
    h_plus : TimeSeries
        A PyCBC TimeSeries vector that contains the plus polarization of the
        gravitational waveform.
    h_cross : TimeSeries
        A PyCBC TimeSeries vector that contains the cross polarization of the
        gravitational waveform.

    Returns
    -------
    GWFrequency : TimeSeries
        A TimeSeries containing the gravitational wave frequency as a function
        of time. 

    Examples
    --------
    >>> from pycbc.waveform import get_td_waveform, phase_from_polarizations
    >>> hp, hc = get_td_waveform(approximant="TaylorT4", mass1=10, mass2=10, 
                         f_lower=30, delta_t=1.0/4096)
    >>> freq = frequency_from_polarizations(hp, hc)

    """
    phase = phase_from_polarizations(h_plus, h_cross)
    freq = numpy.diff(phase) / ( 2 * lal.PI * phase.delta_t )
    start_time = phase.start_time + phase.delta_t / 2
    return TimeSeries(freq, delta_t=phase.delta_t, epoch=start_time)

# map between tapering string in sim_inspiral table or inspiral 
# code option and lalsimulation constants
taper_map = {
    'TAPER_NONE'    : None,
    'TAPER_START'   : sim.SIM_INSPIRAL_TAPER_START,
    'start'         : sim.SIM_INSPIRAL_TAPER_START,
    'TAPER_END'     : sim.SIM_INSPIRAL_TAPER_END,
    'end'           : sim.SIM_INSPIRAL_TAPER_END,
    'TAPER_STARTEND': sim.SIM_INSPIRAL_TAPER_STARTEND,
    'startend'      : sim.SIM_INSPIRAL_TAPER_STARTEND
}

taper_func_map = {
    numpy.dtype(float32): sim.SimInspiralREAL4WaveTaper,
    numpy.dtype(float64): sim.SimInspiralREAL8WaveTaper
}

def taper_timeseries(tsdata, tapermethod=None, return_lal=False):
    """
    Taper either or both ends of a time series using wrapped 
    LALSimulation functions

    Parameters
    ----------
    tsdata : TimeSeries
        Series to be tapered, dtype must be either float32 or float64
    tapermethod : string
        Should be one of ('TAPER_NONE', 'TAPER_START', 'TAPER_END',
        'TAPER_STARTEND', 'start', 'end', 'startend') - NB 'TAPER_NONE' will
        not change the series!
    return_lal : Boolean
        If True, return a wrapped LAL time series object, else return a 
        PyCBC time series.
    """
    if tapermethod is None:
        raise ValueError("Must specify a tapering method (function was called"
                         "with tapermethod=None)")
    if tapermethod not in taper_map.keys():
        raise ValueError("Tapering method %s is not known! Valid methods: "
                         + taper_map.keys() % (tapermethod))
    if not tsdata.dtype in (float32, float64):
        raise TypeError("Strain dtype must be float32 or float64, not "
                    + str(tsdata.dtype))
    taper_func = taper_func_map[tsdata.dtype]
    # make a LAL TimeSeries to pass to the LALSim function
    ts_lal = tsdata.astype(tsdata.dtype).lal()
    if taper_map[tapermethod] is not None:
        taper_func(ts_lal.data, taper_map[tapermethod])
    if return_lal == True:
        return ts_lal
    else:
        return TimeSeries(ts_lal.data.data[:], delta_t=ts_lal.deltaT,
                          epoch=ts_lal.epoch)

