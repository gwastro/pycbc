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
from pycbc.types import TimeSeries, FrequencySeries, Array, float32, float64, complex_same_precision_as
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

def phase_from_frequencyseries(htilde, remove_start_phase=True):
    """Returns the phase from the given frequency-domain waveform. This assumes
    that the waveform has been sampled finely enough that the phase cannot
    change by more than pi radians between each step.

    Parameters
    ----------
    htilde : FrequencySeries
        The waveform to get the phase for; must be a complex frequency series.
    remove_start_phase : {True, bool}
        Subtract the initial phase before returning.

    Returns
    -------
    FrequencySeries
        The phase of the waveform as a function of frequency.
    """
    p = numpy.unwrap(numpy.angle(htilde.data))
    if remove_start_phase:
        p += -p[0]    
    return FrequencySeries(p, delta_f=htilde.delta_f, epoch=htilde.epoch,
        copy=False)

def amplitude_from_frequencyseries(htilde):
    """Returns the amplitude of the given frequency-domain waveform as a
    FrequencySeries.

    Parameters
    ----------
    htilde : FrequencySeries
        The waveform to get the amplitude of.

    Returns
    -------
    FrequencySeries
        The amplitude of the waveform as a function of frequency.
    """
    amp = abs(htilde.data)
    return FrequencySeries(amp, delta_f=htilde.delta_f, epoch=htilde.epoch,
        copy=False)

def time_from_frequencyseries(htilde, sample_frequencies=None):
    """Computes time as a function of frequency from the given
    frequency-domain waveform. This assumes the stationary phase
    approximation. Any frequencies lower than the first non-zero value in
    htilde are assigned the time at the first non-zero value. Times for any
    frequencies above the next-to-last non-zero value in htilde will be
    assigned the time of the next-to-last non-zero value.

    Parameters
    ----------
    htilde : FrequencySeries
        The waveform to get the time evolution of; must be complex.
    sample_frequencies : {None, array}
        The frequencies at which the waveform is sampled. If None, will
        retrieve from ``htilde.sample_frequencies``.

    Returns
    -------
    FrequencySeries
        The time evolution of the waveform as a function of frequency.
    """
    if sample_frequencies is None:
        sample_frequencies = htilde.sample_frequencies.numpy()
    phase = phase_from_frequencyseries(htilde).data
    time = -numpy.diff(phase) / (2.*numpy.pi*numpy.diff(sample_frequencies))
    nzidx = numpy.nonzero(abs(htilde.data))[0]
    kmin, kmax = nzidx[0], nzidx[-2]
    time[:kmin] = time[kmin]
    time[kmax:] = time[kmax]
    return FrequencySeries(time, delta_f=htilde.delta_f, epoch=htilde.epoch,
        copy=False)

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
    p = numpy.unwrap(numpy.arctan2(h_cross.data, h_plus.data))
    if remove_start_phase:
        p += -p[0]    
    return TimeSeries(p, delta_t=h_plus.delta_t, epoch=h_plus.start_time,
        copy=False)

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

def apply_fd_time_shift(htilde, shifttime, fseries=None, copy=True):
    """Shifts a frequency domain waveform in time. The shift applied is
    shiftime - htilde.epoch.

    Parameters
    ----------
    htilde : FrequencySeries
        The waveform frequency series.
    shifttime : float
        The time to shift the frequency series to.
    fseries : {None, numpy array}
        The frequencies of each element in the the FrequencySeries. If None,
        will use htilde.sample_frequencies. Note: providing a frequency series
        can reduce the exectution time of this function by as much as a 1/2.
    copy : {True, bool}
        Make a copy of htilde before applying the time shift. If False, the time
        shift will be applied to htilde's data.

    Returns
    -------
    FrequencySeries
        A frequency series with the waveform shifted to the new time. If makecopy
        is True, will be a new frequency series; if makecopy is False, will be
        the same as htilde.
    """
    dt = float(shifttime - htilde.epoch)
    if fseries is None:
        fseries = htilde.sample_frequencies.numpy()
    shift = Array(numpy.exp(-2j*numpy.pi*dt*fseries), 
                dtype=complex_same_precision_as(htilde))
    if copy:
        htilde = 1. * htilde
    htilde *= shift
    return htilde
