# Copyright (C) 2016 Miriam Cabero Mueller
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
"""Generate ringdown templates in the frequency domain.
"""

import numpy, lal
import lalsimulation as lalsim
from pycbc.types import TimeSeries, FrequencySeries, complex128, zeros
from pycbc.waveform.waveform import get_obj_attrs

default_qnm_args = {'t_0':0, 'phi':0, 'amp':1}
qnm_required_args = ['f_0', 'tau']
lm_required_args = ['Mfinal','Sfinal','l','m','nmodes']

max_freq = 16384.
pi = numpy.pi
two_pi = 2 * numpy.pi
pi_sq = numpy.pi * numpy.pi

# Input parameters ############################################################
def props(obj, required, **kwargs):
    """ Return a dictionary built from the combination of defaults, kwargs,
    and the attributes of the given object.
    """
    # Get the attributes of the template object
    pr = get_obj_attrs(obj)

    # Get the parameters to generate the waveform
    # Note that keyword arguments override values in the template object
    input_params = default_qnm_args.copy()
    input_params.update(pr)
    input_params.update(kwargs)

    # Check if the required arguments are given
    for arg in required:
        if arg not in input_params:
            raise ValueError('Please provide ' + str(arg))

    return input_params

def lm_amps_phases(kwargs)
    """From an input_params dictionary, return the amplitudes and phases of
    each overtone of an lm mode, checking that all of them are given.
    """
    amps, phis = numpy.zeros(nmodes), numpy.zeros(nmodes)
    for n in range(kwargs['nmodes']):
        try:
            amps[n] = kwargs['amp_%d' %n]
        except KeyError:
            raise ValueError('amp_%d is required' %n)
        try:
            phis[n] = kwargs['phi_%d' %n]
        except KeyError:
            raise ValueError('phi_%d is required' %n)

    return amps, phis

# Functions to obtain t_final and f_final #####################################
def qnm_time_decay(tau, decay):
    """Return the time at which the amplitude of the 
    ringdown falls to decay of the peak amplitude.

    Parameters
    ----------
    tau : float
        The damping time of the sinusoid.
    decay: float
        The fraction of the peak amplitude.

    Returns
    -------
    t_decay: float
        The time at which the amplitude of the time-domain
        ringdown falls to decay of the peak amplitude.
    """

    return - tau * numpy.log(decay)

def qnm_freq_decay(f_0, tau, decay):
    """Return the frequency at which the amplitude of the 
    ringdown falls to decay of the peak amplitude.

    Parameters
    ----------
    f_0 : float
        The ringdown-frequency, which gives the peak amplitude.
    tau : float
        The damping time of the sinusoid.
    decay: float
        The fraction of the peak amplitude.

    Returns
    -------
    f_decay: float
        The frequency at which the amplitude of the frequency-domain
        ringdown falls to decay of the peak amplitude.
    """

    q_0 = pi * f_0 * tau
    alpha = 1. / decay
    alpha_sq = 1. / decay / decay

    # Expression obtained analytically under the assumption
    # that 1./alpha_sq, q_0^2 >> 1
    q_sq = (alpha_sq + 4*q_0*q_0 + alpha*numpy.sqrt(alpha_sq + 16*q_0*q_0)) / 4.
    return numpy.sqrt(q_sq) / pi / tau

# Functions to generate ringdown waveforms ####################################
def get_td_qnm(template=None, **kwargs):
    """Return a time domain damped sinusoid.

    Parameters
    ----------
    template: object
        An object that has attached properties. This can be used to substitute
        for keyword arguments. A common example would be a row in an xml table.
    f_0 : float
        The ringdown-frequency.
    tau : float
        The damping time of the sinusoid.
    phi : {0, float}, optional
        The initial phase of the ringdown.
    amp : {1, float}, optional
        The amplitude of the ringdown (constant for now).
    delta_t : {None, float}, optional
        The time step used to generate the ringdown.
        If None, it will be set to the inverse of the frequency at which the
        amplitude is 1/1000 of the peak amplitude.
    t_final : {None, float}, optional
        The ending time of the output time series.
        If None, it will be set to the time at which the amplitude is 
        1/1000 of the peak amplitude.

    Returns
    -------
    hplus: TimeSeries
        The plus phase of the ringdown in time domain.
    hcross: TimeSeries
        The cross phase of the ringdown in time domain.
    """

    input_params = props(template, qnm_required_args, **kwargs)

    # the following args have defaults, and so will be populated
    phi = input_params.pop('phi')
    amp = input_params.pop('amp')
    # the following may not be in input_params
    delta_t = input_params.pop('delta_t', None)
    t_final = input_params.pop('t_final', None)

    if delta_t is None:
        delta_t = 1. / qnm_freq_decay(f_0, tau, 1./1000)
    if t_final is None:
        t_final = qnm_time_decay(tau, 1./1000)
    kmax = int(t_final / delta_t) + 1

    times = numpy.arange(kmax)*delta_t

    hp = amp * numpy.exp(-times/tau) * numpy.cos(two_pi*f_0*times + phi)
    hc = amp * numpy.exp(-times/tau) * numpy.sin(two_pi*f_0*times + phi)

    hplus = TimeSeries(zeros(kmax), delta_t=delta_t)
    hcross = TimeSeries(zeros(kmax), delta_t=delta_t)
    hplus.data[:kmax] = hp
    hcross.data[:kmax] = hc

    return hplus, hcross

def get_fd_qnm(template=None, **kwargs):
    """Return a frequency domain damped sinusoid.

    Parameters
    ----------
    template: object
        An object that has attached properties. This can be used to substitute
        for keyword arguments. A common example would be a row in an xml table.
    f_0 : float
        The ringdown-frequency.
    tau : float
        The damping time of the sinusoid.
    t_0 :  {0, float}, optional
        The starting time of the ringdown.
    phi : {0, float}, optional
        The initial phase of the ringdown.
    amp : {1, float}, optional
        The amplitude of the ringdown (constant for now).
    delta_f : {None, float}, optional
        The frequency step used to generate the ringdown.
        If None, it will be set to the inverse of the time at which the
        amplitude is 1/1000 of the peak amplitude.
    f_lower: {None, float}, optional
        The starting frequency of the output frequency series.
        If None, it will be set to delta_f.
    f_final : {None, float}, optional
        The ending frequency of the output frequency series.
        If None, it will be set to the frequency at which the amplitude is 
        1/1000 of the peak amplitude.

    Returns
    -------
    hplustilde: FrequencySeries
        The plus phase of the ringdown in frequency domain.
    hcrosstilde: FrequencySeries
        The cross phase of the ringdown in frequency domain.
    """

    input_params = props(template, qnm_required_args, **kwargs)

    # get optional args
    # the following have defaults, and so will be populated
    t_0 = input_params.pop('t_0')
    phi = input_params.pop('phi')
    amp = input_params.pop('amp')
    # the following may not be in input_params
    delta_f = input_params.pop('delta_f', None)
    f_lower = input_params.pop('f_lower', None)
    f_final = input_params.pop('f_final', None)

    if delta_f is None:
        delta_f = 1. / qnm_time_decay(tau, 1./1000)
    if f_lower is None:
        f_lower = delta_f
        kmin = 0
    else:
        kmin = int(f_lower / delta_f)
    if f_final is None:
        f_final = qnm_freq_decay(f_0, tau, 1./1000)
    if f_final > max_freq:
            f_final = max_freq
    kmax = int(f_final / delta_f) + 1

    freqs = numpy.arange(kmin, kmax)*delta_f

    denominator = 1 + (4j * pi * freqs * tau)
                    - (4 * pi_sq * ( freqs*freqs - f_0*f_0) * tau*tau)
    norm = amp * tau / denominator
    if t_0 != 0:
        time_shift = numpy.exp(-1j * two_pi * freqs * t_0) 
        norm *= time_shift

    # Anallytical expression for the Fourier transform of the ringdown (damped sinusoid)
    hp_tilde = norm * ( (1 + 2j * pi * freqs * tau) * numpy.cos(phi)
                               - two_pi * f_0 * tau * numpy.sin(phi) )
    hc_tilde = norm * ( (1 + 2j * pi * freqs * tau) * numpy.sin(phi)
                               + two_pi * f_0 * tau * numpy.cos(phi) )

    hplustilde = FrequencySeries(zeros(kmax, dtype=complex128), delta_f=delta_f)
    hcrosstilde = FrequencySeries(zeros(kmax, dtype=complex128), delta_f=delta_f)
    hplustilde.data[kmin:kmax] = hp_tilde
    hcrosstilde.data[kmin:kmax] = hc_tilde

    return hplustilde, hcrosstilde

def get_lm_f0tau(mass, spin, l, m, nmodes):
    """Return the f_0 and the tau of each overtone for a given lm mode 
    """

    qnmfreq = lal.CreateCOMPLEX16Vector(nmodes)
    lalsim.SimIMREOBGenerateQNMFreqV2fromFinal(qnmfreq, Mfinal, Sfinal, l, m, nmodes)

    f_0 = [qnmfreq.data[n].real / (2 * pi) for n in range(nmodes)]
    tau = [1. / qnmfreq.data[n].imag for n in range(nmodes)]

    return f_0, tau

def get_td_lm(template=None, **kwargs):
    """Return frequency domain lm mode with the given number of overtones.
    Parameters
    ----------
    template: object
        An object that has attached properties. This can be used to substitute
        for keyword arguments. A common example would be a row in an xml table.
    Mfinal : float
        Mass of the final black hole.
    Sfinal : float
        Spin of the final black hole.
    l : int
        l mode (lm modes available: 22, 21, 33, 44, 55).
    m : int
        m mode (lm modes available: 22, 21, 33, 44, 55).
    nmodes: int
        Number of overtones desired (maximum n=8)
    amp_n : float
        Amplitude of the n overtone, as many as the number of nmodes.
    phi_n : float
        Phase of the n overtone, as many as the number of nmodes.
    delta_t : {None, float}, optional
        The time step used to generate the ringdown.
        If None, it will be set to the inverse of the frequency at which the
        amplitude is 1/1000 of the peak amplitude (the minimum of all modes).
    t_final : {None, float}, optional
        The ending time of the output time series.
        If None, it will be set to the time at which the amplitude is 
        1/1000 of the peak amplitude (the maximum of all modes).
    Returns
    -------
    hplustilde: FrequencySeries
        The plus phase of a lm mode with overtones (n) in frequency domain.
    hcrosstilde: FrequencySeries
        The cross phase of a lm mode with overtones (n) in frequency domain.
    """

    input_params = props(template, lm_required_args, **kwargs)
    amps, phis = lm_amps_phases(input_params)
    # The following may not be in input_params
    delta_t = input_params.pop('delta_t', None)
    t_final = input_params.pop('t_final', None)

    f_0, tau = get_lm_f0tau(Mfinal, Sfinal, l, m, nmodes)

    if delta_t is None:
        dt = [1. / qnm_freq_decay(f_0[n], tau[n], 1./1000) for n in range(nmodes)]
        delta_t = min(dt)
    if t_final is None:
        t_max = [qnm_time_decay(tau[n], 1./1000) for n in range(nmodes)]
        t_final = max(t_max)
    kmax = int(t_final / delta_t) + 1

    outplus = TimeSeries(zeros(kmax, dtype=complex128), delta_t=delta_t)
    outcross = TimeSeries(zeros(kmax, dtype=complex128), delta_t=delta_t)
    for n in range(nmodes):
        hplus, hcross = get_td_qnm(template=None, f_0=f_0[n], tau=tau[n],
                            phi=phis[n], amp=amps[n], delta_t=delta_t,
                            t_final=t_final)
        outplus.data += hplus.data
        outcross.data += hcross.data

    return outplus, outcross

def get_fd_lm(template=None, **kwargs):
    """Return frequency domain lm mode with a given number of overtones.
    Parameters
    ----------
    template: object
        An object that has attached properties. This can be used to substitute
        for keyword arguments. A common example would be a row in an xml table.
    Mfinal : float
        Mass of the final black hole.
    Sfinal : float
        Spin of the final black hole.
    l : int
        l mode (lm modes available: 22, 21, 33, 44, 55).
    m : int
        m mode (lm modes available: 22, 21, 33, 44, 55).
    nmodes: int
        Number of overtones desired (maximum n=8)
    amp_n : float
        Amplitude of the n overtone, as many as the number of nmodes.
    phi_n : float
        Phase of the n overtone, as many as the number of nmodes.
    delta_f : {None, float}, optional
        The frequency step used to generate the ringdown.
        If None, it will be set to the inverse of the time at which the
        amplitude is 1/1000 of the peak amplitude (the minimum of all modes).
    f_lower: {None, float}, optional
        The starting frequency of the output frequency series.
        If None, it will be set to delta_f.
    f_final : {None, float}, optional
        The ending frequency of the output frequency series.
        If None, it will be set to the frequency at which the amplitude
        is 1/1000 of the peak amplitude (the maximum of all modes).
    Returns
    -------
    hplustilde: FrequencySeries
        The plus phase of a lm mode with overtones (n) in frequency domain.
    hcrosstilde: FrequencySeries
        The cross phase of a lm mode with overtones (n) in frequency domain.
    """

    input_params = props(template, lm_required_args, **kwargs)
    amps, phis = lm_amps_phases(input_params) 
    # the following may not be in input_params
    delta_f = input_params.pop('delta_f', None)
    f_lower = input_params.pop('f_lower', None)
    f_final = input_params.pop('f_final', None)
    # get required amplitudes and phases
#    amps, phis = numpy.zeros(nmodes), numpy.zeros(nmodes)
#    for n in range(nmodes):
#        try:
#            amps[n] = input_params['amp_%d' %n]
#        except KeyError:
#            raise ValueError('amp_%d is required' %n)
#        try:
#            phis[n] = input_params['phi_%d' %n]
#        except KeyError:
#            raise ValueError('phi_%d is required' %n)

    f_0, tau = get_lm_f0tau(Mfinal, Sfinal, l, m, nmodes)

    if delta_f is None:
        df = [1. / qnm_time_decay(tau[n], 1./1000) for n in range(nmodes)]
        delta_f = min(df)
    if f_final is None:
        f_max = [qnm_freq_decay(f_0[n], tau[n], 1./1000) for n in range(nmodes)]
        f_final = max(f_max)
    kmax = int(f_final / delta_f) + 1

    outplus = FrequencySeries(zeros(kmax, dtype=complex128), delta_f=delta_f)
    outcross = FrequencySeries(zeros(kmax, dtype=complex128), delta_f=delta_f)
    for n in range(nmodes):
        hplus, hcross = get_fd_qnm(template=None, f_0=f_0[n], tau=tau[n], 
                            phi=phis[n], amp=amps[n], delta_f=delta_f, 
                            f_lower=f_lower, f_final=f_final)
        outplus.data += hplus.data
        outcross.data += hcross.data

    return outplus, outcross

# Approximant names ###########################################################
ringdown_fd_approximants = {'FdQNM': get_fd_qnm}
ringdown_td_approximants = {'TdQNM': get_td_qnm}
