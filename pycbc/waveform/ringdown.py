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
"""Generate ringdown templates in the time and frequency domain.
"""

import numpy, lal
import lalsimulation as lalsim
from pycbc.types import TimeSeries, FrequencySeries, float64, complex128, zeros
from pycbc.waveform.waveform import get_obj_attrs
from pycbc.conversions import get_lm_f0tau_allmodes

default_qnm_args = {'t_0':0}
qnm_required_args = ['f_0', 'tau', 'amp', 'phi']
lm_required_args = ['freqs','taus','l','m','nmodes']
mass_spin_required_args = ['final_mass','final_spin', 'lmns', 'inclination']
freqtau_required_args = ['lmns']

max_freq = 16384.
min_dt = 1. / (2 * max_freq)
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

def lm_amps_phases(**kwargs):
    """ Take input_params and return dictionaries with amplitudes and phases
    of each overtone of a specific lm mode, checking that all of them are given.
    """
    l, m = kwargs['l'], kwargs['m']
    amps, phis = {}, {}
    # amp220 is always required, because the amplitudes of subdominant modes
    # are given as fractions of amp220.
    try:
        amps['220'] = kwargs['amp220']
    except KeyError:
        raise ValueError('amp220 is always required')

    # Get amplitudes of subdominant modes and all phases
    for n in range(kwargs['nmodes']):
        # If it is the 22 mode, skip 220
        if (l, m, n) != (2, 2, 0):
            try:
                amps['%d%d%d' %(l,m,n)] = kwargs['amp%d%d%d' %(l,m,n)] * amps['220']
            except KeyError:
                raise ValueError('amp%d%d%d is required' %(l,m,n))
        try:
            phis['%d%d%d' %(l,m,n)] = kwargs['phi%d%d%d' %(l,m,n)]
        except KeyError:
            raise ValueError('phi%d%d%d is required' %(l,m,n))

    return amps, phis

def lm_freqs_taus(**kwargs):
    """ Take input_params and return dictionaries with frequencies and damping
    times of each overtone of a specific lm mode, checking that all of them
    are given.
    """
    lmns = kwargs['lmns']
    freqs, taus = {}, {}

    for lmn in lmns:
        l, m, nmodes = int(lmn[0]), int(lmn[1]), int(lmn[2])
        for n in range(nmodes):
            try:
                freqs['%d%d%d' %(l,m,n)] = kwargs['f_%d%d%d' %(l,m,n)]
            except KeyError:
                raise ValueError('f_%d%d%d is required' %(l,m,n))
            try:
                taus['%d%d%d' %(l,m,n)] = kwargs['tau_%d%d%d' %(l,m,n)]
            except KeyError:
                raise ValueError('tau_%d%d%d is required' %(l,m,n))

    return freqs, taus

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

def lm_tfinal(damping_times, modes):
    """Return the maximum t_final of the modes given, with t_final the time
    at which the amplitude falls to 1/1000 of the peak amplitude
    """

    t_max = {}
    for lmn in modes:
        l, m, nmodes = int(lmn[0]), int(lmn[1]), int(lmn[2])
        for n in range(nmodes):
            t_max['%d%d%d' %(l,m,n)] = \
            qnm_time_decay(damping_times['%d%d%d' %(l,m,n)], 1./1000)

    return max(t_max.values())

def lm_deltat(freqs, damping_times, modes):
    """Return the minimum delta_t of all the modes given, with delta_t given by
    the inverse of the frequency at which the amplitude of the ringdown falls to
    1/1000 of the peak amplitude.
    """

    dt = {}
    for lmn in modes:
        l, m, nmodes = int(lmn[0]), int(lmn[1]), int(lmn[2])
        for n in range(nmodes):
            dt['%d%d%d' %(l,m,n)] = 1. / qnm_freq_decay(freqs['%d%d%d' %(l,m,n)],
                                       damping_times['%d%d%d' %(l,m,n)], 1./1000)

    delta_t = min(dt.values())
    if delta_t < min_dt:
        delta_t = min_dt

    return delta_t

def lm_ffinal(freqs, damping_times, modes):
    """Return the maximum f_final of the modes given, with f_final the frequency
    at which the amplitude falls to 1/1000 of the peak amplitude
    """

    f_max = {}
    for lmn in modes:
        l, m, nmodes = int(lmn[0]), int(lmn[1]), int(lmn[2])
        for n in range(nmodes):
            f_max['%d%d%d' %(l,m,n)] = qnm_freq_decay(freqs['%d%d%d' %(l,m,n)],
                                      damping_times['%d%d%d' %(l,m,n)], 1./1000)

    f_final = max(f_max.values())
    if f_final > max_freq:
        f_final = max_freq

    return f_final

def lm_deltaf(damping_times, modes):
    """Return the minimum delta_f of all the modes given, with delta_f given by
    the inverse of the time at which the amplitude of the ringdown falls to
    1/1000 of the peak amplitude.
    """

    df = {}
    for lmn in modes:
        l, m, nmodes = int(lmn[0]), int(lmn[1]), int(lmn[2])
        for n in range(nmodes):
            df['%d%d%d' %(l,m,n)] = \
                1. / qnm_time_decay(damping_times['%d%d%d' %(l,m,n)], 1./1000)

    return min(df.values())

# Spherical harmonics and Kerr factor #########################################

def spher_harms(l, m, inclination):
    """Return spherical harmonic polarizations
    """

    # FIXME: we are using spin -2 weighted spherical harmonics for now,
    # when possible switch to spheroidal harmonics.
    Y_lm = lal.SpinWeightedSphericalHarmonic(inclination, 0., -2, l, m).real
    Y_lminusm = lal.SpinWeightedSphericalHarmonic(inclination, 0., -2, l, -m).real
    Y_plus = Y_lm + (-1)**l * Y_lminusm
    Y_cross = Y_lm - (-1)**l * Y_lminusm

    return Y_plus, Y_cross

def Kerr_factor(final_mass, distance):
    """Return the factor final_mass/distance (in dimensionless units) for Kerr
    ringdowns
    """

    # Convert solar masses to meters
    mass = final_mass * lal.MSUN_SI * lal.G_SI / lal.C_SI**2
    # Convert Mpc to meters
    dist = distance * 1e6 * lal.PC_SI

    return mass / dist

# Functions for tapering ######################################################

def apply_taper(delta_t, taper, f_0, tau, amp, phi, l, m, inclination):
    """Return tapering window.
    """

    # Times of tapering do not include t=0
    taper_times = -numpy.arange(1, int(taper*tau/delta_t))[::-1] * delta_t
    Y_plus, Y_cross = spher_harms(l, m, inclination)
    taper_hp = amp * Y_plus * numpy.exp(10*taper_times/tau) * \
                     numpy.cos(two_pi*f_0*taper_times + phi)
    taper_hc = amp * Y_cross * numpy.exp(10*taper_times/tau) * \
                     numpy.sin(two_pi*f_0*taper_times + phi)

    return taper_hp, taper_hc

def taper_shift(waveform, output):
    """Add waveform to output with waveform shifted accordingly (for tapering
    multi-mode ringdowns)
    """

    if len(waveform) == len(output):
        output.data += waveform.data
    else:
        output.data[len(output)-len(waveform):] += waveform.data

    return output

# Functions to generate ringdown waveforms ####################################

######################################################
#### Basic functions to generate damped sinusoid
######################################################

def get_td_qnm(template=None, taper=None, **kwargs):
    """Return a time domain damped sinusoid.

    Parameters
    ----------
    template: object
        An object that has attached properties. This can be used to substitute
        for keyword arguments. A common example would be a row in an xml table.
    taper: {None, float}, optional
        Tapering at the beginning of the waveform with duration taper * tau.
        This option is recommended with timescales taper=1./2 or 1. for
        time-domain ringdown-only injections.
        The abrupt turn on of the ringdown can cause issues on the waveform
        when doing the fourier transform to the frequency domain. Setting
        taper will add a rapid ringup with timescale tau/10.
    f_0 : float
        The ringdown-frequency.
    tau : float
        The damping time of the sinusoid.
    amp : float
        The amplitude of the ringdown (constant for now).
    phi : float
        The initial phase of the ringdown. Should also include the information
        from the azimuthal angle (phi_0 + m*Phi)
    inclination : {None, float}, optional
        Inclination of the system in radians for the spherical harmonics.
    l : {2, int}, optional
        l mode for the spherical harmonics. Default is l=2.
    m : {2, int}, optional
        m mode for the spherical harmonics. Default is m=2.
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

    f_0 = input_params.pop('f_0')
    tau = input_params.pop('tau')
    amp = input_params.pop('amp')
    phi = input_params.pop('phi')
    # the following may not be in input_params
    inc = input_params.pop('inclination', None)
    l = input_params.pop('l', 2)
    m = input_params.pop('m', 2)
    delta_t = input_params.pop('delta_t', None)
    t_final = input_params.pop('t_final', None)

    if not delta_t:
        delta_t = 1. / qnm_freq_decay(f_0, tau, 1./1000)
        if delta_t < min_dt:
            delta_t = min_dt
    if not t_final:
        t_final = qnm_time_decay(tau, 1./1000)

    kmax = int(t_final / delta_t) + 1
    times = numpy.arange(kmax) * delta_t
    if inc is not None:
        Y_plus, Y_cross = spher_harms(l, m, inc)
    else:
        Y_plus, Y_cross = 1, 1

    hplus = amp * Y_plus * numpy.exp(-times/tau) * \
                                numpy.cos(two_pi*f_0*times + phi)
    hcross = amp * Y_cross * numpy.exp(-times/tau) * \
                                numpy.sin(two_pi*f_0*times + phi)

    if taper and delta_t < taper*tau:
        taper_window = int(taper*tau/delta_t)
        kmax += taper_window

    outplus = TimeSeries(zeros(kmax), delta_t=delta_t)
    outcross = TimeSeries(zeros(kmax), delta_t=delta_t)

    # If size of tapering window is less than delta_t, do not apply taper.
    if not taper or delta_t > taper*tau:
        outplus.data[:kmax] = hplus
        outcross.data[:kmax] = hcross

        return outplus, outcross

    else:
        taper_hp, taper_hc = apply_taper(delta_t, taper, f_0, tau, amp, phi,
                                                                    l, m, inc)
        start = - taper * tau
        outplus.data[:taper_window] = taper_hp
        outplus.data[taper_window:] = hplus
        outcross.data[:taper_window] = taper_hc
        outcross.data[taper_window:] = hcross
        outplus._epoch, outcross._epoch = start, start

        return outplus, outcross

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
    amp : float
        The amplitude of the ringdown (constant for now).
    phi : float
        The initial phase of the ringdown. Should also include the information
        from the azimuthal angle (phi_0 + m*Phi).
    inclination : {None, float}, optional
        Inclination of the system in radians for the spherical harmonics.
    l : {2, int}, optional
        l mode for the spherical harmonics. Default is l=2.
    m : {2, int}, optional
        m mode for the spherical harmonics. Default is m=2.
    t_0 :  {0, float}, optional
        The starting time of the ringdown.
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

    f_0 = input_params.pop('f_0')
    tau = input_params.pop('tau')
    amp = input_params.pop('amp')
    phi = input_params.pop('phi')
    # the following have defaults, and so will be populated
    t_0 = input_params.pop('t_0')
    # the following may not be in input_params
    inc = input_params.pop('inclination', None)
    l = input_params.pop('l', 2)
    m = input_params.pop('m', 2)
    delta_f = input_params.pop('delta_f', None)
    f_lower = input_params.pop('f_lower', None)
    f_final = input_params.pop('f_final', None)

    if not delta_f:
        delta_f = 1. / qnm_time_decay(tau, 1./1000)
    if not f_lower:
        f_lower = delta_f
        kmin = 0
    else:
        kmin = int(f_lower / delta_f)
    if not f_final:
        f_final = qnm_freq_decay(f_0, tau, 1./1000)
    if f_final > max_freq:
        f_final = max_freq
    kmax = int(f_final / delta_f) + 1

    freqs = numpy.arange(kmin, kmax)*delta_f
    if inc is not None:
        Y_plus, Y_cross = spher_harms(l, m, inc)
    else:
        Y_plus, Y_cross = 1, 1

    denominator = 1 + (4j * pi * freqs * tau) - (4 * pi_sq * ( freqs*freqs - f_0*f_0) * tau*tau)
    norm = amp * tau / denominator
    if t_0 != 0:
        time_shift = numpy.exp(-1j * two_pi * freqs * t_0)
        norm *= time_shift

    # Analytical expression for the Fourier transform of the ringdown (damped sinusoid)
    hp_tilde = norm * Y_plus * ( (1 + 2j * pi * freqs * tau) * numpy.cos(phi)
                                        - two_pi * f_0 * tau * numpy.sin(phi) )
    hc_tilde = norm * Y_cross * ( (1 + 2j * pi * freqs * tau) * numpy.sin(phi)
                                        + two_pi * f_0 * tau * numpy.cos(phi) )

    outplus = FrequencySeries(zeros(kmax, dtype=complex128), delta_f=delta_f)
    outcross = FrequencySeries(zeros(kmax, dtype=complex128), delta_f=delta_f)
    outplus.data[kmin:kmax] = hp_tilde
    outcross.data[kmin:kmax] = hc_tilde

    return outplus, outcross

######################################################
#### Single lm mode with overtones
######################################################

def get_td_lm(template=None, taper=None, **kwargs):
    """Return time domain lm mode with the given number of overtones.

    Parameters
    ----------
    template: object
        An object that has attached properties. This can be used to substitute
        for keyword arguments. A common example would be a row in an xml table.
    taper: {None, float}, optional
        Tapering at the beginning of the waveform with duration taper * tau.
        This option is recommended with timescales taper=1./2 or 1. for
        time-domain ringdown-only injections.
        The abrupt turn on of the ringdown can cause issues on the waveform
        when doing the fourier transform to the frequency domain. Setting
        taper will add a rapid ringup with timescale tau/10.
        Each overtone will have a different taper depending on its tau, the
        final taper being the superposition of all the tapers.
    freqs : dict
        {lmn:f_lmn} Dictionary of the central frequencies for each overtone,
        as many as number of modes.
    taus : dict
        {lmn:tau_lmn} Dictionary of the damping times for each overtone,
        as many as number of modes.
    l : int
        l mode (lm modes available: 22, 21, 33, 44, 55).
    m : int
        m mode (lm modes available: 22, 21, 33, 44, 55).
    nmodes: int
        Number of overtones desired (maximum n=8)
    amp220 : float
        Amplitude of the fundamental 220 mode, needed for any lm.
    amplmn : float
        Fraction of the amplitude of the lmn overtone relative to the
        fundamental mode, as many as the number of subdominant modes.
    philmn : float
        Phase of the lmn overtone, as many as the number of modes. Should also
        include the information from the azimuthal angle (phi + m*Phi).
    inclination : {None, float}, optional
        Inclination of the system in radians for the spherical harmonics.
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
    hplus: TimeSeries
        The plus phase of a lm mode with overtones (n) in time domain.
    hcross: TimeSeries
        The cross phase of a lm mode with overtones (n) in time domain.
    """

    input_params = props(template, lm_required_args, **kwargs)

    # Get required args
    amps, phis = lm_amps_phases(**input_params)
    f_0 = input_params.pop('freqs')
    tau = input_params.pop('taus')
    inc = input_params.pop('inclination', None)
    l, m = input_params.pop('l'), input_params.pop('m')
    nmodes = input_params.pop('nmodes')
    if int(nmodes) == 0:
        raise ValueError('Number of overtones (nmodes) must be greater '
                         'than zero.')
    # The following may not be in input_params
    delta_t = input_params.pop('delta_t', None)
    t_final = input_params.pop('t_final', None)

    if not delta_t:
        delta_t = lm_deltat(f_0, tau, ['%d%d%d' %(l,m,nmodes)])
    if not t_final:
        t_final = lm_tfinal(tau, ['%d%d%d' %(l, m, nmodes)])

    kmax = int(t_final / delta_t) + 1
    # Different overtones will have different tapering window-size
    # Find maximum window size to create long enough output vector
    if taper:
        taper_window = int(taper*max(tau.values())/delta_t)
        kmax += taper_window

    outplus = TimeSeries(zeros(kmax, dtype=float64), delta_t=delta_t)
    outcross = TimeSeries(zeros(kmax, dtype=float64), delta_t=delta_t)
    if taper:
        start = - taper * max(tau.values())
        outplus._epoch, outcross._epoch = start, start

    for n in range(nmodes):
        hplus, hcross = get_td_qnm(template=None, taper=taper,
                            f_0=f_0['%d%d%d' %(l,m,n)],
                            tau=tau['%d%d%d' %(l,m,n)],
                            phi=phis['%d%d%d' %(l,m,n)],
                            amp=amps['%d%d%d' %(l,m,n)],
                            inclination=inc, l=l, m=m,
                            delta_t=delta_t, t_final=t_final)
        if not taper:
            outplus.data += hplus.data
            outcross.data += hcross.data
        else:
            outplus = taper_shift(hplus, outplus)
            outcross = taper_shift(hcross, outcross)

    return outplus, outcross

def get_fd_lm(template=None, **kwargs):
    """Return frequency domain lm mode with a given number of overtones.

    Parameters
    ----------
    template: object
        An object that has attached properties. This can be used to substitute
        for keyword arguments. A common example would be a row in an xml table.
    freqs : dict
        {lmn:f_lmn} Dictionary of the central frequencies for each overtone,
        as many as number of modes.
    taus : dict
        {lmn:tau_lmn} Dictionary of the damping times for each overtone,
        as many as number of modes.
    l : int
        l mode (lm modes available: 22, 21, 33, 44, 55).
    m : int
        m mode (lm modes available: 22, 21, 33, 44, 55).
    nmodes: int
        Number of overtones desired (maximum n=8)
    amplmn : float
        Amplitude of the lmn overtone, as many as the number of nmodes.
    philmn : float
        Phase of the lmn overtone, as many as the number of modes. Should also
        include the information from the azimuthal angle (phi + m*Phi).
    inclination : {None, float}, optional
        Inclination of the system in radians for the spherical harmonics.
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
        The plus phase of a lm mode with n overtones in frequency domain.
    hcrosstilde: FrequencySeries
        The cross phase of a lm mode with n overtones in frequency domain.
    """

    input_params = props(template, lm_required_args, **kwargs)

    # Get required args
    amps, phis = lm_amps_phases(**input_params)
    f_0 = input_params.pop('freqs')
    tau = input_params.pop('taus')
    l, m = input_params.pop('l'), input_params.pop('m')
    inc = input_params.pop('inclination', None)
    nmodes = input_params.pop('nmodes')
    if int(nmodes) == 0:
        raise ValueError('Number of overtones (nmodes) must be greater '
                         'than zero.')
    # The following may not be in input_params
    delta_f = input_params.pop('delta_f', None)
    f_lower = input_params.pop('f_lower', None)
    f_final = input_params.pop('f_final', None)

    if not delta_f:
        delta_f = lm_deltaf(tau, ['%d%d%d' %(l,m,nmodes)])
    if not f_final:
        f_final = lm_ffinal(f_0, tau, ['%d%d%d' %(l, m, nmodes)])
    kmax = int(f_final / delta_f) + 1

    outplus = FrequencySeries(zeros(kmax, dtype=complex128), delta_f=delta_f)
    outcross = FrequencySeries(zeros(kmax, dtype=complex128), delta_f=delta_f)

    for n in range(nmodes):
        hplus, hcross = get_fd_qnm(template=None, f_0=f_0['%d%d%d' %(l,m,n)],
                            tau=tau['%d%d%d' %(l,m,n)],
                            amp=amps['%d%d%d' %(l,m,n)],
                            phi=phis['%d%d%d' %(l,m,n)],
                            inclination=inc, l=l, m=m, delta_f=delta_f,
                            f_lower=f_lower, f_final=f_final)
        outplus.data += hplus.data
        outcross.data += hcross.data

    return outplus, outcross

######################################################
#### Approximants
######################################################

def get_td_from_final_mass_spin(template=None, taper=None,
                                distance=None, **kwargs):
    """Return time domain ringdown with all the modes specified.

    Parameters
    ----------
    template: object
        An object that has attached properties. This can be used to substitute
        for keyword arguments. A common example would be a row in an xml table.
    taper: {None, float}, optional
        Tapering at the beginning of the waveform with duration taper * tau.
        This option is recommended with timescales taper=1./2 or 1. for
        time-domain ringdown-only injections.
        The abrupt turn on of the ringdown can cause issues on the waveform
        when doing the fourier transform to the frequency domain. Setting
        taper will add a rapid ringup with timescale tau/10.
        Each mode and overtone will have a different taper depending on its tau,
        the final taper being the superposition of all the tapers.
    distance : {None, float}, optional
        Luminosity distance of the system. If specified, the returned ringdown
        will contain the factor (final_mass/distance).
    final_mass : float
        Mass of the final black hole.
    final_spin : float
        Spin of the final black hole.
    lmns : list
        Desired lmn modes as strings (lm modes available: 22, 21, 33, 44, 55).
        The n specifies the number of overtones desired for the corresponding
        lm pair (maximum n=8).
        Example: lmns = ['223','331'] are the modes 220, 221, 222, and 330
    amp220 : float
        Amplitude of the fundamental 220 mode. Note that if distance is given,
        this parameter will have a completely different order of magnitude.
        See table II in https://arxiv.org/abs/1107.0854 for an estimate.
    amplmn : float
        Fraction of the amplitude of the lmn overtone relative to the
        fundamental mode, as many as the number of subdominant modes.
    philmn : float
        Phase of the lmn overtone, as many as the number of modes. Should also
        include the information from the azimuthal angle (phi + m*Phi).
    inclination : float
        Inclination of the system in radians.
    delta_t : {None, float}, optional
        The time step used to generate the ringdown.
        If None, it will be set to the inverse of the frequency at which the
        amplitude is 1/1000 of the peak amplitude (the minimum of all modes).
    t_final : {None, float}, optional
        The ending time of the output frequency series.
        If None, it will be set to the time at which the amplitude
        is 1/1000 of the peak amplitude (the maximum of all modes).

    Returns
    -------
    hplus: TimeSeries
        The plus phase of a ringdown with the lm modes specified and
        n overtones in time domain.
    hcross: TimeSeries
        The cross phase of a ringdown with the lm modes specified and
        n overtones in time domain.
    """

    input_params = props(template, mass_spin_required_args, **kwargs)

    # Get required args
    final_mass = input_params['final_mass']
    final_spin = input_params['final_spin']
    lmns = input_params['lmns']
    for lmn in lmns:
        if int(lmn[2]) == 0:
            raise ValueError('Number of overtones (nmodes) must be greater '
                             'than zero.')
    # following may not be in input_params
    delta_t = input_params.pop('delta_t', None)
    t_final = input_params.pop('t_final', None)

    f_0, tau = get_lm_f0tau_allmodes(final_mass, final_spin, lmns)

    if not delta_t:
        delta_t = lm_deltat(f_0, tau, lmns)
    if not t_final:
        t_final = lm_tfinal(tau, lmns)

    kmax = int(t_final / delta_t) + 1
    # Different overtones will have different tapering window-size
    # Find maximum window size to create long enough output vector
    if taper:
        taper_window = int(taper*max(tau.values())/delta_t)
        kmax += taper_window

    outplus = TimeSeries(zeros(kmax, dtype=float64), delta_t=delta_t)
    outcross = TimeSeries(zeros(kmax, dtype=float64), delta_t=delta_t)
    if taper:
        start = - taper * max(tau.values())
        outplus._epoch, outcross._epoch = start, start

    for lmn in lmns:
        l, m, nmodes = int(lmn[0]), int(lmn[1]), int(lmn[2])
        hplus, hcross = get_td_lm(taper=taper, freqs=f_0, taus=tau,
                             l=l, m=m, nmodes=nmodes,
                             delta_t=delta_t, t_final=t_final,
                             **input_params)
        if not taper:
            outplus.data += hplus.data
            outcross.data += hcross.data
        else:
            outplus = taper_shift(hplus, outplus)
            outcross = taper_shift(hcross, outcross)

    norm = Kerr_factor(final_mass, distance) if distance else 1.
    return norm*outplus, norm*outcross

def get_fd_from_final_mass_spin(template=None, distance=None, **kwargs):
    """Return frequency domain ringdown with all the modes specified.

    Parameters
    ----------
    template: object
        An object that has attached properties. This can be used to substitute
        for keyword arguments. A common example would be a row in an xml table.
    distance : {None, float}, optional
        Luminosity distance of the system. If specified, the returned ringdown
        will contain the factor (final_mass/distance).
    final_mass : float
        Mass of the final black hole.
    final_spin : float
        Spin of the final black hole.
    lmns : list
        Desired lmn modes as strings (lm modes available: 22, 21, 33, 44, 55).
        The n specifies the number of overtones desired for the corresponding
        lm pair (maximum n=8).
        Example: lmns = ['223','331'] are the modes 220, 221, 222, and 330
    amp220 : float
        Amplitude of the fundamental 220 mode. Note that if distance is given,
        this parameter will have a completely different order of magnitude.
        See table II in https://arxiv.org/abs/1107.0854 for an estimate.
    amplmn : float
        Fraction of the amplitude of the lmn overtone relative to the
        fundamental mode, as many as the number of subdominant modes.
    philmn : float
        Phase of the lmn overtone, as many as the number of modes. Should also
        include the information from the azimuthal angle (phi + m*Phi).
    inclination : float
        Inclination of the system in radians.
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
        The plus phase of a ringdown with the lm modes specified and
        n overtones in frequency domain.
    hcrosstilde: FrequencySeries
        The cross phase of a ringdown with the lm modes specified and
        n overtones in frequency domain.
    """

    input_params = props(template, mass_spin_required_args, **kwargs)

    # Get required args
    final_mass = input_params['final_mass']
    final_spin = input_params['final_spin']
    lmns = input_params['lmns']
    for lmn in lmns:
        if int(lmn[2]) == 0:
            raise ValueError('Number of overtones (nmodes) must be greater '
                             'than zero.')
    # The following may not be in input_params
    delta_f = input_params.pop('delta_f', None)
    f_lower = input_params.pop('f_lower', None)
    f_final = input_params.pop('f_final', None)

    f_0, tau = get_lm_f0tau_allmodes(final_mass, final_spin, lmns)

    if not delta_f:
        delta_f = lm_deltaf(tau, lmns)
    if not f_final:
        f_final = lm_ffinal(f_0, tau, lmns)
    if not f_lower:
        f_lower = delta_f
    kmax = int(f_final / delta_f) + 1

    outplustilde = FrequencySeries(zeros(kmax, dtype=complex128), delta_f=delta_f)
    outcrosstilde = FrequencySeries(zeros(kmax, dtype=complex128), delta_f=delta_f)
    for lmn in lmns:
        l, m, nmodes = int(lmn[0]), int(lmn[1]), int(lmn[2])
        hplustilde, hcrosstilde = get_fd_lm(freqs=f_0, taus=tau,
                                        l=l, m=m, nmodes=nmodes,
                                        delta_f=delta_f, f_lower=f_lower,
                                        f_final=f_final, **input_params)
        outplustilde.data += hplustilde.data
        outcrosstilde.data += hcrosstilde.data

    norm = Kerr_factor(final_mass, distance) if distance else 1.
    return norm*outplustilde, norm*outcrosstilde

def get_td_from_freqtau(template=None, taper=None, **kwargs):
    """Return time domain ringdown with all the modes specified.

    Parameters
    ----------
    template: object
        An object that has attached properties. This can be used to substitute
        for keyword arguments. A common example would be a row in an xml table.
    taper: {None, float}, optional
        Tapering at the beginning of the waveform with duration taper * tau.
        This option is recommended with timescales taper=1./2 or 1. for
        time-domain ringdown-only injections.
        The abrupt turn on of the ringdown can cause issues on the waveform
        when doing the fourier transform to the frequency domain. Setting
        taper will add a rapid ringup with timescale tau/10.
        Each mode and overtone will have a different taper depending on its tau,
        the final taper being the superposition of all the tapers.
    lmns : list
        Desired lmn modes as strings (lm modes available: 22, 21, 33, 44, 55).
        The n specifies the number of overtones desired for the corresponding
        lm pair (maximum n=8).
        Example: lmns = ['223','331'] are the modes 220, 221, 222, and 330
    f_lmn: float
        Central frequency of the lmn overtone, as many as number of modes.
    tau_lmn: float
        Damping time of the lmn overtone, as many as number of modes.
    amp220 : float
        Amplitude of the fundamental 220 mode.
    amplmn : float
        Fraction of the amplitude of the lmn overtone relative to the
        fundamental mode, as many as the number of subdominant modes.
    philmn : float
        Phase of the lmn overtone, as many as the number of modes. Should also
        include the information from the azimuthal angle (phi + m*Phi).
    inclination : {None, float}, optional
        Inclination of the system in radians. If None, the spherical harmonics
        will be set to 1.
    delta_t : {None, float}, optional
        The time step used to generate the ringdown.
        If None, it will be set to the inverse of the frequency at which the
        amplitude is 1/1000 of the peak amplitude (the minimum of all modes).
    t_final : {None, float}, optional
        The ending time of the output frequency series.
        If None, it will be set to the time at which the amplitude
        is 1/1000 of the peak amplitude (the maximum of all modes).

    Returns
    -------
    hplustilde: FrequencySeries
        The plus phase of a ringdown with the lm modes specified and
        n overtones in frequency domain.
    hcrosstilde: FrequencySeries
        The cross phase of a ringdown with the lm modes specified and
        n overtones in frequency domain.
    """

    input_params = props(template, freqtau_required_args, **kwargs)

    # Get required args
    f_0, tau = lm_freqs_taus(**input_params)
    lmns = input_params['lmns']
    for lmn in lmns:
        if int(lmn[2]) == 0:
            raise ValueError('Number of overtones (nmodes) must be greater '
                             'than zero.')
    # following may not be in input_params
    inc = input_params.pop('inclination', None)
    delta_t = input_params.pop('delta_t', None)
    t_final = input_params.pop('t_final', None)

    if not delta_t:
        delta_t = lm_deltat(f_0, tau, lmns)
    if not t_final:
        t_final = lm_tfinal(tau, lmns)

    kmax = int(t_final / delta_t) + 1
    # Different overtones will have different tapering window-size
    # Find maximum window size to create long enough output vector
    if taper:
        taper_window = int(taper*max(tau.values())/delta_t)
        kmax += taper_window

    outplus = TimeSeries(zeros(kmax, dtype=float64), delta_t=delta_t)
    outcross = TimeSeries(zeros(kmax, dtype=float64), delta_t=delta_t)
    if taper:
        start = - taper * max(tau.values())
        outplus._epoch, outcross._epoch = start, start

    for lmn in lmns:
        l, m, nmodes = int(lmn[0]), int(lmn[1]), int(lmn[2])
        hplus, hcross = get_td_lm(freqs=f_0, taus=tau, l=l, m=m, nmodes=nmodes,
                             taper=taper, inclination=inc, delta_t=delta_t,
                             t_final=t_final, **input_params)
        if not taper:
            outplus.data += hplus.data
            outcross.data += hcross.data
        else:
            outplus = taper_shift(hplus, outplus)
            outcross = taper_shift(hcross, outcross)

    return outplus, outcross

def get_fd_from_freqtau(template=None, **kwargs):
    """Return frequency domain ringdown with all the modes specified.

    Parameters
    ----------
    template: object
        An object that has attached properties. This can be used to substitute
        for keyword arguments. A common example would be a row in an xml table.
    lmns : list
        Desired lmn modes as strings (lm modes available: 22, 21, 33, 44, 55).
        The n specifies the number of overtones desired for the corresponding
        lm pair (maximum n=8).
        Example: lmns = ['223','331'] are the modes 220, 221, 222, and 330
    f_lmn: float
        Central frequency of the lmn overtone, as many as number of modes.
    tau_lmn: float
        Damping time of the lmn overtone, as many as number of modes.
    amp220 : float
        Amplitude of the fundamental 220 mode.
    amplmn : float
        Fraction of the amplitude of the lmn overtone relative to the
        fundamental mode, as many as the number of subdominant modes.
    philmn : float
        Phase of the lmn overtone, as many as the number of modes. Should also
        include the information from the azimuthal angle (phi + m*Phi).
    inclination : {None, float}, optional
        Inclination of the system in radians. If None, the spherical harmonics
        will be set to 1.
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
        The plus phase of a ringdown with the lm modes specified and
        n overtones in frequency domain.
    hcrosstilde: FrequencySeries
        The cross phase of a ringdown with the lm modes specified and
        n overtones in frequency domain.
    """

    input_params = props(template, freqtau_required_args, **kwargs)

    # Get required args
    f_0, tau = lm_freqs_taus(**input_params)
    lmns = input_params['lmns']
    for lmn in lmns:
        if int(lmn[2]) == 0:
            raise ValueError('Number of overtones (nmodes) must be greater '
                             'than zero.')
    # The following may not be in input_params
    inc = input_params.pop('inclination', None)
    delta_f = input_params.pop('delta_f', None)
    f_lower = input_params.pop('f_lower', None)
    f_final = input_params.pop('f_final', None)

    if not delta_f:
        delta_f = lm_deltaf(tau, lmns)
    if not f_final:
        f_final = lm_ffinal(f_0, tau, lmns)
    if not f_lower:
        f_lower = delta_f
    kmax = int(f_final / delta_f) + 1

    outplustilde = FrequencySeries(zeros(kmax, dtype=complex128), delta_f=delta_f)
    outcrosstilde = FrequencySeries(zeros(kmax, dtype=complex128), delta_f=delta_f)
    for lmn in lmns:
        l, m, nmodes = int(lmn[0]), int(lmn[1]), int(lmn[2])
        hplustilde, hcrosstilde = get_fd_lm(freqs=f_0, taus=tau,
                                        l=l, m=m, nmodes=nmodes,
                                        inclination=inc,
                                        delta_f=delta_f, f_lower=f_lower,
                                        f_final=f_final, **input_params)
        outplustilde.data += hplustilde.data
        outcrosstilde.data += hcrosstilde.data

    return outplustilde, outcrosstilde

# Approximant names ###########################################################
ringdown_fd_approximants = {'FdQNMfromFinalMassSpin': get_fd_from_final_mass_spin,
                            'FdQNMfromFreqTau': get_fd_from_freqtau}
ringdown_td_approximants = {'TdQNMfromFinalMassSpin': get_td_from_final_mass_spin,
                            'TdQNMfromFreqTau': get_td_from_freqtau}
