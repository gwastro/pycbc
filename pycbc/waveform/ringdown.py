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

import re, numpy, lal
from pycbc.types import TimeSeries, FrequencySeries, float64, complex128, zeros
from pycbc.waveform.waveform import get_obj_attrs
from pycbc.conversions import get_lm_f0tau_allmodes

qnm_required_args = ['f_0', 'tau', 'amp', 'phi']
mass_spin_required_args = ['final_mass','final_spin', 'lmns', 'inclination']
freqtau_required_args = ['lmns']
td_args = {'delta_t':None, 't_final':None, 'taper':False}
fd_args = {'t_0':0, 'delta_f':None, 'f_lower':None, 'f_final':None}

max_freq = 16384/2.
min_dt = 1. / (2 * max_freq)
pi = numpy.pi
two_pi = 2 * numpy.pi
pi_sq = numpy.pi * numpy.pi

# Input parameters ############################################################

def props(obj, required, domain_args, **kwargs):
    """ Return a dictionary built from the combination of defaults, kwargs,
    and the attributes of the given object.
    """
    # Get the attributes of the template object
    pr = get_obj_attrs(obj)

    # Get the parameters to generate the waveform
    # Note that keyword arguments override values in the template object
    input_params = domain_args.copy()
    input_params.update(pr)
    input_params.update(kwargs)

    # Check if the required arguments are given
    for arg in required:
        if arg not in input_params:
            raise ValueError('Please provide ' + str(arg))

    return input_params

def format_lmns(lmns):
    """ Checks if the format of the parameter lmns is correct, returning the
    appropriate format if not, and raise an error if nmodes=0.
    The required format for the ringdown approximants is a list of lmn modes
    as strings, with n the number of overtones desired. For instance,
    lmns = ['223','331'] are the modes 220, 221, 222, and 330.
    The ConfigParser of a workflow might convert that to a single string
    (case 1 below) or a list with a single string (case 2), and this function
    will return the appropriate list of strings. If a different format is
    given, raise an error.
    """

    # Case 1: the list is in a string "['221', '331']"
    # In Python3 this might be "[b'221', b'331']"
    if isinstance(lmns, str):
        # strip off brackets and convert to list
        lmns = lmns.strip('[]').split(',')

    # Case 2: a list with only one string with a list ["221', '331"]
    # In Python3 this might be ["b221', b'331"]
    elif (len(lmns) == 1 and isinstance(lmns[0], str)
          and len(lmns[0]) > 3):
        lmns = lmns[0].split(',')

    out = []
    # Cycle over the lmns to ensure that we get back a list of strings that
    # are three digits long, and that nmodes!=0
    for lmn in lmns:
        # The following line is to be used with Python3 if the lmns are stored
        # as a list of strings in the HDF files and the workflow converts that
        # to a string
        # lmn = lmn.strip(" b'")
        # Try to convert to int and then str, to ensure the right format
        lmn = str(int(lmn))
        if len(lmn) != 3:
            raise ValueError('Format of parameter lmns not recognized. See '
                             'approximant documentation for more info.')
        elif int(lmn[2]) == 0:
            raise ValueError('Number of overtones (nmodes) must be greater '
                             'than zero in lmn={}.'.format(lmn))
        out.append(lmn)

    return out

def lm_amps_phases(**kwargs):
    """ Take input_params and return dictionaries with amplitudes and phases
    of each overtone of a specific lm mode, checking that all of them are given.
    """
    lmns = format_lmns(kwargs['lmns'])
    amps, phis = {}, {}
    # amp220 is always required, because the amplitudes of subdominant modes
    # are given as fractions of amp220.
    try:
        amps['220'] = kwargs['amp220']
    except KeyError:
        raise ValueError('amp220 is always required')

    # Get amplitudes of subdominant modes and all phases
    for lmn in lmns:
        lm, nmodes = lmn[0:-1], int(lmn[2])
        for n in range(nmodes):
            # If it is the 22 mode, skip 220 amplitude
            mode = lm + '{}'.format(n)
            if mode != '220':
                try:
                    amps[mode] = kwargs['amp' + mode] * amps['220']
                except KeyError:
                    raise ValueError('amp{} is required'.format(mode))
            try:
                phis[mode] = kwargs['phi' + mode]
            except KeyError:
                raise ValueError('phi{} is required'.format(mode))

    return amps, phis

def lm_freqs_taus(**kwargs):
    """ Take input_params and return dictionaries with frequencies and damping
    times of each overtone of a specific lm mode, checking that all of them
    are given.
    """
    lmns = format_lmns(kwargs['lmns'])
    freqs, taus = {}, {}

    for lmn in lmns:
        lm, nmodes = lmn[0:-1], int(lmn[2])
        for n in range(nmodes):
            mode = lm + '{}'.format(n)
            try:
                freqs[mode] = kwargs['f_' + mode]
            except KeyError:
                raise ValueError('f_{} is required'.format(mode))
            try:
                taus[mode] = kwargs['tau_' + mode]
            except KeyError:
                raise ValueError('tau_{} is required'.format(mode))

    return freqs, taus

# Functions to obtain t_final, f_final and output vector ######################

def qnm_time_decay(tau, decay):
    """Return the time at which the amplitude of the
    ringdown falls to decay of the peak amplitude.

    Parameters
    ----------
    tau : float
        The damping time of the sinusoid.
    decay : float
        The fraction of the peak amplitude.

    Returns
    -------
    t_decay : float
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
    decay : float
        The fraction of the peak amplitude.

    Returns
    -------
    f_decay : float
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

def lm_tfinal(damping_times):
    """Return the maximum t_final of the modes given, with t_final the time
    at which the amplitude falls to 1/1000 of the peak amplitude
    """

    if isinstance(damping_times, dict):
        t_max = {}
        for lmn in damping_times.keys():
            t_max[lmn] = qnm_time_decay(damping_times[lmn], 1./1000)
        t_final = max(t_max.values())
    else:
        t_final = qnm_time_decay(damping_times, 1./1000)

    return t_final

def lm_deltat(freqs, damping_times):
    """Return the minimum delta_t of all the modes given, with delta_t given by
    the inverse of the frequency at which the amplitude of the ringdown falls to
    1/1000 of the peak amplitude.
    """

    if isinstance(freqs, dict) and isinstance(damping_times, dict):
        dt = {}
        for lmn in freqs.keys():
            dt[lmn] = 1. / qnm_freq_decay(freqs[lmn],
                               damping_times[lmn], 1./1000)
        delta_t = min(dt.values())
    elif isinstance(freqs, dict) and not isinstance(damping_times, dict):
        raise ValueError('Missing damping times.')
    elif isinstance(damping_times, dict) and not isinstance(freqs, dict):
        raise ValueError('Missing frequencies.')
    else:
        delta_t = 1. / qnm_freq_decay(freqs, damping_times, 1./1000)

    if delta_t < min_dt:
        delta_t = min_dt

    return delta_t

def lm_ffinal(freqs, damping_times):
    """Return the maximum f_final of the modes given, with f_final the frequency
    at which the amplitude falls to 1/1000 of the peak amplitude
    """

    if isinstance(freqs, dict) and isinstance(damping_times, dict):
        f_max = {}
        for lmn in freqs.keys():
            f_max[lmn] = qnm_freq_decay(freqs[lmn],
                              damping_times[lmn], 1./1000)
        f_final = max(f_max.values())
    elif isinstance(freqs, dict) and not isinstance(damping_times, dict):
        raise ValueError('Missing damping times.')
    elif isinstance(damping_times, dict) and not isinstance(freqs, dict):
        raise ValueError('Missing frequencies.')
    else:
        f_final = qnm_freq_decay(freqs, damping_times, 1./1000)

    if f_final > max_freq:
        f_final = max_freq

    return f_final

def lm_deltaf(damping_times):
    """Return the minimum delta_f of all the modes given, with delta_f given by
    the inverse of the time at which the amplitude of the ringdown falls to
    1/1000 of the peak amplitude.
    """

    if isinstance(damping_times, dict):
        df = {}
        for lmn in damping_times.keys():
            df[lmn] = 1. / qnm_time_decay(damping_times[lmn], 1./1000)
        delta_f = min(df.values())
    else:
        delta_f = 1. / qnm_time_decay(damping_times, 1./1000)

    return delta_f

def td_output_vector(freqs, damping_times, taper=False,
                     delta_t=None, t_final=None):
    """Return an empty TimeSeries with the appropriate size to fit all
    the quasi-normal modes present in freqs, damping_times
    """

    if not delta_t:
        delta_t = lm_deltat(freqs, damping_times)
    if not t_final:
        t_final = lm_tfinal(damping_times)
    kmax = int(t_final / delta_t) + 1

    # Different modes will have different tapering window-size
    # Find maximum window size to create long enough output vector
    if taper:
        max_tau = max(damping_times.values()) if \
                  isinstance(damping_times, dict) else damping_times
        kmax += int(max_tau/delta_t)

    outplus = TimeSeries(zeros(kmax, dtype=float64), delta_t=delta_t)
    outcross = TimeSeries(zeros(kmax, dtype=float64), delta_t=delta_t)
    if taper:
        # Change epoch of output vector if tapering will be applied
        start = - max_tau
        # To ensure that t=0 is still in output vector
        start -= start % delta_t
        outplus._epoch, outcross._epoch = start, start

    return outplus, outcross

def fd_output_vector(freqs, damping_times, delta_f=None, f_final=None):
    """Return an empty FrequencySeries with the appropriate size to fit all
    the quasi-normal modes present in freqs, damping_times
    """

    if not delta_f:
        delta_f = lm_deltaf(damping_times)
    if not f_final:
        f_final = lm_ffinal(freqs, damping_times)
    kmax = int(f_final / delta_f) + 1

    outplus = FrequencySeries(zeros(kmax, dtype=complex128), delta_f=delta_f)
    outcross = FrequencySeries(zeros(kmax, dtype=complex128), delta_f=delta_f)

    return outplus, outcross

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

def apply_taper(tau, amp, phi, delta_t, l=2, m=2, inclination=None):
    """Return a tapering window of the form exp(10*t/tau).
    """

    taper_times = -numpy.arange(0, int(tau/delta_t))[::-1] * delta_t
    if inclination is not None:
        Y_plus, Y_cross = spher_harms(l, m, inclination)
    else:
        Y_plus, Y_cross = 1, 1
    hp_amp = amp * Y_plus * numpy.cos(phi)
    hc_amp = amp * Y_cross * numpy.sin(phi)

    taper_hp = hp_amp * numpy.exp(10*taper_times/tau)
    taper_hc = hc_amp * numpy.exp(10*taper_times/tau)

    return taper_hp, taper_hc

# Functions to generate ringdown waveforms ####################################

######################################################
#### Basic functions to generate damped sinusoid
######################################################

def td_damped_sinusoid(f_0, tau, amp, phi, delta_t, t_final,
                       l=2, m=2, inclination=None):
    """Return a time domain damped sinusoid (plus and cross polarizations)
    with central frequency f_0, damping time tau, amplitude amp and phase phi.
    The l, m, and inclination parameters are used for the spherical harmonics.
    """

    tlen = int(t_final/delta_t) + 1
    times = numpy.linspace(0, t_final, num=tlen)

    if inclination is not None:
        Y_plus, Y_cross = spher_harms(l, m, inclination)
    else:
        Y_plus, Y_cross = 1, 1

    common_factor = amp * numpy.exp(-times/tau)
    common_angle = (two_pi * f_0 * times) + phi

    hplus = Y_plus * common_factor * numpy.cos(common_angle)
    hcross = Y_cross * common_factor * numpy.sin(common_angle)

    return hplus, hcross

def fd_damped_sinusoid(f_0, tau, amp, phi, delta_f, f_lower, f_final, t_0=0.,
                       l=2, m=2, inclination=None):
    """Return a frequency domain damped sinusoid (plus and cross polarizations)
    with central frequency f_0, damping time tau, amplitude amp and phase phi.
    The l, m, and inclination parameters are used for the spherical harmonics.
    """

    if not f_lower:
        f_lower = delta_f
        kmin = 0
    else:
        kmin = int(f_lower / delta_f)

    # Create output vector with appropriate size
    outplus, outcross = fd_output_vector(f_0, tau, delta_f, f_final)
    freqs = outplus.sample_frequencies[kmin:]

    if inclination is not None:
        Y_plus, Y_cross = spher_harms(l, m, inclination)
    else:
        Y_plus, Y_cross = 1, 1

    denominator = 1 + (4j * pi * freqs * tau) - \
        (4 * pi_sq * (freqs*freqs - f_0*f_0) * tau*tau)
    norm = amp * tau / denominator
    if t_0 != 0:
        time_shift = numpy.exp(-1j * two_pi * freqs * t_0)
        norm *= time_shift
    A1 = (1 + 2j * pi * freqs * tau)
    A2 = two_pi * f_0 * tau

    # Analytical expression for the Fourier transform of the ringdown (damped sinusoid)
    hp_tilde = norm * Y_plus * (A1 * numpy.cos(phi) - A2 * numpy.sin(phi))
    hc_tilde = norm * Y_cross * (A1 * numpy.sin(phi) + A2 * numpy.cos(phi))

    outplus.data[kmin:] = hp_tilde
    outcross.data[kmin:] = hc_tilde

    return outplus, outcross

######################################################
#### Base multi-mode for all approximants
######################################################

def multimode_base(input_params):
    """Return a superposition of damped sinusoids in either time or frequency
    domains with parameters set by input_params.
    """
    input_params['lmns'] = format_lmns(input_params['lmns'])
    amps, phis = lm_amps_phases(**input_params)
    if 'final_mass' in input_params.keys():
        freqs, taus = get_lm_f0tau_allmodes(input_params['final_mass'],
                        input_params['final_spin'], input_params['lmns'])
        norm = Kerr_factor(input_params['final_mass'],
            input_params['distance']) if 'distance' in input_params.keys() \
            else 1.
    else:
        freqs, taus = lm_freqs_taus(**input_params)
        norm = 1.

    if 'delta_t' in input_params.keys():
        outplus, outcross = td_output_vector(freqs, taus,
                            input_params['taper'], input_params['delta_t'],
                            input_params['t_final'])
        for lmn in freqs:
            hplus, hcross = td_damped_sinusoid(freqs[lmn], taus[lmn],
                            amps[lmn], phis[lmn], outplus.delta_t,
                            outplus.sample_times[-1], int(lmn[0]), int(lmn[1]),
                            input_params['inclination'])
            if input_params['taper'] and outplus.delta_t < taus[lmn]:
                taper_hp, taper_hc = apply_taper(taus[lmn], amps[lmn],
                                     phis[lmn], outplus.delta_t, int(lmn[0]),
                                     int(lmn[1]), input_params['inclination'])
                t0 = -int(outplus.start_time * outplus.sample_rate)
                outplus[t0-len(taper_hp):t0].data += taper_hp
                outplus[t0:].data += hplus
                outcross[t0-len(taper_hc):t0].data += taper_hc
                outcross[t0:].data += hcross
            elif input_params['taper'] and outplus.delta_t > taus[lmn]:
                # This mode has taper duration < delta_t, do not apply taper
                t0 = -int(outplus.start_time * outplus.sample_rate)
                outplus[t0:].data += hplus
                outcross[t0:].data += hcross
            else:
                outplus.data += hplus
                outcross.data += hcross
    elif 'delta_f' in input_params.keys():
        outplus, outcross = fd_output_vector(freqs, taus,
                            input_params['delta_f'], input_params['f_final'])
        for lmn in freqs:
            hplus, hcross = fd_damped_sinusoid(freqs[lmn], taus[lmn],
                            amps[lmn], phis[lmn], outplus.delta_f,
                            input_params['f_lower'],
                            outplus.sample_frequencies[-1], int(lmn[0]),
                            int(lmn[1]), input_params['inclination'])
            outplus.data += hplus.data
            outcross.data += hcross.data

    return norm * outplus, norm * outcross

######################################################
#### Approximants
######################################################

def get_td_from_final_mass_spin(template=None, **kwargs):
    """Return time domain ringdown with all the modes specified.

    Parameters
    ----------
    template : object
        An object that has attached properties. This can be used to substitute
        for keyword arguments. A common example would be a row in an xml table.
    taper : {False, bool}, optional
        Add a rapid ringup with timescale tau/10 at the beginning of the
        waveform to avoid the abrupt turn on of the ringdown.
        Each mode and overtone will have a different taper depending on its tau,
        the final taper being the superposition of all the tapers.
    distance : {None, float}, optional
        Luminosity distance of the system. If specified, the returned ringdown
        will include the Kerr factor (final_mass/distance).
    final_mass : float
        Mass of the final black hole in solar masses.
    final_spin : float
        Dimensionless spin of the final black hole.
    lmns : list
        Desired lmn modes as strings (lm modes available: 22, 21, 33, 44, 55).
        The n specifies the number of overtones desired for the corresponding
        lm pair (maximum n=8).
        Example: lmns = ['223','331'] are the modes 220, 221, 222, and 330
    amp220 : float
        Amplitude of the fundamental 220 mode.  Always required, even if 220
        mode has not been selected. Note that if distance is given,
        this parameter will have a completely different order of magnitude.
        See table II in https://arxiv.org/abs/1107.0854 for an estimate.
    amplmn : float
        Fraction of the amplitude of the lmn overtone relative to the
        fundamental mode, i.e. amplmn/amp220. Provide as many as the number
        of selected subdominant modes.
    philmn : float
        Phase of the lmn overtone, as many as the number of modes. Should also
        include the information from the azimuthal angle, philmn=(phi + m*Phi).
    inclination : float
        Inclination of the system in radians (for the spherical harmonics).
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
    hplus : TimeSeries
        The plus phase of a ringdown with the lm modes specified and
        n overtones in time domain.
    hcross : TimeSeries
        The cross phase of a ringdown with the lm modes specified and
        n overtones in time domain.
    """

    input_params = props(template, mass_spin_required_args, td_args, **kwargs)

    return multimode_base(input_params)

def get_fd_from_final_mass_spin(template=None, **kwargs):
    """Return frequency domain ringdown with all the modes specified.

    Parameters
    ----------
    template : object
        An object that has attached properties. This can be used to substitute
        for keyword arguments. A common example would be a row in an xml table.
    distance : {None, float}, optional
        Luminosity distance of the system. If specified, the returned ringdown
        will include the Kerr factor (final_mass/distance).
    final_mass : float
        Mass of the final black hole in solar masses.
    final_spin : float
        Dimensionless spin of the final black hole.
    lmns : list
        Desired lmn modes as strings (lm modes available: 22, 21, 33, 44, 55).
        The n specifies the number of overtones desired for the corresponding
        lm pair (maximum n=8).
        Example: lmns = ['223','331'] are the modes 220, 221, 222, and 330
    amp220 : float
        Amplitude of the fundamental 220 mode.  Always required, even if 220
        mode has not been selected. Note that if distance is given,
        this parameter will have a completely different order of magnitude.
        See table II in https://arxiv.org/abs/1107.0854 for an estimate.
    amplmn : float
        Fraction of the amplitude of the lmn overtone relative to the
        fundamental mode, i.e. amplmn/amp220. Provide as many as the number
        of selected subdominant modes.
    philmn : float
        Phase of the lmn overtone, as many as the number of modes. Should also
        include the information from the azimuthal angle, philmn=(phi + m*Phi).
    inclination : float
        Inclination of the system in radians (for the spherical harmonics).
    delta_f : {None, float}, optional
        The frequency step used to generate the ringdown.
        If None, it will be set to the inverse of the time at which the
        amplitude is 1/1000 of the peak amplitude (the minimum of all modes).
    f_lower : {None, float}, optional
        The starting frequency of the output frequency series.
        If None, it will be set to delta_f.
    f_final : {None, float}, optional
        The ending frequency of the output frequency series.
        If None, it will be set to the frequency at which the amplitude
        is 1/1000 of the peak amplitude (the maximum of all modes).

    Returns
    -------
    hplustile : FrequencySeries
        The plus phase of a ringdown with the lm modes specified and
        n overtones in frequency domain.
    hcrosstilde : FrequencySeries
        The cross phase of a ringdown with the lm modes specified and
        n overtones in frequency domain.
    """

    input_params = props(template, mass_spin_required_args, fd_args, **kwargs)

    return multimode_base(input_params)

def get_td_from_freqtau(template=None, **kwargs):
    """Return time domain ringdown with all the modes specified.

    Parameters
    ----------
    template : object
        An object that has attached properties. This can be used to substitute
        for keyword arguments. A common example would be a row in an xml table.
    taper : {False, bool}, optional
        Add a rapid ringup with timescale tau/10 at the beginning of the
        waveform to avoid the abrupt turn on of the ringdown.
        Each mode and overtone will have a different taper depending on its tau,
        the final taper being the superposition of all the tapers.
    lmns : list
        Desired lmn modes as strings (lm modes available: 22, 21, 33, 44, 55).
        The n specifies the number of overtones desired for the corresponding
        lm pair (maximum n=8).
        Example: lmns = ['223','331'] are the modes 220, 221, 222, and 330
    f_lmn : float
        Central frequency of the lmn overtone, as many as number of modes.
    tau_lmn : float
        Damping time of the lmn overtone, as many as number of modes.
    amp220 : float
        Amplitude of the fundamental 220 mode. Note that if distance is given,
        this parameter will have a completely different order of magnitude.
        See table II in https://arxiv.org/abs/1107.0854 for an estimate.
        Always required, even if 220 mode has not been selected.
    amplmn : float
        Fraction of the amplitude of the lmn overtone relative to the
        fundamental mode, i.e. amplmn/amp220. Provide as many as the number
        of selected subdominant modes.
    philmn : float
        Phase of the lmn overtone, as many as the number of modes. Should also
        include the information from the azimuthal angle, philmn=(phi + m*Phi).
    inclination : float
        Inclination of the system in radians (for the spherical harmonics).
        If None, the spherical harmonics will be set to 1.
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
    hplus : TimeSeries
        The plus phase of a ringdown with the lm modes specified and
        n overtones in time domain.
    hcross : TimeSeries
        The cross phase of a ringdown with the lm modes specified and
        n overtones in time domain.
    """

    input_params = props(template, freqtau_required_args, td_args, **kwargs)

    return multimode_base(input_params)

def get_fd_from_freqtau(template=None, **kwargs):
    """Return frequency domain ringdown with all the modes specified.

    Parameters
    ----------
    template : object
        An object that has attached properties. This can be used to substitute
        for keyword arguments. A common example would be a row in an xml table.
    lmns : list
        Desired lmn modes as strings (lm modes available: 22, 21, 33, 44, 55).
        The n specifies the number of overtones desired for the corresponding
        lm pair (maximum n=8).
        Example: lmns = ['223','331'] are the modes 220, 221, 222, and 330
    f_lmn : float
        Central frequency of the lmn overtone, as many as number of modes.
    tau_lmn : float
        Damping time of the lmn overtone, as many as number of modes.
    amp220 : float
        Amplitude of the fundamental 220 mode. Always required, even if 220
        mode has not been selected.
    amplmn : float
        Fraction of the amplitude of the lmn overtone relative to the
        fundamental mode, i.e. amplmn/amp220. Provide as many as the number
        of selected subdominant modes.
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
    f_lower : {None, float}, optional
        The starting frequency of the output frequency series.
        If None, it will be set to delta_f.
    f_final : {None, float}, optional
        The ending frequency of the output frequency series.
        If None, it will be set to the frequency at which the amplitude
        is 1/1000 of the peak amplitude (the maximum of all modes).

    Returns
    -------
    hplustilde : FrequencySeries
        The plus phase of a ringdown with the lm modes specified and
        n overtones in frequency domain.
    hcrosstilde : FrequencySeries
        The cross phase of a ringdown with the lm modes specified and
        n overtones in frequency domain.
    """

    input_params = props(template, freqtau_required_args, fd_args, **kwargs)

    return multimode_base(input_params)

# Approximant names ###########################################################
ringdown_fd_approximants = {'FdQNMfromFinalMassSpin': get_fd_from_final_mass_spin,
                            'FdQNMfromFreqTau': get_fd_from_freqtau}
ringdown_td_approximants = {'TdQNMfromFinalMassSpin': get_td_from_final_mass_spin,
                            'TdQNMfromFreqTau': get_td_from_freqtau}
