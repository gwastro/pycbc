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
from pycbc.types import FrequencySeries, complex128, zeros

default_ringdown_args = {'t_0':0, 'phi_0':0, 'Amp':1}

def props_ringdown(obj, **kwargs):
    """ DOCUMENT ME !!
    """
    pr = {}
    if obj is not None:
        if hasattr(obj, '__dict__'):
            pr = obj.__dict__
        elif hasattr(obj, '__slots__'):
            for slot in obj.__slots__:
                if hasattr(obj, slot):
                    pr[slot] = getattr(obj, slot)
        else:
            for name in dir(obj):
                try:
                    value = getattr(obj, name)
                    if not name.startswith('__') and not inspect.ismethod(value):
                        pr[name] = value
                except:
                    continue

    # Get the parameters to generate the waveform
    # Note that keyword arguments override values in the template object
    input_params = default_ringdown_args.copy()
    input_params.update(pr)
    input_params.update(kwargs)

    return input_params

def qnm_freq_decay(f_0, tau, decay):
    """Return the frequency at which the amplitude of the 
    ringdown falls to 1/decay of the peak amplitude.

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
        ringdown is 1/decay of the peak amplitude.
    """

    q_0 = numpy.pi * f_0 * tau
    decay_sq = decay*decay

    # Expression obtained analytically under the assumption
    # that decay^2, q_0^2 >> 1
    q_sq = (decay_sq + 4*q_0*q_0 + decay*numpy.sqrt(decay_sq + 16*q_0*q_0)) / 4.
    f_decay = numpy.sqrt(q_sq) / numpy.pi / tau
    return f_decay

def qnm_time_decay(tau, decay):
    """Return the time at which the amplitude of the 
    ringdown falls to 1/decay of the peak amplitude.

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
        ringdown is 1/decay of the peak amplitude.
    """

    t_decay = tau * numpy.log(decay)
    return t_decay

def get_fd_qnm(template=None, delta_f=None, f_lower=None, f_final=None, **kwargs):
    """Return a frequency domain ringdown.

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
    phi_0 : {0, float}, optional
        The initial phase of the ringdown.
    Amp : {1, float}, optional
        The amplitude of the ringdown (constant for now).
    delta_f : {None, float}, optional
        The frequency step used to generate the ringdown.
        If None, it will be set to the inverse of the time at which the
        amplitude is a 1/1000 of the peak amplitude.
    f_lower: {None, float}, optional
        The starting frequency of the output frequency series.
        If None, it will be set to delta_f.
    f_final : {None, float}, optional
        The ending frequency of the output frequency series.
        If None, it will be set to the frequency at which the amplitude is 
        1/100 of the peak amplitude.

    Returns
    -------
    hplustilde: FrequencySeries
        The plus phase of the ringdown in frequency domain.
    hcrosstilde: FrequencySeries
        The cross phase of the ringdown in frequency domain.
    """

    input_params = props_ringdown(template,**kwargs)

    f_0 = input_params['f_0']
    tau = input_params['tau']
    t_0 = input_params['t_0']
    phi_0 = input_params['phi_0']
    Amp = input_params['Amp']
    if delta_f is None:
        delta_f = 1. / qnm_time_decay(tau, 1000.)
    if f_lower is None:
        f_lower = delta_f
        kmin = 0
    else:
        kmin = int(f_lower / delta_f)
    if f_final is None:
        f_final = qnm_freq_decay(f_0, tau, 100.)
    kmax = int(f_final / delta_f)
    n = int(f_final / delta_f) + 1

    pi = numpy.pi
    two_pi = 2 * numpy.pi
    pi_sq = numpy.pi * numpy.pi

    freqs = numpy.arange( f_lower, f_final, delta_f)

    denominator = 1 + ( 4j * pi * freqs * tau ) - ( 4 * pi_sq * ( freqs*freqs - f_0*f_0) * tau*tau )
    time_shift = [ numpy.exp( -1j * two_pi * f * t_0 ) for f in freqs ]

    hp_tilde = Amp * tau * ( ( 1 + 2j * pi * freqs * tau ) * numpy.cos(phi_0) - two_pi * f_0 * tau * numpy.sin(phi_0) ) * time_shift / denominator
    hc_tilde = Amp * tau * ( ( 1 + 2j * pi * freqs * tau ) * numpy.sin(phi_0) + two_pi * f_0 * tau * numpy.cos(phi_0) ) * time_shift / denominator

    hplustilde = FrequencySeries(zeros(n, dtype=complex128), delta_f=delta_f)
    hcrosstilde = FrequencySeries(zeros(n, dtype=complex128), delta_f=delta_f)
    hplustilde.data[kmin:kmax] = hp_tilde
    hcrosstilde.data[kmin:kmax] = hc_tilde

    return hplustilde, hcrosstilde
