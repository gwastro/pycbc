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

def get_fd_ringdown(template=None,**kwargs):
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
    Amp : {1, float}
        The amplitude of the ringdown (constant for now).
    delta_f : float
        The frequency step used to generate the ringdown.
    f_lower: float
        The starting frequency of the output frequency series.
    f_final : float
        The ending frequency of the output frequency series.

        Returns
    -------
    hplustilde: FrequencySeries
        The plus phase of the ringdown in frequency domain.
    hcrosstilde: FrequencySeries
        The cross phase of the ringdown in frequency domain.
    """

    input_params = props_ringdown(template,**kwargs)

    f_0 = input_params['f_0']
    f_lower = input_params['f_lower']
    f_final = input_params['f_final']
    delta_f = input_params['delta_f']
    tau = input_params['tau']
    t_0 = input_params['t_0']
    phi_0 = input_params['phi_0']
    Amp = input_params['Amp']

    pi = numpy.pi
    two_pi = 2 * numpy.pi
    pi_sq = numpy.pi * numpy.pi

    freqs = numpy.arange( f_lower, f_final, delta_f)
    kmin = int(f_lower / delta_f)
    kmax = int(f_final / delta_f)
    n = int(f_final / delta_f) + 1

    denominator = 1 + ( 4j * pi * freqs * tau ) - ( 4 * pi_sq * ( freqs*freqs - f_0*f_0) * tau*tau )
    time_shift = [ numpy.exp( -1j * two_pi * f * t_0 ) for f in freqs ]

    hp_tilde = Amp * tau * ( ( 1 + 2j * pi * freqs * tau ) * numpy.cos(phi_0) - two_pi * f_0 * tau * numpy.sin(phi_0) ) * time_shift / denominator
    hc_tilde = Amp * tau * ( ( 1 + 2j * pi * freqs * tau ) * numpy.sin(phi_0) + two_pi * f_0 * tau * numpy.cos(phi_0) ) * time_shift / denominator

    hplustilde = FrequencySeries(zeros(n, dtype=complex128), delta_f=delta_f)
    hcrosstilde = FrequencySeries(zeros(n, dtype=complex128), delta_f=delta_f)
    hplustilde.data[kmin:kmax] = hp_tilde
    hcrosstilde.data[kmin:kmax] = hc_tilde

    return hplustilde, hcrosstilde
