# Copyright (C) 2020  Collin Capano, Alex Nitz
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

"""Provides functions and utilities for generating waveforms mode-by-mode.
"""

from string import Formatter
import lal
import lalsimulation

from .waveform import props, _check_lal_pars
from . import parameters

from pycbc.types import TimeSeries

def sum_modes(hlms, inclination, phi):
    """Applies spherical harmonics and sums modes to produce a plus and cross
    polarization.

    Parameters
    ----------
    hlms : dict
        Dictionary of ``(l, m)`` -> complex ``hlm``. The ``hlm`` may be a
        complex number or array, complex ``TimeSeries``, or complex
        ``FrequencySeries``. All modes in the dictionary will be summed.
    inclination : float
        The inclination to use.
    phi : float
        The phase to use.

    Returns
    -------
    complex float or array
        The plus and cross polarization as a complex number. The real part
        gives the plus, the imaginary the cross.
    """
    out = None
    for mode in hlms:
        ell, m = mode
        hlm = hlms[ell, m]
        ylm = lal.SpinWeightedSphericalHarmonic(
            inclination, numpy.pi/2 - phi, -2,
            ell, m)
        if out is None:
            out = ylm * hlm
        else:
            out += ylm * hlm
    # return the conjugate, since h = h_+ - ih_x
    return numpy.conj(out)


def get_nrsur_modes(template=None, **kwargs):
    """Generates NRSurrogate waveform mode-by-mode.

    All waveform parameters are should be provided as keyword arguments.
    Recognized parameters are listed below. Unrecognized arguments are ignored.

    Parameters
    ----------
    template: object
        An object that has attached properties. This can be used to substitute
        for keyword arguments. A common example would be a row in an xml table.
    approximant : str
        The approximant to generate. Must be one of the ``NRSur*`` models.
    {delta_t}
    {mass1}
    {mass2}
    {spin1x}
    {spin1y}
    {spin1z}
    {spin2x}
    {spin2y}
    {spin2z}
    {f_lower}
    {f_ref}
    {distance}
    {mode_array}

    Returns
    -------
    dict :
        Dictionary of ``(l, m)`` -> complex ``TimeSeries`` giving each ``hlm``.
    """
    params = props(template, **kwargs)
    laldict = _check_lal_pars(params)
    ret = lalsimulation.SimInspiralPrecessingNRSurModes(
        params['delta_t'],
        params['mass1']*lal.MSUN_SI,
        params['mass2']*lal.MSUN_SI,
        params['spin1x'], params['spin1y'], params['spin1z'],
        params['spin2x'], params['spin2y'], params['spin2z'],
        params['f_lower'], params['f_ref'],
        params['distance']*1e6*lal.PC_SI, laldict,
        getattr(lalsimulation, params['approximant'])
    )
    hlms = {}
    while ret:
        hlms[ret.l, ret.m] =  TimeSeries(ret.mode.data.data,
                                     delta_t=ret.mode.deltaT,
                                     epoch=ret.mode.epoch)
        ret = ret.next
    return hlms

_mode_waveform_td = {'NRSur7dq4':get_nrsur_modes,
                     'NRSur7dq2':get_nrsur_modes,
                     }
_mode_waveform_fd = {}

def get_td_modes(template=None, **kwargs):
    """ Return all modes composing a time domain waveform
    """
    params = props(template, **kwargs)
    return _mode_waveform_td[params['approximant']](**params)

def get_fd_modes(template=None, **kwargs):
    """ Return all modes composing a frequency domain waveform
    """
    params = props(template, **kwargs)
    return _mode_waveform_fd[params['approximant']](**params)

# Collin, this doesn't run, needs fixing
#get_nrsur_modes.__doc__ = get_nrsur_modes.__doc__.format(
#    **{_p[1]: getattr(parameters, _p[1]).docstr(
#        prefix="    ", include_label=False).lstrip(' ')
#       for _p in Formatter().parse(get_nrsur_modes.__doc__)
#       })
