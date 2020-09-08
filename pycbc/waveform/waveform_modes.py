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
import numpy
import lal
import lalsimulation

from .waveform import (props, _check_lal_pars, check_args)
from .waveform import get_td_waveform, get_fd_waveform
from . import parameters
from pycbc.types import TimeSeries


def _formatdocstr(docstr):
    """Utility for formatting docstrings with parameter information.
    """
    return docstr.format(
        **{_p[1]: getattr(parameters, _p[1]).docstr(
            prefix="    ", include_label=False).lstrip(' ')
           for _p in Formatter().parse(docstr) if _p[1] is not None
           })


def sum_modes(hlms, inclination, phi):
    """Applies spherical harmonics and sums modes to produce a plus and cross
    polarization.

    Parameters
    ----------
    hlms : dict
        Dictionary of ``(l, m)`` -> complex ``hlm``. The ``hlm`` may be a
        complex number or array, or complex ``TimeSeries``. All modes in the
        dictionary will be summed.
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
        ylm = lal.SpinWeightedSphericalHarmonic(inclination, phi, -2, ell, m)
        if out is None:
            out = ylm * hlm
        else:
            out += ylm * hlm
    # return the conjugate, since h = h_+ - ih_x
    return numpy.conj(out)


def default_modes(approximant):
    """Returns the default modes for the given approximant.
    """
    # FIXME: this should be replaced to a call to a lalsimulation function,
    # whenever that's added
    if approximant in ['IMRPhenomXPHM', 'IMRPhenomXHM']:
        # according to arXiv:2004.06503
        ma = [(2, 2), (2, 1), (3, 3), (3, 2), (4, 4)]
        # add the -m modes
        ma += [(ell, -m) for ell, m in ma]
        return ma 
    elif approximant in ['IMRPhenomPv3HM', 'IMRPhenomHM']:
        # according to arXiv:1911.06050
        ma = [(2, 2), (2, 1), (3, 3), (3, 2), (4, 4), (4, 3)]
        # add the -m modes
        ma += [(ell, -m) for ell, m in ma]
        return ma 
    elif approximant.startswith('NRSur7dq4'):
        # according to arXiv:1905.09300
        return [(ell, m) for ell in [2, 3, 4] for m in range(-ell, ell+1)]
    else:
        raise ValueError("I don't know what the default modes are for "
                         "approximant {}, sorry!".format(approximant))


def get_nrsur_modes(**params):
    """Generates NRSurrogate waveform mode-by-mode.

    All waveform parameters should be provided as keyword arguments.
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
        Dictionary of ``(l, m)`` -> ``(h_+, h_x)`` ``TimeSeries``.
    """
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
    # the NR surrogates use a different coordinate system than other lal
    # lal waveforms; to make standard, we need to shift the phase by -pi/2
    dphi = -numpy.pi/2
    while ret:
        hlm = TimeSeries(ret.mode.data.data, delta_t=ret.mode.deltaT,
                         epoch=ret.mode.epoch)
        # correct the phase to make LAL standard
        hlm *= numpy.exp(1j * ret.m * dphi)
        # store the conjugate, since h = h_+ - ih_x
        hlms[ret.l, ret.m] = (hlm.real(), -hlm.imag())
        ret = ret.next
    return hlms

get_nrsur_modes.__doc__ = _formatdocstr(get_nrsur_modes.__doc__)


def get_imrphenomx_modes(**params):
    """Generates ``IMRPhenomX[P]HM`` waveforms mode-by-mode.
    """
    raise NotImplementedError("Not currently supported")
    approx = params['approximant']
    if not approx.startswith('IMRPhenomX'):
        raise ValueError("unsupported approximant")
    mode_array = params.pop('mode_array', None)
    if mode_array is None:
        mode_array = default_modes(approx)
    laldict = _check_lal_pars(params)
    if 'f_final' not in params:
        # setting to 0 will default to ringdown frequency
        params['f_final'] = 0. 
    hlms = {}
    for (ell, m) in mode_array:
        hpos, hneg = lalsimulation.SimIMRPhenomXPHMOneMode(
            ell, m, 
            params['mass1']*lal.MSUN_SI,
            params['mass2']*lal.MSUN_SI,
            params['spin1x'], params['spin1y'], params['spin1z'],
            params['spin2x'], params['spin2y'], params['spin2z'],
            params['distance']*1e6*lal.PC_SI, params['coa_phase'],
            params['delta_f'], params['f_lower'], params['f_final'],
            params['f_ref'],
            laldict)
       hlms[ell, m] = hpos, hneg 
    return hlms


def get_imrphenomhm_modes(**kwargs):
    """Generates ``IMRPhenom[Pv3]HM`` waveforms mode-by-mode.
    """
    approx = kwargs['approximant']
    try:
        mode_array = kwargs['mode_array']
    except KeyError:
        mode_array = None
    if mode_array is None:
        mode_array = default_modes(approx)
    hlms = {}
    for mode in mode_array:
        ma = [mode]
        # PhenomPv3HM always needs the 2,2 mode, so we'll generate it then
        # subract it off
        if approx == 'IMRPhenomPv3HM' and mode != (2, 2):
            ma.append((2,2))
        kwargs.update({'mode_array': ma})
        hp, hc = get_fd_waveform(**kwargs)
        hlms[mode] = (hp, hc)
    # subtract off the 2,2 mode for PhenomPv3HM
    if approx == 'IMRPhenomPv3HM':
        try:
            hp22, hc22 = hlms[(2,2)]
        except KeyError:
            # 2,2 mode wasn't requested; generate it separately
            kwargs.update({'mode_array': (2, 2)})
            hp22, hc22 = get_fd_waveform(**kwargs) 
        for mode in hlms:
            if mode != (2, 2):
                hp, hc = hlms[mode]
                hp -= hp22
                hc -= hc22
    return hlms


_mode_waveform_td = {'NRSur7dq4': get_nrsur_modes,
                     'NRSur7dq2': get_nrsur_modes,
                     }


_mode_waveform_fd = {'IMRPhenomHM': get_imrphenomhm_modes,
                     'IMRPhenomPv3HM': get_imrphenomhm_modes,
                     'IMRPhenomXHM': get_imrphenomx_modes,
                     'IMRPhenomXPHM' : get_imrphenomx_modes,
                    }


def get_fd_waveform_modes(template=None, **kwargs):
    """Generates frequency domain waveforms, but does not sum over the modes.
    """
    params = props(template, **kwargs)
    required = parameters.fd_required
    check_args(params, required)
    try:
        return _mode_waveform_fd[params['approximant']](**params)
    except KeyError:
        raise ValueError("I don't support approximant {}, sorry"
                         .format(approx))


def get_td_waveform_modes(template=None, **kwargs):
    """Generates time domain waveforms, but does not sum over the modes.
    """
    params = props(template, **kwargs)
    required = parameters.fd_required
    check_args(params, required)
    try:
        return _mode_waveform_td[params['approximant']](**params)
    except KeyError:
        raise ValueError("I don't support approximant {}, sorry"
                         .format(approx))
