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
from gwsignal.core import waveform as wfm

from pycbc import libutils, pnutils
from pycbc.types import FrequencySeries, TimeSeries

from . import parameters
from .waveform import _check_lal_pars, check_args, props
from .waveform.gwsignal_utils import to_gwsignal_dict

lalsimulation = libutils.import_optional('lalsimulation')


## FIXME: How to actually include extra-GR parameters?
def _lalsim_fd_waveform(**p):
    # convert paramaters to gwsignal parameters
    p_gws = to_gwsignal_dict(p)
    gen = wfm.LALCompactBinaryCoalescenceGenerator(p['approximant'])
    hlm = wfm.GenerateFDModes(python_dict, gen)
    fs = hlm.pop('frequency_array')
    pidx = fs >= 0
    hlms = {}
    for k, v in hlm.items():
        hlms[k[0], k[1]] = v[pidx].to_pycbc()

    return hlms


## FIXME: How to actually include extra-GR parameters?
def _lalsim_td_waveform(**p):
    # convert paramaters to gwsignal parameters
    p_gws = to_gwsignal_dict(p)
    gen = wfm.LALCompactBinaryCoalescenceGenerator(p['approximant'])
    hlm = wfm.GenerateTDModes(python_dict, gen)
    hlms = {}
    for k, v in hlm.items():
        hlms[k[0], k[1]] = v.to_pycbc()
    return hlm


def _formatdocstr(docstr):
    """Utility for formatting docstrings with parameter information.
    """
    return docstr.format(
        **{
            _p[1]: getattr(parameters, _p[1]).docstr(
                prefix="    ", include_label=False).lstrip(' ')
            for _p in Formatter().parse(docstr) if _p[1] is not None
        })


def _formatdocstrlist(docstr, paramlist, skip_params=None):
    """Utility for formatting docstrings with parameter information.
    """
    if skip_params is None:
        skip_params = []
    pl = '\n'.join([
        _p.docstr(prefix="    ", include_label=False) for _p in paramlist
        if _p not in skip_params
    ]).lstrip(' ')
    return docstr.format(params=pl)


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
        gives the plus, the negative imaginary part the cross.
    """
    out = None
    for mode in hlms:
        l, m = mode
        hlm = hlms[l, m]
        ylm = lal.SpinWeightedSphericalHarmonic(inclination, phi, -2, l, m)
        if out is None:
            out = ylm * hlm
        else:
            out += ylm * hlm
    return out


def default_modes(approximant):
    """Returns the default modes for the given approximant.
    """
    # FIXME: this should be replaced to a call to a lalsimulation function,
    # whenever that's added
    if approximant in ['IMRPhenomXPHM', 'IMRPhenomXHM']:
        # according to arXiv:2004.06503
        ma = [(2, 2), (2, 1), (3, 3), (3, 2), (4, 4)]
        # add the -m modes
        ma += [(l, -m) for l, m in ma]
    elif approximant in ['IMRPhenomPv3HM', 'IMRPhenomHM']:
        # according to arXiv:1911.06050
        ma = [(2, 2), (2, 1), (3, 3), (3, 2), (4, 4), (4, 3)]
        # add the -m modes
        ma += [(l, -m) for l, m in ma]
    elif approximant.startswith('NRSur7dq4'):
        # according to arXiv:1905.09300
        ma = [(l, m) for l in [2, 3, 4] for m in range(-l, l + 1)]
    else:
        raise ValueError("I don't know what the default modes are for "
                         "approximant {}, sorry!".format(approximant))
    return ma


def get_glm(l, m, theta):
    r"""The maginitude of the :math:`{}_{-2}Y_{\ell m}`.

    The spin-weighted spherical harmonics can be written as
    :math:`{}_{-2}Y_{\ell m}(\theta, \phi) = g_{\ell m}(\theta)e^{i m \phi}`.
    This returns the `g_{\ell m}(\theta)` part. Note that this is real.

    Parameters
    ----------
    l : int
        The :math:`\ell` index of the spherical harmonic.
    m : int
        The :math:`m` index of the spherical harmonic.
    theta : float
        The polar angle (in radians).

    Returns
    -------
    float :
        The amplitude of the harmonic at the given polar angle.
    """
    return lal.SpinWeightedSphericalHarmonic(theta, 0., -2, l, m).real


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
        Dictionary of ``(l, m)`` -> ``(h_+, -h_x)`` ``TimeSeries``.
    """
    laldict = _check_lal_pars(params)
    ret = lalsimulation.SimInspiralPrecessingNRSurModes(
        params['delta_t'], params['mass1'] * lal.MSUN_SI,
        params['mass2'] * lal.MSUN_SI, params['spin1x'], params['spin1y'],
        params['spin1z'], params['spin2x'], params['spin2y'], params['spin2z'],
        params['f_lower'], params['f_ref'],
        params['distance'] * 1e6 * lal.PC_SI, laldict,
        getattr(lalsimulation, params['approximant']))
    hlms = {}
    while ret:
        hlm = TimeSeries(ret.mode.data.data,
                         delta_t=ret.mode.deltaT,
                         epoch=ret.mode.epoch)
        hlms[ret.l, ret.m] = (hlm.real(), hlm.imag())
        ret = ret.next
    return hlms


get_nrsur_modes.__doc__ = _formatdocstr(get_nrsur_modes.__doc__)


def get_imrphenomxh_modes(**params):
    """Generates ``IMRPhenomXHM`` waveforms mode-by-mode. """
    approx = params['approximant']
    if not approx.startswith('IMRPhenomX'):
        raise ValueError("unsupported approximant")
    mode_array = params.pop('mode_array', None)
    if mode_array is None:
        mode_array = default_modes(approx)
    if 'f_final' not in params:
        # setting to 0 will default to ringdown frequency
        params['f_final'] = 0.
    hlms = {}
    for (l, m) in mode_array:
        params['mode_array'] = [(l, m)]
        laldict = _check_lal_pars(params)
        hpos, hneg = lalsimulation.SimIMRPhenomXPHMOneMode(
            l, m, params['mass1'] * lal.MSUN_SI, params['mass2'] * lal.MSUN_SI,
            params['spin1x'], params['spin1y'], params['spin1z'],
            params['spin2x'], params['spin2y'], params['spin2z'],
            params['distance'] * 1e6 * lal.PC_SI, params['coa_phase'],
            params['delta_f'], params['f_lower'], params['f_final'],
            params['f_ref'], laldict)
        hpos = FrequencySeries(hpos.data.data,
                               delta_f=hpos.deltaF,
                               epoch=hpos.epoch)
        hneg = FrequencySeries(hneg.data.data,
                               delta_f=hneg.deltaF,
                               epoch=hneg.epoch)
        if return_posneg:
            hlms[l, m] = (hpos, hneg)
        else:
            # convert to ulm, vlm
            ulm = 0.5 * (hpos + hneg.conj())
            vlm = 0.5j * (hneg.conj() - hpos)
            hlms[l, m] = (ulm, vlm)
    return hlms


_mode_waveform_td = {
    'NRSur7dq4': get_nrsur_modes,
    'IMRPhenomXHM': _lalsim_fd_waveform,
    'IMRPhenomXPHM': _lalsim_fd_waveform,
    'SEOBNRv4HM_ROM': _lalsim_fd_waveform,
    'SEOBNRv4HM': _lalsim_td_waveform,
    'SEOBNRv4PHM': _lalsim_td_waveform,
    'IMRPhenomTHM': _lalsim_td_waveform,
    'IMRPhenomTPHM': _lalsim_td_waveform,
}

# Remove commented out once IMRPhenomX one mode is fixed
_mode_waveform_fd = {  #'IMRPhenomXHM': get_imrphenomhm_modes,
    #'IMRPhenomXPHM' : get_imrphenomhm_modes,
}

def fd_waveform_mode_approximants():
    """Frequency domain approximants that will return separate modes."""
    return sorted(_mode_waveform_fd.keys())


def td_waveform_mode_approximants():
    """Time domain approximants that will return separate modes."""
    return sorted(_mode_waveform_td.keys())


def get_fd_waveform_modes(template=None, **kwargs):
    r"""Generates frequency domain waveforms, but does not sum over the modes.

    The returned values are the frequency-domain equivalents of the real and
    imaginary parts of the complex :math:`\mathfrak{{h}}_{{\ell m}}(t)` time
    series. In other words, the returned values are equivalent to the Fourier
    Transform of the two time series returned by
    :py:func:`get_td_waveform_modes`; see that function for more details.

    Parameters
    ----------
    template: object
        An object that has attached properties. This can be used to subsitute
        for keyword arguments.
    {params}

    Returns
    -------
    ulm : dict
        Dictionary of mode tuples -> fourier transform of the real part of the
        hlm time series, as a :py:class:`pycbc.types.FrequencySeries`.
    vlm : dict
        Dictionary of mode tuples -> fourier transform of the imaginary part of
        the hlm time series, as a :py:class:`pycbc.types.FrequencySeries`.
    """
    params = props(template, **kwargs)
    required = parameters.fd_required
    check_args(params, required)
    apprx = params['approximant']
    if apprx not in _mode_waveform_fd:
        raise ValueError("I don't support approximant {}, sorry".format(apprx))
    return _mode_waveform_fd[apprx](**params)


get_fd_waveform_modes.__doc__ = _formatdocstrlist(
    get_fd_waveform_modes.__doc__,
    parameters.fd_waveform_params,
    skip_params=['inclination', 'coa_phase'])


def get_td_waveform_modes(template=None, **kwargs):
    r"""Generates time domain waveforms, but does not sum over the modes.

    The returned values are the real and imaginary parts of the complex
    :math:`\mathfrak{{h}}_{{\ell m}}(t)`. These are defined such that the plus
    and cross polarizations :math:`h_{{+,\times}}` are:

    .. math::

       h_{{+,\times}}(\theta, \phi; t) = (\Re, -\Im) \sum_{{\ell m}}
        {{}}_{{-2}}Y_{{\ell m}}(\theta, \phi) \mathfrak{{h}}_{{\ell m}}(t).


    Parameters
    ----------
    template: object
        An object that has attached properties. This can be used to subsitute
        for keyword arguments.
    {params}

    Returns
    -------
    ulm : dict
        Dictionary of mode tuples -> real part of the hlm, as a
        :py:class:`pycbc.types.TimeSeries`.
    vlm : dict
        Dictionary of mode tuples -> imaginary part of the hlm, as a
        :py:class:`pycbc.types.TimeSeries`.
    """
    params = props(template, **kwargs)
    required = parameters.td_required
    check_args(params, required)
    apprx = params['approximant']
    if apprx not in _mode_waveform_td:
        raise ValueError("I don't support approximant {}, sorry".format(apprx))
    return _mode_waveform_td[apprx](**params)


get_td_waveform_modes.__doc__ = _formatdocstrlist(
    get_td_waveform_modes.__doc__,
    parameters.td_waveform_params,
    skip_params=['inclination', 'coa_phase'])
