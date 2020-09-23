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

from pycbc.types import (TimeSeries, FrequencySeries)
from .waveform import (props, _check_lal_pars, check_args)
from . import parameters


def _formatdocstr(docstr):
    """Utility for formatting docstrings with parameter information.
    """
    return docstr.format(
        **{_p[1]: getattr(parameters, _p[1]).docstr(
            prefix="    ", include_label=False).lstrip(' ')
           for _p in Formatter().parse(docstr) if _p[1] is not None
           })


def _formatdocstrlist(docstr, paramlist, skip_params=None):
    """Utility for formatting docstrings with parameter information.
    """
    if skip_params is None:
        skip_params = []
    pl = '\n'.join([_p.docstr(prefix="    ", include_label=False)
                    for _p in paramlist if _p not in skip_params]).lstrip(' ')
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
        ell, m = mode
        hlm = hlms[ell, m]
        ylm = lal.SpinWeightedSphericalHarmonic(inclination, phi, -2, ell, m)
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
        ma += [(ell, -m) for ell, m in ma]
    elif approximant in ['IMRPhenomPv3HM', 'IMRPhenomHM']:
        # according to arXiv:1911.06050
        ma = [(2, 2), (2, 1), (3, 3), (3, 2), (4, 4), (4, 3)]
        # add the -m modes
        ma += [(ell, -m) for ell, m in ma]
    elif approximant.startswith('NRSur7dq4'):
        # according to arXiv:1905.09300
        ma = [(ell, m) for ell in [2, 3, 4] for m in range(-ell, ell+1)]
    else:
        raise ValueError("I don't know what the default modes are for "
                         "approximant {}, sorry!".format(approximant))
    return ma


def get_glm(ell, m, theta):
    """The inclination part of the ylm."""
    return lal.SpinWeightedSphericalHarmonic(theta, 0., -2, ell, m).real


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
        hlm = TimeSeries(ret.mode.data.data, delta_t=ret.mode.deltaT,
                         epoch=ret.mode.epoch)
        hlms[ret.l, ret.m] = (hlm.real(), hlm.imag())
        ret = ret.next
    return hlms


get_nrsur_modes.__doc__ = _formatdocstr(get_nrsur_modes.__doc__)


def get_imrphenomx_modes(return_posneg=False, **params):
    """Generates ``IMRPhenomX[P]HM`` waveforms mode-by-mode.
    """
    # FIXME: raising not implemented error because this currently does not
    # work. The issue is the OneMode function adds the +/- m modes together
    # automatically. Remove once this is fixed in lalsimulation, and/or I
    # figure out a reliable way to separate the +/-m modes.
    raise NotImplementedError("Currently not implemented")
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
    for (ell, m) in mode_array:
        params['mode_array'] = [(ell, m)]
        laldict = _check_lal_pars(params)
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
        hpos = FrequencySeries(hpos.data.data, delta_f=hpos.deltaF,
                               epoch=hpos.epoch)
        hneg = FrequencySeries(hneg.data.data, delta_f=hneg.deltaF,
                               epoch=hneg.epoch)
        if return_posneg:
            hlms[ell, m] = (hpos, hneg)
        else:
            # convert to ulm, vlm
            ulm = 0.5 * (hpos + hneg.conj())
            vlm = 0.5j * (hneg.conj() - hpos)
            hlms[ell, m] = (ulm, vlm)
    return hlms


def phenom_l0frame_to_jframe(approximant, mass1, mass2, f_ref, phiref=0.,
                             inclination=0.,
                             spin1x=0., spin1y=0., spin1z=0.,
                             spin2x=0., spin2y=0., spin2z=0.):
    r"""Converts L0- to J-frame parameters used by IMRPhenomP waveforms.

    Parameters
    ----------
    approximant : str
        Name of the approximant. Must be one of the IMRPhenom approximants.
    {mass1}
    {mass2}
    {f_ref}
    phiref : float
        Reference phase.
    inclination : float
        Angle between the orbital angular momentum at ``f_ref`` and the line
        of sight.
    {spin1x}
    {spin1y}
    {spin1z}
    {spin2x}
    {spin2y}
    {spin2z}

    Returns
    -------
    dict :
        Dictionary of:
        * thetajn : float
            Angle between the line of sight and the total angular momentum.
            This is the thing that goes into the polar angle part of the
            spherical harmonics.
        * alpha0 : float
            Azimuthal angle in the J frame. This is the thing that goes into
            the azimuthal part of the spherical harmonics.
        * phi_aligned : float
            Beats me.
        * zeta_polarization : float
            Another mystery.
        * spin1_l : float
            Component of the larger object's spin that is aligned with the
            orbital angular momentum.
        * spin2_l : float
            Component of the smaller object's spin that is aligned with the
            orbital angular momentum.
        * {chi_p}
    """
    phenomv = approximant.replace('HM', '') + '_V'
    spin1_l, spin2_l, chip, thetajn, alpha0, phi_aligned, zeta_pol = \
        lalsimulation.SimIMRPhenomPCalculateModelParametersFromSourceFrame(
            mass1*lal.MSUN_SI, mass2*lal.MSUN_SI, f_ref, phiref, inclination,
            spin1x, spin1y, spin1z, spin2x, spin2y, spin2z,
            getattr(lalsimulation, phenomv))
    out = {'thetajn': thetajn,
           'alpha0': alpha0,
           'phi_aligned': phi_aligned,
           'zeta_polarization': zeta_pol,
           'spin1_l': spin1_l,
           'spin2_l': spin2_l,
           'chi_p': chip}
    return out


phenom_l0frame_to_jframe.__doc__ = _formatdocstr(
    phenom_l0frame_to_jframe.__doc__)


def l0frame_to_jframe(mass1, mass2, f_ref, phiref=0., inclination=0.,
                      spin1x=0., spin1y=0., spin1z=0.,
                      spin2x=0., spin2y=0., spin2z=0.):
    """Converts L0-frame parameters to J-frame.

    Parameters
    ----------
    {mass1}
    {mass2}
    {f_ref}
    phiref : float
        The orbital phase at ``f_ref``.
    {inclination}
    {spin1x}
    {spin1y}
    {spin1z}
    {spin2x}
    {spin2y}
    {spin2z}

    Returns
    -------
    dict :
        Dictionary of:
        * thetajn : float
            Angle between the line of sight and the total angular momentume J.
        * phijl : float
            Azimuthal angle of L on its cone about J.
        * {spin1_a}
        * {spin2_a}
        * spin1_polar : float
            Angle between L and the spin magnitude of the larger object.
        * spin2_polar : float
            Angle betwen L and the spin magnitude of the smaller object.
        * spin12_deltaphi : float
            Difference between the azimuthal angles of the spin of the larger
            object (S1) and the spin of the smaller object (S2).
    """
    # Note: unlike other LALSimulation functions, this one takes masses in
    # solar masses
    thetajn, phijl, s1pol, s2pol, s12_deltaphi, spin1_a, spin2_a = \
        lalsimulation.SimInspiralTransformPrecessingWvf2PE(
            inclination, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z,
            mass1, mass2, f_ref, phiref)
    out = {'thetajn': thetajn,
           'phijl' : phijl,
           'spin1_polar': s1pol,
           'spin2_polar': s2pol,
           'spin12_deltaphi': s12_deltaphi,
           'spin1_a': spin1_a,
           'spin2_a': spin2_a}
    return out


l0frame_to_jframe.__doc__ = _formatdocstr(l0frame_to_jframe.__doc__)


def jframe_to_l0frame(mass1, mass2, f_ref, phiref=0., thetajn=0., phijl=0.,
                      spin1_a=0., spin2_a=0.,
                      spin1_polar=0., spin2_polar=0.,
                      spin12_deltaphi=0.):
    """Converts J-frame parameters into L0 frame.

    Parameters
    ----------
    {mass1}
    {mass2}
    {f_ref}
    thetajn : float
        Angle between the line of sight and the total angular momentume J.
    phijl : float
        Azimuthal angle of L on its cone about J.
    {spin1_a}
    {spin2_a}
    spin1_polar : float
        Angle between L and the spin magnitude of the larger object.
    spin2_polar : float
        Angle betwen L and the spin magnitude of the smaller object.
    spin12_deltaphi : float
        Difference between the azimuthal angles of the spin of the larger
        object (S1) and the spin of the smaller object (S2).

    Returns
    -------
    dict :
        Dictionary of:
        * {inclination}
        * {spin1x}
        * {spin1y}
        * {spin1z}
        * {spin2x}
        * {spin2y}
        * {spin2z}
    """
    inclination, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z = \
        lalsimulation.SimInspiralTransformPrecessingNewInitialConditions(
            thetajn, phijl, spin1_polar, spin2_polar, spin12_deltaphi,
            spin1_a, spin2_a, mass1*lal.MSUN_SI, mass2*lal.MSUN_SI, f_ref,
            phiref)
    out = {'inclination': inclination,
           'spin1x': spin1x,
           'spin1y': spin1y,
           'spin1z': spin1z,
           'spin2x': spin2x,
           'spin2y': spin2y,
           'spin2z': spin2z}
    return out


jframe_to_l0frame.__doc__ = _formatdocstr(jframe_to_l0frame.__doc__)


_mode_waveform_td = {'NRSur7dq4': get_nrsur_modes,
                     }


# Remove commented out once IMRPhenomX one mode is fixed
_mode_waveform_fd = {#'IMRPhenomXHM': get_imrphenomhm_modes,
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
        raise ValueError("I don't support approximant {}, sorry"
                         .format(apprx))
    return _mode_waveform_fd[apprx](**params)


get_fd_waveform_modes.__doc__ = _formatdocstrlist(
    get_fd_waveform_modes.__doc__, parameters.fd_waveform_params,
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
        Dictionary of mode tuples -> real part of the hlm, as a
        :py:class:`pycbc.types.TimeSeries`.
    """
    params = props(template, **kwargs)
    required = parameters.fd_required
    check_args(params, required)
    apprx = params['approximant']
    if apprx not in _mode_waveform_td:
        raise ValueError("I don't support approximant {}, sorry"
                         .format(apprx))
    return _mode_waveform_td[apprx](**params)


get_td_waveform_modes.__doc__ = _formatdocstrlist(
    get_td_waveform_modes.__doc__, parameters.td_waveform_params,
    skip_params=['inclination', 'coa_phase'])
