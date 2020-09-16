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
from pycbc.types import (TimeSeries, FrequencySeries)


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
            print('returning pos, neg')
            hlms[ell, m] = (hpos, hneg)
        else:
            print('returning ulm, vlm')
            # convert to ulm, vlm
            ulm = 0.5 * (hpos + hneg.conj())
            vlm = 0.5j * (hneg.conj() - hpos)
            hlms[ell, m] = (ulm, vlm)
    return hlms


def l0frame_to_jframe(approximant, mass1, mass2, f_ref, phiref=0.,
                      inclination=0.,
                      spin1x=0., spin1y=0., spin1z=0.,
                      spin2x=0., spin2y=0., spin2z=0.):
    r"""Converts L0 frame parameters to J frame.

    Only works for IMRPhenom approximants.

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
    thetajn : float
        Angle between the line of sight and the total angular momentum. This
        is the thing that goes into the polar angle part of the spherical
        harmonics.
    alpha0 : float
        Azimuthal angle in the J frame. This is the thing that goes into the
        azimuthal part of the spherical harmonics.
    phi_aligned : float
        Beats me.
    zeta_polarization : float
        Another mystery.
    spin1_l : float
        Component of the larger object's spin that is aligned with the orbital
        angular momentum.
    spin2_l : float
        Component of the smaller object's spin that is aligned with the orbital
        angular momentum.
    {chi_p}
    """
    phenomv = approximant.replace('HM', '') + '_V'
    spin1_l, spin2_l, chip, thetajn, alpha0, phi_aligned, zeta_polarization = \
        lalsimulation.SimIMRPhenomPCalculateModelParametersFromSourceFrame(
            mass1*lal.MSUN_SI, mass2*lal.MSUN_SI, f_ref, phiref, inclination,
            spin1x, spin1y, spin1z, spin2x, spin2y, spin2z,
            getattr(lalsimulation, phenomv))
    return thetajn, alpha0, phi_aligned, zeta_pol, spin1_l, spin2_l, chi_p


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
    phiref : float
    thetajn : float

    """
    inclination, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z = \
        lalsimulation.SimInspiralTransformPrecessingNewInitialConditions(
            thetajn, phijl, spin1_polar, spin2_polar, spin12_deltaphi,
            spin1_a, spin2_a, mass1*lal.MSUN_SI, mass2*lal.MSUN_SI, f_ref,
            phiref)
    return inclination, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z


def _shiftparams(theta, phi, **params):
    # XXX: ignoring spins
    fref = params['f_ref']
    if fref == 0:
        fref = params['f_lower']
    jparams = jframe_to_l0frame(
        params['mass1'], params['mass2'], fref,
        phiref=phi, thetajn=theta)
    params.update({'inclination': jparams[0], 'coa_phase': phi})
    return params


def get_imrphenomhm_modes(**params):
    """Generates ``IMRPhenom[Pv3]HM`` waveforms mode-by-mode.
    """
    approx = params['approximant']
    # remove inclination and phase if they were included
    inc = params.pop('inclination', 0.)
    phi = params.pop('coa_phase', 0.)
    if inc != 0.:
        raise ValueError("non zero inclination given")
    if phi != 0.:
        raise ValueError("non zero coalescence phase given")
    mode_array = params.pop('mode_array', None)
    if mode_array is None:
        mode_array = default_modes(approx)
    # make sure everything is specified as +m modes
    mode_array = list(set([(ell, abs(m)) for (ell, m) in mode_array]))
    # Our goal below is to set the polar and azimuthal angles of the
    # spherical harmonics. In the Phenom waveforms, this corresponds to
    # thetaJN and alpha0, respectively. However, since we're using the
    # ChooseFDWaveform interface, we need to put these in he L0 frame (which
    # (the waveform generator will then undo... sigh)
    # XXX: ignoring spins for now

    # we need the in plane spins for some modes
    #spin1perp = (params['spin1x']**2 + params['spin1y']**2.)**0.5
    #spin1az = numpy.arctan2(params['spin1y'], params['spin1x'])
    #spin2perp = (params['spin2x']**2 + params['spin2y']**2.)**0.5
    #spin2az = numpy.arctan2(params['spin2y'], params['spin2x'])

    # We will need to isolate the + and - m modes from each other. How to do
    # that depends on the mode. For:
    # m = 2:
    #   generate once face-on (-m is zero) and once face-off (+m is zero).
    # (l, m) = (2, 1), (4, 4), or (4, 3):
    #   generate at inc = (pi/2). At that angle y_{l, m} = y_{l,-m} != 0. 
    #   we can then ioslate the +/-m by generating once at phi = 0 and once at
    #   phi = (1/m)*(pi/2), then taking combinations of the hplus and hcross.
    # (l, m) = (3, 3):
    #   generate at inc = (pi/2). At that angle, y_{3,3} = -y_{3,-3} != 0.
    #   we can then ioslate the +/-m by generating once at phi = 0 and once at
    #   phi = (1/m)*(pi/2), then taking combinations of the hplus and hcross.
    hlms = {}
    for mode in mode_array:
        ell, m = mode
        ma = [mode]
        if approx == 'IMRPhenomPv3HM' and mode != (2, 2):
            ma.append((2, 2))
        if m == 2:
            # generate once face on and once face off
            ulm, vlm = get_fd_waveform(
                mode_array=ma, **_shiftparams(0., 0., **params))
            # face off
            ulmm, vlmm = get_fd_waveform(
                mode_array=ma, **_shiftparams(numpy.pi, 0., **params))
            # subtract off the 2,2 mode for PhenomPv3HM
            if approx == 'IMRPhenomPv3HM' and mode != (2, 2):
                # +m
                hp22, hc22 = get_fd_waveform(
                    mode_array=[(2, 2)], **_shiftparams(0., 0., **params))
                ulm -= hp22
                vlm -= hc22
                # -m
                hp22, hc22 = get_fd_waveform(
                    mode_array=[(2, 2)], **_shiftparams(numpy.pi, 0.,
                                                        **params))
                ulmm -= hp22
                vlmm -= hc22
            # divide out the inclination part
            glm = get_glm(ell, m, 0.)
            ulm /= glm
            vlm /= glm
            # -m: although same magnitude, may have opposite sign
            glm = get_glm(ell, -m, numpy.pi)
            ulmm /= glm
            vlmm /= glm
        elif mode in [(2, 1), (3, 3), (4, 4), (4, 3)]: 
            inc = numpy.pi/2.
            hp1, hc1 = get_fd_waveform(mode_array=ma,
                                       **_shiftparams(inc, 0., **params))
            phi = (1./m)*(numpy.pi/2.)
            # rotate the spins to accomodate the change in phase
            #shifted_ps = params.copy()
            #shifted_ps['spin1x'] = spin1perp * numpy.cos(spin1az-phi)
            #shifted_ps['spin1y'] = spin1perp * numpy.sin(spin1az-phi)
            #shifted_ps['spin2x'] = spin2perp * numpy.cos(spin2az-phi)
            #shifted_ps['spin2y'] = spin2perp * numpy.sin(spin2az-phi)
            hp2, hc2 = get_fd_waveform(mode_array=ma,
                                       **_shiftparams(inc, phi, **params))
            if approx == 'IMRPhenomPv3HM':
                # generate the 2,2 and subtract it off
                hp22, hc22 = get_fd_waveform(mode_array=[(2, 2)],
                                             **_shiftparams(inc, 0., **params))
                hp1 -= hp22
                hc1 -= hc22
                hp22, hc22 = get_fd_waveform(mode_array=[(2, 2)],
                                             **_shiftparams(inc, phi, **params))
                hp2 -= hp22
                hc2 -= hc22
            # divide out the inclination part
            glm = get_glm(ell, m, inc)
            hp1 /= glm
            hc1 /= glm
            hp2 /= glm
            hc2 /= glm
            # if gl(-m) = -glm, we pick up a negative sign in the -m
            sgn = numpy.sign(get_glm(ell, -m, inc)/glm)
            # now convert to the ulm and vlm
            ulm = (hp1 - hc2)/2.
            vlm = -(hp2 + hc1)/2.
            ulmm = sgn*(hp1 + hc2)/2.
            vlmm = sgn*(hp2 - hc1)/2.
        else:
            raise ValueError("I don't know what to do with mode {}"
                             .format(mode))
        # XXX: tests showed that the vlms have the wrong sign, so adding a
        # negative here. Not sure what the cause is, but doing so seems to
        # work.
        hlms[ell, m] = (ulm, -vlm)
        hlms[ell, -m] = (ulmm, -vlmm)
    return hlms


_mode_waveform_td = {'NRSur7dq4': get_nrsur_modes,
                     'NRSur7dq2': get_nrsur_modes,
                     }


_mode_waveform_fd = {'IMRPhenomHM': get_imrphenomhm_modes,
                     'IMRPhenomPv3HM': get_imrphenomhm_modes,
                     'IMRPhenomXHM': get_imrphenomhm_modes,
                     'IMRPhenomXPHM' : get_imrphenomhm_modes,
                    }


def get_fd_waveform_modes(template=None, **kwargs):
    """Generates frequency domain waveforms, but does not sum over the modes.
    """
    params = props(template, **kwargs)
    required = parameters.fd_required
    check_args(params, required)
    apprx = params['approximant']
    if apprx not in _mode_waveform_fd:
        raise ValueError("I don't support approximant {}, sorry"
                         .format(apprx))
    return _mode_waveform_fd[apprx](**params)


def get_td_waveform_modes(template=None, **kwargs):
    """Generates time domain waveforms, but does not sum over the modes.
    """
    params = props(template, **kwargs)
    required = parameters.fd_required
    check_args(params, required)
    apprx = params['approximant']
    if apprx not in _mode_waveform_td:
        raise ValueError("I don't support approximant {}, sorry"
                         .format(apprx))
    return _mode_waveform_td[apprx](**params)
