# Copyright (C) 2022  Shichao Wu, Alex Nitz
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
"""
This module provides analytical PSDs and sensitivity curves for space borne
detectors, such as LISA. Based on LISA technical note <LISA-LCST-SGS-TN-001>,
LDC manual <LISA-LCST-SGS-MAN-001>, and paper <10.1088/1361-6382/ab1101>.
"""

import numpy as np
from astropy import constants


def psd_lisa_acc_noise(f, acc_noise_level):
    r""" The PSD of LISA's acceleration noise.
    Parameters
    ----------
    f : float or numpy.array
        The frequency or frequency range, in the unit of "Hz".
    acc_noise_level : float
        The level of acceleration noise.
    Returns
    -------
    s_acc_nu : float or numpy.array
        The PSD value or array for acceleration noise.
    Notes
    -----
        Pease see Eq.(11-13) in <LISA-LCST-SGS-TN-001> for more details.
    """
    s_acc = acc_noise_level**2 * (1+(4e-4/f)**2)*(1+(f/8e-3)**4)
    s_acc_d = s_acc * (2*np.pi*f)**(-4)
    s_acc_nu = s_acc_d * (2*np.pi*f/constants.c.value)**2

    return s_acc_nu

def psd_lisa_oms_noise(f, oms_noise_level):
    r""" The PSD of LISA's OMS noise.
    Parameters
    ----------
    f : float or numpy.array
        The frequency or frequency range, in the unit of "Hz".
    oms_noise_level : float
        The level of OMS noise.
    Returns
    -------
    s_oms_nu : float or numpy.array
        The PSD value or array for OMS noise.
    Notes
    -----
        Pease see Eq.(9-10) in <LISA-LCST-SGS-TN-001> for more details.
    """
    s_oms_d = oms_noise_level**2 * (1+(2e-3/f)**4)
    s_oms_nu = s_oms_d * (2*np.pi*f/constants.c.value)**2

    return s_oms_nu

def lisa_psd_components(f, acc_noise_level, oms_noise_level):
    r""" The PSD of LISA's acceleration and OMS noise.
    Parameters
    ----------
    f : float or numpy.array
        The frequency or frequency range, in the unit of "Hz".
    acc_noise_level : float
        The level of acceleration noise.
    oms_noise_level : float
        The level of OMS noise.
    Returns
    -------
    [low_freq_component, high_freq_component] : list
        The PSD value or array for acceleration and OMS noise.
    """
    low_freq_component = psd_lisa_acc_noise(f, acc_noise_level)
    high_freq_component = psd_lisa_oms_noise(f, oms_noise_level)

    return [low_freq_component, high_freq_component]

def analytical_psd_lisa_tdi_1p5(f, len_arm, acc_noise_level, oms_noise_level):
    r""" The TDI-1.5 analytical PSD for LISA.
    Parameters
    ----------
    f : float or numpy.array
        The frequency or frequency range, in the unit of "Hz".
    len_arm : float
        The arm length of LISA, in the unit of "m".
    acc_noise_level : float
        The level of acceleration noise.
    oms_noise_level : float
        The level of OMS noise.
    Returns
    -------
    psd : float or numpy.array
        The TDI-1.5 PSD for LISA.
    Notes
    -----
        Pease see Eq.(19) in <LISA-LCST-SGS-TN-001> for more details.
    """
    [s_acc_nu, s_oms_nu] = lisa_psd_components(
                              f, acc_noise_level, oms_noise_level)
    omega_len = 2*np.pi*f * len_arm/constants.c.value
    psd = 16*(np.sin(omega_len))**2 * (s_oms_nu+s_acc_nu*(3+np.cos(omega_len)))

    return psd

def analytical_psd_lisa_tdi_2p0(f, len_arm, acc_noise_level, oms_noise_level):
    r""" The TDI-2.0 analytical PSD for LISA.
    Parameters
    ----------
    f : float or numpy.array
        The frequency or frequency range, in the unit of "Hz".
    len_arm : float
        The arm length of LISA, in the unit of "m".
    acc_noise_level : float
        The level of acceleration noise.
    oms_noise_level : float
        The level of OMS noise.
    Returns
    -------
    psd : float or numpy.array
        The TDI-2.0 PSD for LISA.
    Notes
    -----
        Pease see Eq.(20) in <LISA-LCST-SGS-TN-001> for more details.
    """
    [s_acc_nu, s_oms_nu] = lisa_psd_components(
                              f, acc_noise_level, oms_noise_level)
    omega_len = 2*np.pi*f * len_arm/constants.c.value
    psd = 64*(np.sin(omega_len))**2 * (np.sin(2*omega_len))**2 * \
          (s_oms_nu+s_acc_nu*(3+np.cos(2*omega_len)))

    return psd

def averaged_lisa_fplus_sq_approx(f, len_arm):
    r""" An approximant for LISA's squared antenna response function,
    averaged over sky and polarization.
    Parameters
    ----------
    f : float or numpy.array
        The frequency or frequency range, in the unit of "Hz".
    len_arm : float
        The arm length of LISA, in the unit of "m".
    Returns
    -------
    fp_sq_approx : float or numpy.array
        The sky and polarization averaged squared antenna response.
    Notes
    -----
        Pease see Eq.(36) in <LISA-LCST-SGS-TN-001> for more details.
    """
    omega_len = 2*np.pi*f * len_arm/constants.c.value
    fp_sq_approx = 3/20/(1+0.6*omega_len**2)

    return fp_sq_approx

def averaged_response_lisa_tdi_1p5(f, len_arm):
    r""" LISA's TDI-1.5 response function to GW,
    averaged over sky and polarization.
    Parameters
    ----------
    f : float or numpy.array
        The frequency or frequency range, in the unit of "Hz".
    len_arm : float
        The arm length of LISA, in the unit of "m".
    Returns
    -------
    response_tdi_1p5 : float or numpy.array
        The sky and polarization averaged TDI-1.5 response function to GW.
    Notes
    -----
        Pease see Eq.(39) in <LISA-LCST-SGS-TN-001> for more details.
    """
    omega_len = 2*np.pi*f * len_arm/constants.c.value
    ave_fp2 = averaged_lisa_fplus_sq_approx(f, len_arm)
    response_tdi_1p5 = (4*omega_len)**2 * np.sin(omega_len)**2 * ave_fp2

    return response_tdi_1p5

def averaged_response_lisa_tdi_2p0(f, len_arm):
    r""" LISA's TDI-2.0 response function to GW,
    averaged over sky and polarization.
    Parameters
    ----------
    f : float or numpy.array
        The frequency or frequency range, in the unit of "Hz".
    len_arm : float
        The arm length of LISA, in the unit of "m".
    Returns
    -------
    response_tdi_2p0 : float or numpy.array
        The sky and polarization averaged TDI-2.0 response function to GW.
    Notes
    -----
        Pease see Eq.(40) in <LISA-LCST-SGS-TN-001> for more details.
    """
    omega_len = 2*np.pi*f * len_arm/constants.c.value
    response_tdi_1p5 = averaged_response_lisa_tdi_1p5(f, len_arm)
    response_tdi_2p0 = response_tdi_1p5 * (2*np.sin(2*omega_len))**2
    
    return response_tdi_2p0

def semi_sensitivity_curve_lisa(f, len_arm, acc_noise_level, oms_noise_level):
    r""" The semi-analytical LISA's sensitivity curve,
    averaged over sky and polarization.
    Parameters
    ----------
    f : float or numpy.array
        The frequency or frequency range, in the unit of "Hz".
    len_arm : float
        The arm length of LISA, in the unit of "m".
    acc_noise_level : float
        The level of acceleration noise.
    oms_noise_level : float
        The level of OMS noise.
    Returns
    -------
    sense_curve : float or numpy.array
        The sky and polarization averaged semi-analytical
        LISA's sensitivity curve.
    Notes
    -----
        Pease see Eq.(42) in <LISA-LCST-SGS-TN-001> for more details.
    """
    psd = analytical_psd_lisa_tdi_2p0(
               f, len_arm, acc_noise_level, oms_noise_level)
    response = averaged_response_lisa_tdi_2p0(f, len_arm)
    sense_curve = psd / response

    return sense_curve

def analytical_csd_lisa_tdi_1p5(f, len_arm, acc_noise_level, oms_noise_level):
    r""" The cross-spectrum density between LISA's TDI channel X and Y.
    Parameters
    ----------
    f : float or numpy.array
        The frequency or frequency range, in the unit of "Hz".
    len_arm : float
        The arm length of LISA, in the unit of "m".
    acc_noise_level : float
        The level of acceleration noise.
    oms_noise_level : float
        The level of OMS noise.
    Returns
    -------
    csd : float or numpy.array
        The CSD (value or array) between LISA's TDI channel X and Y.
    Notes
    -----
        Pease see Eq.(56) in <LISA-LCST-SGS-MAN-001> for more details.
    """
    omega_len = 2*np.pi*f * len_arm/constants.c.value
    [s_acc_nu, s_oms_nu] = lisa_psd_components(
                              f, acc_noise_level, oms_noise_level)
    csd = -8*np.sin(omega_len)**2 * np.cos(omega_len) * (s_oms_nu+4*s_acc_nu)

    return csd

def psd_lisa_tdi_A_E_1p5(f, len_arm, acc_noise_level, oms_noise_level):
    r""" The PSD of LISA's TDI-1.5 channel A and E.
    Parameters
    ----------
    f : float or numpy.array
        The frequency or frequency range, in the unit of "Hz".
    len_arm : float
        The arm length of LISA, in the unit of "m".
    acc_noise_level : float
        The level of acceleration noise.
    oms_noise_level : float
        The level of OMS noise.
    Returns
    -------
    psd_A_E : float or numpy.array
        The PSD (value or array) of LISA's TDI-1.5 channel A and E.
    Notes
    -----
        Pease see Eq.(58) in <LISA-LCST-SGS-MAN-001> for more details.
    """
    psd_X = analytical_psd_lisa_tdi_1p5(
               f, len_arm, acc_noise_level, oms_noise_level)
    csd_XY = analytical_csd_lisa_tdi_1p5(
               f, len_arm, acc_noise_level, oms_noise_level)
    psd_A_E = psd_A - csd_XY

    return psd_A_E

def psd_lisa_tdi_T_1p5(f, len_arm, acc_noise_level, oms_noise_level):
    r""" The PSD of LISA's TDI-1.5 channel T.
    Parameters
    ----------
    f : float or numpy.array
        The frequency or frequency range, in the unit of "Hz".
    len_arm : float
        The arm length of LISA, in the unit of "m".
    acc_noise_level : float
        The level of acceleration noise.
    oms_noise_level : float
        The level of OMS noise.
    Returns
    -------
    psd_T : float or numpy.array
        The PSD (value or array) of LISA's TDI-1.5 channel T.
    Notes
    -----
        Pease see Eq.(59) in <LISA-LCST-SGS-MAN-001> for more details.
    """
    psd_X = analytical_psd_lisa_tdi_1p5(
               f, len_arm, acc_noise_level, oms_noise_level)
    csd_XY = analytical_csd_lisa_tdi_1p5(
               f, len_arm, acc_noise_level, oms_noise_level)
    psd_T = psd_X + 2*csd_XY

    return psd_T

def sensitivity_curve_lisa_SciRD(f):
    r""" The analytical LISA's sensitivity curve in SciRD,
    averaged over sky and polarization.
    Parameters
    ----------
    f : float or numpy.array
        The frequency or frequency range, in the unit of "Hz".
    Returns
    -------
    sense_curve : float or numpy.array
        The sky and polarization averaged analytical
        LISA's sensitivity curve in SciRD.
    Notes
    -----
        Pease see Eq.(114) in <LISA-LCST-SGS-TN-001> for more details.
    """
    s_I = 5.76e-48 * (1+(4e-4/f)**2)
    s_II = 3.6e-41
    R = 1 + (f/2.5e-2)**2
    sense_curve = 10/3 * (s_I/(2*np.pi*f)**4+s_II) * R

    return sense_curve
