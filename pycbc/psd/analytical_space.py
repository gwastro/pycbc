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
This module provides (semi-)analytical PSDs and sensitivity curves for space
borne detectors, such as LISA, Taiji, and TianQin. Based on LISA technical note
<LISA-LCST-SGS-TN-001>, LDC manual <LISA-LCST-SGS-MAN-001>,
paper <10.1088/1361-6382/ab1101>, <10.1088/0264-9381/33/3/035010>,
<10.1103/PhysRevD.102.063021>, <10.1103/PhysRevD.100.043003>,
and <10.1103/PhysRevD.107.064021>.
"""

import numpy as np
from scipy.interpolate import interp1d
from astropy.constants import c
from pycbc.psd.read import from_numpy_arrays


def _psd_acc_noise(f, acc_noise_level=None):
    """ The PSD of TDI-based space-borne GW
    detectors' acceleration noise. Note that
    this is suitable for LISA and Taiji, TianQin
    has a different form.

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
        Please see Eq.(11-13) in <LISA-LCST-SGS-TN-001> for more details.
    """
    s_acc = acc_noise_level**2 * (1+(4e-4/f)**2)*(1+(f/8e-3)**4)
    s_acc_d = s_acc * (2*np.pi*f)**(-4)
    s_acc_nu = (2*np.pi*f/c.value)**2 * s_acc_d

    return s_acc_nu


def psd_lisa_acc_noise(f, acc_noise_level=3e-15):
    """ The PSD of LISA's acceleration noise.

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
        Please see Eq.(11-13) in <LISA-LCST-SGS-TN-001> for more details.
    """
    s_acc_nu = _psd_acc_noise(f, acc_noise_level)

    return s_acc_nu


def psd_tianqin_acc_noise(f, acc_noise_level=1e-15):
    """ The PSD of TianQin's acceleration noise.

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
        Please see Table(1) in <10.1088/0264-9381/33/3/035010>
        and that paper for more details.
    """
    s_acc_d = acc_noise_level**2 * (2*np.pi*f)**(-4) * (1+1e-4/f)
    s_acc_nu = (2*np.pi*f/c.value)**2 * s_acc_d

    return s_acc_nu


def psd_taiji_acc_noise(f, acc_noise_level=3e-15):
    """ The PSD of Taiji's acceleration noise.

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
        Please see Eq.(2) in <10.1103/PhysRevD.107.064021> for more details.
    """
    s_acc_nu = _psd_acc_noise(f, acc_noise_level)

    return s_acc_nu


def _psd_oms_noise(f, oms_noise_level=None):
    """ The PSD of TDI-based space-borne GW detectors' OMS noise.
    Note that this is suitable for LISA and Taiji, TianQin
    has a different form.

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
        Please see Eq.(9-10) in <LISA-LCST-SGS-TN-001> for more details.
    """
    s_oms_d = oms_noise_level**2 * (1+(2e-3/f)**4)
    s_oms_nu = s_oms_d * (2*np.pi*f/c.value)**2

    return s_oms_nu


def psd_lisa_oms_noise(f, oms_noise_level=15e-12):
    """ The PSD of LISA's OMS noise.

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
        Please see Eq.(9-10) in <LISA-LCST-SGS-TN-001> for more details.
    """
    s_oms_nu = _psd_oms_noise(f, oms_noise_level)

    return s_oms_nu


def psd_tianqin_oms_noise(f, oms_noise_level=1e-12):
    """ The PSD of TianQin's OMS noise.

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
        Please see Table(1) in <10.1088/0264-9381/33/3/035010>
        and that paper for more details.
    """
    s_oms_d = oms_noise_level**2
    s_oms_nu = s_oms_d * (2*np.pi*f/c.value)**2

    return s_oms_nu


def psd_taiji_oms_noise(f, oms_noise_level=8e-12):
    """ The PSD of Taiji's OMS noise.

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
        Please see Eq.(1) in <10.1103/PhysRevD.107.064021> for more details.
    """
    s_oms_nu = _psd_oms_noise(f, oms_noise_level)

    return s_oms_nu


def lisa_psd_components(f, acc_noise_level=3e-15, oms_noise_level=15e-12):
    """ The PSD of LISA's acceleration and OMS noise.

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
    low_freq_component, high_freq_component :
        The PSD value or array for acceleration and OMS noise.
    """
    acc_noise_level = np.float64(acc_noise_level)
    oms_noise_level = np.float64(oms_noise_level)
    low_freq_component = psd_lisa_acc_noise(f, acc_noise_level)
    high_freq_component = psd_lisa_oms_noise(f, oms_noise_level)

    return low_freq_component, high_freq_component


def tianqin_psd_components(f, acc_noise_level=1e-15, oms_noise_level=1e-12):
    """ The PSD of TianQin's acceleration and OMS noise.

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
    low_freq_component, high_freq_component :
        The PSD value or array for acceleration and OMS noise.
    """
    acc_noise_level = np.float64(acc_noise_level)
    oms_noise_level = np.float64(oms_noise_level)
    low_freq_component = psd_tianqin_acc_noise(f, acc_noise_level)
    high_freq_component = psd_tianqin_oms_noise(f, oms_noise_level)

    return low_freq_component, high_freq_component


def taiji_psd_components(f, acc_noise_level=3e-15, oms_noise_level=8e-12):
    """ The PSD of Taiji's acceleration and OMS noise.

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
    low_freq_component, high_freq_component :
        The PSD value or array for acceleration and OMS noise.
    """
    acc_noise_level = np.float64(acc_noise_level)
    oms_noise_level = np.float64(oms_noise_level)
    low_freq_component = psd_taiji_acc_noise(f, acc_noise_level)
    high_freq_component = psd_taiji_oms_noise(f, oms_noise_level)

    return low_freq_component, high_freq_component


def omega_length(f, len_arm=None):
    """ The function to calculate 2*pi*f*arm_length.

    Parameters
    ----------
    f : float or numpy.array
        The frequency or frequency range, in the unit of "Hz".
    len_arm : float
        The arm length of LISA, TianQin, or Taiji. The default
        value here is None.

    Returns
    -------
    omega_len : float or numpy.array
        The value of 2*pi*f*arm_length.
    """
    omega_len = 2*np.pi*f * len_arm/c.value

    return omega_len


def _analytical_psd_tdi_XYZ(length, delta_f, low_freq_cutoff,
                            len_arm=None, psd_components=None, tdi=None):
    """ The TDI-1.5/2.0 analytical PSD (X,Y,Z channel) for TDI-based
    space-borne GW detectors.

    Parameters
    ----------
    length : int
        Length of output Frequencyseries.
    delta_f : float
        Frequency step for output FrequencySeries.
    low_freq_cutoff : float
        Low-frequency cutoff for output FrequencySeries.
    len_arm : float
        The arm length of the detector, in the unit of "m".
    psd_components : numpy.array
        The PSDs of acceleration noise and OMS noise.
    tdi : str
        The version of TDI. Choose from "1.5" or "2.0".

    Returns
    -------
    fseries : FrequencySeries
        The TDI-1.5/2.0 PSD (X,Y,Z channel).
    Notes
    -----
        Please see Eq.(19-20) in <LISA-LCST-SGS-TN-001> for more details.
    """
    len_arm = np.float64(len_arm)
    fr = np.linspace(low_freq_cutoff, (length-1)*2*delta_f, length)
    s_acc_nu, s_oms_nu = psd_components
    omega_len = omega_length(fr, len_arm)
    psd = 16*(np.sin(omega_len))**2 * (s_oms_nu +
                                       s_acc_nu*(3+np.cos(2*omega_len)))
    if str(tdi) not in ["1.5", "2.0"]:
        raise ValueError("The version of TDI, currently only for 1.5 or 2.0.")
    if str(tdi) == "2.0":
        tdi2_factor = 4*(np.sin(2*omega_len))**2
        psd *= tdi2_factor
    fseries = from_numpy_arrays(fr, psd, length, delta_f, low_freq_cutoff)

    return fseries


def analytical_psd_lisa_tdi_XYZ(length, delta_f, low_freq_cutoff,
                                len_arm=2.5e9, acc_noise_level=3e-15,
                                oms_noise_level=15e-12, tdi=None):
    """ The TDI-1.5/2.0 analytical PSD (X,Y,Z channel) for LISA.

    Parameters
    ----------
    length : int
        Length of output Frequencyseries.
    delta_f : float
        Frequency step for output FrequencySeries.
    low_freq_cutoff : float
        Low-frequency cutoff for output FrequencySeries.
    len_arm : float
        The arm length of LISA, in the unit of "m".
    acc_noise_level : float
        The level of acceleration noise.
    oms_noise_level : float
        The level of OMS noise.
    tdi : str
        The version of TDI. Choose from "1.5" or "2.0".

    Returns
    -------
    fseries : FrequencySeries
        The TDI-1.5/2.0 PSD (X,Y,Z channel) for LISA.
    Notes
    -----
        Please see Eq.(19-20) in <LISA-LCST-SGS-TN-001> for more details.
    """
    fr = np.linspace(low_freq_cutoff, (length-1)*2*delta_f, length)
    psd_components = np.array(lisa_psd_components(
        fr, acc_noise_level, oms_noise_level))
    fseries = _analytical_psd_tdi_XYZ(length, delta_f, low_freq_cutoff,
                                      len_arm, psd_components, tdi)

    return fseries


def analytical_psd_tianqin_tdi_XYZ(length, delta_f, low_freq_cutoff,
                                   len_arm=np.sqrt(3)*1e8,
                                   acc_noise_level=1e-15,
                                   oms_noise_level=1e-12, tdi=None):
    """ The TDI-1.5/2.0 analytical PSD (X,Y,Z channel) for TianQin.

    Parameters
    ----------
    length : int
        Length of output Frequencyseries.
    delta_f : float
        Frequency step for output FrequencySeries.
    low_freq_cutoff : float
        Low-frequency cutoff for output FrequencySeries.
    len_arm : float
        The arm length of TianQin, in the unit of "m".
    acc_noise_level : float
        The level of acceleration noise.
    oms_noise_level : float
        The level of OMS noise.
    tdi : str
        The version of TDI. Choose from "1.5" or "2.0".

    Returns
    -------
    fseries : FrequencySeries
        The TDI-1.5/2.0 PSD (X,Y,Z channel) for TianQin.
    Notes
    -----
        Please see Table(1) in <10.1088/0264-9381/33/3/035010>
        for more details.
    """
    fr = np.linspace(low_freq_cutoff, (length-1)*2*delta_f, length)
    psd_components = np.array(tianqin_psd_components(
        fr, acc_noise_level, oms_noise_level))
    fseries = _analytical_psd_tdi_XYZ(length, delta_f, low_freq_cutoff,
                                      len_arm, psd_components, tdi)

    return fseries


def analytical_psd_taiji_tdi_XYZ(length, delta_f, low_freq_cutoff,
                                 len_arm=3e9, acc_noise_level=3e-15,
                                 oms_noise_level=8e-12, tdi=None):
    """ The TDI-1.5/2.0 analytical PSD (X,Y,Z channel) for Taiji.

    Parameters
    ----------
    length : int
        Length of output Frequencyseries.
    delta_f : float
        Frequency step for output FrequencySeries.
    low_freq_cutoff : float
        Low-frequency cutoff for output FrequencySeries.
    len_arm : float
        The arm length of Taiji, in the unit of "m".
    acc_noise_level : float
        The level of acceleration noise.
    oms_noise_level : float
        The level of OMS noise.
    tdi : str
        The version of TDI. Choose from "1.5" or "2.0".

    Returns
    -------
    fseries : FrequencySeries
        The TDI-1.5/2.0 PSD (X,Y,Z channel) for Taiji.
    Notes
    -----
        Please see <10.1103/PhysRevD.107.064021> for more details.
    """
    fr = np.linspace(low_freq_cutoff, (length-1)*2*delta_f, length)
    psd_components = np.array(taiji_psd_components(
        fr, acc_noise_level, oms_noise_level))
    fseries = _analytical_psd_tdi_XYZ(length, delta_f, low_freq_cutoff,
                                      len_arm, psd_components, tdi)

    return fseries


def _analytical_csd_tdi_XY(length, delta_f, low_freq_cutoff,
                           len_arm=None, psd_components=None, tdi=None):
    """ The cross-spectrum density between TDI channel X and Y.

    Parameters
    ----------
    length : int
        Length of output Frequencyseries.
    delta_f : float
        Frequency step for output FrequencySeries.
    low_freq_cutoff : float
        Low-frequency cutoff for output FrequencySeries.
    len_arm : float
        The arm length of the detector, in the unit of "m".
    psd_components : numpy.array
        The PSDs of acceleration noise and OMS noise.
    tdi : str
        The version of TDI. Choose from "1.5" or "2.0".

    Returns
    -------
    fseries : FrequencySeries
        The CSD between TDI-1.5/2.0 channel X and Y.
    Notes
    -----
        Please see Eq.(56) in <LISA-LCST-SGS-MAN-001(Radler)> for more details.
    """
    len_arm = np.float64(len_arm)
    fr = np.linspace(low_freq_cutoff, (length-1)*2*delta_f, length)
    s_acc_nu, s_oms_nu = psd_components
    omega_len = omega_length(fr, len_arm)
    csd = (-8*np.sin(omega_len)**2 * np.cos(omega_len) *
           (s_oms_nu+4*s_acc_nu))
    if str(tdi) not in ["1.5", "2.0"]:
        raise ValueError("The version of TDI, currently only for 1.5 or 2.0.")
    if str(tdi) == "2.0":
        tdi2_factor = 4*(np.sin(2*omega_len))**2
        csd *= tdi2_factor
    fseries = from_numpy_arrays(fr, csd, length, delta_f, low_freq_cutoff)

    return fseries


def analytical_csd_lisa_tdi_XY(length, delta_f, low_freq_cutoff,
                               len_arm=2.5e9, acc_noise_level=3e-15,
                               oms_noise_level=15e-12, tdi=None):
    """ The cross-spectrum density between LISA's TDI channel X and Y.

    Parameters
    ----------
    length : int
        Length of output Frequencyseries.
    delta_f : float
        Frequency step for output FrequencySeries.
    low_freq_cutoff : float
        Low-frequency cutoff for output FrequencySeries.
    len_arm : float
        The arm length of LISA, in the unit of "m".
    acc_noise_level : float
        The level of acceleration noise.
    oms_noise_level : float
        The level of OMS noise.
    tdi : str
        The version of TDI. Choose from "1.5" or "2.0".

    Returns
    -------
    fseries : FrequencySeries
        The CSD between LISA's TDI-1.5/2.0 channel X and Y.
    Notes
    -----
        Please see Eq.(56) in <LISA-LCST-SGS-MAN-001(Radler)> for more details.
    """
    fr = np.linspace(low_freq_cutoff, (length-1)*2*delta_f, length)
    psd_components = np.array(lisa_psd_components(
        fr, acc_noise_level, oms_noise_level))
    fseries = _analytical_csd_tdi_XY(length, delta_f, low_freq_cutoff,
                                     len_arm, psd_components, tdi)

    return fseries


def _analytical_psd_tdi_AE(length, delta_f, low_freq_cutoff,
                           len_arm=None, psd_components=None, tdi=None):
    """ The PSD of TDI-1.5/2.0 channel A and E.

    Parameters
    ----------
    length : int
        Length of output Frequencyseries.
    delta_f : float
        Frequency step for output FrequencySeries.
    low_freq_cutoff : float
        Low-frequency cutoff for output FrequencySeries.
    len_arm : float
        The arm length of the detector, in the unit of "m".
    psd_components : numpy.array
        The PSDs of acceleration noise and OMS noise.
    tdi : str
        The version of TDI. Choose from "1.5" or "2.0".

    Returns
    -------
    fseries : FrequencySeries
        The PSD of TDI-1.5/2.0 channel A and E.
    Notes
    -----
        Please see Eq.(58) in <LISA-LCST-SGS-MAN-001(Radler)> for more details.
    """
    len_arm = np.float64(len_arm)
    fr = np.linspace(low_freq_cutoff, (length-1)*2*delta_f, length)
    s_acc_nu, s_oms_nu = psd_components
    omega_len = omega_length(fr, len_arm)
    psd = (8*(np.sin(omega_len))**2 *
           (4*(1+np.cos(omega_len)+np.cos(omega_len)**2)*s_acc_nu +
           (2+np.cos(omega_len))*s_oms_nu))
    if str(tdi) not in ["1.5", "2.0"]:
        raise ValueError("The version of TDI, currently only for 1.5 or 2.0.")
    if str(tdi) == "2.0":
        tdi2_factor = 4*(np.sin(2*omega_len))**2
        psd *= tdi2_factor
    fseries = from_numpy_arrays(fr, psd, length, delta_f, low_freq_cutoff)

    return fseries


def analytical_psd_lisa_tdi_AE(length, delta_f, low_freq_cutoff,
                               len_arm=2.5e9, acc_noise_level=3e-15,
                               oms_noise_level=15e-12, tdi=None):
    """ The PSD of LISA's TDI-1.5/2.0 channel A and E.

    Parameters
    ----------
    length : int
        Length of output Frequencyseries.
    delta_f : float
        Frequency step for output FrequencySeries.
    low_freq_cutoff : float
        Low-frequency cutoff for output FrequencySeries.
    len_arm : float
        The arm length of LISA, in the unit of "m".
    acc_noise_level : float
        The level of acceleration noise.
    oms_noise_level : float
        The level of OMS noise.
    tdi : str
        The version of TDI. Choose from "1.5" or "2.0".

    Returns
    -------
    fseries : FrequencySeries
        The PSD of LISA's TDI-1.5/2.0 channel A and E.
    Notes
    -----
        Please see Eq.(58) in <LISA-LCST-SGS-MAN-001(Radler)> for more details.
    """
    fr = np.linspace(low_freq_cutoff, (length-1)*2*delta_f, length)
    psd_components = np.array(lisa_psd_components(
        fr, acc_noise_level, oms_noise_level))
    fseries = _analytical_psd_tdi_AE(length, delta_f, low_freq_cutoff,
                                     len_arm, psd_components, tdi)

    return fseries


def analytical_psd_tianqin_tdi_AE(length, delta_f, low_freq_cutoff,
                                  len_arm=np.sqrt(3)*1e8,
                                  acc_noise_level=1e-15,
                                  oms_noise_level=1e-12, tdi=None):
    """ The PSD of TianQin's TDI-1.5/2.0 channel A and E.

    Parameters
    ----------
    length : int
        Length of output Frequencyseries.
    delta_f : float
        Frequency step for output FrequencySeries.
    low_freq_cutoff : float
        Low-frequency cutoff for output FrequencySeries.
    len_arm : float
        The arm length of TianQin, in the unit of "m".
    acc_noise_level : float
        The level of acceleration noise.
    oms_noise_level : float
        The level of OMS noise.
    tdi : str
        The version of TDI. Choose from "1.5" or "2.0".

    Returns
    -------
    fseries : FrequencySeries
        The PSD of TianQin's TDI-1.5/2.0 channel A and E.
    Notes
    -----
        Please see Table(1) in <10.1088/0264-9381/33/3/035010>
        for more details.
    """
    fr = np.linspace(low_freq_cutoff, (length-1)*2*delta_f, length)
    psd_components = np.array(tianqin_psd_components(
        fr, acc_noise_level, oms_noise_level))
    fseries = _analytical_psd_tdi_AE(length, delta_f, low_freq_cutoff,
                                     len_arm, psd_components, tdi)

    return fseries


def analytical_psd_taiji_tdi_AE(length, delta_f, low_freq_cutoff,
                                len_arm=3e9, acc_noise_level=3e-15,
                                oms_noise_level=8e-12, tdi=None):
    """ The PSD of Taiji's TDI-1.5/2.0 channel A and E.

    Parameters
    ----------
    length : int
        Length of output Frequencyseries.
    delta_f : float
        Frequency step for output FrequencySeries.
    low_freq_cutoff : float
        Low-frequency cutoff for output FrequencySeries.
    len_arm : float
        The arm length of Taiji, in the unit of "m".
    acc_noise_level : float
        The level of acceleration noise.
    oms_noise_level : float
        The level of OMS noise.
    tdi : str
        The version of TDI. Choose from "1.5" or "2.0".

    Returns
    -------
    fseries : FrequencySeries
        The PSD of Taiji's TDI-1.5/2.0 channel A and E.
    Notes
    -----
        Please see <10.1103/PhysRevD.107.064021> for more details.
    """
    fr = np.linspace(low_freq_cutoff, (length-1)*2*delta_f, length)
    psd_components = np.array(taiji_psd_components(
        fr, acc_noise_level, oms_noise_level))
    fseries = _analytical_psd_tdi_AE(length, delta_f, low_freq_cutoff,
                                     len_arm, psd_components, tdi)

    return fseries


def _analytical_psd_tdi_T(length, delta_f, low_freq_cutoff,
                          len_arm=None, psd_components=None, tdi=None):
    """ The PSD of TDI-1.5/2.0 channel T.

    Parameters
    ----------
    length : int
        Length of output Frequencyseries.
    delta_f : float
        Frequency step for output FrequencySeries.
    low_freq_cutoff : float
        Low-frequency cutoff for output FrequencySeries.
    len_arm : float
        The arm length of the detector, in the unit of "m".
    psd_components : numpy.array
        The PSDs of acceleration noise and OMS noise.
    tdi : str
        The version of TDI. Choose from "1.5" or "2.0".

    Returns
    -------
    fseries : FrequencySeries
        The PSD of TDI-1.5/2.0 channel T.
    Notes
    -----
        Please see Eq.(59) in <LISA-LCST-SGS-MAN-001(Radler)> for more details.
    """
    len_arm = np.float64(len_arm)
    fr = np.linspace(low_freq_cutoff, (length-1)*2*delta_f, length)
    s_acc_nu, s_oms_nu = psd_components
    omega_len = omega_length(fr, len_arm)
    psd = (32*np.sin(omega_len)**2 * np.sin(omega_len/2)**2 *
           (4*s_acc_nu*np.sin(omega_len/2)**2 + s_oms_nu))
    if str(tdi) not in ["1.5", "2.0"]:
        raise ValueError("The version of TDI, currently only for 1.5 or 2.0.")
    if str(tdi) == "2.0":
        tdi2_factor = 4*(np.sin(2*omega_len))**2
        psd *= tdi2_factor
    fseries = from_numpy_arrays(fr, psd, length, delta_f, low_freq_cutoff)

    return fseries


def analytical_psd_lisa_tdi_T(length, delta_f, low_freq_cutoff,
                              len_arm=2.5e9, acc_noise_level=3e-15,
                              oms_noise_level=15e-12, tdi=None):
    """ The PSD of LISA's TDI-1.5/2.0 channel T.

    Parameters
    ----------
    length : int
        Length of output Frequencyseries.
    delta_f : float
        Frequency step for output FrequencySeries.
    low_freq_cutoff : float
        Low-frequency cutoff for output FrequencySeries.
    len_arm : float
        The arm length of LISA, in the unit of "m".
    acc_noise_level : float
        The level of acceleration noise.
    oms_noise_level : float
        The level of OMS noise.
    tdi : str
        The version of TDI. Choose from "1.5" or "2.0".

    Returns
    -------
    fseries : FrequencySeries
        The PSD of LISA's TDI-1.5/2.0 channel T.
    Notes
    -----
        Please see Eq.(59) in <LISA-LCST-SGS-MAN-001(Radler)> for more details.
    """
    fr = np.linspace(low_freq_cutoff, (length-1)*2*delta_f, length)
    psd_components = np.array(lisa_psd_components(
        fr, acc_noise_level, oms_noise_level))
    fseries = _analytical_psd_tdi_T(length, delta_f, low_freq_cutoff,
                                    len_arm, psd_components, tdi)

    return fseries


def analytical_psd_tianqin_tdi_T(length, delta_f, low_freq_cutoff,
                                 len_arm=np.sqrt(3)*1e8,
                                 acc_noise_level=1e-15,
                                 oms_noise_level=1e-12, tdi=None):
    """ The PSD of TianQin's TDI-1.5/2.0 channel T.

    Parameters
    ----------
    length : int
        Length of output Frequencyseries.
    delta_f : float
        Frequency step for output FrequencySeries.
    low_freq_cutoff : float
        Low-frequency cutoff for output FrequencySeries.
    len_arm : float
        The arm length of TianQin, in the unit of "m".
    acc_noise_level : float
        The level of acceleration noise.
    oms_noise_level : float
        The level of OMS noise.
    tdi : str
        The version of TDI. Choose from "1.5" or "2.0".

    Returns
    -------
    fseries : FrequencySeries
        The PSD of TianQin's TDI-1.5/2.0 channel T.
    Notes
    -----
        Please see Table(1) in <10.1088/0264-9381/33/3/035010>
        for more details.
    """
    fr = np.linspace(low_freq_cutoff, (length-1)*2*delta_f, length)
    psd_components = np.array(tianqin_psd_components(
        fr, acc_noise_level, oms_noise_level))
    fseries = _analytical_psd_tdi_T(length, delta_f, low_freq_cutoff,
                                    len_arm, psd_components, tdi)

    return fseries


def analytical_psd_taiji_tdi_T(length, delta_f, low_freq_cutoff,
                               len_arm=3e9, acc_noise_level=3e-15,
                               oms_noise_level=8e-12, tdi=None):
    """ The PSD of Taiji's TDI-1.5/2.0 channel T.

    Parameters
    ----------
    length : int
        Length of output Frequencyseries.
    delta_f : float
        Frequency step for output FrequencySeries.
    low_freq_cutoff : float
        Low-frequency cutoff for output FrequencySeries.
    len_arm : float
        The arm length of Taiji, in the unit of "m".
    acc_noise_level : float
        The level of acceleration noise.
    oms_noise_level : float
        The level of OMS noise.
    tdi : str
        The version of TDI. Choose from "1.5" or "2.0".

    Returns
    -------
    fseries : FrequencySeries
        The PSD of Taiji's TDI-1.5/2.0 channel T.
    Notes
    -----
        Please see <10.1103/PhysRevD.107.064021> for more details.
    """
    fr = np.linspace(low_freq_cutoff, (length-1)*2*delta_f, length)
    psd_components = np.array(taiji_psd_components(
        fr, acc_noise_level, oms_noise_level))
    fseries = _analytical_psd_tdi_T(length, delta_f, low_freq_cutoff,
                                    len_arm, psd_components, tdi)

    return fseries


def averaged_lisa_fplus_sq_numerical(f, len_arm=2.5e9):
    """ A numerical fit for LISA's squared antenna response function,
    averaged over sky and polarization angle.

    Parameters
    ----------
    f : float or numpy.array
        The frequency or frequency range, in the unit of "Hz".
    len_arm : float
        The arm length of LISA, in the unit of "m".

    Returns
    -------
    fp_sq_numerical : float or numpy.array
        The sky and polarization angle averaged squared antenna response.
    Notes
    -----
        Please see Eq.(36) in <LISA-LCST-SGS-TN-001> for more details.
    """
    from astropy.utils.data import download_file

    if len_arm != 2.5e9:
        raise ValueError("Currently only support 'len_arm=2.5e9'.")
    # Download the numerical LISA averaged response.
    url = "https://zenodo.org/record/7497853/files/AvFXp2_Raw.npy"
    file_path = download_file(url, cache=True)
    freqs, fp_sq = np.load(file_path)
    # Padding the end.
    freqs = np.append(freqs, 2)
    fp_sq = np.append(fp_sq, 0.0012712348970728724)
    fp_sq_interp = interp1d(freqs, fp_sq, kind='linear',
                            fill_value="extrapolate")
    fp_sq_numerical = fp_sq_interp(f)/16

    return fp_sq_numerical


def averaged_fplus_sq_approximated(f, len_arm=None):
    r""" A simplified fit for TDI-based space-borne GW detectors'
    squared antenna response function, averaged over sky and
    polarization angle.

    .. math::
        <\left(F_{X}^{+}\right)^{2}>\approx \frac{3}{20} \frac{1}{
            1+0.6(\omega L)^{2}}

    Parameters
    ----------
    f : float or numpy.array
        The frequency or frequency range, in the unit of "Hz".
    len_arm : float
        The arm length of the detector, in the unit of "m".

    Returns
    -------
    fp_sq_approx : float or numpy.array
        The sky and polarization angle averaged squared antenna response.
    Notes
    -----
        Please see Eq.(9) in <10.1088/1361-6382/ab1101> for more details.
    """
    fp_sq_approx = (3./20.)*(1./(1.+0.6*omega_length(f, len_arm)**2))

    return fp_sq_approx


def averaged_tianqin_fplus_sq_numerical(f, len_arm=np.sqrt(3)*1e8):
    """ A numerical fit for TianQin's squared antenna response function,
    averaged over sky and polarization angle.

    Parameters
    ----------
    f : float or numpy.array
        The frequency or frequency range, in the unit of "Hz".
    len_arm : float
        The arm length of TianQin, in the unit of "m".

    Returns
    -------
    fp_sq_numerical : float or numpy.array
        The sky and polarization angle averaged squared antenna response.
    Notes
    -----
        Please see Eq.(15-16) in <10.1103/PhysRevD.100.043003>
        for more details.
    """
    base = averaged_fplus_sq_approximated(f, len_arm)
    a = [1, 1e-4, 2639e-4, 231/5*1e-4, -2093/1.25*1e-4, 2173e-5,
         2101e-6, 3027/2*1e-5, -42373/5*1e-6, 176087e-8,
         -8023/5*1e-7, 5169e-9]
    omega_len = omega_length(f, len_arm)
    omega_len_low_f = omega_len[omega_len < 4.1]
    omega_len_high_f = omega_len[omega_len >= 4.1]
    base_low_f = base[:len(omega_len_low_f)]
    base_high_f = base[len(omega_len_low_f):]
    low_f_modulation = np.polyval(a[::-1], omega_len_low_f)
    high_f_modulation = np.exp(
        -0.322 * np.sin(2*omega_len_high_f-4.712) + 0.078
    )
    low_f_result = np.multiply(base_low_f, low_f_modulation)
    high_f_result = np.multiply(base_high_f, high_f_modulation)
    fp_sq_numerical = np.concatenate((low_f_result, high_f_result))

    return fp_sq_numerical


def averaged_response_lisa_tdi(f, len_arm=2.5e9, tdi=None):
    """ LISA's TDI-1.5/2.0 response function to GW,
    averaged over sky and polarization angle.

    Parameters
    ----------
    f : float or numpy.array
        The frequency or frequency range, in the unit of "Hz".
    len_arm : float
        The arm length of LISA, in the unit of "m".
    tdi : str
        The version of TDI. Choose from "1.5" or "2.0".

    Returns
    -------
    response_tdi : float or numpy.array
        The sky and polarization angle averaged TDI-1.5/2.0 response to GW.
    Notes
    -----
        Please see Eq.(39-40) in <LISA-LCST-SGS-TN-001> for more details.
    """
    omega_len = omega_length(f, len_arm)
    ave_fp2 = averaged_lisa_fplus_sq_numerical(f, len_arm)
    response_tdi = (4*omega_len)**2 * np.sin(omega_len)**2 * ave_fp2
    if str(tdi) not in ["1.5", "2.0"]:
        raise ValueError("The version of TDI, currently only for 1.5 or 2.0.")
    if str(tdi) == "2.0":
        tdi2_factor = 4*(np.sin(2*omega_len))**2
        response_tdi *= tdi2_factor

    return response_tdi


def averaged_response_tianqin_tdi(f, len_arm=np.sqrt(3)*1e8, tdi=None):
    """ TianQin's TDI-1.5/2.0 response function to GW,
    averaged over sky and polarization angle.

    Parameters
    ----------
    f : float or numpy.array
        The frequency or frequency range, in the unit of "Hz".
    len_arm : float
        The arm length of TianQin, in the unit of "m".
    tdi : str
        The version of TDI. Choose from "1.5" or "2.0".

    Returns
    -------
    response_tdi : float or numpy.array
        The sky and polarization angle averaged TDI-1.5/2.0 response to GW.
    """
    omega_len = omega_length(f, len_arm)
    ave_fp2 = averaged_tianqin_fplus_sq_numerical(f, len_arm)
    response_tdi = (4*omega_len)**2 * np.sin(omega_len)**2 * ave_fp2
    if str(tdi) not in ["1.5", "2.0"]:
        raise ValueError("The version of TDI, currently only for 1.5 or 2.0.")
    if str(tdi) == "2.0":
        tdi2_factor = 4*(np.sin(2*omega_len))**2
        response_tdi *= tdi2_factor

    return response_tdi


def averaged_response_taiji_tdi(f, len_arm=3e9, tdi=None):
    """ Taiji's TDI-1.5/2.0 response function to GW,
    averaged over sky and polarization angle.

    Parameters
    ----------
    f : float or numpy.array
        The frequency or frequency range, in the unit of "Hz".
    len_arm : float
        The arm length of Taiji, in the unit of "m".
    tdi : str
        The version of TDI. Choose from "1.5" or "2.0".

    Returns
    -------
    response_tdi : float or numpy.array
        The sky and polarization angle averaged TDI-1.5/2.0 response to GW.
    """
    omega_len = omega_length(f, len_arm)
    ave_fp2 = averaged_fplus_sq_approximated(f, len_arm)
    response_tdi = (4*omega_len)**2 * np.sin(omega_len)**2 * ave_fp2
    if str(tdi) not in ["1.5", "2.0"]:
        raise ValueError("The version of TDI, currently only for 1.5 or 2.0.")
    if str(tdi) == "2.0":
        tdi2_factor = 4*(np.sin(2*omega_len))**2
        response_tdi *= tdi2_factor

    return response_tdi


def sensitivity_curve_lisa_semi_analytical(length, delta_f, low_freq_cutoff,
                                           len_arm=2.5e9,
                                           acc_noise_level=3e-15,
                                           oms_noise_level=15e-12):
    """ The semi-analytical LISA's sensitivity curve (6-links),
    averaged over sky and polarization angle.

    Parameters
    ----------
    length : int
        Length of output Frequencyseries.
    delta_f : float
        Frequency step for output FrequencySeries.
    low_freq_cutoff : float
        Low-frequency cutoff for output FrequencySeries.
    len_arm : float
        The arm length of LISA, in the unit of "m".
    acc_noise_level : float
        The level of acceleration noise.
    oms_noise_level : float
        The level of OMS noise.

    Returns
    -------
    fseries : FrequencySeries
        The sky and polarization angle averaged semi-analytical
        LISA's sensitivity curve (6-links).
    Notes
    -----
        Please see Eq.(42-43) in <LISA-LCST-SGS-TN-001> for more details.
    """
    len_arm = np.float64(len_arm)
    acc_noise_level = np.float64(acc_noise_level)
    oms_noise_level = np.float64(oms_noise_level)
    fr = np.linspace(low_freq_cutoff, (length-1)*2*delta_f, length)
    fp_sq = averaged_lisa_fplus_sq_numerical(fr, len_arm)
    s_acc_nu, s_oms_nu = lisa_psd_components(
                            fr, acc_noise_level, oms_noise_level)
    omega_len = omega_length(fr, len_arm)
    sense_curve = ((s_oms_nu + s_acc_nu*(3+np.cos(2*omega_len))) /
                   (omega_len**2*fp_sq))
    fseries = from_numpy_arrays(fr, sense_curve/2,
                                length, delta_f, low_freq_cutoff)

    return fseries


def sensitivity_curve_tianqin_analytical(length, delta_f, low_freq_cutoff,
                                         len_arm=np.sqrt(3)*1e8,
                                         acc_noise_level=1e-15,
                                         oms_noise_level=1e-12):
    """ The analytical TianQin's sensitivity curve (6-links),
    averaged over sky and polarization angle.

    Parameters
    ----------
    length : int
        Length of output Frequencyseries.
    delta_f : float
        Frequency step for output FrequencySeries.
    low_freq_cutoff : float
        Low-frequency cutoff for output FrequencySeries.
    len_arm : float
        The arm length of TianQin, in the unit of "m".
    acc_noise_level : float
        The level of acceleration noise.
    oms_noise_level : float
        The level of OMS noise.

    Returns
    -------
    fseries : FrequencySeries
        The sky and polarization angle averaged analytical
        TianQin's sensitivity curve (6-links).
    """
    len_arm = np.float64(len_arm)
    acc_noise_level = np.float64(acc_noise_level)
    oms_noise_level = np.float64(oms_noise_level)
    fr = np.linspace(low_freq_cutoff, (length-1)*2*delta_f, length)
    fp_sq = averaged_tianqin_fplus_sq_numerical(fr, len_arm)
    s_acc_nu, s_oms_nu = tianqin_psd_components(
                            fr, acc_noise_level, oms_noise_level)
    omega_len = omega_length(fr, len_arm)
    sense_curve = ((s_oms_nu + s_acc_nu*(3+np.cos(2*omega_len))) /
                   (omega_len**2*fp_sq))
    fseries = from_numpy_arrays(fr, sense_curve/2,
                                length, delta_f, low_freq_cutoff)

    return fseries


def sensitivity_curve_taiji_analytical(length, delta_f, low_freq_cutoff,
                                       len_arm=3e9, acc_noise_level=3e-15,
                                       oms_noise_level=8e-12):
    """ The analytical Taiji's sensitivity curve (6-links),
    averaged over sky and polarization angle.

    Parameters
    ----------
    length : int
        Length of output Frequencyseries.
    delta_f : float
        Frequency step for output FrequencySeries.
    low_freq_cutoff : float
        Low-frequency cutoff for output FrequencySeries.
    len_arm : float
        The arm length of Taiji, in the unit of "m".
    acc_noise_level : float
        The level of acceleration noise.
    oms_noise_level : float
        The level of OMS noise.

    Returns
    -------
    fseries : FrequencySeries
        The sky and polarization angle averaged analytical
        Taiji's sensitivity curve (6-links).
    """
    len_arm = np.float64(len_arm)
    acc_noise_level = np.float64(acc_noise_level)
    oms_noise_level = np.float64(oms_noise_level)
    fr = np.linspace(low_freq_cutoff, (length-1)*2*delta_f, length)
    fp_sq = averaged_fplus_sq_approximated(fr, len_arm)
    s_acc_nu, s_oms_nu = taiji_psd_components(
                            fr, acc_noise_level, oms_noise_level)
    omega_len = omega_length(fr, len_arm)
    sense_curve = ((s_oms_nu + s_acc_nu*(3+np.cos(2*omega_len))) /
                   (omega_len**2*fp_sq))
    fseries = from_numpy_arrays(fr, sense_curve/2,
                                length, delta_f, low_freq_cutoff)

    return fseries


def sensitivity_curve_lisa_SciRD(length, delta_f, low_freq_cutoff):
    """ The analytical LISA's sensitivity curve in SciRD,
    averaged over sky and polarization angle.

    Parameters
    ----------
    length : int
        Length of output Frequencyseries.
    delta_f : float
        Frequency step for output FrequencySeries.
    low_freq_cutoff : float
        Low-frequency cutoff for output FrequencySeries.

    Returns
    -------
    fseries : FrequencySeries
        The sky and polarization angle averaged analytical
        LISA's sensitivity curve in SciRD.
    Notes
    -----
        Please see Eq.(114) in <LISA-LCST-SGS-TN-001> for more details.
    """
    fr = np.linspace(low_freq_cutoff, (length-1)*2*delta_f, length)
    s_I = 5.76e-48 * (1+(4e-4/fr)**2)
    s_II = 3.6e-41
    R = 1 + (fr/2.5e-2)**2
    sense_curve = 10/3 * (s_I/(2*np.pi*fr)**4+s_II) * R
    fseries = from_numpy_arrays(fr, sense_curve, length,
                                delta_f, low_freq_cutoff)

    return fseries


def confusion_fit_lisa(length, delta_f, low_freq_cutoff, duration=1.0):
    """ The LISA's sensitivity curve for Galactic confusion noise,
    averaged over sky and polarization angle. No instrumental noise.

    Parameters
    ----------
    length : int
        Length of output Frequencyseries.
    delta_f : float
        Frequency step for output FrequencySeries.
    low_freq_cutoff : float
        Low-frequency cutoff for output FrequencySeries.
    duration : float
        The duration of observation, between 0 and 10, in the unit of years.

    Returns
    -------
    fseries : FrequencySeries
        The sky and polarization angle averaged
        LISA's sensitivity curve for Galactic confusion noise.
        No instrumental noise.
    Notes
    -----
        Please see Eq.(85-86) in <LISA-LCST-SGS-TN-001> for more details.
    """
    fr = np.linspace(low_freq_cutoff, (length-1)*2*delta_f, length)
    f1 = 10**(-0.25*np.log10(duration)-2.7)
    fk = 10**(-0.27*np.log10(duration)-2.47)
    sh_confusion = (0.5*1.14e-44*fr**(-7/3)*np.exp(-(fr/f1)**1.8) *
                    (1.0+np.tanh((fk-fr)/0.31e-3)))
    fseries = from_numpy_arrays(fr, sh_confusion, length, delta_f,
                                low_freq_cutoff)

    return fseries


def confusion_fit_tianqin(length, delta_f, low_freq_cutoff, duration=1.0):
    """ The TianQin's sensitivity curve for Galactic confusion noise,
    averaged over sky and polarization angle. No instrumental noise.
    Only valid for 0.5 mHz < f < 10 mHz. Note that the results between
    0.5, 1, 2, 4, and 5 years are extrapolated, might be non-physical.

    Parameters
    ----------
    length : int
        Length of output Frequencyseries.
    delta_f : float
        Frequency step for output FrequencySeries.
    low_freq_cutoff : float
        Low-frequency cutoff for output FrequencySeries.
    duration : float
        The duration of observation, between 0 and 5, in the unit of years.

    Returns
    -------
    fseries : FrequencySeries
        The sky and polarization angle averaged
        TianQin's sensitivity curve for Galactic confusion noise.
        No instrumental noise.
    Notes
    -----
        Please see Table(II) in <10.1103/PhysRevD.102.063021>
        for more details.
    """
    fr = np.linspace(low_freq_cutoff, (length-1)*2*delta_f, length)
    t_obs = [0.5, 1, 2, 4, 5]
    a0 = [-18.6, -18.6, -18.6, -18.6, -18.6]
    a1 = [-1.22, -1.13, -1.45, -1.43, -1.51]
    a2 = [0.009, -0.945, 0.315, -0.687, -0.710]
    a3 = [-1.87, -1.02, -1.19, 0.24, -1.13]
    a4 = [0.65, 4.05, -4.48, -0.15, -0.83]
    a5 = [3.6, -4.5, 10.8, -1.8, 13.2]
    a6 = [-4.6, -0.5, -9.4, -3.2, -19.1]
    fit_a0 = interp1d(t_obs, a0, kind='cubic', fill_value="extrapolate")
    fit_a1 = interp1d(t_obs, a1, kind='cubic', fill_value="extrapolate")
    fit_a2 = interp1d(t_obs, a2, kind='cubic', fill_value="extrapolate")
    fit_a3 = interp1d(t_obs, a3, kind='cubic', fill_value="extrapolate")
    fit_a4 = interp1d(t_obs, a4, kind='cubic', fill_value="extrapolate")
    fit_a5 = interp1d(t_obs, a5, kind='cubic', fill_value="extrapolate")
    fit_a6 = interp1d(t_obs, a6, kind='cubic', fill_value="extrapolate")
    if duration not in t_obs:
        raise Warning("Note that the results between " +
                      "0.5, 1, 2, 4, and 5 years are extrapolated, " +
                      "might be non-physical.")
    # 10/3 is the factor for sky-average, the original fit in the paper
    # is not sky-averaged.
    sh_confusion = 10./3 * np.power(
        10,
        fit_a0(duration) +
        fit_a1(duration) * np.log10(fr*1e3) +
        fit_a2(duration) * np.log10(fr*1e3)**2 +
        fit_a3(duration) * np.log10(fr*1e3)**3 +
        fit_a4(duration) * np.log10(fr*1e3)**4 +
        fit_a5(duration) * np.log10(fr*1e3)**5 +
        fit_a6(duration) * np.log10(fr*1e3)**6
    )**2
    # avoid the jump of values
    sh_confusion[(fr > 3e-4) & (fr < 5e-4)] = \
        sh_confusion[(np.abs(fr - 5e-4)).argmin()]
    sh_confusion[(fr < 3e-4) | (fr > 1e-2)] = 0
    fseries = from_numpy_arrays(fr, sh_confusion, length, delta_f,
                                low_freq_cutoff)

    return fseries


def confusion_fit_taiji(length, delta_f, low_freq_cutoff, duration=1.0):
    """ The Taiji's sensitivity curve for Galactic confusion noise,
    averaged over sky and polarization angle. No instrumental noise.
    Only valid for 0.1 mHz < f < 10 mHz. Note that the results between
    0.5, 1, 2, and 4 years are extrapolated, might be non-physical.

    Parameters
    ----------
    length : int
        Length of output Frequencyseries.
    delta_f : float
        Frequency step for output FrequencySeries.
    low_freq_cutoff : float
        Low-frequency cutoff for output FrequencySeries.
    duration : float
        The duration of observation, between 0 and 4, in the unit of years.

    Returns
    -------
    fseries : FrequencySeries
        The sky and polarization angle averaged
        Taiji's sensitivity curve for Galactic confusion noise.
        No instrumental noise.
    Notes
    -----
        Please see Eq.(6) and Table(I) in <10.1103/PhysRevD.107.064021>
        for more details.
    """
    fr = np.linspace(low_freq_cutoff, (length-1)*2*delta_f, length)
    t_obs = [0.5, 1, 2, 4]
    a0 = [-85.3498, -85.4336, -85.3919, -85.5448]
    a1 = [-2.64899, -2.46276, -2.69735, -3.23671]
    a2 = [-0.0699707, -0.183175, -0.749294, -1.64187]
    a3 = [-0.478447, -0.884147, -1.15302, -1.14711]
    a4 = [-0.334821, -0.427176, -0.302761, 0.0325887]
    a5 = [0.0658353, 0.128666, 0.175521, 0.187854]
    fit_a0 = interp1d(t_obs, a0, kind='cubic', fill_value="extrapolate")
    fit_a1 = interp1d(t_obs, a1, kind='cubic', fill_value="extrapolate")
    fit_a2 = interp1d(t_obs, a2, kind='cubic', fill_value="extrapolate")
    fit_a3 = interp1d(t_obs, a3, kind='cubic', fill_value="extrapolate")
    fit_a4 = interp1d(t_obs, a4, kind='cubic', fill_value="extrapolate")
    fit_a5 = interp1d(t_obs, a5, kind='cubic', fill_value="extrapolate")
    if duration not in t_obs:
        raise Warning("Note that the results between " +
                      "0.5, 1, 2, and 4 years are extrapolated, " +
                      "might be non-physical.")
    sh_confusion = np.exp(
        fit_a0(duration) +
        fit_a1(duration) * np.log(fr*1e3) +
        fit_a2(duration) * np.log(fr*1e3)**2 +
        fit_a3(duration) * np.log(fr*1e3)**3 +
        fit_a4(duration) * np.log(fr*1e3)**4 +
        fit_a5(duration) * np.log(fr*1e3)**5
    )
    sh_confusion[(fr < 1e-4) | (fr > 1e-2)] = 0
    fseries = from_numpy_arrays(fr, sh_confusion, length, delta_f,
                                low_freq_cutoff)

    return fseries


def sensitivity_curve_lisa_confusion(length, delta_f, low_freq_cutoff,
                                     len_arm=2.5e9, acc_noise_level=3e-15,
                                     oms_noise_level=15e-12,
                                     base_model="semi", duration=1.0):
    """ The LISA's sensitivity curve with Galactic confusion noise,
    averaged over sky and polarization angle.

    Parameters
    ----------
    length : int
        Length of output Frequencyseries.
    delta_f : float
        Frequency step for output FrequencySeries.
    low_freq_cutoff : float
        Low-frequency cutoff for output FrequencySeries.
    len_arm : float
        The arm length of LISA, in the unit of "m".
    acc_noise_level : float
        The level of acceleration noise.
    oms_noise_level : float
        The level of OMS noise.
    base_model : string
        The base model of sensitivity curve, chosen from "semi" or "SciRD".
    duration : float
        The duration of observation, between 0 and 10, in the unit of years.

    Returns
    -------
    fseries : FrequencySeries
        The sky and polarization angle averaged
        LISA's sensitivity curve with Galactic confusion noise.
    Notes
    -----
        Please see Eq.(85-86) in <LISA-LCST-SGS-TN-001> for more details.
    """
    if base_model == "semi":
        base_curve = sensitivity_curve_lisa_semi_analytical(
            length, delta_f, low_freq_cutoff,
            len_arm, acc_noise_level, oms_noise_level)
    elif base_model == "SciRD":
        base_curve = sensitivity_curve_lisa_SciRD(
            length, delta_f, low_freq_cutoff)
    else:
        raise ValueError("Must choose from 'semi' or 'SciRD'.")
    if duration < 0 or duration > 10:
        raise ValueError("Must between 0 and 10.")
    fseries_confusion = confusion_fit_lisa(
        length, delta_f, low_freq_cutoff, duration)
    fseries = from_numpy_arrays(base_curve.sample_frequencies,
                                base_curve+fseries_confusion,
                                length, delta_f, low_freq_cutoff)

    return fseries


def sensitivity_curve_tianqin_confusion(length, delta_f, low_freq_cutoff,
                                        len_arm=np.sqrt(3)*1e8,
                                        acc_noise_level=1e-15,
                                        oms_noise_level=1e-12, duration=1.0):
    """ The TianQin's sensitivity curve with Galactic confusion noise,
    averaged over sky and polarization angle.

    Parameters
    ----------
    length : int
        Length of output Frequencyseries.
    delta_f : float
        Frequency step for output FrequencySeries.
    low_freq_cutoff : float
        Low-frequency cutoff for output FrequencySeries.
    len_arm : float
        The arm length of TianQin, in the unit of "m".
    acc_noise_level : float
        The level of acceleration noise.
    oms_noise_level : float
        The level of OMS noise.
    duration : float
        The duration of observation, between 0 and 5, in the unit of years.

    Returns
    -------
    fseries : FrequencySeries
        The sky and polarization angle averaged
        TianQin's sensitivity curve with Galactic confusion noise.
    """
    base_curve = sensitivity_curve_tianqin_analytical(
        length, delta_f, low_freq_cutoff,
        len_arm, acc_noise_level, oms_noise_level)
    if duration < 0 or duration > 5:
        raise ValueError("Must between 0 and 5.")
    fseries_confusion = confusion_fit_tianqin(
        length, delta_f, low_freq_cutoff, duration)
    fseries = from_numpy_arrays(base_curve.sample_frequencies,
                                base_curve+fseries_confusion,
                                length, delta_f, low_freq_cutoff)

    return fseries


def sensitivity_curve_taiji_confusion(length, delta_f, low_freq_cutoff,
                                      len_arm=3e9, acc_noise_level=3e-15,
                                      oms_noise_level=8e-12, duration=1.0):
    """ The Taiji's sensitivity curve with Galactic confusion noise,
    averaged over sky and polarization angle.

    Parameters
    ----------
    length : int
        Length of output Frequencyseries.
    delta_f : float
        Frequency step for output FrequencySeries.
    low_freq_cutoff : float
        Low-frequency cutoff for output FrequencySeries.
    len_arm : float
        The arm length of Taiji, in the unit of "m".
    acc_noise_level : float
        The level of acceleration noise.
    oms_noise_level : float
        The level of OMS noise.
    duration : float
        The duration of observation, between 0 and 4, in the unit of years.

    Returns
    -------
    fseries : FrequencySeries
        The sky and polarization angle averaged
        Taiji's sensitivity curve with Galactic confusion noise.
    """
    base_curve = sensitivity_curve_taiji_analytical(
        length, delta_f, low_freq_cutoff,
        len_arm, acc_noise_level, oms_noise_level)
    if duration < 0 or duration > 4:
        raise ValueError("Must between 0 and 4.")
    fseries_confusion = confusion_fit_taiji(
        length, delta_f, low_freq_cutoff, duration)
    fseries = from_numpy_arrays(base_curve.sample_frequencies,
                                base_curve+fseries_confusion,
                                length, delta_f, low_freq_cutoff)

    return fseries


def sh_transformed_psd_lisa_tdi_XYZ(length, delta_f, low_freq_cutoff,
                                    len_arm=2.5e9, acc_noise_level=3e-15,
                                    oms_noise_level=15e-12,
                                    base_model="semi", duration=1.0,
                                    tdi=None):
    """ The TDI-1.5/2.0 PSD (X,Y,Z channel) for LISA
    with Galactic confusion noise, transformed from LISA sensitivity curve.

    Parameters
    ----------
    length : int
        Length of output Frequencyseries.
    delta_f : float
        Frequency step for output FrequencySeries.
    low_freq_cutoff : float
        Low-frequency cutoff for output FrequencySeries.
    len_arm : float
        The arm length of LISA, in the unit of "m".
    acc_noise_level : float
        The level of acceleration noise.
    oms_noise_level : float
        The level of OMS noise.
    base_model : string
        The base model of sensitivity curve, chosen from "semi" or "SciRD".
    duration : float
        The duration of observation, between 0 and 10, in the unit of years.
    tdi : string
        The version of TDI, currently only for 1.5 or 2.0.

    Returns
    -------
    fseries : FrequencySeries
        The TDI-1.5/2.0 PSD (X,Y,Z channel) for LISA with Galactic confusion
        noise, transformed from LISA sensitivity curve.
    Notes
    -----
        Please see Eq.(7,41-43) in <LISA-LCST-SGS-TN-001> for more details.
    """
    fr = np.linspace(low_freq_cutoff, (length-1)*2*delta_f, length)
    if str(tdi) in ["1.5", "2.0"]:
        response = averaged_response_lisa_tdi(fr, len_arm, tdi)
    else:
        raise ValueError("The version of TDI, currently only for 1.5 or 2.0.")
    fseries_response = from_numpy_arrays(fr, np.array(response),
                                         length, delta_f, low_freq_cutoff)
    sh = sensitivity_curve_lisa_confusion(length, delta_f, low_freq_cutoff,
                                          len_arm, acc_noise_level,
                                          oms_noise_level, base_model,
                                          duration)
    psd = 2*sh.data * fseries_response.data
    fseries = from_numpy_arrays(sh.sample_frequencies, psd,
                                length, delta_f, low_freq_cutoff)

    return fseries


def semi_analytical_psd_lisa_confusion_noise(length, delta_f, low_freq_cutoff,
                                             len_arm=2.5e9, duration=1.0,
                                             tdi=None):
    """ The TDI-1.5/2.0 PSD (X,Y,Z channel) for LISA Galactic confusion noise,
    no instrumental noise.

    Parameters
    ----------
    length : int
        Length of output Frequencyseries.
    delta_f : float
        Frequency step for output FrequencySeries.
    low_freq_cutoff : float
        Low-frequency cutoff for output FrequencySeries.
    len_arm : float
        The arm length of LISA, in the unit of "m".
    duration : float
        The duration of observation, between 0 and 10, in the unit of years.
    tdi : string
        The version of TDI, currently only for 1.5 or 2.0.

    Returns
    -------
    fseries : FrequencySeries
        The TDI-1.5/2.0 PSD (X,Y,Z channel) for LISA Galactic confusion
        noise, no instrumental noise.
    """
    fr = np.linspace(low_freq_cutoff, (length-1)*2*delta_f, length)
    if str(tdi) in ["1.5", "2.0"]:
        response = averaged_response_lisa_tdi(fr, len_arm, tdi)
    else:
        raise ValueError("The version of TDI, currently only for 1.5 or 2.0.")
    fseries_response = from_numpy_arrays(fr, np.array(response),
                                         length, delta_f, low_freq_cutoff)
    fseries_confusion = confusion_fit_lisa(
        length, delta_f, low_freq_cutoff, duration)
    psd_confusion = 2*fseries_confusion.data * fseries_response.data
    fseries = from_numpy_arrays(fseries_confusion.sample_frequencies,
                                psd_confusion, length, delta_f,
                                low_freq_cutoff)

    return fseries


def analytical_psd_tianqin_confusion_noise(length, delta_f, low_freq_cutoff,
                                           len_arm=np.sqrt(3)*1e8,
                                           duration=1.0, tdi=None):
    """ The TDI-1.5/2.0 PSD (X,Y,Z channel) for TianQin Galactic confusion
    noise, no instrumental noise.

    Parameters
    ----------
    length : int
        Length of output Frequencyseries.
    delta_f : float
        Frequency step for output FrequencySeries.
    low_freq_cutoff : float
        Low-frequency cutoff for output FrequencySeries.
    len_arm : float
        The arm length of TianQin, in the unit of "m".
    duration : float
        The duration of observation, between 0 and 5, in the unit of years.
    tdi : string
        The version of TDI. Choose from "1.5" or "2.0".

    Returns
    -------
    fseries : FrequencySeries
        The TDI-1.5/2.0 PSD (X,Y,Z channel) for TianQin Galactic confusion
        noise, no instrumental noise.
    """
    fr = np.linspace(low_freq_cutoff, (length-1)*2*delta_f, length)
    if str(tdi) in ["1.5", "2.0"]:
        response = averaged_response_tianqin_tdi(fr, len_arm, tdi)
    else:
        raise ValueError("The version of TDI, currently only for 1.5 or 2.0.")
    fseries_response = from_numpy_arrays(fr, np.array(response),
                                         length, delta_f, low_freq_cutoff)
    fseries_confusion = confusion_fit_tianqin(
        length, delta_f, low_freq_cutoff, duration)
    psd_confusion = 2*fseries_confusion.data * fseries_response.data
    fseries = from_numpy_arrays(fseries_confusion.sample_frequencies,
                                psd_confusion, length, delta_f,
                                low_freq_cutoff)

    return fseries


def analytical_psd_taiji_confusion_noise(length, delta_f, low_freq_cutoff,
                                         len_arm=3e9, duration=1.0,
                                         tdi=None):
    """ The TDI-1.5/2.0 PSD (X,Y,Z channel) for Taiji Galactic confusion
    noise, no instrumental noise.

    Parameters
    ----------
    length : int
        Length of output Frequencyseries.
    delta_f : float
        Frequency step for output FrequencySeries.
    low_freq_cutoff : float
        Low-frequency cutoff for output FrequencySeries.
    len_arm : float
        The arm length of Taiji, in the unit of "m".
    duration : float
        The duration of observation, between 0 and 4, in the unit of years.
    tdi : string
        The version of TDI. Choose from "1.5" or "2.0".

    Returns
    -------
    fseries : FrequencySeries
        The TDI-1.5/2.0 PSD (X,Y,Z channel) for Taiji Galactic confusion
        noise, no instrumental noise.
    """
    fr = np.linspace(low_freq_cutoff, (length-1)*2*delta_f, length)
    if str(tdi) in ["1.5", "2.0"]:
        response = averaged_response_taiji_tdi(fr, len_arm, tdi)
    else:
        raise ValueError("The version of TDI, currently only for 1.5 or 2.0.")
    fseries_response = from_numpy_arrays(fr, np.array(response),
                                         length, delta_f, low_freq_cutoff)
    fseries_confusion = confusion_fit_taiji(
        length, delta_f, low_freq_cutoff, duration)
    psd_confusion = 2*fseries_confusion.data * fseries_response.data
    fseries = from_numpy_arrays(fseries_confusion.sample_frequencies,
                                psd_confusion, length, delta_f,
                                low_freq_cutoff)

    return fseries


def analytical_psd_lisa_tdi_AE_confusion(length, delta_f, low_freq_cutoff,
                                         len_arm=2.5e9, acc_noise_level=3e-15,
                                         oms_noise_level=15e-12,
                                         duration=1.0, tdi=None):
    """ The TDI-1.5/2.0 PSD (A,E channel) for LISA
    with Galactic confusion noise.

    Parameters
    ----------
    length : int
        Length of output Frequencyseries.
    delta_f : float
        Frequency step for output FrequencySeries.
    low_freq_cutoff : float
        Low-frequency cutoff for output FrequencySeries.
    len_arm : float
        The arm length of LISA, in the unit of "m".
    acc_noise_level : float
        The level of acceleration noise.
    oms_noise_level : float
        The level of OMS noise.
    duration : float
        The duration of observation, between 0 and 10, in the unit of years.
    tdi : string
        The version of TDI. Choose from "1.5" or "2.0".

    Returns
    -------
    fseries : FrequencySeries
        The TDI-1.5/2.0 PSD (A,E channel) for LISA with Galactic confusion
        noise.
    """
    psd_AE = analytical_psd_lisa_tdi_AE(length, delta_f, low_freq_cutoff,
                                        len_arm, acc_noise_level,
                                        oms_noise_level, tdi)
    psd_X_confusion = semi_analytical_psd_lisa_confusion_noise(
                        length, delta_f, low_freq_cutoff,
                        len_arm, duration, tdi)
    # S_A = S_E = S_X - S_XY, confusion noise's contribution to
    # S_XY is -0.5 * psd_X_confusion, while for S_X is psd_X_confusion.
    # S_T = S_X + 2*S_XY, so S_T keeps the same.
    fseries = psd_AE + 1.5 * psd_X_confusion

    return fseries


def analytical_psd_tianqin_tdi_AE_confusion(length, delta_f, low_freq_cutoff,
                                            len_arm=np.sqrt(3)*1e8,
                                            acc_noise_level=1e-15,
                                            oms_noise_level=1e-12,
                                            duration=1.0, tdi=None):
    """ The TDI-1.5/2.0 PSD (A,E channel) for TianQin
    with Galactic confusion noise.

    Parameters
    ----------
    length : int
        Length of output Frequencyseries.
    delta_f : float
        Frequency step for output FrequencySeries.
    low_freq_cutoff : float
        Low-frequency cutoff for output FrequencySeries.
    len_arm : float
        The arm length of TianQin, in the unit of "m".
    acc_noise_level : float
        The level of acceleration noise.
    oms_noise_level : float
        The level of OMS noise.
    duration : float
        The duration of observation, between 0 and 5, in the unit of years.
    tdi : string
        The version of TDI. Choose from "1.5" or "2.0".

    Returns
    -------
    fseries : FrequencySeries
        The TDI-1.5/2.0 PSD (A,E channel) for TianQin with Galactic confusion
        noise.
    """
    psd_AE = analytical_psd_tianqin_tdi_AE(length, delta_f,
                                           low_freq_cutoff,
                                           len_arm, acc_noise_level,
                                           oms_noise_level, tdi)
    psd_X_confusion = analytical_psd_tianqin_confusion_noise(
                        length, delta_f, low_freq_cutoff,
                        len_arm, duration, tdi)
    # S_A = S_E = S_X - S_XY, confusion noise's contribution to
    # S_XY is -0.5 * psd_X_confusion, while for S_X is psd_X_confusion.
    # S_T = S_X + 2*S_XY, so S_T keeps the same.
    fseries = psd_AE + 1.5 * psd_X_confusion

    return fseries


def analytical_psd_taiji_tdi_AE_confusion(length, delta_f, low_freq_cutoff,
                                          len_arm=3e9, acc_noise_level=3e-15,
                                          oms_noise_level=8e-12,
                                          duration=1.0, tdi=None):
    """ The TDI-1.5/2.0 PSD (A,E channel) for Taiji
    with Galactic confusion noise.

    Parameters
    ----------
    length : int
        Length of output Frequencyseries.
    delta_f : float
        Frequency step for output FrequencySeries.
    low_freq_cutoff : float
        Low-frequency cutoff for output FrequencySeries.
    len_arm : float
        The arm length of Taiji, in the unit of "m".
    acc_noise_level : float
        The level of acceleration noise.
    oms_noise_level : float
        The level of OMS noise.
    duration : float
        The duration of observation, between 0 and 4, in the unit of years.
    tdi : string
        The version of TDI. Choose from "1.5" or "2.0".

    Returns
    -------
    fseries : FrequencySeries
        The TDI-1.5/2.0 PSD (A,E channel) for Taiji with Galactic confusion
        noise.
    """
    psd_AE = analytical_psd_taiji_tdi_AE(length, delta_f, low_freq_cutoff,
                                         len_arm, acc_noise_level,
                                         oms_noise_level, tdi)
    psd_X_confusion = analytical_psd_taiji_confusion_noise(
                        length, delta_f, low_freq_cutoff,
                        len_arm, duration, tdi)
    # S_A = S_E = S_X - S_XY, confusion noise's contribution to
    # S_XY is -0.5 * psd_X_confusion, while for S_X is psd_X_confusion.
    # S_T = S_X + 2*S_XY, so S_T keeps the same.
    fseries = psd_AE + 1.5 * psd_X_confusion

    return fseries
