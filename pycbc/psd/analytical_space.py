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
borne detectors, such as LISA. Based on LISA technical note
<LISA-LCST-SGS-TN-001>, LDC manual <LISA-LCST-SGS-MAN-001>,
and paper <10.1088/1361-6382/ab1101>.
"""

import numpy as np
from astropy import constants
from pycbc.psd.read import from_numpy_arrays


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
        Pease see Eq.(11-13) in <LISA-LCST-SGS-TN-001> for more details.
    """
    s_acc = acc_noise_level**2 * (1+(4e-4/f)**2)*(1+(f/8e-3)**4)
    s_acc_d = s_acc * (2*np.pi*f)**(-4)
    s_acc_nu = (2*np.pi*f/constants.c.value)**2 * s_acc_d

    return s_acc_nu


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
        Pease see Eq.(9-10) in <LISA-LCST-SGS-TN-001> for more details.
    """
    s_oms_d = oms_noise_level**2 * (1+(2e-3/f)**4)
    s_oms_nu = s_oms_d * (2*np.pi*f/constants.c.value)**2

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
    [low_freq_component, high_freq_component] : list
        The PSD value or array for acceleration and OMS noise.
    """
    low_freq_component = psd_lisa_acc_noise(f, acc_noise_level)
    high_freq_component = psd_lisa_oms_noise(f, oms_noise_level)

    return [low_freq_component, high_freq_component]


def omega_length(f, len_arm=2.5e9):
    """ The function to calculate 2*pi*f*LISA_arm_length.
    Parameters
    ----------
    f : float or numpy.array
        The frequency or frequency range, in the unit of "Hz".
    len_arm : float
        The arm length of LISA.
    Returns
    -------
    omega_len : float or numpy.array
        The value of 2*pi*f*LISA_arm_length.
    """
    omega_len = 2*np.pi*f * len_arm/constants.c.value

    return omega_len


def analytical_psd_lisa_tdi_1p5_XYZ(length, delta_f, low_freq_cutoff,
                                    len_arm=2.5e9, acc_noise_level=3e-15,
                                    oms_noise_level=15e-12):
    """ The TDI-1.5 analytical PSD (X,Y,Z channel) for LISA.
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
        The TDI-1.5 PSD (X,Y,Z channel) for LISA.
    Notes
    -----
        Pease see Eq.(19) in <LISA-LCST-SGS-TN-001> for more details.
    """
    len_arm = np.float64(len_arm)
    acc_noise_level = np.float64(acc_noise_level)
    oms_noise_level = np.float64(oms_noise_level)
    psd = []
    fr = np.linspace(low_freq_cutoff, (length-1)*2*delta_f, length)
    for f in fr:
        [s_acc_nu, s_oms_nu] = lisa_psd_components(
                                f, acc_noise_level, oms_noise_level)
        omega_len = omega_length(f, len_arm)
        psd.append(16*(np.sin(omega_len))**2 *
                   (s_oms_nu+s_acc_nu*(3+np.cos(omega_len))))
    fseries = from_numpy_arrays(fr, np.array(psd),
                                length, delta_f, low_freq_cutoff)

    return fseries


def analytical_psd_lisa_tdi_2p0_XYZ(length, delta_f, low_freq_cutoff,
                                    len_arm=2.5e9, acc_noise_level=3e-15,
                                    oms_noise_level=15e-12):
    """ The TDI-2.0 analytical PSD (X,Y,Z channel) for LISA.
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
        The TDI-2.0 PSD (X,Y,Z channel) for LISA.
    Notes
    -----
        Pease see Eq.(20) in <LISA-LCST-SGS-TN-001> for more details.
    """
    len_arm = np.float64(len_arm)
    acc_noise_level = np.float64(acc_noise_level)
    oms_noise_level = np.float64(oms_noise_level)
    psd = []
    fr = np.linspace(low_freq_cutoff, (length-1)*2*delta_f, length)
    for f in fr:
        [s_acc_nu, s_oms_nu] = lisa_psd_components(
                                f, acc_noise_level, oms_noise_level)
        omega_len = omega_length(f, len_arm)
        psd.append(64*(np.sin(omega_len))**2 * (np.sin(2*omega_len))**2 *
                   (s_oms_nu+s_acc_nu*(3+np.cos(2*omega_len))))
    fseries = from_numpy_arrays(fr, np.array(psd),
                                length, delta_f, low_freq_cutoff)

    return fseries


def analytical_csd_lisa_tdi_1p5_XY(length, delta_f, low_freq_cutoff,
                                   len_arm=2.5e9, acc_noise_level=3e-15,
                                   oms_noise_level=15e-12):
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
    Returns
    -------
    fseries : FrequencySeries
        The CSD between LISA's TDI-1.5 channel X and Y.
    Notes
    -----
        Pease see Eq.(56) in <LISA-LCST-SGS-MAN-001(Radler)> for more details.
    """
    len_arm = np.float64(len_arm)
    acc_noise_level = np.float64(acc_noise_level)
    oms_noise_level = np.float64(oms_noise_level)
    csd = []
    fr = np.linspace(low_freq_cutoff, (length-1)*2*delta_f, length)
    for f in fr:
        omega_len = omega_length(f, len_arm)
        [s_acc_nu, s_oms_nu] = lisa_psd_components(
                                f, acc_noise_level, oms_noise_level)
        csd.append(-8*np.sin(omega_len)**2 * np.cos(omega_len) *
                   (s_oms_nu+4*s_acc_nu))
    fseries = from_numpy_arrays(fr, np.array(csd),
                                length, delta_f, low_freq_cutoff)
    return fseries


def analytical_psd_lisa_tdi_1p5_AE(length, delta_f, low_freq_cutoff,
                                   len_arm=2.5e9, acc_noise_level=3e-15,
                                   oms_noise_level=15e-12):
    """ The PSD of LISA's TDI-1.5 channel A and E.
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
        The PSD of LISA's TDI-1.5 channel A and E.
    Notes
    -----
        Pease see Eq.(58) in <LISA-LCST-SGS-MAN-001(Radler)> for more details.
    """
    len_arm = np.float64(len_arm)
    acc_noise_level = np.float64(acc_noise_level)
    oms_noise_level = np.float64(oms_noise_level)
    psd = []
    fr = np.linspace(low_freq_cutoff, (length-1)*2*delta_f, length)
    for f in fr:
        [s_acc_nu, s_oms_nu] = lisa_psd_components(
                                f, acc_noise_level, oms_noise_level)
        omega_len = omega_length(f, len_arm)
        psd.append(8*(np.sin(omega_len))**2 *
                   (4*(1+np.cos(omega_len)+np.cos(omega_len)**2)*s_acc_nu +
                   (2+np.cos(omega_len))*s_oms_nu))
    fseries = from_numpy_arrays(fr, np.array(psd),
                                length, delta_f, low_freq_cutoff)

    return fseries


def analytical_psd_lisa_tdi_1p5_T(length, delta_f, low_freq_cutoff,
                                  len_arm=2.5e9, acc_noise_level=3e-15,
                                  oms_noise_level=15e-12):
    """ The PSD of LISA's TDI-1.5 channel T.
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
        The PSD of LISA's TDI-1.5 channel T.
    Notes
    -----
        Pease see Eq.(59) in <LISA-LCST-SGS-MAN-001(Radler)> for more details.
    """
    len_arm = np.float64(len_arm)
    acc_noise_level = np.float64(acc_noise_level)
    oms_noise_level = np.float64(oms_noise_level)
    psd = []
    fr = np.linspace(low_freq_cutoff, (length-1)*2*delta_f, length)
    for f in fr:
        [s_acc_nu, s_oms_nu] = lisa_psd_components(
                                f, acc_noise_level, oms_noise_level)
        omega_len = omega_length(f, len_arm)
        psd.append(32*np.sin(omega_len)**2 * np.sin(omega_len/2)**2 *
                   (4*s_acc_nu*np.sin(omega_len/2)**2 + s_oms_nu))
    fseries = from_numpy_arrays(fr, np.array(psd),
                                length, delta_f, low_freq_cutoff)

    return fseries


def averaged_lisa_fplus_sq_approx(f, len_arm=2.5e9):
    """ An approximant for LISA's squared antenna response function,
    averaged over sky and polarization angle.
    Parameters
    ----------
    f : float or numpy.array
        The frequency or frequency range, in the unit of "Hz".
    len_arm : float
        The arm length of LISA, in the unit of "m".
    Returns
    -------
    fp_sq_approx : float or numpy.array
        The sky and polarization angle averaged squared antenna response.
    Notes
    -----
        Pease see Eq.(36) in <LISA-LCST-SGS-TN-001> for more details.
    """
    from os import getcwd, path
    from urllib import request
    from scipy.interpolate import interp1d

    if len_arm != 2.5e9:
        raise Exception("Currently only support 'len_arm=2.5e9'.")
    cwd = getcwd()
    if path.exists(cwd+"/AvFXp2_Raw.npy") is False:
        url = "https://zenodo.org/record/7497853/files/AvFXp2_Raw.npy"
        request.urlretrieve(url, cwd+"/AvFXp2_Raw.npy")
    freqs, fp_sq = np.load(cwd+"/AvFXp2_Raw.npy")
    # Padding the end.
    freqs = np.append(freqs, 2)
    fp_sq = np.append(fp_sq, 0.0012712348970728724)
    fp_sq_interp = interp1d(freqs, fp_sq, kind='linear',
                            fill_value="extrapolate")
    fp_sq_approx = fp_sq_interp(f)/16

    return fp_sq_approx


def averaged_response_lisa_tdi_1p5(f, len_arm=2.5e9):
    """ LISA's TDI-1.5 response function to GW,
    averaged over sky and polarization angle.
    Parameters
    ----------
    f : float or numpy.array
        The frequency or frequency range, in the unit of "Hz".
    len_arm : float
        The arm length of LISA, in the unit of "m".
    Returns
    -------
    response_tdi_1p5 : float or numpy.array
        The sky and polarization angle averaged TDI-1.5 response to GW.
    Notes
    -----
        Pease see Eq.(39) in <LISA-LCST-SGS-TN-001> for more details.
    """
    omega_len = omega_length(f, len_arm)
    ave_fp2 = averaged_lisa_fplus_sq_approx(f, len_arm)
    response_tdi_1p5 = (4*omega_len)**2 * np.sin(omega_len)**2 * ave_fp2

    return response_tdi_1p5


def averaged_response_lisa_tdi_2p0(f, len_arm=2.5e9):
    """ LISA's TDI-2.0 response function to GW,
    averaged over sky and polarization angle.
    Parameters
    ----------
    f : float or numpy.array
        The frequency or frequency range, in the unit of "Hz".
    len_arm : float
        The arm length of LISA, in the unit of "m".
    Returns
    -------
    response_tdi_2p0 : float or numpy.array
        The sky and polarization angle averaged TDI-2.0 response to GW.
    Notes
    -----
        Pease see Eq.(40) in <LISA-LCST-SGS-TN-001> for more details.
    """
    omega_len = omega_length(f, len_arm)
    response_tdi_1p5 = averaged_response_lisa_tdi_1p5(f, len_arm)
    response_tdi_2p0 = response_tdi_1p5 * (2*np.sin(2*omega_len))**2

    return response_tdi_2p0


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
        Pease see Eq.(42-43) in <LISA-LCST-SGS-TN-001> for more details.
    """
    sense_curve = []
    len_arm = np.float64(len_arm)
    acc_noise_level = np.float64(acc_noise_level)
    oms_noise_level = np.float64(oms_noise_level)
    fr = np.linspace(low_freq_cutoff, (length-1)*2*delta_f, length)
    fp_sq = averaged_lisa_fplus_sq_approx(fr, len_arm)
    for i in range(len(fr)):
        [s_acc_nu, s_oms_nu] = lisa_psd_components(
                                fr[i], acc_noise_level, oms_noise_level)
        omega_len = 2*np.pi*fr[i] * len_arm/constants.c.value
        sense_curve.append((s_oms_nu + s_acc_nu*(3+np.cos(2*omega_len))) /
                           (omega_len**2*fp_sq[i]))
    fseries = from_numpy_arrays(fr, np.array(sense_curve)/2,
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
        Pease see Eq.(114) in <LISA-LCST-SGS-TN-001> for more details.
    """
    sense_curve = []
    fr = np.linspace(low_freq_cutoff, (length-1)*2*delta_f, length)
    for f in fr:
        s_I = 5.76e-48 * (1+(4e-4/f)**2)
        s_II = 3.6e-41
        R = 1 + (f/2.5e-2)**2
        sense_curve.append(10/3 * (s_I/(2*np.pi*f)**4+s_II) * R)
    fseries = from_numpy_arrays(fr, sense_curve,
                                length, delta_f, low_freq_cutoff)

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
        Pease see Eq.(85-86) in <LISA-LCST-SGS-TN-001> for more details.
    """
    if base_model == "semi":
        base_curve = sensitivity_curve_lisa_semi_analytical(
            length, delta_f, low_freq_cutoff,
            len_arm, acc_noise_level, oms_noise_level)
    elif base_curve == "SciRD":
        base_curve = sensitivity_curve_lisa_SciRD(
            length, delta_f, low_freq_cutoff)
    else:
        raise Exception("Must choose from 'semi' or 'SciRD'.")
    if duration < 0 or duration > 10:
        raise Exception("Must between 0 and 10.")
    fr = np.linspace(low_freq_cutoff, (length-1)*2*delta_f, length)
    sh_confusion = []
    f1 = 10**(-0.25*np.log10(duration)-2.7)
    fk = 10**(-0.27*np.log10(duration)-2.47)
    for f in fr:
        sh_confusion.append(0.5*1.14e-44*f**(-7/3)*np.exp(-(f/f1)**1.8) *
                            (1.0+np.tanh((fk-f)/(0.31e-3))))
    fseries_confusion = from_numpy_arrays(fr, np.array(sh_confusion),
                                          length, delta_f, low_freq_cutoff)
    fseries = from_numpy_arrays(base_curve.sample_frequencies,
                                base_curve+fseries_confusion,
                                length, delta_f, low_freq_cutoff)

    return fseries


def sh_transformed_psd_lisa_tdi_XYZ(length, delta_f, low_freq_cutoff,
                                    len_arm=2.5e9, acc_noise_level=3e-15,
                                    oms_noise_level=15e-12,
                                    base_model="semi", duration=1.0,
                                    tdi="1.5"):
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
        Pease see Eq.(7,41-43) in <LISA-LCST-SGS-TN-001> for more details.
    """
    fr = np.linspace(low_freq_cutoff, (length-1)*2*delta_f, length)
    if tdi == "1.5":
        response = averaged_response_lisa_tdi_1p5(fr, len_arm)
    elif tdi == "2.0":
        response = averaged_response_lisa_tdi_2p0(fr, len_arm)
    else:
        raise Exception("The version of TDI, currently only for 1.5 or 2.0.")
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
