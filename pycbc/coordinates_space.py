# Copyright (C) 2023  Shichao Wu, Alex Nitz
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
This module provides coordinate transformations related to space-borne
detectors, such as coordinate transformations between space-borne detectors
and ground-based detectors.
"""

import numpy as np
from scipy.spatial.transform import Rotation
from scipy.optimize import fsolve
from astropy.constants import c, au
from astropy.coordinates import BarycentricMeanEcliptic, PrecessedGeocentric


YRSID_SI = 31558149.763545603  # same as BBHx, but not for lal.YRJUL_SI
OMEGA_0 = 1.99098659277e-7


def localization_to_propagation_vector(lamda, beta):
    """ Converting the sky localization to the corresponding
    propagation unit vector of a GW signal.

    Parameters
    ----------
    lamda : float
        The right ascension or ecliptic longitude, in the unit of 'deg'.
    beta : float
        The declination or ecliptic latitude, in the unit of 'deg'.

    Returns
    -------
    [[x], [y], [z]] : numpy.array
        The propagation unit vector of that GW signal.
    """
    x = -np.cos(beta) * np.cos(lamda)
    y = -np.cos(beta) * np.sin(lamda)
    z = -np.sin(beta)

    return np.array([[x], [y], [z]])


def rotation_matrix_SSBtoLISA(alpha):
    """ The rotation matrix (of frame basis) from SSB frame to LISA frame.

    Parameters
    ----------
    alpha : float
        The angular displacement of LISA in SSB frame, in the unit of 'deg'.

    Returns
    -------
    r_total : numpy.array
        A 3x3 rotation matrix from SSB frame to LISA frame.
    """
    r = Rotation.from_rotvec([
        [0, 0, alpha],
        [0, -np.pi/3, 0],
        [0, 0, -alpha]
    ]).as_matrix()
    r_total = np.array(r[0]) @ np.array(r[1]) @ np.array(r[2])

    return r_total


def propagation_vector_to_localization(k):
    """ Converting the propagation unit vector to the corresponding
    sky localization of a GW signal.

    Parameters
    ----------
    k : numpy.array
        The propagation unit vector of a GW signal.

    Returns
    -------
    [lamda, beta] : numpy.array
        The sky localization of that GW signal.
    """
    # beta already within [-pi/2, pi/2]
    beta = np.float64(np.arcsin(-k[2]))
    lamda = np.float64(np.arctan2(-k[1]/np.cos(beta), -k[0]/np.cos(beta)))
    # lamda should within [0, 2*pi]
    lamda = np.mod(lamda, 2*np.pi)

    return np.array([lamda, beta])


def polarization_newframe(psi, k, rotation_matrix):
    """ Converting a polarization angle from a frame to a new frame
    by using rotation matrix method.

    Parameters
    ----------
    psi : float
        The polarization angle in the old frame, in the unit of 'deg'.
    k : numpy.array
        The propagation unit vector of a GW signal in the old frame.
    rotation_matrix : numpy.array
        The rotation matrix (of frame basis) from the old frame to
        the new frame.

    Returns
    -------
    psi_newframe : float
        The polarization angle in the new frame of that GW signal.
    """
    lamda, beta = propagation_vector_to_localization(k)
    u = np.array([[np.sin(lamda)], [-np.cos(lamda)], [0]])
    rotation_vector = psi * k
    rotation_psi = Rotation.from_rotvec(rotation_vector.T[0])
    p = rotation_psi.apply(u.T[0]).reshape(3, 1)
    p_newframe = rotation_matrix.T @ p
    k_newframe = rotation_matrix.T @ k
    lamda_newframe, beta_newframe = \
        propagation_vector_to_localization(k_newframe)
    u_newframe = np.array([[np.sin(lamda_newframe)],
                           [-np.cos(lamda_newframe)], [0]])
    v_newframe = np.array([[-np.sin(beta_newframe) * np.cos(lamda_newframe)],
                           [-np.sin(beta_newframe) * np.sin(lamda_newframe)],
                           [np.cos(beta_newframe)]])
    p_dot_u_newframe = np.vdot(p_newframe, u_newframe)
    p_dot_v_newframe = np.vdot(p_newframe, v_newframe)
    psi_newframe = np.arctan2(p_dot_v_newframe, p_dot_u_newframe)
    psi_newframe = np.mod(psi_newframe, np.pi)

    return psi_newframe


def tL_from_SSB(tSSB, lamdaSSB, betaSSB, t0=0.0):
    """ Calculating the time when a GW signal arrives at the barycenter
    of LISA, by using the time and sky localization in SSB frame.

    Parameters
    ----------
    tSSB : float
        The time when a GW signal arrives at the origin of SSB frame.
        In the unit of 's'.
    lamdaSSB : float
        The ecliptic longitude of a GW signal in SSB frame.
        In the unit of 'deg'.
    betaSSB : float
        The ecliptic latitude of a GW signal in SSB frame.
        In the unit of 'deg'.
    t0 : float
        The initial time offset, in the unit of 'year', default is 0.

    Returns
    -------
    tL : float
        The time when a GW signal arrives at the origin of LISA frame.
    """
    R_ORBIT = au.value
    OMEGA_L = 2 * np.pi / YRSID_SI
    t0 *= YRSID_SI
    alpha = (tSSB + t0) * OMEGA_L
    k = localization_to_propagation_vector(lamdaSSB, betaSSB)
    v_LISA = np.array([[R_ORBIT * np.cos(alpha)],
                       [R_ORBIT * np.sin(alpha)], [0]])
    tL = tSSB + np.vdot(k, v_LISA) / c.value

    return tL


def tSSB_from_tL(tL, lamdaSSB, betaSSB, t0=0.0):
    """ Calculating the time when a GW signal arrives at the barycenter
    of SSB, by using the time in LISA frame and sky localization in SSB frame.

    Parameters
    ----------
    tL : float
        The time when a GW signal arrives at the origin of LISA frame.
        In the unit of 's'.
    lamdaSSB : float
        The ecliptic longitude of a GW signal in SSB frame.
        In the unit of 'deg'.
    betaSSB : float
        The ecliptic latitude of a GW signal in SSB frame.
        In the unit of 'deg'.
    t0 : float
        The initial time offset, in the unit of 'year', default is 0.

    Returns
    -------
    tSSB : float
        The time when a GW signal arrives at the origin of SSB frame.
    """
    t0 *= YRSID_SI
    R_ORBIT = au.value
    OMEGA_L = 2 * np.pi / YRSID_SI

    def equation(tSSB):
        return tL - tSSB + R_ORBIT * np.cos(betaSSB) * (
            np.cos(lamdaSSB) * np.cos(OMEGA_L * (tSSB + t0)) +
            np.sin(lamdaSSB) * np.sin(OMEGA_L * (tSSB + t0))
        ) / c.value

    return fsolve(equation, tL)[0]


def SSB_to_LISA(tSSB, lamdaSSB, betaSSB, psiSSB, t0):
    """ Converting the arrive time, the sky localization, and the polarization
    from the SSB frame to the LISA frame.

    Parameters
    ----------
    tSSB : float
        The time when a GW signal arrives at the origin of SSB frame.
        In the unit of 's'.
    lamdaSSB : float
        The ecliptic longitude of a GW signal in SSB frame.
        In the unit of 'deg'.
    betaSSB : float
        The ecliptic latitude of a GW signal in SSB frame.
        In the unit of 'deg'.
    psiSSB : float
        The polarization angle of a GW signal in SSB frame.
        In the unit of 'deg'.
    t0 : float
        The initial time offset, in the unit of 'year', default is 0.

    Returns
    -------
    (tL, lamdaL, betaL, psiL) : tuple
    tL : float
        The time when a GW signal arrives at the origin of LISA frame.
        In the unit of 's'.
    lamdaL : float
        The longitude of a GW signal in LISA frame, in the unit of 'deg'.
    betaL : float
        The latitude of a GW signal in LISA frame, in the unit of 'deg'.
    psiL : float
        The polarization angle of a GW signal in LISA frame.
        In the unit of 'deg'.
    """
    tL = tL_from_SSB(tSSB, lamdaSSB, betaSSB, t0)
    kSSB = localization_to_propagation_vector(lamdaSSB, betaSSB)
    alpha = OMEGA_0 * (tSSB + t0*YRSID_SI)
    rotation_matrix_L = rotation_matrix_SSBtoLISA(alpha)
    kL = rotation_matrix_L.T @ kSSB
    lamdaL, betaL = propagation_vector_to_localization(kL)
    psiL = polarization_newframe(psiSSB, kSSB, rotation_matrix_L)

    return (tL, lamdaL, betaL, psiL)


def LISA_to_SSB(tL, lamdaL, betaL, psiL, t0):
    """ Converting the arrive time, the sky localization, and the polarization
    from the LISA frame to the SSB frame.

    Parameters
    ----------
    tL : float
        The time when a GW signal arrives at the origin of LISA frame.
        In the unit of 's'.
    lamdaL : float
        The longitude of a GW signal in LISA frame, in the unit of 'deg'.
    betaL : float
        The latitude of a GW signal in LISA frame, in the unit of 'deg'.
    psiL : float
        The polarization angle of a GW signal in LISA frame.
        In the unit of 'deg'.
    t0 : float
        The initial time offset, in the unit of 'year', default is 0.

    Returns
    -------
    (tSSB_approx, lamdaSSB_approx, betaSSB_approx, psiSSB) : tuple
    tSSB_approx : float
        The time when a GW signal arrives at the origin of SSB frame.
        In the unit of 's'.
    lamdaSSB_approx : float
        The ecliptic longitude of a GW signal in SSB frame.
        In the unit of 'deg'.
    betaSSB_approx : float
        The ecliptic latitude of a GW signal in SSB frame.
        In the unit of 'deg'.
    psiSSB : float
        The polarization angle of a GW signal in SSB frame.
        In the unit of 'deg'.
    """
    lamdaSSB_approx = 0.0
    betaSSB_approx = 0.0
    tSSB_approx = tL
    kL = localization_to_propagation_vector(lamdaL, betaL)
    for i in range(3):
        alpha = OMEGA_0 * (tSSB_approx + t0*YRSID_SI)
        rotation_matrix_L = rotation_matrix_SSBtoLISA(alpha)
        kSSB_approx = rotation_matrix_L @ kL
        lamdaSSB_approx, betaSSB_approx = \
            propagation_vector_to_localization(kSSB_approx)
        tSSB_approx = tSSB_from_tL(tL, lamdaSSB_approx, betaSSB_approx, t0)
    psiSSB = polarization_newframe(psiL, kL, rotation_matrix_L.T)

    return (tSSB_approx, lamdaSSB_approx, betaSSB_approx, psiSSB)


def rotation_matrix_SSBtoGEO(epsilon=np.deg2rad(23.44)):
    """ The rotation matrix (of frame basis) from SSB frame to
    geocentric frame.

    Parameters
    ----------
    epsilon : float
        The Earth's axial tilt (obliquity), in the unit of 'deg'.

    Returns
    -------
    r : numpy.array
        A 3x3 rotation matrix from SSB frame to geocentric frame.
    """
    r = Rotation.from_rotvec([
        [-epsilon, 0, 0]
    ]).as_matrix()

    return np.array(r[0])


def tG_from_SSB(tSSB, lamdaSSB, betaSSB, t0=0.0):
    """ Calculating the time when a GW signal arrives at the barycenter
    of the Earth, by using the time and sky localization in SSB frame.

    Parameters
    ----------
    tSSB : float
        The time when a GW signal arrives at the origin of SSB frame.
        In the unit of 's'.
    lamdaSSB : float
        The ecliptic longitude of a GW signal in SSB frame.
        In the unit of 'deg'.
    betaSSB : float
        The ecliptic latitude of a GW signal in SSB frame.
        In the unit of 'deg'.
    t0 : float
        The initial time offset, in the unit of 'year', default is 0.

    Returns
    -------
    tG : float
        The time when a GW signal arrives at the origin of geocentric frame.
    """
    R_ORBIT = au.value
    OMEGA_G = 2 * np.pi / YRSID_SI
    t0 *= YRSID_SI
    alpha = (tSSB + t0) * OMEGA_G
    k = localization_to_propagation_vector(lamdaSSB, betaSSB)
    v_GEO = np.array([[R_ORBIT * np.cos(alpha)],
                      [R_ORBIT * np.sin(alpha)], [0]])
    tG = tSSB + np.vdot(k, v_GEO) / c.value

    return tG


def tSSB_from_tG(tG, lamdaSSB, betaSSB, t0=0.0):
    """ Calculating the time when a GW signal arrives at the barycenter
    of SSB, by using the time in geocentric frame and sky localization
    in SSB frame.

    Parameters
    ----------
    tG : float
        The time when a GW signal arrives at the origin of geocentric frame.
        In the unit of 's'.
    lamdaSSB : float
        The ecliptic longitude of a GW signal in SSB frame.
        In the unit of 'deg'.
    betaSSB : float
        The ecliptic latitude of a GW signal in SSB frame.
        In the unit of 'deg'.
    t0 : float
        The initial time offset, in the unit of 'year', default is 0.

    Returns
    -------
    tSSB : float
        The time when a GW signal arrives at the origin of SSB frame.
    """
    t0 *= YRSID_SI
    R_ORBIT = au.value
    OMEGA_G = 2 * np.pi / YRSID_SI

    def equation(tSSB):
        return tG - tSSB + R_ORBIT * np.cos(betaSSB) * (
            np.cos(lamdaSSB) * np.cos(OMEGA_G * (tSSB + t0)) +
            np.sin(lamdaSSB) * np.sin(OMEGA_G * (tSSB + t0))
        ) / c.value

    return fsolve(equation, tG)[0]


def SSB_to_GEO(tSSB, lamdaSSB, betaSSB, psiSSB, t0, useAstropy=True):
    """ Converting the arrive time, the sky localization, and the polarization
    from the SSB frame to the geocentric frame.

    Parameters
    ----------
    tSSB : float
        The time when a GW signal arrives at the origin of SSB frame.
        In the unit of 's'.
    lamdaSSB : float
        The ecliptic longitude of a GW signal in SSB frame.
        In the unit of 'deg'.
    betaSSB : float
        The ecliptic latitude of a GW signal in SSB frame.
        In the unit of 'deg'.
    psiSSB : float
        The polarization angle of a GW signal in SSB frame.
        In the unit of 'deg'.
    t0 : float
        The initial time offset, in the unit of 'year', default is 0.
    useAstropy : bool
        Using Astropy to calculate the sky localization or not.
        Default is True.

    Returns
    -------
    (tG, lamdaG, betaG, psiG) : tuple
    tG : float
        The time when a GW signal arrives at the origin of geocentric frame.
        In the unit of 's'.
    lamdaG : float
        The longitude of a GW signal in geocentric frame.
        In the unit of 'deg'.
    betaG : float
        The latitude of a GW signal in geocentric frame.
        In the unit of 'deg'.
    psiG : float
        The polarization angle of a GW signal in geocentric frame.
        In the unit of 'deg'.
    """
    tG = tG_from_SSB(tSSB, lamdaSSB, betaSSB, t0)
    kSSB = localization_to_propagation_vector(lamdaSSB, betaSSB)
    rotation_matrix_G = rotation_matrix_SSBtoGEO()
    if useAstropy:
        c = BarycentricMeanEcliptic(lon=lamdaSSB, lat=betaSSB,
                                    equinox='J2000', unit='deg')
        geo_sky = c.transform_to(PrecessedGeocentric)
        lamdaG = geo_sky.ra.rad
        betaG = geo_sky.dec.rad
    else:
        kG = rotation_matrix_G.T @ kSSB
        lamdaG, betaG = propagation_vector_to_localization(kG)
    psiG = polarization_newframe(psiSSB, kSSB, rotation_matrix_G)

    return (tG, lamdaG, betaG, psiG)


def GEO_to_SSB(tG, lamdaG, betaG, psiG, t0, useAstropy=True):
    """ Converting the arrive time, the sky localization, and the polarization
    from the geocentric frame to the SSB frame.

    Parameters
    ----------
    tG : float
        The time when a GW signal arrives at the origin of geocentric frame.
        In the unit of 's'.
    lamdaG : float
        The longitude of a GW signal in geocentric frame.
        In the unit of 'deg'.
    betaG : float
        The latitude of a GW signal in geocentric frame.
        In the unit of 'deg'.
    psiG : float
        The polarization angle of a GW signal in geocentric frame.
        In the unit of 'deg'.
    t0 : float
        The initial time offset, in the unit of 'year', default is 0.
    useAstropy : bool
        Using Astropy to calculate the sky localization or not.
        Default is True.

    Returns
    -------
    (tSSB_approx, lamdaSSB_approx, betaSSB_approx, psiSSB) : tuple
    tSSB_approx : float
        The time when a GW signal arrives at the origin of SSB frame.
        In the unit of 's'.
    lamdaSSB_approx : float
        The ecliptic longitude of a GW signal in SSB frame.
        In the unit of 'deg'.
    betaSSB_approx : float
        The ecliptic latitude of a GW signal in SSB frame.
        In the unit of 'deg'.
    psiSSB : float
        The polarization angle of a GW signal in SSB frame.
        In the unit of 'deg'.
    """
    lamdaSSB_approx = 0.0
    betaSSB_approx = 0.0
    tSSB_approx = tG
    kG = localization_to_propagation_vector(lamdaG, betaG)
    for i in range(3):
        rotation_matrix_G = rotation_matrix_SSBtoGEO()
        kSSB_approx = rotation_matrix_G @ kG
        lamdaSSB_approx, betaSSB_approx = \
            propagation_vector_to_localization(kSSB_approx)
        tSSB_approx = tSSB_from_tG(tG, lamdaSSB_approx, betaSSB_approx, t0)
    psiSSB = polarization_newframe(psiG, kG, rotation_matrix_G.T)
    if useAstropy:
        c = PrecessedGeocentric(ra=lamdaG, dec=betaG,
                                equinox='J2000', unit='deg')
        ssb_sky = c.transform_to(BarycentricMeanEcliptic)
        lamdaSSB_approx = ssb_sky.lat.rad
        betaSSB_approx = ssb_sky.lon.rad

    return (tSSB_approx, lamdaSSB_approx, betaSSB_approx, psiSSB)


def LISA_to_GEO(tL, lamdaL, betaL, psiL, t0, useAstropy=True):
    """ Converting the arrive time, the sky localization, and the polarization
    from the LISA frame to the geocentric frame.

    Parameters
    ----------
    tL : float
        The time when a GW signal arrives at the origin of LISA frame.
        In the unit of 's'.
    lamdaL : float
        The longitude of a GW signal in LISA frame, in the unit of 'deg'.
    betaL : float
        The latitude of a GW signal in LISA frame, in the unit of 'deg'.
    psiL : float
        The polarization angle of a GW signal in LISA frame.
        In the unit of 'deg'.
    t0 : float
        The initial time offset, in the unit of 'year', default is 0.
    useAstropy : bool
        Using Astropy to calculate the sky localization or not.
        Default is True.

    Returns
    -------
    (tG, lamdaG, betaG, psiG) : tuple
    tG : float
        The time when a GW signal arrives at the origin of geocentric frame.
        In the unit of 's'.
    lamdaG : float
        The ecliptic longitude of a GW signal in geocentric frame.
        In the unit of 'deg'.
    betaG : float
        The ecliptic latitude of a GW signal in geocentric frame.
        In the unit of 'deg'.
    psiG : float
        The polarization angle of a GW signal in geocentric frame.
        In the unit of 'deg'.
    """
    tSSB, lamdaSSB, betaSSB, psiSSB = LISA_to_SSB(tL, lamdaL, betaL, psiL, t0)
    tG, lamdaG, betaG, psiG = SSB_to_GEO(
        tSSB, lamdaSSB, betaSSB, psiSSB, t0, useAstropy=True)

    return (tG, lamdaG, betaG, psiG)


def GEO_to_LISA(tG, lamdaG, betaG, psiG, t0, useAstropy=True):
    """ Converting the arrive time, the sky localization, and the polarization
    from the geocentric frame to the LISA frame.

    Parameters
    ----------
    tG : float
        The time when a GW signal arrives at the origin of geocentric frame.
        In the unit of 's'.
    lamdaG : float
        The longitude of a GW signal in geocentric frame.
        In the unit of 'deg'.
    betaG : float
        The latitude of a GW signal in geocentric frame.
        In the unit of 'deg'.
    psiG : float
        The polarization angle of a GW signal in geocentric frame.
        In the unit of 'deg'.
    t0 : float
        The initial time offset, in the unit of 'year', default is 0.
    useAstropy : bool
        Using Astropy to calculate the sky localization or not.
        Default is True.

    Returns
    -------
    (tL, lamdaL, betaL, psiL) : tuple
    tL : float
        The time when a GW signal arrives at the origin of LISA frame.
        In the unit of 's'.
    lamdaL : float
        The longitude of a GW signal in LISA frame, in the unit of 'deg'.
    betaL : float
        The latitude of a GW signal in LISA frame, in the unit of 'deg'.
    psiG : float
        The polarization angle of a GW signal in LISA frame.
        In the unit of 'deg'.
    """
    tSSB, lamdaSSB, betaSSB, psiSSB = GEO_to_SSB(
        tG, lamdaG, betaG, psiG, t0, useAstropy=True)
    tL, lamdaL, betaL, psiL = SSB_to_LISA(tSSB, lamdaSSB, betaSSB, psiSSB, t0)

    return (tL, lamdaL, betaL, psiL)


__all__ = ['localization_to_propagation_vector', 'rotation_matrix_SSBtoLISA',
           'propagation_vector_to_localization', 'polarization_newframe',
           'tL_from_SSB', 'tSSB_from_tL', 'SSB_to_LISA', 'LISA_to_SSB',
           'rotation_matrix_SSBtoGEO', 'tG_from_SSB', 'tSSB_from_tG',
           'SSB_to_GEO', 'GEO_to_SSB', 'LISA_to_GEO', 'GEO_to_LISA',
           ]
