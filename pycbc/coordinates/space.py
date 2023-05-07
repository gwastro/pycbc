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
from astropy import units
from astropy.constants import c, au
from astropy.time import Time
from astropy.coordinates import BarycentricMeanEcliptic, PrecessedGeocentric
from astropy.coordinates import get_body_barycentric
from astropy.coordinates import SkyCoord


# This constant makes sure LISA is behind the Earth by 19-23 degrees.
# Making this a stand-alone constant will also make it readable by
# the waveform plugin. In the unit of 's'.
TIME_OFFSET = 7365189.431698299

# "rotation_matrix_ssb_to_lisa" and "lisa_position_ssb" should be
# more general for other detectors in the near future.


def rotation_matrix_ssb_to_lisa(alpha):
    """ The rotation matrix (of frame basis) from SSB frame to LISA frame.
    This function assumes the angle between LISA plane and the ecliptic
    is 60 degrees, and the period of LISA's self-rotation and orbital
    revolution is both one year.

    Parameters
    ----------
    alpha : float
        The angular displacement of LISA in SSB frame.
        In the unit of 'radian'.

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


def lisa_position_ssb(t_lisa, t0=TIME_OFFSET):
    """ Calculating the position vector and angular displacement of LISA
    in the SSB frame, at a given time. This function assumes LISA's barycenter
    is orbiting around a circular orbit within the ecliptic behind the Earth.
    The period of it is one year.

    Parameters
    ----------
    t_lisa : float
        The time when a GW signal arrives at the origin of LISA frame,
        or any other time you want.
    t0 : float
        The initial time offset of LISA, in the unit of 's',
        default is 7365189.431698299. This makes sure LISA is behind
        the Earth by 19-23 degrees.

    Returns
    -------
    (p, alpha) : tuple
    p : numpy.array
        The position vector of LISA in the SSB frame. In the unit of 'm'.
    alpha : float
        The angular displacement of LISA in the SSB frame.
        In the unit of 'radian'.
    """
    OMEGA_0 = 1.99098659277e-7
    R_ORBIT = au.value
    alpha = np.mod(OMEGA_0 * (t_lisa + t0), 2*np.pi)
    p = np.array([[R_ORBIT * np.cos(alpha)],
                  [R_ORBIT * np.sin(alpha)],
                  [0]], dtype=object)
    return (p, alpha)


def localization_to_propagation_vector(lamda, beta):
    """ Converting the sky localization to the corresponding
    propagation unit vector of a GW signal.

    Parameters
    ----------
    lamda : float
        The right ascension or ecliptic longitude, in the unit of 'radian'.
    beta : float
        The declination or ecliptic latitude, in the unit of 'radian'.

    Returns
    -------
    [[x], [y], [z]] : numpy.array
        The propagation unit vector of that GW signal.
    """
    x = -np.cos(beta) * np.cos(lamda)
    y = -np.cos(beta) * np.sin(lamda)
    z = -np.sin(beta)

    return np.array([[x], [y], [z]])


def propagation_vector_to_localization(k):
    """ Converting the propagation unit vector to the corresponding
    sky localization of a GW signal.

    Parameters
    ----------
    k : numpy.array
        The propagation unit vector of a GW signal.

    Returns
    -------
    (lamda, beta) : tuple
        The sky localization of that GW signal.
    """
    # beta already within [-pi/2, pi/2]
    beta = np.float64(np.arcsin(-k[2]))
    lamda = np.float64(np.arctan2(-k[1]/np.cos(beta), -k[0]/np.cos(beta)))
    # lamda should within [0, 2*pi)
    lamda = np.mod(lamda, 2*np.pi)

    return (lamda, beta)


def polarization_newframe(psi, k, rotation_matrix):
    """ Converting a polarization angle from a frame to a new frame
    by using rotation matrix method.

    Parameters
    ----------
    psi : float
        The polarization angle in the old frame, in the unit of 'radian'.
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
    lamda, _ = propagation_vector_to_localization(k)
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
    psi_newframe = np.mod(psi_newframe, 2*np.pi)

    return psi_newframe


def t_lisa_from_ssb(t_ssb, lamda_ssb, beta_ssb, t0=TIME_OFFSET):
    """ Calculating the time when a GW signal arrives at the barycenter
    of LISA, by using the time and sky localization in SSB frame.

    Parameters
    ----------
    t_ssb : float
        The time when a GW signal arrives at the origin of SSB frame.
        In the unit of 's'.
    lamda_ssb : float
        The ecliptic longitude of a GW signal in SSB frame.
        In the unit of 'radian'.
    beta_ssb : float
        The ecliptic latitude of a GW signal in SSB frame.
        In the unit of 'radian'.
    t0 : float
        The initial time offset of LISA, in the unit of 's',
        default is 7365189.431698299. This makes sure LISA is behind
        the Earth by 19-23 degrees.

    Returns
    -------
    t_lisa : float
        The time when a GW signal arrives at the origin of LISA frame.
    """
    k = localization_to_propagation_vector(lamda_ssb, beta_ssb)

    def equation(t_lisa):
        # LISA is moving, when GW arrives at LISA center,
        # time is t_lisa, not t_ssb.
        p = lisa_position_ssb(t_lisa, t0)[0]
        return t_lisa - t_ssb - np.vdot(k, p) / c.value

    return fsolve(equation, t_ssb)[0]


def t_ssb_from_t_lisa(t_lisa, lamda_ssb, beta_ssb, t0=TIME_OFFSET):
    """ Calculating the time when a GW signal arrives at the barycenter
    of SSB, by using the time in LISA frame and sky localization in SSB frame.

    Parameters
    ----------
    t_lisa : float
        The time when a GW signal arrives at the origin of LISA frame.
        In the unit of 's'.
    lamda_ssb : float
        The ecliptic longitude of a GW signal in SSB frame.
        In the unit of 'radian'.
    beta_ssb : float
        The ecliptic latitude of a GW signal in SSB frame.
        In the unit of 'radian'.
    t0 : float
        The initial time offset of LISA, in the unit of 's',
        default is 7365189.431698299. This makes sure LISA is behind
        the Earth by 19-23 degrees.

    Returns
    -------
    t_ssb : float
        The time when a GW signal arrives at the origin of SSB frame.
    """
    k = localization_to_propagation_vector(lamda_ssb, beta_ssb)
    # LISA is moving, when GW arrives at LISA center,
    # time is t_lisa, not t_ssb.
    p = lisa_position_ssb(t_lisa, t0)[0]

    def equation(t_ssb):
        return t_lisa - t_ssb - np.vdot(k, p) / c.value

    return fsolve(equation, t_lisa)[0]


def ssb_to_lisa(t_ssb, lamda_ssb, beta_ssb, psi_ssb, t0=TIME_OFFSET):
    """ Converting the arrive time, the sky localization, and the polarization
    from the SSB frame to the LISA frame.

    Parameters
    ----------
    t_ssb : float
        The time when a GW signal arrives at the origin of SSB frame.
        In the unit of 's'.
    lamda_ssb : float
        The ecliptic longitude of a GW signal in SSB frame.
        In the unit of 'radian'.
    beta_ssb : float
        The ecliptic latitude of a GW signal in SSB frame.
        In the unit of 'radian'.
    psi_ssb : float
        The polarization angle of a GW signal in SSB frame.
        In the unit of 'radian'.
    t0 : float
        The initial time offset of LISA, in the unit of 's',
        default is 7365189.431698299. This makes sure LISA is behind
        the Earth by 19-23 degrees.

    Returns
    -------
    (t_lisa, lamda_lisa, beta_lisa, psi_lisa) : tuple
    t_lisa : float
        The time when a GW signal arrives at the origin of LISA frame.
        In the unit of 's'.
    lamda_lisa : float
        The longitude of a GW signal in LISA frame, in the unit of 'radian'.
    beta_lisa : float
        The latitude of a GW signal in LISA frame, in the unit of 'radian'.
    psi_lisa : float
        The polarization angle of a GW signal in LISA frame.
        In the unit of 'radian'.
    """
    if isinstance(t_ssb, np.ndarray):
        t_ssb = t_ssb.item()
    if isinstance(lamda_ssb, np.ndarray):
        lamda_ssb = lamda_ssb.item()
    if isinstance(beta_ssb, np.ndarray):
        beta_ssb = beta_ssb.item()
    if isinstance(psi_ssb, np.ndarray):
        psi_ssb = psi_ssb.item()
    if lamda_ssb < 0 or lamda_ssb >= 2*np.pi:
        raise ValueError("Longitude should within [0, 2*pi).")
    if beta_ssb < -np.pi/2 or beta_ssb > np.pi/2:
        raise ValueError("Latitude should within [-pi/2, pi/2].")
    if psi_ssb < 0 or psi_ssb >= 2*np.pi:
        raise ValueError("Polarization angle should within [0, 2*pi).")
    t_lisa = t_lisa_from_ssb(t_ssb, lamda_ssb, beta_ssb, t0)
    k_ssb = localization_to_propagation_vector(lamda_ssb, beta_ssb)
    # Although t_lisa calculated above using the corrected LISA position vector by
    # adding t0, it corresponds to the true t_ssb, not t_ssb+t0,
    # we need to include t0 again to correct LISA position.
    alpha = lisa_position_ssb(t_lisa, t0)[1]
    rotation_matrix_lisa = rotation_matrix_ssb_to_lisa(alpha)
    k_lisa = rotation_matrix_lisa.T @ k_ssb
    lamda_lisa, beta_lisa = propagation_vector_to_localization(k_lisa)
    psi_lisa = polarization_newframe(psi_ssb, k_ssb, rotation_matrix_lisa)

    return (t_lisa, lamda_lisa, beta_lisa, psi_lisa)


def lisa_to_ssb(t_lisa, lamda_lisa, beta_lisa, psi_lisa, t0=TIME_OFFSET):
    """ Converting the arrive time, the sky localization, and the polarization
    from the LISA frame to the SSB frame.

    Parameters
    ----------
    t_lisa : float
        The time when a GW signal arrives at the origin of LISA frame.
        In the unit of 's'.
    lamda_lisa : float
        The longitude of a GW signal in LISA frame, in the unit of 'radian'.
    beta_lisa : float
        The latitude of a GW signal in LISA frame, in the unit of 'radian'.
    psi_lisa : float
        The polarization angle of a GW signal in LISA frame.
        In the unit of 'radian'.
    t0 : float
        The initial time offset of LISA, in the unit of 's',
        default is 7365189.431698299. This makes sure LISA is behind
        the Earth by 19-23 degrees.

    Returns
    -------
    (t_ssb, lamda_ssb, beta_ssb, psi_ssb) : tuple
    t_ssb : float
        The time when a GW signal arrives at the origin of SSB frame.
        In the unit of 's'.
    lamda_ssb : float
        The ecliptic longitude of a GW signal in SSB frame.
        In the unit of 'radian'.
    beta_ssb : float
        The ecliptic latitude of a GW signal in SSB frame.
        In the unit of 'radian'.
    psi_ssb : float
        The polarization angle of a GW signal in SSB frame.
        In the unit of 'radian'.
    """
    if isinstance(t_lisa, np.ndarray):
        t_lisa = t_lisa.item()
    if isinstance(lamda_lisa, np.ndarray):
        lamda_lisa = lamda_lisa.item()
    if isinstance(beta_lisa, np.ndarray):
        beta_lisa = beta_lisa.item()
    if isinstance(psi_lisa, np.ndarray):
        psi_lisa = psi_lisa.item()
    if lamda_lisa < 0 or lamda_lisa >= 2*np.pi:
        raise ValueError("Longitude should within [0, 2*pi).")
    if beta_lisa < -np.pi/2 or beta_lisa > np.pi/2:
        raise ValueError("Latitude should within [-pi/2, pi/2].")
    if psi_lisa < 0 or psi_lisa >= 2*np.pi:
        raise ValueError("Polarization angle should within [0, 2*pi).")
    k_lisa = localization_to_propagation_vector(lamda_lisa, beta_lisa)
    alpha = lisa_position_ssb(t_lisa, t0)[1]
    rotation_matrix_lisa = rotation_matrix_ssb_to_lisa(alpha)
    k_ssb = rotation_matrix_lisa @ k_lisa
    lamda_ssb, beta_ssb = propagation_vector_to_localization(k_ssb)
    t_ssb = t_ssb_from_t_lisa(t_lisa, lamda_ssb, beta_ssb, t0)
    psi_ssb = polarization_newframe(psi_lisa, k_lisa, rotation_matrix_lisa.T)

    return (t_ssb, lamda_ssb, beta_ssb, psi_ssb)


def rotation_matrix_ssb_to_geo(epsilon=np.deg2rad(23.439281)):
    """ The rotation matrix (of frame basis) from SSB frame to
    geocentric frame.

    Parameters
    ----------
    epsilon : float
        The Earth's axial tilt (obliquity), in the unit of 'radian'.

    Returns
    -------
    r : numpy.array
        A 3x3 rotation matrix from SSB frame to geocentric frame.
    """
    r = Rotation.from_rotvec([
        [-epsilon, 0, 0]
    ]).as_matrix()

    return np.array(r[0])


def earth_position_ssb(t_geo):
    """ Calculating the position vector and angular displacement of the Earth
    in the SSB frame, at a given time. By using Astropy.

    Parameters
    ----------
    t_geo : float
        The time when a GW signal arrives at the origin of geocentric frame,
        or any other time you want.

    Returns
    -------
    (p, alpha) : tuple
    p : numpy.array
        The position vector of the Earth in the SSB frame. In the unit of 'm'.
    alpha : float
        The angular displacement of the Earth in the SSB frame.
        In the unit of 'radian'.
    """
    t = Time(t_geo, format='gps')
    pos = get_body_barycentric('earth', t)
    icrs_coord = SkyCoord(pos.x, pos.y, pos.z, frame='icrs',
                          representation_type='cartesian')
    bme_coord = icrs_coord.transform_to(
                    BarycentricMeanEcliptic(equinox='J2000'))
    x = bme_coord.cartesian.x.to(units.m).value
    y = bme_coord.cartesian.y.to(units.m).value
    z = bme_coord.cartesian.z.to(units.m).value
    p = np.array([[x], [y], [z]])
    alpha = bme_coord.lon.rad

    return (p, alpha)


def t_geo_from_ssb(t_ssb, lamda_ssb, beta_ssb):
    """ Calculating the time when a GW signal arrives at the barycenter
    of the Earth, by using the time and sky localization in SSB frame.

    Parameters
    ----------
    t_ssb : float
        The time when a GW signal arrives at the origin of SSB frame.
        In the unit of 's'.
    lamda_ssb : float
        The ecliptic longitude of a GW signal in SSB frame.
        In the unit of 'radian'.
    beta_ssb : float
        The ecliptic latitude of a GW signal in SSB frame.
        In the unit of 'radian'.

    Returns
    -------
    t_geo : float
        The time when a GW signal arrives at the origin of geocentric frame.
    """
    k = localization_to_propagation_vector(lamda_ssb, beta_ssb)

    def equation(t_geo):
        # Earth is moving, when GW arrives at Earth center,
        # time is t_geo, not t_ssb.
        p = earth_position_ssb(t_geo)[0]
        return t_geo - t_ssb - np.vdot(k, p) / c.value

    return fsolve(equation, t_ssb)[0]


def t_ssb_from_t_geo(t_geo, lamda_ssb, beta_ssb):
    """ Calculating the time when a GW signal arrives at the barycenter
    of SSB, by using the time in geocentric frame and sky localization
    in SSB frame.

    Parameters
    ----------
    t_geo : float
        The time when a GW signal arrives at the origin of geocentric frame.
        In the unit of 's'.
    lamda_ssb : float
        The ecliptic longitude of a GW signal in SSB frame.
        In the unit of 'radian'.
    beta_ssb : float
        The ecliptic latitude of a GW signal in SSB frame.
        In the unit of 'radian'.

    Returns
    -------
    t_ssb : float
        The time when a GW signal arrives at the origin of SSB frame.
    """
    k = localization_to_propagation_vector(lamda_ssb, beta_ssb)
    # Earth is moving, when GW arrives at Earth center,
    # time is t_geo, not t_ssb.
    p = earth_position_ssb(t_geo)[0]

    def equation(t_ssb):
        return t_geo - t_ssb - np.vdot(k, p) / c.value

    return fsolve(equation, t_geo)[0]


def ssb_to_geo(t_ssb, lamda_ssb, beta_ssb, psi_ssb, use_astropy=True):
    """ Converting the arrive time, the sky localization, and the polarization
    from the SSB frame to the geocentric frame.

    Parameters
    ----------
    t_ssb : float
        The time when a GW signal arrives at the origin of SSB frame.
        In the unit of 's'.
    lamda_ssb : float
        The ecliptic longitude of a GW signal in SSB frame.
        In the unit of 'radian'.
    beta_ssb : float
        The ecliptic latitude of a GW signal in SSB frame.
        In the unit of 'radian'.
    psi_ssb : float
        The polarization angle of a GW signal in SSB frame.
        In the unit of 'radian'.
    use_astropy : bool
        Using Astropy to calculate the sky localization or not.
        Default is True.

    Returns
    -------
    (t_geo, lamda_geo, beta_geo, psi_geo) : tuple
    t_geo : float
        The time when a GW signal arrives at the origin of geocentric frame.
        In the unit of 's'.
    lamda_geo : float
        The longitude of a GW signal in geocentric frame.
        In the unit of 'radian'.
    beta_geo : float
        The latitude of a GW signal in geocentric frame.
        In the unit of 'radian'.
    psi_geo : float
        The polarization angle of a GW signal in geocentric frame.
        In the unit of 'radian'.
    """
    if isinstance(t_ssb, np.ndarray):
        t_ssb = t_ssb.item()
    if isinstance(lamda_ssb, np.ndarray):
        lamda_ssb = lamda_ssb.item()
    if isinstance(beta_ssb, np.ndarray):
        beta_ssb = beta_ssb.item()
    if isinstance(psi_ssb, np.ndarray):
        psi_ssb = psi_ssb.item()
    if lamda_ssb < 0 or lamda_ssb >= 2*np.pi:
        raise ValueError("Longitude should within [0, 2*pi).")
    if beta_ssb < -np.pi/2 or beta_ssb > np.pi/2:
        raise ValueError("Latitude should within [-pi/2, pi/2].")
    if psi_ssb < 0 or psi_ssb >= 2*np.pi:
        raise ValueError("Polarization angle should within [0, 2*pi).")
    t_geo = t_geo_from_ssb(t_ssb, lamda_ssb, beta_ssb)
    k_ssb = localization_to_propagation_vector(lamda_ssb, beta_ssb)
    rotation_matrix_geo = rotation_matrix_ssb_to_geo()
    if use_astropy:
        bme_coord = BarycentricMeanEcliptic(lon=lamda_ssb*units.radian,
                                            lat=beta_ssb*units.radian,
                                            equinox='J2000')
        geo_sky = bme_coord.transform_to(PrecessedGeocentric(equinox='J2000'))
        lamda_geo = geo_sky.ra.rad
        beta_geo = geo_sky.dec.rad
    else:
        k_geo = rotation_matrix_geo.T @ k_ssb
        lamda_geo, beta_geo = propagation_vector_to_localization(k_geo)
    psi_geo = polarization_newframe(psi_ssb, k_ssb, rotation_matrix_geo)

    return (t_geo, lamda_geo, beta_geo, psi_geo)


def geo_to_ssb(t_geo, lamda_geo, beta_geo, psi_geo, use_astropy=True):
    """ Converting the arrive time, the sky localization, and the polarization
    from the geocentric frame to the SSB frame.

    Parameters
    ----------
    t_geo : float
        The time when a GW signal arrives at the origin of geocentric frame.
        In the unit of 's'.
    lamda_geo : float
        The longitude of a GW signal in geocentric frame.
        In the unit of 'radian'.
    beta_geo : float
        The latitude of a GW signal in geocentric frame.
        In the unit of 'radian'.
    psi_geo : float
        The polarization angle of a GW signal in geocentric frame.
        In the unit of 'radian'.
    use_astropy : bool
        Using Astropy to calculate the sky localization or not.
        Default is True.

    Returns
    -------
    (t_ssb, lamda_ssb, beta_ssb, psi_ssb) : tuple
    t_ssb : float
        The time when a GW signal arrives at the origin of SSB frame.
        In the unit of 's'.
    lamda_ssb : float
        The ecliptic longitude of a GW signal in SSB frame.
        In the unit of 'radian'.
    beta_ssb : float
        The ecliptic latitude of a GW signal in SSB frame.
        In the unit of 'radian'.
    psi_ssb : float
        The polarization angle of a GW signal in SSB frame.
        In the unit of 'radian'.
    """
    if isinstance(t_geo, np.ndarray):
        t_geo = t_geo.item()
    if isinstance(lamda_geo, np.ndarray):
        lamda_geo = lamda_geo.item()
    if isinstance(beta_geo, np.ndarray):
        beta_geo = beta_geo.item()
    if isinstance(psi_geo, np.ndarray):
        psi_geo = psi_geo.item()
    if lamda_geo < 0 or lamda_geo >= 2*np.pi:
        raise ValueError("Longitude should within [0, 2*pi).")
    if beta_geo < -np.pi/2 or beta_geo > np.pi/2:
        raise ValueError("Latitude should within [-pi/2, pi/2].")
    if psi_geo < 0 or psi_geo >= 2*np.pi:
        raise ValueError("Polarization angle should within [0, 2*pi).")
    k_geo = localization_to_propagation_vector(lamda_geo, beta_geo)
    rotation_matrix_geo = rotation_matrix_ssb_to_geo()
    k_ssb = rotation_matrix_geo @ k_geo
    lamda_ssb, beta_ssb = propagation_vector_to_localization(k_ssb)
    psi_ssb = polarization_newframe(psi_geo, k_geo, rotation_matrix_geo.T)
    if use_astropy:
        bme_coord = PrecessedGeocentric(ra=lamda_geo*units.radian,
                                        dec=beta_geo*units.radian,
                                        equinox='J2000')
        ssb_sky = bme_coord.transform_to(
                    BarycentricMeanEcliptic(equinox='J2000'))
        lamda_ssb = ssb_sky.lon.rad
        beta_ssb = ssb_sky.lat.rad
    t_ssb = t_ssb_from_t_geo(t_geo, lamda_ssb, beta_ssb)

    return (t_ssb, lamda_ssb, beta_ssb, psi_ssb)


def lisa_to_geo(t_lisa, lamda_lisa, beta_lisa, psi_lisa,
                t0=TIME_OFFSET, use_astropy=True):
    """ Converting the arrive time, the sky localization, and the polarization
    from the LISA frame to the geocentric frame.

    Parameters
    ----------
    t_lisa : float
        The time when a GW signal arrives at the origin of LISA frame.
        In the unit of 's'.
    lamda_lisa : float
        The longitude of a GW signal in LISA frame, in the unit of 'radian'.
    beta_lisa : float
        The latitude of a GW signal in LISA frame, in the unit of 'radian'.
    psi_lisa : float
        The polarization angle of a GW signal in LISA frame.
        In the unit of 'radian'.
    t0 : float
        The initial time offset of LISA, in the unit of 's',
        default is 7365189.431698299. This makes sure LISA is behind
        the Earth by 19-23 degrees.
    use_astropy : bool
        Using Astropy to calculate the sky localization or not.
        Default is True.

    Returns
    -------
    (t_geo, lamda_geo, beta_geo, psi_geo) : tuple
    t_geo : float
        The time when a GW signal arrives at the origin of geocentric frame.
        In the unit of 's'.
    lamda_geo : float
        The ecliptic longitude of a GW signal in geocentric frame.
        In the unit of 'radian'.
    beta_geo : float
        The ecliptic latitude of a GW signal in geocentric frame.
        In the unit of 'radian'.
    psi_geo : float
        The polarization angle of a GW signal in geocentric frame.
        In the unit of 'radian'.
    """
    t_ssb, lamda_ssb, beta_ssb, psi_ssb = lisa_to_ssb(
        t_lisa, lamda_lisa, beta_lisa, psi_lisa, t0)
    t_geo, lamda_geo, beta_geo, psi_geo = ssb_to_geo(
        t_ssb, lamda_ssb, beta_ssb, psi_ssb, use_astropy)

    return (t_geo, lamda_geo, beta_geo, psi_geo)


def geo_to_lisa(t_geo, lamda_geo, beta_geo, psi_geo,
                t0=TIME_OFFSET, use_astropy=True):
    """ Converting the arrive time, the sky localization, and the polarization
    from the geocentric frame to the LISA frame.

    Parameters
    ----------
    t_geo : float
        The time when a GW signal arrives at the origin of geocentric frame.
        In the unit of 's'.
    lamda_geo : float
        The longitude of a GW signal in geocentric frame.
        In the unit of 'radian'.
    beta_geo : float
        The latitude of a GW signal in geocentric frame.
        In the unit of 'radian'.
    psi_geo : float
        The polarization angle of a GW signal in geocentric frame.
        In the unit of 'radian'.
    t0 : float
        The initial time offset of LISA, in the unit of 's',
        default is 7365189.431698299. This makes sure LISA is behind
        the Earth by 19-23 degrees.
    use_astropy : bool
        Using Astropy to calculate the sky localization or not.
        Default is True.

    Returns
    -------
    (t_lisa, lamda_lisa, beta_lisa, psi_lisa) : tuple
    t_lisa : float
        The time when a GW signal arrives at the origin of LISA frame.
        In the unit of 's'.
    lamda_lisa : float
        The longitude of a GW signal in LISA frame, in the unit of 'radian'.
    beta_lisa : float
        The latitude of a GW signal in LISA frame, in the unit of 'radian'.
    psi_geo : float
        The polarization angle of a GW signal in LISA frame.
        In the unit of 'radian'.
    """
    t_ssb, lamda_ssb, beta_ssb, psi_ssb = geo_to_ssb(
        t_geo, lamda_geo, beta_geo, psi_geo, use_astropy)
    t_lisa, lamda_lisa, beta_lisa, psi_lisa = ssb_to_lisa(
        t_ssb, lamda_ssb, beta_ssb, psi_ssb, t0)

    return (t_lisa, lamda_lisa, beta_lisa, psi_lisa)
