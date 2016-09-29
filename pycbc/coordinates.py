# Copyright (C) 2016 Christopher M. Biwer
#
#
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
""" Coordinate transformations.
"""
import numpy


def spherical_rho(x, y, z):
    """ Calculates rho from Cartesian coordinates.

    Parameters
    ----------
    x : {numpy.array, float}
        X-coordinate.
    y : {numpy.array, float}
        Y-coordinate.
    z : {numpy.array, float}
        Z-coordinate.

    Returns
    -------
    rho : {numpy.array, float}
        The radial amplitude.
    """
    return numpy.sqrt(x**2 + y**2 + z**2)


def spherical_phi(x, y):
    """ Calculates phi from Cartesian coordinates. phi is in [0,2*pi].

    Parameters
    ----------
    x : {numpy.array, float}
        X-coordinate.
    y : {numpy.array, float}
        Y-coordinate.

    Returns
    -------
    phi : {numpy.array, float}
        The azmuthal angle.
    """
    if type(y) is int:
        y = float(y)
    return numpy.arctan(y / x)


def spherical_theta(x, y, z):
    """ Calculates theta from Cartesian coordinates. theta is in [0,pi].

    Parameters
    ----------
    x : {numpy.array, float}
        X-coordinate.
    y : {numpy.array, float}
        Y-coordinate.
    z : {numpy.array, float}
        Z-coordinate.

    Returns
    -------
    theta : {numpy.array, float}
        The polar angle.
    """
    return numpy.arccos(z / spherical_rho(x, y, z))


def cartesian_to_spherical(x, y, z):
    """ Maps cartesian coordinates (x,y,z) to spherical coordinates
    (rho,phi,theta) where phi is in [0,2*pi] and theta is in [0,pi].

    Parameters
    ----------
    x : {numpy.array, float}
        X-coordinate.
    y : {numpy.array, float}
        Y-coordinate.
    z : {numpy.array, float}
        Z-coordinate.

    Returns
    -------
    rho : {numpy.array, float}
        The radial amplitude.
    phi : {numpy.array, float}
        The azimuthal angle.
    theta : {numpy.array, float}
        The polar angle.
    """
    rho = spherical_rho(x, y, z)
    phi = spherical_phi(x, y)
    theta = spherical_theta(x, y, z)
    return rho, phi, theta


def spherical_to_cartesian(rho, phi, theta):
    """ Maps spherical coordinates (rho,phi,theta) to cartesian coordinates
    (x,y,z) where phi is in [0,2*pi] and theta is in [0,pi].

    Parameters
    ----------
    rho : {numpy.array, float}
        The radial amplitude.
    phi : {numpy.array, float}
        The azimuthal angle.
    theta : {numpy.array, float}
        The polar angle.

    Returns
    -------
    x : {numpy.array, float}
        X-coordinate.
    y : {numpy.array, float}
        Y-coordinate.
    z : {numpy.array, float}
        Z-coordinate.
    """
    x = rho * numpy.cos(phi) * numpy.sin(theta)
    y = rho * numpy.sin(phi) * numpy.sin(theta)
    z = rho * numpy.cos(theta)
    return x, y, z
