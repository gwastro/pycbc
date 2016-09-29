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
    if type(y) is int:
        y = float(y)
    rho = numpy.sqrt(x**2 + y**2 + z**2)
    phi = numpy.arctan(y / x)
    theta = numpy.arccos(z / rho)
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
