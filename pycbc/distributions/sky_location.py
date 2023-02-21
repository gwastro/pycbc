# Copyright (C) 2016  Collin Capano
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
"""This modules provides classes for evaluating sky distributions in
right acension and declination.
"""


import numpy
from pycbc.distributions import angular
from pycbc.transforms import new_z_to_euler, rotate_euler
from pycbc.transforms import decra2polaz, polaz2radec

class UniformSky(angular.UniformSolidAngle):
    """A distribution that is uniform on the sky. This is the same as
    UniformSolidAngle, except that the polar angle varies from pi/2 (the north
    pole) to -pi/2 (the south pole) instead of 0 to pi. Also, the default
    names are "dec" (declination) for the polar angle and "ra" (right
    ascension) for the azimuthal angle, instead of "theta" and "phi".
    """
    name = 'uniform_sky'
    _polardistcls = angular.CosAngle
    _default_polar_angle = 'dec'
    _default_azimuthal_angle = 'ra'

class FisherDist():
    """A distribution that returns a random (ra, dec) angle drawn from the
    Fisher distribution. Assume that the concentration parameter (kappa)
    is large so that we can use a Rayleigh distribution about the north
    pole and rotate it to be centered at the (ra, dec) coordinate mu.
    Assume kappa = 1 / sigma**2
    As in UniformSky, the declination (dec) varies from pi/2 to-pi/2
    and right ascension (ra) varies from 0 to 2pi. And the angles
    should be provided in (ra,dec) format in radians (mu_radians=True),
    rather than factors of pi, or in degrees (mu_radians=False).
    References:
      * http://en.wikipedia.org/wiki/Von_Mises-Fisher_distribution
      * http://arxiv.org/pdf/0902.0737v1 (states the Rayleigh limit)
    """
    name = 'fisher_dist'
    _params=['ra','dec']

    def __init__(self,ra,dec):
        self.ra = ra
        self.dec = dec
        print('params ok')
    
    @property
    def params(self):
        return self._params
    
    @classmethod
    def from_config(cls,cp, section,variable_args):
        tag = variable_args
        variable_args = variable_args.split(VARARGS_DELIM)
        if not set(variable_args) == set(cls._params):
            raise ValueError("Not all parameters used by this distribution "
                             "included in tag portion of section name")
        ra = get_param_bounds_from_config(cp, section, tag, 'ra')
        dec = get_param_bounds_from_config(cp, section, tag, 'dec')
        print('In from_config()')
        return cls(ra=ra, dec=dec)

    def rvs(self,size):
        mu_values =numpy.array([self.ra,self.dec])
        print(mu_values)
        kappa=900
        size=size
        arr=numpy.array([
            numpy.random.rayleigh(scale=1./numpy.sqrt(kappa),
                                  size=size),
            numpy.random.uniform(low=0,
                                 high=(2*numpy.pi),
                                 size=size)]).reshape((2, size)).T
        alpha, beta = new_z_to_euler(mu_values)
        print(alpha, beta)
        print('rvs ok')
        print(self.ra, self.dec)
        return rotate_euler(arr, alpha, beta, 0)


__all__ = ['UniformSky', 'Fisher']
