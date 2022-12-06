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


from pycbc.distributions import angular
import numpy
from pycbc.transforms import new_z_to_euler, rotate_euler

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
    """A distribution that returns a random (ra, dec) angle drawn from the Fisher
    distribution. Assume that the concentration parameter (kappa) is large
    so that we can use a Rayleigh distribution about the north pole and
    rotate it to be centered at the (ra, dec) coordinate mu.

    Assume kappa = 1 / sigma**2
   
    As in UniformSky, the polar angle varies from pi/2 (the north pole) to
    -pi/2 (the south pole) and dec varies from 0 to 2pi.

    References:
      * http://en.wikipedia.org/wiki/Von_Mises-Fisher_distribution
      * http://arxiv.org/pdf/0902.0737v1 (states the Rayleigh limit)
    """
    name = 'fisher_dist'
    def __init__(self, mu, kappa, size, if_radians=True):
        self.kappa = kappa
        if if_radians is True:
            self.mu= (mu[1],mu[0])
        elif if_radians is False :
            self.mu=numpy.array(numpy.deg2rad([mu[1],mu[0]]))
        else:
            raise ValueError("You can choose to not give any option if angles are in radians or you should give either 'True' or 'False' ")  
            
    def rvs(self,size):
        arr=numpy.array([numpy.random.rayleigh(scale=1. / numpy.sqrt(self.kappa), size=size),
                         numpy.random.uniform(low=0, high=2*numpy.pi, size=size)]).reshape((2, size)).T
        a, b = new_z_to_euler(self.mu)
        return rotate_euler(arr, b, ((numpy.pi)/2) -a, 0)

__all__ = ['UniformSky', 'FisherDist']
