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
from pycbc import VARARGS_DELIM
from pycbc.io import FieldArray

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

class FisherSky():
    """A distribution that returns a random (ra, dec) angle drawn from the
    Fisher distribution. Assume that the concentration parameter (kappa)
    is large so that we can use a Rayleigh distribution about the north
    pole and rotate it to be centered at the (ra, dec) coordinate mu.
    Assume kappa = 1 / (0.66*sigma)**2
    As in UniformSky, the declination (dec) varies from pi/2 to-pi/2
    and right ascension (ra) varies from 0 to 2pi. The angles
    should be provided in (ra,dec) format in radians (the variable
    angle_unit='rad' by default), rather than factors of pi, or
    in degrees (in this case set angle_unit='deg').
    References:
      * http://en.wikipedia.org/wiki/Von_Mises-Fisher_distribution
      * http://arxiv.org/pdf/0902.0737v1 (states the Rayleigh limit)
    """
    name = 'fisher_sky'
    _params = ['ra', 'dec']

    def __init__(self, angle_unit='rad', **params):
        self.sigma = params['sigma']
        if angle_unit == 'rad':
            self.mu_values = numpy.array([params['mean_ra'],
                                          params['mean_dec']])
            self.kappa = 1./((0.66*params['sigma']))**2
        elif angle_unit == 'deg':
            self.mu_values = numpy.deg2rad([params['mean_ra'],
                                            params['mean_dec']])
            self.kappa = 1./((0.66*numpy.deg2rad(params['sigma'])))**2
        else:
            raise ValueError("Only deg or rad is allowed as unit")
        self.alpha, self.beta = new_z_to_euler(self.mu_values)

    @property
    def params(self):
        return self._params

    @classmethod
    def from_config(cls, cp, section, variable_args):
        tag = variable_args
        variable_args = variable_args.split(VARARGS_DELIM)
        if not set(variable_args) == set(cls._params):
            raise ValueError("Not all parameters used by this distribution "
                             "included in tag portion of section name")
        mean_ra = float(cp.get_opt_tag(section, 'mean_ra', tag))
        mean_dec = float(cp.get_opt_tag(section, 'mean_dec', tag))
        sigma = float(cp.get_opt_tag(section, 'sigma', tag))
        angle_unit = cp.get_opt_tag(section, 'angle_unit', tag)
        return cls(mean_ra=mean_ra, mean_dec=mean_dec, sigma=sigma,
                   angle_unit=angle_unit)

    def rvs(self, size):
        arr = numpy.array([
            numpy.random.rayleigh(scale=1./numpy.sqrt(self.kappa),
                                  size=size),
            numpy.random.uniform(low=0,
                                 high=(2*numpy.pi),
                                 size=size)]).T
        euler = rotate_euler(arr, self.alpha, self.beta, 0)
        rot_euler = FieldArray(size, dtype=[('ra', '<f8'), ('dec', '<f8')])
        rot_euler['ra'], rot_euler['dec'] = euler[:, 0], euler[:, 1]
        return rot_euler

__all__ = ['UniformSky', 'FisherSky']
