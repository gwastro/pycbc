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
right ascension and declination.
"""

import logging
import numpy


from scipy.spatial.transform import Rotation

from pycbc.distributions import angular
from pycbc import VARARGS_DELIM
from pycbc.io import FieldArray
from pycbc.types import angle_as_radians

logger = logging.getLogger('pycbc.distributions.sky_location')


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


class FisherSky:
    """A distribution that returns a random angle drawn from an approximate
    `Von_Mises-Fisher distribution`_. Assumes that the Fisher concentration
    parameter is large, so that we can draw the samples from a simple
    rotationally-invariant distribution centered at the North Pole (which
    factors as a uniform distribution for the right ascension, and a Rayleigh
    distribution for the declination, as described in
    `Fabrycky and Winn 2009 ApJ 696 1230`) and then rotate the samples to be
    centered around the specified mean position. As in UniformSky, the
    declination varies from π/2 to -π/2 and the right ascension varies from
    0 to 2π.

    .. _Von_Mises-Fisher distribution:
        http://en.wikipedia.org/wiki/Von_Mises-Fisher_distribution

    .. _Fabrycky and Winn 2009 ApJ 696 1230:
        https://doi.org/10.1088/0004-637X/696/2/1230

    .. _Briggs et al 1999 ApJS 122 503:
        https://doi.org/10.1086/313221

    Parameters
    ----------
    mean_ra: float or str
        RA of the center of the distribution. Use the rad or deg suffix to
        specify units, otherwise radians are assumed.
    mean_dec: float or str
        Declination of the center of the distribution. Use the rad or deg 
        suffix to specify units, otherwise radians are assumed.
    sigma: float or str
        Spread of the distribution. For the precise interpretation, see Eq 8
        of `Briggs et al 1999 ApJS 122 503`_. This should be smaller than
        about 20 deg for the approximation to be valid. Use the rad or deg 
        suffix to specify units, otherwise radians are assumed.

    """

    name = 'fisher_sky'
    _params = ['ra', 'dec']

    def __init__(self, **params):
        mean_ra = angle_as_radians(params['mean_ra'])
        mean_dec = angle_as_radians(params['mean_dec'])
        sigma = angle_as_radians(params['sigma'])
        if mean_ra < 0 or mean_ra > 2 * numpy.pi:
            raise ValueError(
                f'The mean RA must be between 0 and 2π, {mean_ra} rad given'
            )
        if mean_dec < -numpy.pi / 2 or mean_dec > numpy.pi / 2:
            raise ValueError(
                'The mean declination must be between '
                f'-π/2 and π/2, {mean_dec} rad given'
            )
        if sigma <= 0 or sigma > 2 * numpy.pi:
            raise ValueError(
                'Sigma must be positive and smaller than 2π '
                '(preferably much smaller)'
            )
        if sigma > 0.35:
            logger.warning(
                'Warning: sigma = %s rad is probably too large for the '
                'Fisher approximation to be valid',
                sigma,
            )
        self.rayleigh_scale = 0.66 * sigma
        # Prepare a rotation that puts the North Pole at the mean position
        self.rotation = Rotation.from_euler(
            'yz', [numpy.pi / 2 - mean_dec, mean_ra]
        )

    @property
    def params(self):
        return self._params

    @classmethod
    def from_config(cls, cp, section, variable_args):
        tag = variable_args
        variable_args = variable_args.split(VARARGS_DELIM)
        if set(variable_args) != set(cls._params):
            raise ValueError(
                "Not all parameters used by this distribution "
                "included in tag portion of section name"
            )
        mean_ra = cp.get_opt_tag(section, 'mean_ra', tag)
        mean_dec = cp.get_opt_tag(section, 'mean_dec', tag)
        sigma = cp.get_opt_tag(section, 'sigma', tag)
        return cls(
            mean_ra=mean_ra,
            mean_dec=mean_dec,
            sigma=sigma,
        )

    def rvs(self, size):
        # Draw samples from a distribution centered on the North pole
        np_ra = numpy.random.uniform(low=0, high=(2 * numpy.pi), size=size)
        np_dec = numpy.random.rayleigh(scale=self.rayleigh_scale, size=size)

        # Convert the samples to intermediate cartesian representation
        np_cart = numpy.empty(shape=(size, 3))
        np_cart[:, 0] = numpy.cos(np_ra) * numpy.sin(np_dec)
        np_cart[:, 1] = numpy.sin(np_ra) * numpy.sin(np_dec)
        np_cart[:, 2] = numpy.cos(np_dec)

        # Rotate the samples according to our pre-built rotation
        rot_cart = self.rotation.apply(np_cart)

        # Convert the samples back to spherical coordinates.
        # Some unpleasant conditional operations are needed
        # to get the correct angle convention.
        rot_radec = FieldArray(size, dtype=[('ra', '<f8'), ('dec', '<f8')])
        rot_radec['ra'] = numpy.arctan2(rot_cart[:, 1], rot_cart[:, 0])
        neg_mask = rot_radec['ra'] < 0
        rot_radec['ra'][neg_mask] += 2 * numpy.pi
        rot_radec['dec'] = numpy.arcsin(rot_cart[:, 2])
        return rot_radec


class HealpixSky:
    """Sample the distribution given by a HEALPix map by using the rejection
    sampling method.
    To have an acceptable acceptance rate, a rectangle in (alpha,delta) is
    first found, that encloses a given amount of probability (specified by
    the `coverage` parameter). In order to do this the map is rasterized to a
    resolution specified by `rasterization_nside`.

    The declination (delta) varies from π/2 to -π/2 and the right ascension
    (alpha) varies from 0 to 2π.

    Parameters
    ----------
    healpix_file : str
        Path to a fits file containing probability distribution encoded in a
        HEALPix scheme.
    coverage : float
        Fraction of the map covered by the method. Must be between 0 and 1.
    rasterization_nside : int
        Nside of the rasterized map used to determine the
        boundaries of the input map. Must be a power of 2,
        preferably between 64 and 256.
    """

    name = 'healpix_sky'
    _params = ['ra', 'dec']

    def __init__(self, **params):
        import mhealpy

        def boundaries(healpix_map, nside, coverage):
            """Boundaries of the part of the celestial sphere which we are
            looking to do the distribution.

            Parameters
            ----------
            healpix_map : HealpixMap instance

            nside : int
                nside of the rasterized map used to determine the boundaries of
                the input map.
                must be a power of 2
            coverage : float
                percentage of the map covered by the method
                must be betwen 0 and 1

            Returns
            -------
            delta_min : float
                minimum declination of the map in radians.
            delta_max : float
                maximum declination of the map in radians.
            alpha_min : float
                minimum right ascention of the map in radians.
            alpha_max : float
                maximum right ascention of the map in radians.

            """

            nside = min(nside, healpix_map.nside)

            delta_max = -numpy.pi / 2
            delta_min = numpy.pi / 2
            alpha_max = 0
            alpha_min = 2 * numpy.pi
            # FIXME this only works with 'NESTED', why?
            rasterized_map = healpix_map.rasterize(
                scheme='NESTED', nside=nside
            )
            data = rasterized_map.data
            renormalization_constant = data.sum()

            normalized_data = data / renormalization_constant
            sorter = numpy.argsort(-normalized_data)

            map_coverage = 0

            for i in range(len(data)):
                j = sorter[i]
                pix = rasterized_map.uniq[j] - 4 * nside * nside
                delta, alpha = rasterized_map.pix2ang(pix, lonlat=False)
                # converts colatitude to latitude
                delta = numpy.pi / 2 - delta

                map_coverage += normalized_data[j]
                delta_max, delta_min = max(delta_max, delta), min(
                    delta_min, delta
                )
                alpha_max, alpha_min = max(alpha_max, alpha), min(
                    alpha_min, alpha
                )

                if map_coverage > coverage:
                    break

            # A safety margin is added to ensure that the pixels at the edges
            # of the distribution are fully accounted for.
            # width of one pixel : < π/4nside-1

            margin = 2 * numpy.pi / (4 * nside - 1)

            delta_max = min(delta_max + margin, numpy.pi / 2)
            delta_min = max(delta_min - margin, -numpy.pi / 2)
            alpha_max = min(alpha_max + margin, 2 * numpy.pi)
            alpha_min = max(alpha_min - margin, 0)

            return delta_min, delta_max, alpha_min, alpha_max

        file_name = params['healpix_file']

        coverage = params['coverage']

        if coverage > 1 or coverage < 0:
            raise ValueError(
                f'Coverage must be between 0 and 1, {coverage} is not correct'
            )

        rasterization_nside = params['rasterization_nside']

        if bin(rasterization_nside).count('1') != 1:
            raise ValueError(
                f'Rasterization_nside must be a power of 2,'
                f'{rasterization_nside} is not correct'
            )
        self.healpix_map = mhealpy.HealpixMap.read_map(file_name)
        self.boundaries = boundaries(
            self.healpix_map, rasterization_nside, coverage
        )
        logging.info('HealpixSky boundary is %s', self.boundaries)

    @property
    def params(self):
        return self._params

    @classmethod
    def from_config(cls, cp, section, variable_args):
        tag = variable_args
        variable_args = variable_args.split(VARARGS_DELIM)
        if set(variable_args) != set(cls._params):
            raise ValueError(
                "Not all parameters used by this distribution "
                "included in tag portion of section name"
            )
        healpix_file = str(cp.get_opt_tag(section, 'healpix_file', tag))
        coverage = 0.9999
        if cp.has_option_tag(section, 'coverage', tag):
            coverage = float(cp.get_opt_tag(section, 'coverage', tag))

        rasterization_nside = 64
        if cp.has_option_tag(section, 'rasterization_nside', tag):
            rasterization_nside = int(
                cp.get_opt_tag(section, 'rasterization_nside', tag)
            )
        return cls(
            healpix_file=healpix_file,
            coverage=coverage,
            rasterization_nside=rasterization_nside,
        )

    def rvs(self, size):
        def simple_rejection_sampling(healpix_map, size, boundaries):
            """Start from a uniform distribution of points, and accepts those
             whose values on the map are greater than a random value
             following a uniform law

             Parameters
             ----------
             healpix_map  : HealpixMap instance
             size : int
                 number of points tested by the method
             boundaries : list -> tuple of 4 floats
                 delta_min,delta_max,alpha_min,alpha_max = boundaries
                 delta is the declination in radians [-pi/2, pi/2]
                 alpha is the right ascention in radians [0,2pi]

             Returns
             -------
            coordinates of the accepted points,
            following the mhealpy conventions

            """
            # The angles are in radians and follow the radec convention
            delta_min, delta_max, alpha_min, alpha_max = boundaries

            # draw points uniformly distributed inside the region delimited by
            # boundaries
            u = numpy.random.uniform(0, 1, size)
            delta = numpy.arcsin(
                (numpy.sin(delta_max) - numpy.sin(delta_min)) * u
                + numpy.sin(delta_min)
            )
            alpha = numpy.random.uniform(alpha_min, alpha_max, size)
            # a conversion is required to use get_interp_val
            theta, phi = numpy.pi / 2 - delta, alpha

            data = healpix_map.data
            random_data = numpy.random.uniform(0, data.max(), size)

            # the version of mhealpy 0.3.4 or later is needed to run
            # get_interp_val
            # get_interp_val might return an astropy object depending on the
            # units of the column of the fits file,
            # hence the need to transform into an array
            d_data = numpy.array(
                healpix_map.get_interp_val(theta, phi, lonlat=False)
            )

            dist_theta = theta[d_data > random_data]
            dist_phi = phi[d_data > random_data]

            return (dist_theta, dist_phi)

        # Sampling method to generate the desired number of points
        theta, phi = numpy.array([]), numpy.array([])
        while len(theta) < size:
            new_theta, new_phi = simple_rejection_sampling(
                self.healpix_map, size, self.boundaries
            )
            theta = numpy.concatenate((theta, new_theta), axis=0)
            phi = numpy.concatenate((phi, new_phi), axis=0)

        if len(theta) > size:
            theta = theta[:size]
            phi = phi[:size]

        # convert back to the radec convention
        radec = FieldArray(size, dtype=[('ra', '<f8'), ('dec', '<f8')])
        radec['ra'] = phi
        radec['dec'] = numpy.pi / 2 - theta
        return radec


__all__ = ['UniformSky', 'FisherSky', 'HealpixSky']
