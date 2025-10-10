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

import copy
import logging
import warnings
import numpy
from scipy.spatial.transform import Rotation

from pycbc.distributions import angular
from pycbc import VARARGS_DELIM
from pycbc.io import FieldArray
from pycbc.types import angle_as_radians
from pycbc.libutils import import_optional


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

    def to_uniform_patch(self, coverage):
        if coverage < 1:
            logging.warning(
                'Attempt to convert UniformSky to a '
                'uniform patch assumes 100% coverage'
            )
        return self


class UniformDiskSky:
    """A distribution that represents a uniform disk on the sky. The declination
    varies from π/2 to -π/2 and the right ascension varies from 0 to 2π.

    Parameters
    ----------
    mean_ra: float or str
        RA of the center of the distribution. Use the rad or deg suffix to
        specify units, otherwise radians are assumed.
    mean_dec: float or str
        Declination of the center of the distribution. Use the rad or deg
        suffix to specify units, otherwise radians are assumed.
    radius: float or str
        Radius of the disk. Use the rad or deg suffix to specify units,
        otherwise radians are assumed.
    """
    name = 'uniform_disk_sky'
    _params = ['ra', 'dec']

    def __init__(self, **params):
        mean_ra = angle_as_radians(params['mean_ra'])
        mean_dec = angle_as_radians(params['mean_dec'])
        radius = angle_as_radians(params['radius'])
        if mean_ra < 0 or mean_ra > 2 * numpy.pi:
            raise ValueError(
                f'The mean RA must be between 0 and 2π, {mean_ra} rad given'
            )
        if mean_dec < -numpy.pi / 2 or mean_dec > numpy.pi / 2:
            raise ValueError(
                'The mean declination must be between '
                f'-π/2 and π/2, {mean_dec} rad given'
            )
        if radius < 0 or radius > 2 * numpy.pi:
            raise ValueError(
                'Radius must be non-negative and smaller than 2π'
            )
        # Prepare a rotation that puts the North Pole at the mean position
        self.rotation = Rotation.from_euler(
            'yz', [numpy.pi / 2 - mean_dec, mean_ra]
        )
        self.mean_ra, self.mean_dec, self.radius = mean_ra, mean_dec, radius

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
        radius = cp.get_opt_tag(section, 'radius', tag)
        return cls(
            mean_ra=mean_ra,
            mean_dec=mean_dec,
            radius=radius,
        )

    def get_max_prob_point(self):
        return (self.mean_ra, self.mean_dec)

    def rvs(self, size):
        # Draw samples from a distribution centered on the North pole
        np_ra = numpy.random.uniform(low=0, high=(2 * numpy.pi), size=size)
        np_dec = angular.SinAngle(
            polar_bounds=(0,self.radius)).rvs(
            size=size).astype(numpy.float64)

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

    def to_uniform_patch(self, coverage):
        if coverage < 1:
            logging.warning(
                'Attempt to convert UniformDiskSky to a '
                'uniform patch assumes 100% coverage'
            )
        return self


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
        if sigma < 0 or sigma > 2 * numpy.pi:
            raise ValueError(
                'Sigma must be non-negative and smaller than 2π '
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
        # storing center position for `to_uniform_patch()`
        self.mean_ra, self.mean_dec = mean_ra, mean_dec

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

    def get_max_prob_point(self):
        return (self.mean_ra, self.mean_dec)

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

    def to_uniform_patch(self, coverage):
        radius = numpy.sqrt(self.rayleigh_scale**2 * (-2*numpy.log(1-coverage)))
        return UniformDiskSky(mean_ra=self.mean_ra, mean_dec=self.mean_dec, radius=radius)


class HealpixSky:
    """Sky-location distribution given by a HEALPix map from an external file.

    The declination (delta) varies from π/2 to -π/2 and the right ascension
    (alpha) varies from 0 to 2π.

    Parameters
    ----------
    healpix_file : str
        Path to a FITS file containing a probability distribution encoded in a
        HEALPix scheme.
    """

    name = 'healpix_sky'
    _params = ['ra', 'dec']

    def __init__(self, **params):
        # Read the map file.
        file_name = params['healpix_file']
        mhealpy = import_optional('mhealpy')
        self.healpix_map = mhealpy.HealpixMap.read_map(file_name)

        # Get the probabilities at each pixel.
        if self.healpix_map.density():
            # If the map stores a probability density, then we must convert to
            # a probability per pixel.  We assume that the angular variation
            # scale of the density function is much larger than the pixel area,
            # so that we can just scale the density by the pixel area.
            pix_areas = self.healpix_map.pixarea(
                numpy.arange(self.healpix_map.npix)
            )
            self.pix_probs = (self.healpix_map * pix_areas).data
        else:
            # The map stores directly a probability per pixel.
            self.pix_probs = self.healpix_map.data

        # Sanity-check the probabilities, and ensure they are normalized
        # correctly (sum to one).
        assert type(self.pix_probs) == numpy.ndarray
        sum_pix_probs = sum(self.pix_probs)
        if not numpy.isclose(sum_pix_probs, 1):
            warnings.warn(
                f'Sum of probs in HEALPix map is {sum(self.pix_probs)}, '
                'far from 1. Something might be wrong with that map'
            )
        self.pix_probs /= sum_pix_probs

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
        return cls(healpix_file=healpix_file)

    def get_max_prob_point(self):
        coords = self.healpix_map.pix2ang(
                numpy.where(self.pix_probs==self.pix_probs.max())[0],
                lonlat=True
        )
        return (numpy.deg2rad(coords[0][0]), numpy.deg2rad(coords[1][0]))

    def pixel_corners(self, indices):
        """Return the Cartesian vectors corresponding to the corners of one or
        more HEALPix pixels. Dimension 0 is the pixel index, 1 is the Cartesian
        coordinate, 2 is the corner index.
        """
        return numpy.array(
            [self.healpix_map.boundaries(pi, step=1) for pi in indices]
        )

    def normalize_azimuth(self, phi):
        """Helper function to ensure the azimuthal coordinate stays within
        [0,2π).
        """
        phi[phi < 0] += 2 * numpy.pi
        phi[phi >= 2 * numpy.pi] -= 2 * numpy.pi
        return phi

    def to_uniform_patch(self, coverage):
        """Return a new HealpixSky object that represents a patch of uniform probability,
        defined as the region that encloses a given probability in the original map.
        """
        uniform_patch = copy.deepcopy(self)
        non_zero_ind = numpy.flatnonzero(uniform_patch.pix_probs)
        dtype = numpy.dtype([('index', numpy.ndarray), ('prob', numpy.float64)])
        prob = uniform_patch.pix_probs[non_zero_ind]
        ind_prob = numpy.array(list(zip(non_zero_ind, prob)), dtype=dtype)
        ind_prob_sorted = numpy.sort(ind_prob, order='prob')
        prob = numpy.sort(prob)
        ind_prob_rev_sorted = ind_prob_sorted[::-1]
        list_ind = [(ind_prob_rev_sorted[i][0]) for i in range(len(ind_prob_sorted))]
        cum_sum_prob = numpy.cumsum(prob[::-1])
        covered_sky = cum_sum_prob[cum_sum_prob<coverage]
        uniform_patch.pix_probs[list_ind[:len(covered_sky)]] = 1/len(covered_sky)
        uniform_patch.pix_probs[uniform_patch.pix_probs!= 1/len(covered_sky)] = 0
        return uniform_patch

    def rvs(self, size):
        # First of all, draw a random sample of pixel indices following the
        # given probabilities. This is the trivial part of the algorithm.
        pix_indices = numpy.random.choice(
            self.healpix_map.npix, p=self.pix_probs, size=size
        )

        # Then comes the nasty part of the algorithm. For each selected pixel,
        # we need to draw a random point uniformly distributed on the patch of
        # the sky defined by that pixel. We cannot skip points, every pixel
        # (index) needs precisely one point. We do this via rejection sampling:
        # we find the (theta,phi) box that encloses the corners of a pixel,
        # then draw a uniform point inside that box and repeat until the point
        # falls inside the pixel. The acceptance rate is always around 50% due
        # to the orientation of the HEALPix pixels. The result seems correct
        # visually and by checking the acceptance rates. This method is
        # relatively slow due to the calls to `boundaries()`, which is slow
        # especially for MOC maps. For MOC maps, I obtain around one million
        # samples per minute.

        boundaries_vec = self.pixel_corners(pix_indices)

        boundaries_z_min = boundaries_vec[:,2,:].min(axis=1)
        boundaries_z_max = boundaries_vec[:,2,:].max(axis=1)

        # Find ranges of azimuthal angle. `arctan2()` produces angles in the
        # range [-π,π] which is not the right convention, so fix that right
        # away.
        boundaries_phi = numpy.arctan2(
            boundaries_vec[:,1,:], boundaries_vec[:,0,:]
        )
        del boundaries_vec
        boundaries_phi[boundaries_phi < 0] += 2 * numpy.pi
        boundaries_phi_min = boundaries_phi.min(axis=1)
        boundaries_phi_max = boundaries_phi.max(axis=1)
        # Pixels that straddle phi = 0 need a special treatment :( Shift them
        # by 90 deg so that they are not straddlers for the time being.
        # We will shift the corresponding samples back after drawing them.
        straddler_mask = (boundaries_phi_max - boundaries_phi_min) > numpy.pi / 2
        boundaries_phi[straddler_mask,:] -= numpy.pi / 2
        boundaries_phi[straddler_mask,:] = self.normalize_azimuth(
            boundaries_phi[straddler_mask,:]
        )
        boundaries_phi_min[straddler_mask] = boundaries_phi[straddler_mask,:].min(axis=1)
        boundaries_phi_max[straddler_mask] = boundaries_phi[straddler_mask,:].max(axis=1)
        del boundaries_phi

        final_thetas = numpy.array([])
        final_phis = numpy.array([])

        # Here comes the rejection sampling.
        while len(pix_indices) > 0:
            random_thetas = numpy.arccos(
                numpy.random.uniform(boundaries_z_min, boundaries_z_max)
            )
            random_phis = numpy.random.uniform(
                boundaries_phi_min, boundaries_phi_max
            )

            # Now we can undo the shift for the straddlers. Ugly!
            random_phis[straddler_mask] += numpy.pi / 2
            random_phis[straddler_mask] = self.normalize_azimuth(
                random_phis[straddler_mask]
            )

            # Find which points fall inside their pixels.
            sampled_indices = self.healpix_map.ang2pix(
                random_thetas, random_phis
            )
            acceptance_mask = sampled_indices == pix_indices
            final_thetas = numpy.concatenate(
                (final_thetas, random_thetas[acceptance_mask])
            )
            final_phis = numpy.concatenate(
                (final_phis, random_phis[acceptance_mask])
            )

            # Iterate if we are still missing some pixels.
            rej_mask = numpy.logical_not(acceptance_mask)
            pix_indices = pix_indices[rej_mask]
            boundaries_z_min = boundaries_z_min[rej_mask]
            boundaries_z_max = boundaries_z_max[rej_mask]
            boundaries_phi_min = boundaries_phi_min[rej_mask]
            boundaries_phi_max = boundaries_phi_max[rej_mask]
            straddler_mask = straddler_mask[rej_mask]

        # Convert back to the radec convention
        radec = FieldArray(size, dtype=[('ra', '<f8'), ('dec', '<f8')])
        radec['ra'] = final_phis
        radec['dec'] = numpy.pi / 2 - final_thetas
        return radec


__all__ = ['UniformSky', 'UniformDiskSky', 'FisherSky', 'HealpixSky']
