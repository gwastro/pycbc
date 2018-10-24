# Copyright (C) 2017  Collin Capano
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
This modules provides functions for computing cosmological quantities, such as
redshift. This is mostly a wrapper around ``astropy.cosmology``.

Note: in all functions, ``distance`` is short hand for ``luminosity_distance``.
Any other distance measure is explicitly named; e.g., ``comoving_distance``.
"""

import astropy.cosmology
from astropy import units

import numpy
import lal
from scipy import interpolate
from pycbc.conversions import ensurearray, formatreturn

DEFAULT_COSMOLOGY='Planck15'

def get_cosmology(cosmology=DEFAULT_COSMOLOGY):
    """Gets an astropy cosmology class based on the given string.

    Parameters
    ----------
    cosmology : str, optional
        The name of the cosmology to use. For the list of options, see
        :py:attr:`astropy.cosmology.parameters.available`. Default is
        ``'Planck15'``.

    Returns
    -------
    astropy.cosmology instance
        The cosmology to use from astropy. (Currently, all are instances of
        :py:class:`astropy.cosmology.FlatLambdaCDM`.)
    """
    if cosmology not in astropy.cosmology.parameters.available:
        raise ValueError("unrecognized cosmology {}".format(cosmology))
    return getattr(astropy.cosmology, cosmology)


def z_at_value(func, fval, unit, **kwargs):
    """Wrapper around astropy.cosmology.z_at_value to handle numpy arrays.
    
    The unit must be specified as a separate argument.
    """
    fval, input_is_array = ensurearray(fval)
    # make sure fval is atleast 1D
    if fval.size == 1 and fval.ndim == 0:
        fval = fval.reshape(1)
    zs = [astropy.cosmology.z_at_value(func, val*unit, **kwargs)
          for val in fval]
    return formatreturn(numpy.array(zs), input_is_array)


def redshift(distance, cosmology=DEFAULT_COSMOLOGY):
    """Returns the redshift associated with the given luminosity distance.

    Parameters
    ----------
    distance : float
        The luminosity distance, in Mpc.
    cosmology : str, optional
        The name of the cosmology to use. For the list of options, see
        :py:attr:`astropy.cosmology.parameters.available`. Default is
        ``'Planck15'``.

    Returns
    -------
    float :
        The redshift corresponding to the given luminosity distance.
    """
    cosmology = get_cosmology(cosmology)
    return z_at_value(cosmology.luminosity_distance, distance, units.Mpc)


def distance_from_comoving_volume(vc, cosmology=DEFAULT_COSMOLOGY):
    """Returns the luminosity distance from the given comoving volume.

    Parameters
    ----------
    vc : float
        The comoving volume, in units of cubed Mpc.
    cosmology : str, optional
        The name of the cosmology to use. For the list of options, see
        :py:attr:`astropy.cosmology.parameters.available`. Default is
        ``'Planck15'``.

    Returns
    -------
    float :
        The luminosity distance at the given comoving volume.
    """
    cosmology = get_cosmology(cosmology)
    # first get the redshift associated with the given comoving volume
    z = z_at_value(cosmology.comoving_volume, vc, units.Mpc**3)
    # now convert redshift to luminosity distance
    return cosmology.luminosity_distance(z).value

class _DistToZ(object):
    """Class to convert luminosity distance to redshift using the given
    cosmology.

    Cosmological constants `h, om, ol, w0, w1, w2` may optionally be provided
    (see below). Otherwise, standard LambdaCDM will be used by default; see
    `LALCosmologyCalculator` for details.

    LAL provides functions to convert redshift to luminosity distance,
    assuming a cosmology. This class works by setting up a dense grid of
    redshifts, then using linear interpolation to find the inverse function.
    The interpolation uses a grid linear in z for z < 1, and log in z for z >
    1.  For speed, a default interpolater is setup on initialization for
    redshifts out to `default_maxz`. If a requested distance is found to have
    a redshift greater than `default_maxz`, a new set of grid points is
    created from `0.1*default_maxz` to `1e5*default_maxz` and re-interpolated.
    This continues until a redshift can be found.

    Instances of this class can be called like a function on luminosity
    distances, which will return the corresponding redshifts.

    Parameters
    ----------
    h : {None, float}
        The normalized Hubble constant to use. If None, will default to
        standard cosmology.
    om : {None, float}
        The matter energy density. If None, will default to standard cosmology.
    ol : {None, float}
        The cosmological-constant density. If None, will default to standard
        cosmology.
    w0 : {None, float}
        The 0th-order dark energy equation of state parameter. If None, will
        default to standard cosmology.
    w1 : {None, float}
        The first-order dark energy equation of state parameter. If None, will
        default to standard cosmology.
    w2 : {None, float}
        The second-order dark energy equation of state parameter. If None, will
        default to standard cosmology.
    default_maxz : {2., float}
        The default maximum redshift to use for the interpolation.
    numpoints : {5e4, int}
        The number of points to use between for the interpolation.
    """
    def __init__(self, h=None, om=None, ol=None, w0=None, w1=None, w2=None,
                 default_maxz=1e5, numpoints=1e4):
        default_omega = lal.CreateCosmologicalParameters(0., 0., 0., 0., 0.,
                                                         0.)
        lal.SetCosmologicalParametersDefaultValue(default_omega)
        # if any cosmological parameter was specified, create a custom
        # cosmology
        if any([p is not None for p in [h, om, ol, w0, w1, w2]]):
            h = h if h is not None else default_omega.h
            om = om if om is not None else default_omega.om
            ol = ol if ol is not None else default_omega.ol
            w0 = w0 if w0 is not None else default_omega.w0
            w1 = w1 if w1 is not None else default_omega.w1
            w2 = w2 if w2 is not None else default_omega.w2
            self.omega = lal.CreateCosmologicalParameters(h, om, ol,
                                                          w0, w1, w2)
        else:
            self.omega = default_omega
        self.default_maxz = default_maxz
        self.numpoints = int(numpoints)
        self.z2d = numpy.vectorize(lal.LuminosityDistance)
        # for computing nearby (z < 1) redshifts
        zs = numpy.linspace(0., 1., num=self.numpoints)
        ds = self.z2d(self.omega, zs)
        self.nearby_d2z = interpolate.interp1d(ds, zs, kind='linear',
                                                bounds_error=False)
        # for computing far away (z > 1) redshifts
        zs = numpy.logspace(0, numpy.log10(default_maxz), num=self.numpoints)
        ds = self.z2d(self.omega, zs)
        self.faraway_d2z = interpolate.interp1d(ds, zs, kind='linear',
                                                 bounds_error=False)
        # store the default maximum distance
        self.default_maxdist = ds.max()

    def get_redshift(self, dist):
        """Returns the redshift for the given distance.
        """
        dist, input_is_array = ensurearray(dist)
        zs = self.nearby_d2z(dist)
        # if any points had red shifts beyond the nearby, will have nans;
        # replace using the faraway interpolation
        replacemask = numpy.isnan(zs)
        if replacemask.any():
            zs[replacemask] = self.faraway_d2z(dist[replacemask])
            replacemask = numpy.isnan(zs)
            # if we still have nans, means that some distances are beyond our
            # furthest default; replace
            if replacemask.any():
                # well... check that the distance is positive and finite first
                if not (dist > 0.).all() and numpy.isfinite(dist).all():
                    raise ValueError("distance must be finite and > 0")
                minz = numpy.log10(self.default_maxz) - 1.
                maxz = minz + 5
                while replacemask.any():
                    gridzs = numpy.logspace(minz, maxz, num=self.numpoints)
                    ds = self.z2d(self.omega, gridzs)
                    d2z = interpolate.interp1d(ds, gridzs, kind='linear',
                                                bounds_error=False)
                    zs[replacemask] = d2z(dist[replacemask])
                    replacemask = numpy.isnan(zs)
                    minz = maxz - 1
                    maxz = minz + 5
        return formatreturn(zs, input_is_array)

    def __call__(self, dist):
        return self.get_redshift(dist)

# we'll use the default cosmology for computing red shifts.
_d2z = _DistToZ()

def lalredshift(distance, h=None, om=None, ol=None, w0=None, w1=None, w2=None):
    """Returns the redshift associated with the given distance.

    Cosmological constants `h, om, ol, w0, w1, w2` may optionally be provided
    (see below). Otherwise, standard LambdaCDM will be used by default; see
    `LALCosmologyCalculator` for details.

    Parameters
    ----------
    distance : float
        The luminosity distance, in Mpc.
    h : {None, float}
        The normalized Hubble constant to use. If None, will default to
        standard cosmology.
    om : {None, float}
        The matter energy density. If None, will default to standard cosmology.
    ol : {None, float}
        The cosmological-constant density. If None, will default to standard
        cosmology.
    w0 : {None, float}
        The 0th-order dark energy equation of state parameter. If None, will
        default to standard cosmology.
    w1 : {None, float}
        The first-order dark energy equation of state parameter. If None, will
        default to standard cosmology.
    w2 : {None, float}
        The second-order dark energy equation of state parameter. If None, will
        default to standard cosmology.

    Returns
    -------
    float :
        The redshift corresponding to the given distance.
    """
    if all([p is None for p in [h, om, ol, w0, w1, w2]]):
        # just use the default converter
        d2z = _d2z
    else:
        d2z = _DistToZ(h=h, om=om, ol=ol, w0=w0, w1=w1, w2=w2)
    return d2z(distance)


__all__ = ['redshift', 'distance_from_comoving_volume', 'lalredshift']
