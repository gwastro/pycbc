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

import logging
import numpy
from scipy import interpolate
import astropy.cosmology
from astropy import units
from astropy.cosmology.core import CosmologyError
import pycbc.conversions

DEFAULT_COSMOLOGY = 'Planck15'


def get_cosmology(cosmology=None, **kwargs):
    r"""Gets an astropy cosmology class.

    Parameters
    ----------
    cosmology : str or astropy.cosmology.FlatLambdaCDM, optional
        The name of the cosmology to use. For the list of options, see
        :py:attr:`astropy.cosmology.parameters.available`. If None, and no
        other keyword arguments are provided, will default to
        :py:attr:`DEFAULT_COSMOLOGY`. If an instance of
        :py:class:`astropy.cosmology.FlatLambdaCDM`, will just return that.
    \**kwargs :
        If any other keyword arguments are provided they will be passed to
        :py:attr:`astropy.cosmology.FlatLambdaCDM` to create a custom
        cosmology.

    Returns
    -------
    astropy.cosmology.FlatLambdaCDM
        The cosmology to use.

    Examples
    --------
    Use the default:

    >>> from pycbc.cosmology import get_cosmology
    >>> get_cosmology()
    FlatLambdaCDM(name="Planck15", H0=67.7 km / (Mpc s), Om0=0.307,
                  Tcmb0=2.725 K, Neff=3.05, m_nu=[0.   0.   0.06] eV,
                  Ob0=0.0486)

    Use properties measured by WMAP instead:

    >>> get_cosmology("WMAP9")
    FlatLambdaCDM(name="WMAP9", H0=69.3 km / (Mpc s), Om0=0.286, Tcmb0=2.725 K,
                  Neff=3.04, m_nu=[0. 0. 0.] eV, Ob0=0.0463)

    Create your own cosmology (see :py:class:`astropy.cosmology.FlatLambdaCDM`
    for details on the default values used):

    >>> get_cosmology(H0=70., Om0=0.3)
    FlatLambdaCDM(H0=70 km / (Mpc s), Om0=0.3, Tcmb0=0 K, Neff=3.04, m_nu=None,
                  Ob0=None)

    """
    if kwargs and cosmology is not None:
        raise ValueError("if providing custom cosmological parameters, do "
                         "not provide a `cosmology` argument")
    if isinstance(cosmology, astropy.cosmology.FlatLambdaCDM):
        # just return
        return cosmology
    if kwargs:
        cosmology = astropy.cosmology.FlatLambdaCDM(**kwargs)
    else:
        if cosmology is None:
            cosmology = DEFAULT_COSMOLOGY
        if cosmology not in astropy.cosmology.parameters.available:
            raise ValueError("unrecognized cosmology {}".format(cosmology))
        cosmology = getattr(astropy.cosmology, cosmology)
    return cosmology


def z_at_value(func, fval, unit, zmax=1000., **kwargs):
    r"""Wrapper around astropy.cosmology.z_at_value to handle numpy arrays.

    Getting a z for a cosmological quantity involves numerically inverting
    ``func``. The ``zmax`` argument sets how large of a z to guess (see
    :py:func:`astropy.cosmology.z_at_value` for details). If a z is larger than
    ``zmax``, this will try a larger zmax up to ``zmax * 10**5``. If that still
    is not large enough, will just return ``numpy.inf``.

    Parameters
    ----------
    func : function or method
        A function that takes redshift as input.
    fval : float
        The value of ``func(z)``.
    unit : astropy.unit
        The unit of ``fval``.
    zmax : float, optional
        The initial maximum search limit for ``z``. Default is 1000.
    \**kwargs :
        All other keyword arguments are passed to
        :py:func:``astropy.cosmology.z_at_value``.

    Returns
    -------
    float
        The redshift at the requested values.
    """
    fval, input_is_array = pycbc.conversions.ensurearray(fval)
    # make sure fval is atleast 1D
    if fval.size == 1 and fval.ndim == 0:
        fval = fval.reshape(1)
    zs = numpy.zeros(fval.shape, dtype=float)  # the output array
    for (ii, val) in enumerate(fval):
        try:
            zs[ii] = astropy.cosmology.z_at_value(func, val*unit, zmax=zmax,
                                                  **kwargs)
        except CosmologyError:
            # we'll get this if the z was larger than zmax; in that case we'll
            # try bumping up zmax later to get a value
            zs[ii] = numpy.inf
    # check if there were any zs > zmax
    replacemask = numpy.isinf(zs)
    # try bumping up zmax to get a result
    if replacemask.any():
        # we'll keep bumping up the maxz until we can get a result
        counter = 0  # to prevent running forever
        while replacemask.any():
            kwargs['zmin'] = zmax
            zmax = 10 * zmax
            idx = numpy.where(replacemask)
            for ii in idx:
                val = fval[ii]
                try:
                    zs[ii] = astropy.cosmology.z_at_value(
                        func, val*unit, zmax=zmax, **kwargs)
                    replacemask[ii] = False
                except CosmologyError:
                    # didn't work, try on next loop
                    pass
            counter += 1
            if counter == 5:
                # give up and warn the user
                logging.warning("One or more values correspond to a "
                                "redshift > {0:.1e}. The redshift for these "
                                "have been set to inf. If you would like "
                                "better precision, call God.".format(zmax))
                break
    return pycbc.conversions.formatreturn(zs, input_is_array)


def _redshift(distance, **kwargs):
    r"""Uses astropy to get redshift from the given luminosity distance.

    Parameters
    ----------
    distance : float
        The luminosity distance, in Mpc.
    \**kwargs :
        All other keyword args are passed to :py:func:`get_cosmology` to
        select a cosmology. If none provided, will use
        :py:attr:`DEFAULT_COSMOLOGY`.

    Returns
    -------
    float :
        The redshift corresponding to the given luminosity distance.
    """
    cosmology = get_cosmology(**kwargs)
    return z_at_value(cosmology.luminosity_distance, distance, units.Mpc)


class DistToZ(object):
    r"""Interpolates luminosity distance as a function of redshift to allow for
    fast conversion.

    The :mod:`astropy.cosmology` module provides methods for converting any
    cosmological parameter (like luminosity distance) to redshift. This can be
    very slow when operating on a large array, as it involves numerically
    inverting :math:`z(D)` (where :math:`D` is the luminosity distance). This
    class speeds that up by pre-interpolating :math:`D(z)`. It works by setting
    up a dense grid of redshifts, then using linear interpolation to find the
    inverse function.  The interpolation uses a grid linear in z for z < 1, and
    log in z for ``default_maxz`` > z > 1. This interpolater is setup the first
    time `get_redshift` is called.  If a distance is requested that results in
    a z > ``default_maxz``, the class falls back to calling astropy directly.

    Instances of this class can be called like a function on luminosity
    distances, which will return the corresponding redshifts.

    Parameters
    ----------
    default_maxz : float, optional
        The maximum z to interpolate up to before falling back to calling
        astropy directly. Default is 1000.
    numpoints : int, optional
        The number of points to use in the linear interpolation between 0 to 1
        and 1 to ``default_maxz``. Default is 10000.
    \**kwargs :
        All other keyword args are passed to :py:func:`get_cosmology` to
        select a cosmology. If none provided, will use
        :py:attr:`DEFAULT_COSMOLOGY`.
    """
    def __init__(self, default_maxz=1000., numpoints=10000, **kwargs):
        self.numpoints = int(numpoints)
        self.default_maxz = default_maxz
        self.cosmology = get_cosmology(**kwargs)
        # the interpolating functions; we'll set them to None for now, then set
        # them up when get_redshift is first called
        self.nearby_d2z = None
        self.faraway_d2z = None
        self.default_maxdist = None

    def setup_interpolant(self):
        """Initializes the z(d) interpolation."""
        # for computing nearby (z < 1) redshifts
        zs = numpy.linspace(0., 1., num=self.numpoints)
        ds = self.cosmology.luminosity_distance(zs).value
        self.nearby_d2z = interpolate.interp1d(ds, zs, kind='linear',
                                                bounds_error=False)
        # for computing far away (z > 1) redshifts
        zs = numpy.logspace(0, numpy.log10(self.default_maxz),
                            num=self.numpoints)
        ds = self.cosmology.luminosity_distance(zs).value
        self.faraway_d2z = interpolate.interp1d(ds, zs, kind='linear',
                                                 bounds_error=False)
        # store the default maximum distance
        self.default_maxdist = ds.max()

    def get_redshift(self, dist):
        """Returns the redshift for the given distance.
        """
        dist, input_is_array = pycbc.conversions.ensurearray(dist)
        try:
            zs = self.nearby_d2z(dist)
        except TypeError:
            # interpolant hasn't been setup yet
            self.setup_interpolant()
            zs = self.nearby_d2z(dist)
        # if any points had red shifts beyond the nearby, will have nans;
        # replace using the faraway interpolation
        replacemask = numpy.isnan(zs)
        if replacemask.any():
            zs[replacemask] = self.faraway_d2z(dist[replacemask])
            replacemask = numpy.isnan(zs)
        # if we still have nans, means that some distances are beyond our
        # furthest default; fall back to using astropy
        if replacemask.any():
            # well... check that the distance is positive and finite first
            if not (dist > 0.).all() and numpy.isfinite(dist).all():
                raise ValueError("distance must be finite and > 0")
            zs[replacemask] = _redshift(dist[replacemask],
                                        cosmology=self.cosmology)
        return pycbc.conversions.formatreturn(zs, input_is_array)

    def __call__(self, dist):
        return self.get_redshift(dist)


# set up D(z) interpolating classes for the standard cosmologies
_d2zs = {_c: DistToZ(cosmology=_c)
         for _c in astropy.cosmology.parameters.available}


def redshift(distance, **kwargs):
    r"""Returns the redshift associated with the given luminosity distance.

    If the requested cosmology is one of the pre-defined ones in
    :py:attr:`astropy.cosmology.parameters.available`, :py:class:`DistToZ` is
    used to provide a fast interpolation. This takes a few seconds to setup
    on the first call.

    Parameters
    ----------
    distance : float
        The luminosity distance, in Mpc.
    \**kwargs :
        All other keyword args are passed to :py:func:`get_cosmology` to
        select a cosmology. If none provided, will use
        :py:attr:`DEFAULT_COSMOLOGY`.

    Returns
    -------
    float :
        The redshift corresponding to the given distance.
    """
    cosmology = get_cosmology(**kwargs)
    try:
        z = _d2zs[cosmology.name](distance)
    except KeyError:
        # not a standard cosmology, call the redshift function
        z = _redshift(distance, cosmology=cosmology)
    return z


def redshift_from_comoving_volume(vc, **kwargs):
    r"""Returns the redshift from the given comoving volume.

    Parameters
    ----------
    vc : float
        The comoving volume, in units of cubed Mpc.
    \**kwargs :
        All other keyword args are passed to :py:func:`get_cosmology` to
        select a cosmology. If none provided, will use
        :py:attr:`DEFAULT_COSMOLOGY`.

    Returns
    -------
    float :
        The redshift at the given comoving volume.
    """
    cosmology = get_cosmology(**kwargs)
    return z_at_value(cosmology.comoving_volume, vc, units.Mpc**3)


def distance_from_comoving_volume(vc, **kwargs):
    r"""Returns the luminosity distance from the given comoving volume.

    Parameters
    ----------
    vc : float
        The comoving volume, in units of cubed Mpc.
    \**kwargs :
        All other keyword args are passed to :py:func:`get_cosmology` to
        select a cosmology. If none provided, will use
        :py:attr:`DEFAULT_COSMOLOGY`.

    Returns
    -------
    float :
        The luminosity distance at the given comoving volume.
    """
    cosmology = get_cosmology(**kwargs)
    z = redshift_from_comoving_volume(vc, cosmology=cosmology)
    return cosmology.luminosity_distance(z).value


def cosmological_quantity_from_redshift(z, quantity, strip_unit=True,
                                        **kwargs):
    r"""Returns the value of a cosmological quantity (e.g., age) at a redshift.

    Parameters
    ----------
    z : float
        The redshift.
    quantity : str
        The name of the quantity to get. The name may be any attribute of
        :py:class:`astropy.cosmology.FlatLambdaCDM`.
    strip_unit : bool, optional
        Just return the value of the quantity, sans units. Default is True.
    \**kwargs :
        All other keyword args are passed to :py:func:`get_cosmology` to
        select a cosmology. If none provided, will use
        :py:attr:`DEFAULT_COSMOLOGY`.

    Returns
    -------
    float or astropy.units.quantity :
        The value of the quantity at the requested value. If ``strip_unit`` is
        ``True``, will return the value. Otherwise, will return the value with
        units.
    """
    cosmology = get_cosmology(**kwargs)
    val = getattr(cosmology, quantity)(z)
    if strip_unit:
        val = val.value
    return val


__all__ = ['redshift', 'redshift_from_comoving_volume',
           'distance_from_comoving_volume',
           'cosmological_quantity_from_redshift',
           ]
