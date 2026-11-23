# Copyright (C) 2025  Sumit kumar, Shichao Wu
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

import warnings
from abc import ABC, abstractmethod

import astropy.units as u
import numpy as np
import scipy.integrate as scipy_integrate
from pycbc.cosmology import get_cosmology
from scipy.interpolate import CubicSpline

# Base redshift evolution class
class BaseRedshiftEvolution(ABC):
    """
    Abstract Base class for redshift evolution models with Astropy unit
    support.

    This class provides a framework for modeling the redshift dependence of
    astrophysical event rates, including differential comoving volume,
    spacetime volume, and redshift probability distributions. It supports
    various cosmological models and integrates
    with time delay distributions.

    The core functionality includes:
    - Differential comoving volume (3-volume): dVc/dz
    - Differential spacetime volume (4-volume): dVT/dz = psi(z)/(1+z) * dVc/dz
    - Redshift probability density functions (PDFs)
    - Normalization and total spacetime volume computation


    Notes
    -----
    - Redshift `z` is dimensionless.
    - Comoving volumes are returned in units of Gpc^3.
    - Spacetime volumes are returned in Gpc^3 (with implicit time integration).
    - The class uses numerical integration over a redshift grid for performance.

    Parameters
    ----------
    z_max : float, optional
        Maximum redshift for integration and normalization. Default is 2.0.
    num_zbins : int, optional
        Number of redshift bins for numerical integration. Must be at least 2.
        Default is 1000.
    cosmology : astropy.cosmology.Cosmology, optional
        Cosmological model to use. If None, defaults to
        `pycbc.cosmology.get_cosmology()`.
    z_grid : array_like, optional
        Custom redshift grid (must be strictly increasing, start at 0,
                              end at z_max).
        If None, a uniform grid from 0 to z_max is used.

    Attributes
    ----------
    z_max : float
        Maximum redshift used for integration.
    num_zbins : int
        Number of redshift bins.
    cosmology : astropy.cosmology.Cosmology
        Cosmological model used for volume calculations.
    _z_grid : numpy.ndarray
        Internal redshift grid (dimensionless).
    _dvc_dz_grid : astropy.units.Quantity
        Differential comoving volume per steradian on the grid (Gpc³/sr).
    _cached_z : numpy.ndarray or None
        Cached redshift values for interpolation.
    _cached_dvc_dz : astropy.units.Quantity or None
        Cached interpolated dVc/dz values.
    td_min : float
        Minimum time delay (in Gyr).
    td_max : float
        Maximum time delay (in Gyr), derived from cosmology and z_formation_max.

    See Also
    --------
    PowerLawRedshift, GRB2008SFR, MadauDickinson2014SFR, MadauFragos2017SFR,
    SFRTimeDelayRedshift : Concrete implementations of redshift evolution
        models.
    """

    def __init__(
        self,
        z_max: float = 2.0,
        num_zbins: int = 1000,
        cosmology=None,
        z_grid = None,
        z_max_cosmo: float = 5000.0,
        z_grid_cosmo = None,
        ):
        """
        Initialize the redshift evolution model.

        Parameters
        ----------
        z_max : float
            Maximum redshift for normalization and integration.
            Must be positive.
        num_zbins : int
            Number of redshift bins for numerical integration.
            Must be at least 2.
        cosmology : astropy.cosmology.Cosmology, optional
            Cosmological model to use. If None, defaults to `get_cosmology()`.
        z_grid : array_like, optional
            Custom redshift grid
            (must be strictly increasing, start at 0, end at z_max).
            If None, a uniform grid from 0 to z_max is used.
        z_max_cosmo : float
            Maximum redshift for grid for cosmological interpolation.
            Must be positive. It is different from z_max which is required
            for normalization.
        z_grid_cosmo : array_like, optional
            Custom redshift grid for cosmological functions
            (must be strictly increasing, start at 0, end at z_max).
            If None, a uniform grid from 0 to z_max is used.
        """

        # Validate z_max, z_max_cosmo, num_bin,
        if not np.isfinite(z_max) or z_max <= 0:
            raise ValueError("z_max must be finite and positive")

        if not isinstance(num_zbins, (int, np.integer)) or num_zbins < 2:
            raise ValueError("num_zbins must be an integer >= 2")

        if not np.isfinite(z_max_cosmo) or z_max_cosmo <= 0:
            raise ValueError("z_max_cosmo must be finite and positive")

        if z_max_cosmo < z_max:
            raise ValueError("z_max_cosmo must be >= z_max")

        self.z_max = float(z_max)
        self.z_max_cosmo = float(z_max_cosmo)
        self.num_zbins = int(num_zbins)

        # Validate cosmology
        if cosmology is not None:
            self.cosmology = cosmology
        else:
            self.cosmology = get_cosmology()

        # Redshift grid (dimensionless) validation
        if z_grid is None:
            self._z_grid = np.linspace(0.0, self.z_max, num_zbins)
        else:
            self._z_grid = np.asarray(z_grid, dtype=float)
            if self._z_grid.ndim != 1:
                raise ValueError("z_grid must be one-dimensional")

            if np.any(self._z_grid < 0):
                raise ValueError("z_grid must be non-negative")

            if np.any(np.diff(self._z_grid) <= 0):
                raise ValueError("z_grid must be strictly increasing")

            # Validation of end points on the grid
            if not np.isclose(self._z_grid[0], 0):
                raise ValueError("z_grid should start at z=0")

            if not np.isclose(self._z_grid[-1], self.z_max):
                raise ValueError("z_grid must end at z_max")

        # Differential comoving volume on grid
        # Full-sky differential comoving volume dVc/dz = 4pi dVc/dz
        self._dvc_dz_grid = (
            4.0
            * np.pi
            * self.cosmology.differential_comoving_volume(self._z_grid)
            .to(u.Gpc**3 / u.sr)
        )*u.sr

        # Cache for interpolated values
        self._cached_z = None
        self._cached_dvc_dz = None

        # An additional grid for the cosmological quantities
        # that need not depend on the z_max and normalization
        # is not needed.
        self._z_grid_cosmo = \
                np.linspace(0.0, self.z_max_cosmo, int(10*self.num_zbins))
        self._dvc_dz_grid_cosmo = (
            4.0
            * np.pi
            * self.cosmology.differential_comoving_volume(self._z_grid_cosmo)
            .to(u.Gpc**3 / u.sr)
        )*u.sr

        # Cache for interpolated values
        #self._cached_z = None
        #self._cached_dvc_dz = None

    # ------------------------------------------------------------------
    # Utilities; it makes it compatible with scalar/numpy.arrays
    # ------------------------------------------------------------------

    @staticmethod
    def _to_1d_array(redshift):
        """
        Convert scalar or array-like redshift to 1D ndarray.

        This utility ensures consistent handling of scalar and array inputs.

        Parameters
        ----------
        redshift : scalar or array_like
            Dimensionless redshift(s).

        Returns
        -------
        z_array : numpy.ndarray
            1D array of redshift values (float).
        is_scalar : bool
            True if input was a scalar, False otherwise.

        Raises
        ------
        ValueError
            If any redshift value is negative.
        """
        if np.isscalar(redshift):
            z =  np.array([redshift], dtype=float)
            is_scalar = True
        else:
            z = np.asarray(redshift, dtype=float).ravel()
            is_scalar = False

        # Validate for positive and finite redshift values
        if np.any(z < 0):
            raise ValueError("Redshift must be non-negative")

        if not np.all(np.isfinite(z)):
            raise ValueError("Redshift must contain only finite values")

        return z, is_scalar

    # ------------------------------------------------------------------
    # Cosmology utilities
    # ------------------------------------------------------------------

    def fullsky_differential_comoving_volume(self, redshift):
        r"""
        Differential comoving 3-volume per unit redshift.

        This is the full-sky differential comoving volume:
            dVc/dz = 4 \pi dVc/dz/d_\Omega.

        Parameters
        ----------
        redshift : array_like
            Dimensionless redshift(s).

        Returns
        -------
        astropy.units.Quantity
            Differential comoving volume in Gpc^3.

        """

        zz, is_scalar = self._to_1d_array(redshift)

        dvc = (
            4.0
            * np.pi
            * self.cosmology.differential_comoving_volume(zz)
            .to(u.Gpc**3 / u.sr)
        )*u.sr

        return dvc[0] if is_scalar else dvc

    def fullsky_differential_spacetime_volume(self, redshift):
        r"""
        Differential comoving 4-volume per unit redshift.

        This is the full-sky differential comoving volume:
            dVc/dz = 4 \pi 1/(1+z) dVc/dz/d_\Omega.

        Parameters
        ----------
        redshift : array_like
            Dimensionless redshift(s).

        Returns
        -------
        astropy.units.Quantity
            Differential comoving volume in Gpc^3.

        """

        zz, is_scalar = self._to_1d_array(redshift)

        dvc = (
            4.0
            * np.pi/(1 + zz)
            * self.cosmology.differential_comoving_volume(zz)
            .to(u.Gpc**3 / u.sr)
        )*u.sr

        return dvc[0] if is_scalar else dvc

    def comoving_volume(self, redshift):
        """
        Comoving volume enclosed within a given redshift.

        Parameters
        ----------
        redshift : array_like

        Returns
        -------
        astropy.units.Quantity
            Comoving volume in Gpc^3.
        """

        zz, is_scalar = self._to_1d_array(redshift)
        vol = self.cosmology.comoving_volume(zz).to(u.Gpc**3)

        return vol[0] if is_scalar else vol


    def comoving_spacetime_volume(self, redshift):
        """
        Cumulative full-sky comoving volume-time factor up to redshift z.

        Computes:
            V_ct(z) = integral_0^z [dVc/dz' / (1 + z')] dz'

        Parameters
        ----------
        redshift : float or array_like
            Redshift value(s). Must be finite and non-negative.

        Returns
        -------
        astropy.units.Quantity
            Cumulative volume-time factor in Gpc^3 yr.
            Returns a scalar Quantity for scalar input and an array
            of Quantities for array-like input.

        Notes
        -----
        Redshifts above self.z_max are assigned zero. Values within
        [0, self.z_max] are evaluated by interpolation of the cumulative
        integral on the internal redshift grid.
        """
        z, is_scalar = self._to_1d_array(redshift)

        # Full-sky dVc/dz in Gpc^3
        dvc_dz = self._dvc_dz_grid_cosmo.to_value(u.Gpc**3)

        # Observer-frame time-dilation factor
        integrand = dvc_dz / (1.0 + self._z_grid_cosmo)

        # Cumulative integral from z = 0 to each grid point
        cumulative_volume = np.zeros_like(self._z_grid_cosmo, dtype=float)
        dz = np.diff(self._z_grid_cosmo)

        cumulative_volume[1:] = np.cumsum(
            0.5 * (integrand[1:] + integrand[:-1]) * dz
        )

        # Interpolate to requested redshifts.
        # Below zero: zero; above z_max: saturate at the final value.
        result = np.interp(
            z,
            self._z_grid_cosmo,
            cumulative_volume,
            left=0.0,
            right=cumulative_volume[-1],
        ) * u.Gpc**3 * u.yr

        return result[0] if is_scalar else result



    # ------------------------------------------------------------------
    # Redshift evolution model
    # ------------------------------------------------------------------

    @abstractmethod
    def psi_z(self, redshift, **parameters):
        """
        Redshift evolution function psi(z).

        This function defines the intrinsic redshift dependence of the event
        rate.

        Parameters
        ----------
        redshift : array_like
            Dimensionless redshift(s).
        **parameters
            Model-specific parameters (e.g., k, gamma, kappa, z_peak, etc.).

        Returns
        -------
        numpy.ndarray
            Dimensionless evolution factor psi(z) at each redshift.

        Raises
        ------
        NotImplementedError
            Must be implemented by subclasses.
        """
        raise NotImplementedError

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _cache_dvc_dz(self, redshift):
        """
        Cache interpolated differential comoving volume for faster repeated
        evaluation.

        This method caches the interpolated dVc/dz values on a given redshift
        grid.

        Parameters
        ----------
        redshift : numpy.ndarray
            Dimensionless redshift array.

        Returns
        -------
        astropy.units.Quantity
            Cached dVc/dz values in Gpc³.

        Notes
        -----
        - Uses linear interpolation with zero padding at boundaries.
        - The cache is updated only if the input redshift differs from the
        cached one.
        """
        redshift = np.asarray(redshift)
        self._cached_z = redshift
        self._cached_dvc_dz = np.interp(
            redshift,
            self._z_grid,
            self._dvc_dz_grid.value,
            left=0.0,
            right=0.0,
        ) * u.Gpc**3
        return self._cached_dvc_dz

    # ------------------------------------------------------------------
    # 4-volume and normalization
    # ------------------------------------------------------------------

    def dVT_dz(self, redshift, **parameters):
        """
        Differential spacetime (4-volume) element.

        Defined as:
            dVT/dz = psi(z) / (1 + z) * dVc/dz

        This represents the effective rate of events per unit redshift,
        accounting for cosmological volume and redshift evolution.

        Parameters
        ----------
        redshift : array_like
            Dimensionless redshift(s).
        **parameters
            Parameters passed to `psi_z`.

        Returns
        -------
        astropy.units.Quantity
            Differential spacetime volume with units of Gpc^3.
        """

        z, is_scalar = self._to_1d_array(redshift)
        psi = self.psi_z(z, **parameters)

        # Warn if any redshift exceeds z_max
        invalid = z > self.z_max
        if np.any(invalid):
            warnings.warn(
                f"Some redshift values exceed z_max={self.z_max}. "
                "These values will be assigned zero spacetime volume.",
                UserWarning
            )

        # Use cache for array inputs
        if not is_scalar:
            if (self._cached_z is None) or (not np.array_equal(z, self._cached_z)):
                dvc_dz = self._cache_dvc_dz(z)
            else:
                dvc_dz = self._cached_dvc_dz
        else:
            # For scalar, interpolate on-the-fly
            dvc_dz = np.interp(z, self._z_grid, self._dvc_dz_grid.value) * u.Gpc**3

        result = psi * dvc_dz / (1.0 + z)
        return result[0] if is_scalar else result

    def normalize(self, **parameters):
        """
        Normalization constant for the redshift probability distribution.

        This is the integral of dVT/dz over redshift, used to normalize the PDF.

        Parameters
        ----------
        parameters : dict
            Parameters passed to `psi_z`.

        Returns
        -------
        astropy.units.Quantity
            Normalization constant with units of Gpc^3.

        Notes
        -----
        - Uses trapezoidal rule integration over the internal redshift grid.
        - The result is used in `prob_redshift` to compute the normalized PDF.
        """
        psi = self.psi_z(self._z_grid, **parameters)

        integrand = (
            psi
            * self._dvc_dz_grid
            / (1.0 + self._z_grid)
        )

        return np.trapz(integrand.value, self._z_grid) * u.Gpc**3

    def total_4volume(self, analysis_time: u.Quantity, **parameters):
        """
        Total spacetime volume over an observation time.

        This is the product of the normalization constant and the observation
        time.

        Parameters
        ----------
        analysis_time : astropy.units.Quantity
            Observation time (e.g. years).
        parameters : dict
            Parameters passed to `psi_z`.

        Returns
        -------
        astropy.units.Quantity
            Total spacetime volume (Gpc^3 * time).

        Raises
        ------
        astropy.units.UnitTypeError
            If analysis_time is not a time quantity.

        See Also
        --------
        normalize : Compute normalization constant.
        dVT_dz : Differential spacetime volume.
        """
        if not analysis_time.unit.is_equivalent(u.s):
            raise u.UnitTypeError("analysis_time must be a time quantity")

        return self.normalize(**parameters) * analysis_time

    def prob_redshift(self, redshift, **parameters):
        """
        Normalized redshift probability density function.

        Parameters
        ----------
        redshift : array_like
            Dimensionless redshift.
        parameters : dict
            Parameters passed to `psi_z`.

        Returns
        -------
        ndarray
            Dimensionless probability density p(z).
        """
        z, is_scalar = self._to_1d_array(redshift)
        pdf = np.zeros_like(z, dtype=float)

        # Identify redshifts beyond z_max
        # This is for the normalization
        invalid = z > self.z_max
        if np.any(invalid):
            warnings.warn(
                f"Some redshift values exceed z_max={self.z_max}. "
                "These values will be assigned zero probability.",
                UserWarning
            )

        # Only valid redshifts contribute
        valid = z <= self.z_max
        norm = self.normalize(**parameters)
        pdf[valid] = \
                (self.dVT_dz(z[valid], **parameters) / norm).decompose().value
        return pdf[0] if is_scalar else pdf



class PowerLawRedshift(BaseRedshiftEvolution):
    """
    PPower-law redshift distribution model.

    The redshift evolution function is defined as:
        psi(z) = (1 + z)^k

    This model is commonly used for simple, analytic redshift evolution.

    Parameters
    ----------
    z_max : float, optional
        Maximum redshift. Default is 2.0.
    num_zbins : int, optional
        Number of redshift bins. Default is 1000.
    cosmology : astropy.cosmology.Cosmology, optional
        Cosmological model. Default is `get_cosmology()`.
    z_grid : array_like, optional
        Custom redshift grid.

    Attributes
    ----------
    name : str
        Model name: "power_law"
    param_names : tuple
        Parameter names: ("k",)
    """
    name = "power_law"
    param_names = ("k",)

    def __call__(self, redshift, **parameters):
        return self.prob_redshift(redshift, **parameters)

    def psi_z(self, redshift, *, k: float):
        """
        Redshift evolution function psi(z):
            psi(z) = (1 + z)^k

        Parameters
        ----------
        redshift : array_like
            Dimensionless redshift.
        k : float
            Power-law index.

        Returns
        -------
        ndarray
            Dimensionless evolution factor.
        """
        return (1.0 + np.asarray(redshift)) ** k

# Default instance
power_law_redshift = PowerLawRedshift()

class GRB2008SFR(BaseRedshiftEvolution):
    """
    Star formation rate (SFR) calibrated by high-redshift gamma-ray bursts
    (GRBs).

    Based on the model from GRB 2008 data, this model captures the observed
    SFR evolution with a flexible functional form.

    The SFR is modeled as:
        psi(z) =
        rho_local *
        [(1+z)^(3.4*eta)+((1+z)/5000)^(-0.3*eta)+((1+z)/9)^(-3.5*eta)]^(1/eta)

    with eta = -10 and rho_local = 0.02 Msun/yr/Mpc^3.

    Parameters
    ----------
    z_max : float, optional
        Maximum redshift. Default is 2.0.
    num_zbins : int, optional
        Number of redshift bins. Default is 1000.
    cosmology : astropy.cosmology.Cosmology, optional
        Cosmological model. Default is `get_cosmology()`.
    z_grid : array_like, optional
        Custom redshift grid.
    """
    name = "sfr_grb_2008"
    param_names = ()

    def __call__(self, redshift, **parameters):
        return self.prob_redshift(redshift, **parameters)

    def psi_z(self, redshift, **parameters):
        redshift = np.asarray(redshift)
        rho_local = 0.02  # Msolar/yr/Mpc^3
        eta = -10
        psiz = ((1+redshift)**(3.4*eta) + ((1+redshift)/5000)**(-0.3*eta) +\
                ((1+redshift)/9)**(-3.5*eta))**(1./eta)
        psiz *= rho_local
        return psiz


sfr_grb_2008_redshift = GRB2008SFR()


class MadauDickinson2014SFR(BaseRedshiftEvolution):
    """
    Star formation rate (SFR) from Madau & Dickinson (2014).

    The model is:
        psi(z) = 0.015 * (1+z)^gamma / [1 + ((1+z)/(1+z_peak))^kappa]

    This model captures the rise and decline of SFR with redshift.

    Parameters
    ----------
    z_max : float, optional
        Maximum redshift. Default is 2.0.
    num_zbins : int, optional
        Number of redshift bins. Default is 1000.
    cosmology : astropy.cosmology.Cosmology, optional
        Cosmological model. Default is `get_cosmology()`.
    z_grid : array_like, optional
        Custom redshift grid.
    """
    name = "sfr_madau_dickinson_2014"
    param_names = ("gamma", "kappa", "z_peak")

    def __call__(self, redshift, **parameters):
        return self.prob_redshift(redshift, **parameters)

    def psi_z(self, redshift, *, gamma=2.7, kappa=5.6, z_peak=1.9):
        redshift = np.asarray(redshift)
        return 0.015 * (1+redshift)**gamma / (1 + ((1+redshift)/(1+z_peak))**kappa)

sfr_madau_dickinson_2014_redshift = MadauDickinson2014SFR()


class MadauFragos2017SFR(BaseRedshiftEvolution):
    """
    Star formation rate (SFR) from Madau & Fragos (2017).

    The model is:
        psi(z) = k_imf * 0.015 * (1+z)^a / [1 + ((1+z)/b)^c]

    with parameters depending on the `mode` ('high' or 'low').

    Parameters
    ----------
    z_max : float, optional
        Maximum redshift. Default is 2.0.
    num_zbins : int, optional
        Number of redshift bins. Default is 1000.
    cosmology : astropy.cosmology.Cosmology, optional
        Cosmological model. Default is `get_cosmology()`.
    z_grid : array_like, optional
        Custom redshift grid.
    """
    name = "sfr_madau_fragos_2017"
    param_names = ("k_imf", "mode")

    _MODE_PARAMETERS = {
        "low": (2.6, 3.2, 6.2),
        "high": (2.7, 3.0, 5.35),
        }

    def __call__(self, redshift, **parameters):
        return self.prob_redshift(redshift, **parameters)

    def psi_z(self, redshift, *, k_imf=0.66, mode="high"):
        redshift = np.asarray(redshift)

        if mode in ["low", "high"]:
            factor_a, factor_b, factor_c = self._MODE_PARAMETERS[mode]
        else:
            raise ValueError("'mode' must choose from 'high' or 'low'.")
        return k_imf * 0.015 * (1+redshift)**factor_a / (1 + ((1+redshift)/factor_b)**factor_c)

sfr_madau_fragos_2017_redshift = MadauFragos2017SFR()

# ---------------------------------------------
# Classes related to delay-time distributions
# ---------------------------------------------

class BaseDelayTimeDistribution(ABC):
    """Probability distribution of astrophysical delay times.

    All numerical times are expressed in Gyr.
    The PDF is expressed in Gyr^-1.
    """

    @property
    @abstractmethod
    def min_delay(self) -> float:
        """Minimum delay time in Gyr."""
        pass

    @property
    @abstractmethod
    def max_delay(self) -> float:
        """Maximum delay time in Gyr."""
        pass

    @abstractmethod
    def pdf(self, tau: np.ndarray) -> np.ndarray:
        """Return the normalized delay-time PDF.

        Parameters
        ----------
        tau : array_like
            Delay times in Gyr.

        Returns
        -------
        ndarray
            Probability density at tau (dimensionless, in Gyr^-1).
        """
        pass

    def cdf(self, tau: np.ndarray) -> np.ndarray:
        """Compute the cumulative distribution function.

        Parameters
        ----------
        tau : array_like
            Delay times in Gyr.

        Returns
        -------
        ndarray
            Cumulative probability at tau.
        """
        tau = np.asarray(tau, dtype=float)
        # Use numerical integration for CDF
        tau_vals = np.linspace(self.min_delay, self.max_delay, 1000)
        pdf_vals = self.pdf(tau_vals)
        cdf_vals = np.cumsum(pdf_vals) * (tau_vals[1] - tau_vals[0])
        cdf_vals = cdf_vals / cdf_vals[-1]  # Normalize to 1

        # Interpolate to get CDF at requested tau values
        return np.interp(tau, tau_vals, cdf_vals, left=0.0, right=1.0)

    def sample(self, size: int = 1, random_state: np.random.RandomState = None) -> np.ndarray:
        """Sample delay times from the distribution.

        Parameters
        ----------
        size : int, optional
            Number of samples to generate. Default is 1.
        random_state : np.random.RandomState, optional
            Random state for reproducibility.

        Returns
        -------
        ndarray
            Sampled delay times in Gyr.
        """
        if random_state is None:
            random_state = np.random.RandomState()

        # Use inverse transform sampling
        u = random_state.uniform(0, 1, size)
        tau_vals = np.linspace(self.min_delay, self.max_delay, 1000)
        cdf_vals = self.cdf(tau_vals)

        return np.interp(u, cdf_vals, tau_vals)


class PowerLawDelayDistribution(BaseDelayTimeDistribution):
    """Power-law delay time distribution.

    The PDF is defined as:
        p(tau) ∝ tau^(-alpha) for tau_min <= tau <= tau_max
    """

    def __init__(self, td_min: float = 0.02, td_max: float = 13.0, alpha: float = 0.81):
        """
        Parameters
        ----------
        td_min : float, optional
            Minimum delay time in Gyr. Must be positive. Default is 0.02.
        td_max : float, optional
            Maximum delay time in Gyr. Must be greater than td_min. Default is 13.0.
        alpha : float, optional
            Power-law index. Default is 0.81.
        """
        if not np.isfinite(td_min) or td_min <= 0:
            raise ValueError("td_min must be finite and positive")
        if not np.isfinite(td_max) or td_max <= td_min:
            raise ValueError("td_max must be greater than td_min")
        if not np.isfinite(alpha):
            raise ValueError("alpha must be finite")

        self._td_min = float(td_min)
        self._td_max = float(td_max)
        self._alpha = float(alpha)

        # Calculate normalization constant
        if np.isclose(self._alpha, 1.0):
            self._norm = np.log(self._td_max / self._td_min)
        else:
            self._norm = (self._td_max**(1-self._alpha) - self._td_min**(1-self._alpha)) / (1-self._alpha)

    @property
    def min_delay(self) -> float:
        return self._td_min

    @property
    def max_delay(self) -> float:
        return self._td_max

    def pdf(self, tau: np.ndarray) -> np.ndarray:
        tau = np.asarray(tau, dtype=float)
        result = np.zeros_like(tau, dtype=float)

        valid = (tau >= self.min_delay) & (tau <= self.max_delay)
        result[valid] = tau[valid]**(-self._alpha) / self._norm

        return result


class InverseDelayDistribution(BaseDelayTimeDistribution):
    """Inverse power-law delay time distribution (special case of PowerLawDelayDistribution with alpha=1).

    The PDF is defined as:
        p(tau) = 1 / (tau * ln(tau_max/tau_min)) for tau_min <= tau <= tau_max
    """

    def __init__(self, td_min: float = 0.02, td_max: float = 13.0):
        """
        Parameters
        ----------
        td_min : float, optional
            Minimum delay time in Gyr. Must be positive. Default is 0.02.
        td_max : float, optional
            Maximum delay time in Gyr. Must be greater than td_min. Default is 13.0.
        """
        if not np.isfinite(td_min) or td_min <= 0:
            raise ValueError("td_min must be finite and positive")
        if not np.isfinite(td_max) or td_max <= td_min:
            raise ValueError("td_max must be greater than td_min")

        self._td_min = float(td_min)
        self._td_max = float(td_max)
        self._norm = np.log(td_max / td_min)

    @property
    def min_delay(self) -> float:
        return self._td_min

    @property
    def max_delay(self) -> float:
        return self._td_max

    def pdf(self, tau: np.ndarray) -> np.ndarray:
        tau = np.asarray(tau, dtype=float)
        result = np.zeros_like(tau, dtype=float)

        valid = (tau >= self.min_delay) & (tau <= self.max_delay)
        result[valid] = 1.0 / (tau[valid] * self._norm)

        return result


class LogNormalDelayDistribution(BaseDelayTimeDistribution):
    """Log-normal delay time distribution.

    The PDF is defined as:
        p(tau) = 1 / (tau * sigma * sqrt(2*pi)) * exp(-(ln(tau) - ln(t_peak))^2 / (2*sigma^2))
    """

    def __init__(self, td_min: float = 0.02, td_max: float = 13.0, t_peak: float = 2.9, sigma: float = 0.2):
        """
        Parameters
        ----------
        td_min : float, optional
            Minimum delay time in Gyr. Must be positive. Default is 0.02.
        td_max : float, optional
            Maximum delay time in Gyr. Must be greater than td_min. Default is 13.0.
        t_peak : float, optional
            Peak of the distribution in Gyr. Default is 2.9.
        sigma : float, optional
            Standard deviation of the log-normal distribution. Default is 0.2.
        """
        if not np.isfinite(td_min) or td_min <= 0:
            raise ValueError("td_min must be finite and positive")
        if not np.isfinite(td_max) or td_max <= td_min:
            raise ValueError("td_max must be greater than td_min")
        if not np.isfinite(t_peak) or t_peak <= 0:
            raise ValueError("t_peak must be finite and positive")
        if not np.isfinite(sigma) or sigma <= 0:
            raise ValueError("sigma must be finite and positive")

        self._td_min = float(td_min)
        self._td_max = float(td_max)
        self._t_peak = float(t_peak)
        self._sigma = float(sigma)

        # Calculate normalization constant using numerical integration
        tau_vals = np.logspace(np.log10(td_min), np.log10(td_max), 1000)
        unnormalized_pdf = 1.0 / (tau_vals * self._sigma * np.sqrt(2*np.pi)) * \
                          np.exp(-(np.log(tau_vals) - np.log(self._t_peak))**2 / (2*self._sigma**2))
        self._norm = np.trapz(unnormalized_pdf, tau_vals)

    @property
    def min_delay(self) -> float:
        return self._td_min

    @property
    def max_delay(self) -> float:
        return self._td_max

    def pdf(self, tau: np.ndarray) -> np.ndarray:
        tau = np.asarray(tau, dtype=float)
        result = np.zeros_like(tau, dtype=float)

        valid = (tau >= self.min_delay) & (tau <= self.max_delay)
        result[valid] = (1.0 / (tau[valid] * self._sigma * np.sqrt(2*np.pi)) *
                        np.exp(-(np.log(tau[valid]) - np.log(self._t_peak))**2 / (2*self._sigma**2))) / self._norm

        return result


class GaussianDelayDistribution(BaseDelayTimeDistribution):
    """Gaussian delay time distribution.

    The PDF is defined as:
        p(tau) = 1 / (sigma * sqrt(2*pi)) * exp(-(tau - t_mean)^2 / (2*sigma^2))
    """

    def __init__(self, td_min: float = 0.02, td_max: float = 13.0, t_mean: float = 2.0, sigma: float = 0.3):
        """
        Parameters
        ----------
        td_min : float, optional
            Minimum delay time in Gyr. Must be positive. Default is 0.02.
        td_max : float, optional
            Maximum delay time in Gyr. Must be greater than td_min. Default is 13.0.
        t_mean : float, optional
            Mean of the Gaussian distribution in Gyr. Default is 2.0.
        sigma : float, optional
            Standard deviation of the Gaussian distribution. Default is 0.3.
        """
        if not np.isfinite(td_min) or td_min <= 0:
            raise ValueError("td_min must be finite and positive")
        if not np.isfinite(td_max) or td_max <= td_min:
            raise ValueError("td_max must be greater than td_min")
        if not np.isfinite(t_mean) or t_mean <= 0:
            raise ValueError("t_mean must be finite and positive")
        if not np.isfinite(sigma) or sigma <= 0:
            raise ValueError("sigma must be finite and positive")

        self._td_min = float(td_min)
        self._td_max = float(td_max)
        self._t_mean = float(t_mean)
        self._sigma = float(sigma)

        # Calculate normalization constant using numerical integration
        tau_vals = np.linspace(td_min, td_max, 1000)
        unnormalized_pdf = np.exp(-(tau_vals - self._t_mean)**2 / (2*self._sigma**2))
        self._norm = np.trapz(unnormalized_pdf, tau_vals)

    @property
    def min_delay(self) -> float:
        return self._td_min

    @property
    def max_delay(self) -> float:
        return self._td_max

    def pdf(self, tau: np.ndarray) -> np.ndarray:
        tau = np.asarray(tau, dtype=float)
        result = np.zeros_like(tau, dtype=float)

        valid = (tau >= self.min_delay) & (tau <= self.max_delay)
        result[valid] = np.exp(-(tau[valid] - self._t_mean)**2 / (2*self._sigma**2)) / (self._sigma * np.sqrt(2*np.pi) * self._norm)

        return result

class SFRTimeDelayRedshift(BaseRedshiftEvolution):
    r"""
    Redshift evolution from convolving a star formation rate (SFR) with a time
    delay distribution.

    This model accounts for the fact that events (e.g., compact binary mergers)
    occur after star formation, with a time delay distribution.

    The evolution is computed as:
        psi(z) = \int_{z}^{z_max} SFR(z_f; θ_sfr) * p(t_d; θ_td) * dt_f/dz_f dz_f

    where t_d = t_lookback(z_f) - t_lookback(z), and θ_sfr, θ_td are parameters
    for the SFR and time delay models respectively.

    Parameters
    ----------
    sfr_model : BaseRedshiftEvolution
        SFR model (e.g., MadauDickinson2014SFR).
    time_delay_model : DelayTimeDistribution
        Time delay distribution (e.g., PowerLawDelayDistribution).
    z_max : float, optional
        Maximum redshift for integration. Default is 10.0.
    num_zbins : int, optional
        Number of redshift bins. Default is 1000.
    cosmology : astropy.cosmology.Cosmology, optional
        Cosmological model. Default is `get_cosmology()`.
    z_grid : array_like, optional
        Custom redshift grid.
    z_formation_max : float, optional
        Maximum formation redshift for SFR grid. Default is 20.0.
    """
    name = "sfr_time_delay"
    param_names = ()  # Will be populated dynamically based on sub-models

    def __init__(
        self,
        sfr_model,
        time_delay_model,
        z_max=10.0,
        num_zbins=1000,
        cosmology=None,
        z_grid=None,
        z_formation_max=20.0,
    ):
        # Validate inputs
        if not isinstance(sfr_model, BaseRedshiftEvolution):
            raise TypeError("sfr_model must be an instance of BaseRedshiftEvolution")
        if not isinstance(time_delay_model, BaseDelayTimeDistribution):
            raise TypeError("time_delay_model must be an instance of DelayTimeDistribution")

        super().__init__(z_max, num_zbins, cosmology, z_grid=z_grid)
        self.sfr_model = sfr_model
        self.time_delay_model = time_delay_model

        # Validate z_formation_max
        if z_formation_max <= z_max:
            raise ValueError("z_formation_max must be greater than z_max")
        self.z_formation_max = z_formation_max

        # Precompute lookback times for formation redshift grid
        self._zf_grid = np.linspace(0, self.z_formation_max, 5000)
        self._tf_grid = self.cosmology.lookback_time(self._zf_grid).to(u.Gyr).value

        # dt/dz = 1 / (H(z) * (1+z))
        H_z = self.cosmology.H(self._zf_grid).to(1/u.Gyr).value
        self._dt_dz_f = 1.0 / (H_z * (1.0 + self._zf_grid))

        # Extract parameter names from sub-models
        self._sfr_param_names = getattr(sfr_model, 'param_names', ())
        self._td_param_names = ()  # Delay distributions typically don't have parameters

        # Combine parameter names
        self.param_names = self._sfr_param_names + self._td_param_names

        # Cache for precomputed values
        self._psi_z_cache = {}
        self._last_sfr_params = None
        self._last_td_params = None

    def _update_psi_z_grid(self, sfr_params=None, td_params=None):
        """
        Evaluate and cache the convolution over the redshift grid.
        Uses numerical integration for all delay time models.

        Parameters
        ----------
        sfr_params : dict, optional
            Parameters for the SFR model.
        td_params : dict, optional
            Parameters for the time delay model (if supported).
        """
        # Use default parameters if not provided
        if sfr_params is None:
            sfr_params = {}
        if td_params is None:
            td_params = {}

        # Check if we need to recompute (cache invalidation)
        current_sfr_params = tuple(sorted(sfr_params.items())) if sfr_params else ()
        current_td_params = tuple(sorted(td_params.items())) if td_params else ()

        cache_key = (current_sfr_params, current_td_params)
        if cache_key in self._psi_z_cache:
            self._psi_z_grid = self._psi_z_cache[cache_key]
            return

        self._psi_z_grid = np.zeros(len(self._z_grid))

        # Pre-compute lookback times for z_grid to vectorize delay calculation
        tm_grid = self.cosmology.lookback_time(self._z_grid).to(u.Gyr).value

        # Create splines for efficient interpolation
        # Evaluate SFR on the grid with current parameters
        sfr_values = self.sfr_model.psi_z(self._zf_grid, **sfr_params)
        sfr_spline = CubicSpline(self._zf_grid, sfr_values, extrapolate=True)
        dt_dz_spline = CubicSpline(self._zf_grid, self._dt_dz_f, extrapolate=True)
        tf_spline = CubicSpline(self._zf_grid, self._tf_grid, extrapolate=True)

        for i, (zm, tm) in enumerate(zip(self._z_grid, tm_grid, strict=True)):
            # Define the integrand function
            def integrand(zf):
                td = tf_spline(zf) - tm
                # Time delay must be positive and within the model's range
                if td <= 0 or td < self.time_delay_model.min_delay or td > self.time_delay_model.max_delay:
                    return 0.0
                # Evaluate time delay PDF with parameters if supported
                try:
                    p_td = self.time_delay_model.pdf(td, **td_params)
                except TypeError:
                    # If pdf doesn't accept **kwargs, call without parameters
                    p_td = self.time_delay_model.pdf(td)
                return sfr_spline(zf) * p_td * dt_dz_spline(zf)

            # Set integration limits
            z_start = zm
            # Find the redshift where td = td_max
            try:
                # Use root finding to find zf where tf(zf) = tm + td_max
                from scipy.optimize import brentq
                def td_func(zf):
                    return tf_spline(zf) - tm - self.time_delay_model.max_delay
                # Check if root exists within our grid
                if td_func(self.z_formation_max) < 0:
                    # td_max is larger than available lookback time
                    z_end = self.z_formation_max
                else:
                    z_end = brentq(td_func, zm, self.z_formation_max)
            except (ValueError, RuntimeError):
                # If root finding fails, use the maximum redshift
                z_end = self.z_formation_max

            # Perform numerical integration
            if z_end > z_start:
                # Add points near the lower limit for better resolution
                points = [z_start + 1e-8, z_start + 1e-6, z_start + 1e-4, z_start + 1e-2]
                points = [p for p in points if p < z_end]

                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    result = scipy_integrate.quad(
                        integrand, z_start, z_end, points=points, limit=1000
                    )[0]
                    self._psi_z_grid[i] = result
            else:
                self._psi_z_grid[i] = 0.0

        # Cache the result
        self._psi_z_cache[cache_key] = self._psi_z_grid.copy()

    def psi_z(self, redshift, **parameters):
        """
        Redshift evolution function from convolution.

        Parameters
        ----------
        redshift : array_like
            Dimensionless redshift.
        **parameters
            Parameters for the SFR model and time delay model.
            Parameters are passed to the respective sub-models based on their names.

        Returns
        -------
        ndarray
            Dimensionless evolution factor.
        """
        redshift = np.asarray(redshift)

        # Separate SFR and time delay parameters
        sfr_params = {}
        td_params = {}

        for key, value in parameters.items():
            if key in self._sfr_param_names:
                sfr_params[key] = value
            elif key in self._td_param_names:
                td_params[key] = value
            else:
                # Ignore unknown parameters with a warning
                warnings.warn(f"Unknown parameter '{key}' ignored in SFRTimeDelayRedshift.psi_z")

        # Update the psi(z) grid if parameters have changed
        current_sfr_params = tuple(sorted(sfr_params.items())) if sfr_params else ()
        current_td_params = tuple(sorted(td_params.items())) if td_params else ()

        if (current_sfr_params != self._last_sfr_params or
            current_td_params != self._last_td_params):
            self._update_psi_z_grid(sfr_params, td_params)
            self._last_sfr_params = current_sfr_params
            self._last_td_params = current_td_params

        return np.interp(redshift, self._z_grid, self._psi_z_grid, right=0.0)

    def __call__(self, redshift, **parameters):
        return self.prob_redshift(redshift, **parameters)

__all__ = ['PowerLawRedshift', 'power_law_redshift', 'GRB2008SFR',
'sfr_grb_2008_redshift', 'MadauDickinson2014SFR', 'sfr_madau_dickinson_2014_redshift',
'MadauFragos2017SFR', 'sfr_madau_fragos_2017_redshift', 'SFRTimeDelayRedshift'
]



