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

import numpy as np
import warnings
import astropy.units as u
from pycbc.cosmology import get_cosmology
from abc import ABC, abstractmethod

class BaseRedshiftEvolution(ABC):
    """
    Abstract Base class for redshift evolution models with Astropy unit support.

    This class provides a framework for modeling the redshift dependence of astrophysical
    event rates, including differential comoving volume, spacetime volume, and redshift
    probability distributions. It supports various cosmological models and integrates
    with time delay distributions.

    The core functionality includes:
    - Differential comoving volume (3-volume): dVc/dz
    - Differential spacetime volume (4-volume): dVT/dz = psi(z)/(1+z) * dVc/dz
    - Redshift probability density functions (PDFs)
    - Normalization and total spacetime volume computation


    Notes
    -----
    - Redshift `z` is dimensionless.
    - Comoving volumes are returned in units of Gpc³.
    - Spacetime volumes are returned in Gpc³ (with implicit time integration).
    - The class uses numerical integration over a redshift grid for performance.

    Parameters
    ----------
    z_max : float, optional
        Maximum redshift for integration and normalization. Default is 2.0.
    num_zbins : int, optional
        Number of redshift bins for numerical integration. Must be at least 2. Default is 1000.
    cosmology : astropy.cosmology.Cosmology, optional
        Cosmological model to use. If None, defaults to `pycbc.cosmology.get_cosmology()`.
    z_grid : array_like, optional
        Custom redshift grid (must be strictly increasing, start at 0, end at z_max).
        If None, a uniform grid from 0 to z_max is used.

    Attributes
    ----------
    zmax : float
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
    SFRTimeDelayRedshift : Concrete implementations of redshift evolution models.
    """

    def __init__(self, z_max: float = 2.0, num_zbins: int = 1000, 
                 cosmology=None, z_grid = None):
        """
        Initialize the redshift evolution model.

        Parameters
        ----------
        z_max : float
            Maximum redshift for normalization and integration. Must be positive.
        num_zbins : int
            Number of redshift bins for numerical integration. Must be at least 2.
        cosmology : astropy.cosmology.Cosmology, optional
            Cosmological model to use. If None, defaults to `get_cosmology()`.
        z_grid : array_like, optional
            Custom redshift grid (must be strictly increasing, start at 0, end at z_max).
            If None, a uniform grid from 0 to z_max is used.
        """
        self.zmax = float(zmax)
        self.num_zbins = int(num_zbins)

        # Validate z_max, num_bin, and cosmology
        if z_max <= 0:
            raise ValueError("z_max must be positive")
        self.z_max = float(z_max)
        
        if num_zbins < 2:
            raise ValueError("num_zbins must be at least 2")
        self.num_zbins = int(num_zbins)

        if cosmology is not None:
            self.cosmology = cosmology
        else:
            self.cosmology = get_cosmology()

        # Redshift grid (dimensionless) validation
        if z_grid is None:
            #self._z_grid = np.linspace(1e-3, self.z_max, num_zbins)
            self._z_grid = np.linspace(0.0, self.z_max, num_zbins)
        else:
            self._z_grid = np.asarray(z_grid)
            if self._z_grid.ndim != 1:
                raise ValueError("z_grid must be one-dimensional")

            if np.any(self._z_grid < 0):
                raise ValueError("z_grid must be non-negative")

            if np.any(np.diff(self._z_grid) <= 0):
                raise ValueError("z_grid must be strictly increasing")

            if self._z_grid[0] != 0:
                raise ValueError("z_grid should start at z=0")

            if self._z_grid[-1] != self.z_max:
                raise ValueError("z_grid must end at z_max")

        # Differential comoving volume on grid
        # Full-sky differential comoving volume dVc/dz = 4pi dVc/dz/dOmega
        self._dvc_dz_grid = (
            4.0
            * np.pi
            * self.cosmology.differential_comoving_volume(self._z_grid)
            .to(u.Gpc**3 / u.sr)
        )

        # Cache for interpolated values
        self._cached_z = None
        self._cached_dvc_dz = None

        # Time delay boundaries (in Gyr)
        self.td_max = self.cosmology.age(0).to(u.Gyr).value
        self.td_min = 0.02

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

        # Validate for positive redshift values
        if np.any(z < 0):
            raise ValueError("Redshift must be non-negative")

        return z, is_scalar

    # ------------------------------------------------------------------
    # Cosmology utilities
    # ------------------------------------------------------------------

    def dVc_dz(self, redshift):
        """
        Differential comoving 3-volume per unit redshift.

        This is the full-sky differential comoving volume: dVc/dz = 4 \pi dVc/dz/d_\Omega.

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
        )

        return dvc[0] if is_scalar else dvc

    def comoving_3volume(self, redshift):
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

        z, is_scalar = self._to_1d_array(redshift)
        vol = self.cosmology.comoving_volume(z).to(u.Gpc**3)

        return vol[0] if is_scalar else vol

    # ------------------------------------------------------------------
    # Redshift evolution model
    # ------------------------------------------------------------------

    @abstractmethod
    def psi_z(self, redshift, **parameters):
        """
        Redshift evolution function psi(z).

        This function defines the intrinsic redshift dependence of the event rate.

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

    def time_delay_prob(self, tau, td_model):
        """
         Compute the probability density of a time delay tau.

        This method supports several standard time delay models used in astrophysics.

        Parameters
        ----------
        tau : array_like
            Time delay(s) in Gyr.
        td_model : str
            Name of the time delay model. Must be one of:
            - 'log_normal': Log-normal distribution
            - 'gaussian': Gaussian distribution
            - 'power_law': Power-law distribution
            - 'inverse': Inverse power-law with bounded support

        Returns
        -------
        numpy.ndarray
            Probability density at tau (dimensionless).

        Raises
        ------
        ValueError
            If td_model is not one of the supported types.

        Notes
        -----
        - For 'power_law' and 'inverse', normalization is approximate and may require
          careful handling near zero.
        - 'inverse' model uses a bounded support from td_min to td_max.
        """
        tau = np.asarray(tau, dtype=float)
        p_t = np.zeros_like(tau)
        
        if td_model == "log_normal":
            t_ln = 2.9  # Gyr
            sigma_ln = 0.2
            p_t = np.exp(-(np.log(tau)-np.log(t_ln))**2/(2*sigma_ln**2)) / (np.sqrt(2*np.pi)*sigma_ln)
        elif td_model == "gaussian":
            t_g = 2  # Gyr
            sigma_g = 0.3
            p_t = np.exp(-(tau-t_g)**2/(2*sigma_g**2)) / (np.sqrt(2*np.pi)*sigma_g)
        elif td_model == "power_law":
            alpha_t = 0.81
            valid = (tau > 0)
            p_t[valid] = tau[valid]**(-alpha_t)
        elif td_model == "inverse":
            norm_const = 1/np.log(self.td_max/self.td_min)
            valid = (tau >= self.td_min) & (tau <= self.td_max)
            p_t[valid] = norm_const * tau[valid]**(-0.999)
        else:
            raise ValueError(f"Currently, 'td_model' must choose from: "
                             f"['log_normal', 'gaussian', 'power_law', 'inverse'].")
        return p_t

    def _cache_dvc_dz(self, redshift):
        """
        Cache interpolated differential comoving volume for faster repeated evaluation.

        This method caches the interpolated dVc/dz values on a given redshift grid.

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
        - The cache is updated only if the input redshift differs from the cached one.
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

        PThis represents the effective rate of events per unit redshift, accounting for
        cosmological volume and redshift evolution.

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

        This is the product of the normalization constant and the observation time.

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
        pdf[valid] = (self.dVT_dz(z[valid], **parameters) / norm).decompose().value
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
    Star formation rate (SFR) calibrated by high-redshift gamma-ray bursts (GRBs).

    Based on the model from GRB 2008 data, this model captures the observed SFR evolution
    with a flexible functional form.

    The SFR is modeled as:
        psi(z) = rho_local * [ (1+z)^(3.4*eta) + ((1+z)/5000)^(-0.3*eta) + ((1+z)/9)^(-3.5*eta) ]^(1/eta)

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
        return rho_local*((1+redshift)**(3.4*eta) + ((1+redshift)/5000)**(-0.3*eta) +
                       ((1+redshift)/9)**(-3.5*eta))**(1./eta)

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

    def __call__(self, redshift, **parameters):
        return self.prob_redshift(redshift, **parameters)

    def psi_z(self, redshift, *, k_imf=0.66, mode='high'):
        redshift = np.asarray(redshift)
        if mode == 'low':
            factor_a = 2.6
            factor_b = 3.2
            factor_c = 6.2
        elif mode == 'high':
            factor_a = 2.7
            factor_b = 3.0
            factor_c = 5.35
        else:
            raise ValueError("'mode' must choose from 'high' or 'low'.")
        return k_imf * 0.015 * (1+redshift)**factor_a / (1 + ((1+redshift)/factor_b)**factor_c)

sfr_madau_fragos_2017_redshift = MadauFragos2017SFR()


class SFRTimeDelayRedshift(BaseRedshiftEvolution):
    """
    Redshift evolution from convolving a star formation rate (SFR) with a time delay distribution.

    This model accounts for the fact that events (e.g., compact binary mergers) occur after
    star formation, with a time delay distribution.

    The evolution is computed as:
        psi(z) = \int_{z}^{z_max} SFR(z_f) * p(t_d) * dt_f/dz_f dz_f

    where t_d = t_lookback(z_f) - t_lookback(z).

    Parameters
    ----------
    sfr_model : BaseRedshiftEvolution
        SFR model (e.g., MadauDickinson2014SFR).
    td_model : str
        Time delay model: 'log_normal', 'gaussian', 'power_law', 'inverse'.
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
    **kwargs
        Additional arguments:
        - td_min : float, optional (default: 0.02 Gyr)
        - td_max : float, optional (default: lookback time at z_formation_max)
    """
    name = "sfr_time_delay"
    param_names = ()

    def __init__(self, sfr_model, td_model, z_max=10.0, num_zbins=1000, 
                 cosmology=None, z_grid=None, z_formation_max=20.0, **kwargs):
        super().__init__(z_max, num_zbins, cosmology, z_grid)
        self.sfr_model = sfr_model
        self.td_model = td_model
        self.z_formation_max = z_formation_max

        #from astropy.cosmology import Planck18
        #import astropy.units as u
        #self.cosmology = cosmology if cosmology is not None else Planck18
        
        # Define boundaries for time delay
        self.td_min = kwargs.get('td_min', 0.02)  # Gyr (20 Myr)
        self.td_max = kwargs.get('td_max', self.cosmology.lookback_time(self.z_formation_max).to(u.Gyr).value)
        
        #from astropy.cosmology import Planck18
        #import astropy.units as u
        #self.cosmology = cosmology if cosmology is not None else Planck18
        
        # Define boundaries for time delay
        self.td_min = kwargs.get('td_min', 0.02)  # Gyr (20 Myr)
        self.td_max = kwargs.get('td_max', self.cosmology.lookback_time(self.z_formation_max).to(u.Gyr).value)
        
        # Precompute lookback times for formation redshift grid
        # 5000 points is enough because adaptive quad handles sharp features perfectly
        self._zf_grid = np.linspace(0, self.z_formation_max, 5000)
        self._tf_grid = self.cosmology.lookback_time(self._zf_grid).to(u.Gyr).value
        
        # dt/dz = 1 / (H(z) * (1+z))
        H_z = self.cosmology.H(self._zf_grid).to(1/u.Gyr).value
        self._dt_dz_f = 1.0 / (H_z * (1.0 + self._zf_grid))
        
        # Evaluate SFR on the grid
        self._sfr_f = self.sfr_model(self._zf_grid)
        
        self._update_psi_z_grid()

    def _update_psi_z_grid(self):
        """
        Evaluate and cache the convolution over the redshift grid.
        Uses fast numeric grid integration for bounded models, and exact scipy quad 
        for models with integrable singularities (power_law) or discontinuities (inverse).
        """
        self._psi_z_grid = np.zeros(len(self._z_grid))
        
        # Pre-compute time grid for z_grid to vectorize delay calculation
        tm_grid = self.cosmology.lookback_time(self._z_grid).to(u.Gyr).value
        
        if self.td_model in ["power_law", "inverse"]:
            import scipy.integrate as scipy_integrate
            from scipy.interpolate import CubicSpline
            import warnings
            
            # Use CubicSpline to ensure C2 continuous derivatives for quad convergence
            sfr_spline = CubicSpline(self._zf_grid, self._sfr_f, extrapolate=True)
            dt_dz_spline = CubicSpline(self._zf_grid, self._dt_dz_f, extrapolate=True)
            tf_spline = CubicSpline(self._zf_grid, self._tf_grid, extrapolate=True)
            
            if self.td_model == "inverse":
                # Ensure strictly monotonic for inverse mapping
                valid_idx = np.argsort(self._tf_grid)
                z_of_t_spline = CubicSpline(self._tf_grid[valid_idx], self._zf_grid[valid_idx], extrapolate=True)
            
            for i, zm in enumerate(self._z_grid):
                # Pin the singularity mathematically perfectly to the zm boundary
                tm_local = tf_spline(zm)
                
                z_start = zm
                if self.td_model == "inverse":
                    t_start = tm_local + self.td_min
                    if t_start >= self._tf_grid[-1]:
                        self._psi_z_grid[i] = 0.0
                        continue
                    z_start = max(zm, float(z_of_t_spline(t_start)))
                
                def integrand(zf):
                    td = tf_spline(zf) - tm_local
                    if td <= 0:
                        return 0.0
                    p_td = float(self.time_delay_prob(td, self.td_model))
                    return sfr_spline(zf) * p_td * dt_dz_spline(zf)
                
                # Guide quad to densely sample the incredibly narrow peak near z_start
                # when td_min is extremely small (e.g. 1e-10) to prevent missing the peak entirely.
                pts = [z_start + 1e-8, z_start + 1e-6, z_start + 1e-4, z_start + 1e-2]
                pts = [p for p in pts if p < self.z_formation_max]
                
                # Catch any minor roundoff warnings from quad to keep terminal clean
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    self._psi_z_grid[i] = scipy_integrate.quad(
                        integrand, z_start, self.z_formation_max, points=pts, limit=1000
                    )[0]
        else:
            for i, (zm, tm) in enumerate(zip(self._z_grid, tm_grid)):
                valid = self._zf_grid >= zm
                zf_valid = self._zf_grid[valid]
                tf_valid = self._tf_grid[valid]
                td = tf_valid - tm
                
                p_td = self.time_delay_prob(td, self.td_model)
                integrand = self._sfr_f[valid] * p_td * self._dt_dz_f[valid]
                self._psi_z_grid[i] = np.trapz(integrand, zf_valid)

    def psi_z(self, redshift, **parameters):
        """
        Redshift evolution function from convolution.
        
        Parameters
        ----------
        redshift : array_like
            Dimensionless redshift.

        Returns
        -------
        ndarray
            Dimensionless evolution factor.
        """
        redshift = np.asarray(redshift)
        return np.interp(redshift, self._z_grid, self._psi_z_grid, right=0.0)

    def __call__(self, redshift, **parameters):
        return self.prob_redshift(redshift, **parameters)

__all__ = ['PowerLawRedshift', 'power_law_redshift', 'GRB2008SFR',
'sfr_grb_2008_redshift', 'MadauDickinson2014SFR', 'sfr_madau_dickinson_2014_redshift',
'MadauFragos2017SFR', 'sfr_madau_fragos_2017_redshift', 'SFRTimeDelayRedshift'
]

 

