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
    Base class for redshift evolution models with Astropy unit support.

    This class provides:
    - Differential comoving volume (3-volume)
    - Differential spacetime volume (4-volume)
    - Redshift probability distributions

    Cosmology is fixed to Planck15.

    Notes
    -----
    - Redshift `z` is dimensionless.
    - Comoving volumes are returned in units of Gpc^3.
    - Spacetime volumes are returned in Gpc^3.
    """

    def __init__(self, zmax: float = 2.0, num_zbins: int = 1000, 
                 cosmology=None, z_grid = None):
        """
        Parameters
        ----------
        zmax : float
            Maximum redshift for normalization and integration.
        num_zbins : int
            Number of redshift bins for numerical integration.
        """
        self.zmax = float(zmax)
        self.num_zbins = int(num_zbins)
        if cosmology is not None:
            self.cosmology = cosmology
        else:
            self.cosmology = get_cosmology()

        # Redshift grid (dimensionless)
        if z_grid is None:
            self._z_grid = np.linspace(1e-3, self.zmax, num_zbins)
        else:
            self._z_grid = np.asarray(z_grid)

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

        Returns
        -------
        z_array : ndarray
        is_scalar : bool
        """
        if np.isscalar(redshift):
            return np.array([redshift], dtype=float), True
        z = np.asarray(redshift, dtype=float)
        if np.any(z < 0):
            raise ValueError("Redshift must be non-negative")
        return z.ravel(), False

    # ------------------------------------------------------------------
    # Cosmology utilities
    # ------------------------------------------------------------------

    def dV3cdz(self, redshift):
        """
        Differential comoving 3-volume per unit redshift.

        Parameters
        ----------
        redshift : array_like
            Dimensionless redshift(s).

        Returns
        -------
        astropy.units.Quantity
            dVc/dz with units of Gpc^3.
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
        pass

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def time_delay_prob(self, tau, td_model):
        """
        Pure numpy array evaluation of time delay probability distributions.
        
        Parameters
        ----------
        tau : array_like
            Time delay in Gyr.
        td_model : str
            The name of time delay model.
            
        Returns
        -------
        ndarray
            Probability density at tau.
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
            raise ValueError(f"'td_model' must choose from "
                             f"['log_normal', 'gaussian', 'power_law', 'inverse'].")
        return p_t

    def _cache_dvc_dz(self, redshift):
        """
        Cache interpolated differential comoving volume.

        Parameters
        ----------
        redshift : ndarray
            Dimensionless redshift array.
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

    def dV4cdz(self, redshift, **parameters):
        """
        Differential spacetime (4-volume) element.

        Defined as:
            dV4/dz = psi(z) / (1 + z) * dVc/dz

        Parameters
        ----------
        redshift : array_like
            Dimensionless redshift.
        **parameters
            Parameters passed to `psi_z`.

        Returns
        -------
        astropy.units.Quantity
            Differential spacetime volume with units of Gpc^3.
        """

        z, is_scalar = self._to_1d_array(redshift)
        psi = self.psi_z(z, **parameters)

        # Warn if any redshift exceeds zmax
        invalid = z > self.zmax
        if np.any(invalid):
            warnings.warn(
                f"Some redshift values exceed zmax={self.zmax}. "
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

        Parameters
        ----------
        parameters : dict
            Parameters passed to `psi_z`.

        Returns
        -------
        astropy.units.Quantity
            Normalization constant with units of Gpc^3.
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
        """
        if not analysis_time.unit.is_equivalent(u.s):
            raise u.UnitTypeError("analysis_time must be a time quantity")

        return self.normalize(parameters) * analysis_time

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

        # Identify redshifts beyond zmax
        # This is for the normalization
        invalid = z > self.zmax
        if np.any(invalid):
            warnings.warn(
                f"Some redshift values exceed zmax={self.zmax}. "
                "These values will be assigned zero probability.",
                UserWarning
            )

        # Only valid redshifts contribute
        valid = z <= self.zmax
        norm = self.normalize(**parameters)
        pdf[valid] = (self.dV4cdz(z[valid], **parameters) / norm).decompose().value
        return pdf[0] if is_scalar else pdf



class PowerLawRedshift(BaseRedshiftEvolution):
    """
    Power-law redshift distribution model.

    psi(z) = (1 + z)^k
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


class GRB2008Redshift(BaseRedshiftEvolution):
    """
    The star formation rate (SFR) calibrated by high-z GRBs data.
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

sfr_grb_2008_redshift = GRB2008Redshift()


class MadauDickinson2014Redshift(BaseRedshiftEvolution):
    """
    The madau-dickinson 2014 star formation rate (SFR).
    """
    name = "sfr_madau_dickinson_2014"
    param_names = ("gamma", "kappa", "z_peak")

    def __call__(self, redshift, **parameters):
        return self.prob_redshift(redshift, **parameters)

    def psi_z(self, redshift, *, gamma=2.7, kappa=5.6, z_peak=1.9):
        redshift = np.asarray(redshift)
        return 0.015 * (1+redshift)**gamma / (1 + ((1+redshift)/(1+z_peak))**kappa)

sfr_madau_dickinson_2014_redshift = MadauDickinson2014Redshift()


class MadauFragos2017Redshift(BaseRedshiftEvolution):
    """
    The madau-fragos 2017 star formation rate (SFR).
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

sfr_madau_fragos_2017_redshift = MadauFragos2017Redshift()


class SFRTimeDelayRedshift(BaseRedshiftEvolution):
    """
    Redshift evolution from convoluting a star formation rate (SFR) with a time delay distribution.
    """
    name = "sfr_time_delay"
    param_names = ()

    def __init__(self, sfr_model, td_model, zmax=10.0, num_zbins=1000, 
                 cosmology=None, z_grid=None, z_formation_max=20.0, **kwargs):
        super().__init__(zmax, num_zbins, cosmology, z_grid)
        self.sfr_model = sfr_model
        self.td_model = td_model
        self.z_formation_max = z_formation_max
        
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

