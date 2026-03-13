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
import astropy.units as u
from astropy.cosmology import Planck15, Planck18
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
                 cosmology=Planck15, z_grid = None):
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
        self.cosmology = cosmology

        # Redshift grid (dimensionless)
        if z_grid is None:
            self._z_grid = np.linspace(1e-3, zmax, num_zbins)
        else:
            self._z_grid = np.asarray(z_grid)
        #self._z_grid = np.linspace(1e-3, self.zmax, self.num_zbins)

        # Differential comoving volume on grid
        # Full-sky differential comoving volume dVc/dz = 4pi dVc/dz/dOmega
        self._dvc_dz_grid = (
            4.0
            * np.pi
            * Planck15.differential_comoving_volume(self._z_grid)
            .to(u.Gpc**3 / u.sr)
        )

        # Cache for interpolated values
        self._cached_z = None
        self._cached_dvc_dz = None

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



class PowerLawRedshift(RedshiftEvolution):
    """
    Power-law redshift distribution model.

    ψ(z) = (1 + z)^k
    """

    def __call__(self, redshift, parameters):
        return self.prob_redshift(redshift, parameters)


# Default instance
power_law_redshift = PowerLawRedshift()

