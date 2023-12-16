# Copyright (C) 2023  Shichao Wu, Alex Nitz
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
"""
This modules provides functions for base coordinate transformations, and
more advanced transformations between ground-based detectors and space-borne
detectors.
"""

from pycbc.coordinates.base import *
from pycbc.coordinates.space import *


__all__ = ['cartesian_to_spherical_rho', 'cartesian_to_spherical_azimuthal',
           'cartesian_to_spherical_polar', 'cartesian_to_spherical',
           'spherical_to_cartesian',
           'TIME_OFFSET_20_DEGREES',
           'localization_to_propagation_vector',
           'propagation_vector_to_localization', 'polarization_newframe',
           't_lisa_from_ssb', 't_ssb_from_t_lisa',
           'ssb_to_lisa', 'lisa_to_ssb',
           'rotation_matrix_ssb_to_lisa', 'rotation_matrix_ssb_to_geo',
           'lisa_position_ssb', 'earth_position_ssb',
           't_geo_from_ssb', 't_ssb_from_t_geo', 'ssb_to_geo', 'geo_to_ssb',
           'lisa_to_geo', 'geo_to_lisa',
           ]
