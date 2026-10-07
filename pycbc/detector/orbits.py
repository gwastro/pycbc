# Copyright (C) 2026 Shichao Wu, Alex Nitz, Alex Correia
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

"""Base class for spacecraft orbit models."""

from abc import ABC, abstractmethod

import numpy

__all__ = ["BaseOrbit"]


class BaseOrbit(ABC):
    """Interface for native PyCBC spacecraft orbit providers.

    Positions use solar-system-barycentric Cartesian coordinates with fixed
    J2000 ecliptic axes. Positions, velocities and accelerations use SI units.
    Times are SSB coordinate times in seconds; a subclass must document its
    epoch and time scale. Evaluation times and interpolation knots must use
    the same convention. This interface does not convert epochs or frames.

    Each ``compute_*`` method accepts scalar or one-dimensional times and
    one-based spacecraft labels. It returns an array of shape
    ``(number of times, number of selected spacecraft, 3)``. Singleton axes
    are retained, and input order (including repeats) is preserved. ``sc=None``
    selects all spacecraft. Subclasses must document their supported time
    interval and whether they permit extrapolation.

    Subclasses implement orbital kinematics and should call
    ``super().__init__`` in their constructors. They can use the input helpers
    below, but remain responsible for output shapes and units. Compatible
    third-party providers need not inherit from this class.

    Parameters
    ----------
    t_interp : array-like or None, optional
        Finite, strictly increasing sample times for an interpolated orbit.
        At least two samples are required. ``None`` denotes an analytic orbit
        without an interpolation grid. The base class stores an independent
        copy; interpolation order and sample data belong to the subclass.
    num_sc : int, optional
        Number of spacecraft, default 3. Labels run from 1 to ``num_sc``.
    """

    def __init__(self, t_interp=None, num_sc=3):
        if (
            isinstance(num_sc, (bool, numpy.bool_))
            or not isinstance(num_sc, (int, numpy.integer))
            or num_sc < 1
        ):
            raise ValueError("num_sc must be a positive integer")
        self.num_sc = int(num_sc)
        self.t_interp = None
        if t_interp is not None:
            times = self._prepare_times(t_interp)
            if len(times) < 2 or numpy.any(numpy.diff(times) <= 0):
                raise ValueError(
                    "t_interp needs at least two strictly increasing times"
                )
            self.t_interp = times.copy()

    @staticmethod
    def _prepare_times(t):
        """Normalize evaluation times without sorting or removing repeats."""
        times = numpy.asarray(t)
        if times.ndim > 1 or times.dtype.kind not in "iuf":
            raise ValueError("times must be real scalars or a 1D array")
        times = numpy.atleast_1d(numpy.asarray(times, dtype=float))
        if not numpy.all(numpy.isfinite(times)):
            raise ValueError("times must be finite")
        return times

    def _sc_indices(self, sc):
        """Map one-based spacecraft labels to indices, preserving order."""
        if sc is None:
            return numpy.arange(self.num_sc)
        labels = numpy.atleast_1d(numpy.asarray(sc))
        if labels.ndim != 1 or labels.dtype.kind not in "iu":
            raise ValueError("spacecraft labels must be integers")
        if numpy.any(labels < 1) or numpy.any(labels > self.num_sc):
            raise ValueError(f"spacecraft labels must be 1..{self.num_sc}")
        return labels.astype(numpy.intp) - 1

    @abstractmethod
    def compute_position(self, t, sc=None):
        """Return spacecraft positions in metres.

        Parameters
        ----------
        t : float or one-dimensional array-like
            Evaluation times in seconds, in the orbit's time convention.
            Times need not be increasing or distinct.
        sc : int or one-dimensional array-like of int or None, optional
            One-based spacecraft labels; ``None`` selects all spacecraft.

        Returns
        -------
        numpy.ndarray
            Positions with shape ``(N_time, N_selected, 3)``.

        Raises
        ------
        ValueError
            For invalid times or spacecraft labels. Subclasses define their
            own policy for times outside their supported interval.
        """

    @abstractmethod
    def compute_velocity(self, t, sc=None):
        """Return velocities in metres/second.

        Inputs, frame and output shape follow :meth:`compute_position`.
        """

    @abstractmethod
    def compute_acceleration(self, t, sc=None):
        """Return accelerations in metres/second squared.

        Inputs, frame and output shape follow :meth:`compute_position`.
        """
