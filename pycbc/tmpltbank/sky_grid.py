"""Functionality for handling grids of points in the sky for coherent SNR
calculation via `pycbc_multi_inspiral`. The main operation to be performed on
these points is calculating the antenna pattern functions and time delays from
the Earth center for a network of detectors.
"""

import numpy as np
import h5py

from pycbc.detector import Detector
from pycbc.conversions import ensurearray


class SkyGrid:
    def __init__(self, ra, dec, detectors, ref_gps_time):
        """Initialize a sky grid from a list of RA/dec coordinates.

        Parameters
        ----------
        ra: iterable of floats
            Right ascensions for each point in radians, in the interval [0,2π).
        dec: iterable of floats
            Declination for each point in radians, where π/2 is the North pole,
            -π/2 is the South pole, and 0 is the celestial equator.
        detectors: iterable of str
            List of detector names associated with the sky grid, typically the
            detectors for which the grid has been placed or is to be placed.
            The detectors will be used when calculating the antenna pattern
            functions and time delays from Earth center.
        ref_gps_time: float
            Reference GPS time associated with the sky grid. This will be used
            when calculating the antenna pattern functions and time delays from
            Earth center.
        """
        # We store the points in a 2D array internally, first dimension runs
        # over the list of points, second dimension is RA/dec.
        # Question: should we use Astropy sky positions instead?
        ra, dec, _ = ensurearray(ra,dec)
        if (ra < 0).any() or (ra > 2 * np.pi).any():
            raise ValueError('RA must be in the range [0,2π]')
        if (dec < -np.pi/2).any() or (dec > np.pi/2).any():
            raise ValueError('DEC must be in the range [-π/2, π/2]')
        self.positions = np.vstack([ra, dec]).T
        self.detectors = sorted(detectors)
        self.ref_gps_time = ref_gps_time

    def __len__(self):
        """Returns the number of points in the sky grid."""
        return self.positions.shape[0]

    def __getitem__(self, index):
        """Returns the coordinates of a single point in the grid."""
        return self.positions[index]

    @property
    def ras(self):
        """Returns all right ascensions in radians, in the interval [0,2π)."""
        return self.positions[:, 0]

    @property
    def decs(self):
        """Returns all declinations in radians, where π/2 is the North pole,
        -π/2 is the South pole, and 0 is the celestial equator."""
        return self.positions[:, 1]

    @classmethod
    def from_cli(cls, cli_parser, cli_args):
        """Initialize a sky grid from command-line interface, via argparse
        objects.
        """
        if cli_args.sky_grid is not None:
            if cli_args.ra is not None or cli_args.dec is not None:
                cli_parser.error(
                    'Please provide either a sky grid via --sky-grid or a '
                    'single sky position via --ra and --dec, not both'
                )
            return cls.read_from_file(cli_args.sky_grid)
        if cli_args.ra is not None and cli_args.dec is not None:
            return cls(
                [cli_args.ra],
                [cli_args.dec],
                cli_args.instruments,
                cli_args.trigger_time
            )
        cli_parser.error(
            'Please specify a sky grid via --sky-grid or a single sky '
            'position via --ra and --dec'
        )

    @classmethod
    def read_from_file(cls, path):
        """Initialize a sky grid from a given HDF5 file."""
        with h5py.File(path, 'r') as hf:
            ra = hf['ra'][:]
            dec = hf['dec'][:]
            detectors = hf.attrs['detectors']
            ref_gps_time = hf.attrs['ref_gps_time']
        return cls(ra, dec, detectors, ref_gps_time)

    def write_to_file(self, path, extra_attrs=None, extra_datasets=None):
        """Writes a sky grid to an HDF5 file."""
        with h5py.File(path, 'w') as hf:
            hf['ra'] = self.ras
            hf['dec'] = self.decs
            hf.attrs['detectors'] = self.detectors
            hf.attrs['ref_gps_time'] = self.ref_gps_time
            for attribute in (extra_attrs or {}):
                hf.attrs[attribute] = extra_attrs[attribute]
            for dataset in (extra_datasets or {}):
                hf[dataset] = extra_datasets[dataset]

    def calculate_antenna_patterns(self):
        """Calculate the antenna pattern functions at each point in the grid
        for the list of GW detectors specified at instantiation. Return a dict,
        keyed by detector name, whose items are 2-dimensional Numpy arrays.
        The first dimension of these arrays runs over the sky grid, and the
        second dimension runs over the plus and cross polarizations.
        """
        result = {}
        for det_name in self.detectors:
            det = Detector(det_name)
            result[det_name] = np.empty((len(self), 2))
            for i, (ra, dec) in enumerate(self):
                result[det_name][i] = det.antenna_pattern(
                    ra, dec, 0, t_gps=self.ref_gps_time
                )
        return result

    def calculate_time_delays(self):
        """Calculate the time delays from the Earth center to each GW detector
        specified at instantiation, for each point in the grid. Return a dict,
        keyed by detector name, whose items are 1-dimensional Numpy arrays
        containing the time delays for each sky point.
        """
        result = {}
        for det_name in self.detectors:
            det = Detector(det_name)
            result[det_name] = np.empty(len(self))
            for i, (ra, dec) in enumerate(self):
                result[det_name][i] = det.time_delay_from_earth_center(
                    ra, dec, t_gps=self.ref_gps_time
                )
        return result


__all__ = ['SkyGrid']
