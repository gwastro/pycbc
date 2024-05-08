"""Functionality for handling grids of points in the sky for coherent SNR
calculation via `pycbc_multi_inspiral`. The main operation to be performed on
these points is calculating the antenna pattern functions and time delays from
the Earth center for a network of detectors.
"""

import numpy as np
import h5py

from pycbc.detector import Detector


class SkyGrid:
    def __init__(self, ra, dec):
        """Initialize a sky grid from a list of RA/dec coordinates.

        Parameters
        ----------
        ra: iterable of floats
            Right ascensions for each point UNITS AND ANGULAR CONVENTION.
        dec: iterable of floats
            Declination for each point UNITS AND ANGULAR CONVENTION.
        """
        # Question: should the list of detectors and reference times be part of
        # the SkyGrid object? Or should they be given each time they are needed,
        # independently? This decision will affect the `calculate_*()` methods
        # down below.
        # We store the points in a 2D array internally, first dimension runs
        # over the list of points, second dimension is RA/dec.
        # Question: should we use Astropy sky positions instead?
        self.positions = np.vstack([ra, dec]).T

    def __len__(self):
        """Returns the number of points in the sky grid.
        """
        return self.positions.shape[0]

    def __getitem__(self, index):
        """Returns the coordinates of a single point in the grid.
        """
        return self.positions[index]

    @property
    def ras(self):
        """Returns all right ascensions.
        """
        return self.positions[:,0]

    @property
    def decs(self):
        """Returns all declinations.
        """
        return self.positions[:,1]

    @classmethod
    def from_cli(cls, cli_parser, cli_args):
        """Initialize a sky grid from command-line interface, via argparse
        objects.
        """
        if cli_args.sky_grid is not None:
            if cli_args.ra is not None or cli_args.dec is not None:
                cli_parser.error('Use --sky-grid or --ra & --dec, not both')
            return cls.read_from_file(cli_args.sky_grid)
        if cli_args.ra is not None and cli_args.dec is not None:
            return cls([cli_args.ra], [cli_args.dec])
        cli_parser.error('--ra and --dec must be used together')

    @classmethod
    def read_from_file(cls, path):
        """Initialize a sky grid from a given HDF5 file.
        """
        with h5py.File(path, 'r') as hf:
            ra = hf['ra'][:]
            dec = hf['dec'][:]
        return cls(ra, dec)

    def write_to_file(self, path):
        """Writes a sky grid to an HDF5 file.
        """
        with h5py.File(path, 'w') as hf:
            hf['ra'] = self.ras
            hf['dec'] = self.decs

    def calculate_antenna_patterns(self, detector_names, gps_time):
        """Calculate the antenna pattern functions at each point in the grid
        for a list of GW detectors.
        """
        result = {}
        for det_name in detector_names:
            det = Detector(det_name)
            result[det_name] = np.empty((len(self), 2))
            for i, (ra, dec) in enumerate(self):
                result[det_name][i] = det.antenna_pattern(
                    ra, dec, 0, t_gps=gps_time
                )
        return result

    def calculate_time_delays(self, detectors, gps_time):
        """Calculate the time delays from the Earth center to a list of GW
        detectors at each point in the grid.
        """
        result = {}
        for det_name in detector_names:
            det = Detector(det_name)
            result[det_name] = np.empty(len(self))
            for i, (ra, dec) in enumerate(self):
                result[det_name][i] = det.time_delay_from_earth_center(
                    ra, dec, t_gps=gps_time
                )
        return result


__all__ = ['SkyGrid']