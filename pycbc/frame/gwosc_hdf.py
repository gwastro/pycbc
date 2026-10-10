# Copyright (C) 2026 Alex Nitz
#
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
Read GWOSC HDF5 strain files.

GWOSC (https://gwosc.org) distributes the same data as GWF and HDF5 files.
In the HDF5 files the strain is the dataset ``strain/Strain`` and the 1 Hz
data quality and injection bitmasks are ``quality/simple/DQmask`` and
``quality/injections/Injmask``. Each has the attributes ``Xstart`` (GPS time
of the first sample) and ``Xspacing`` (sample spacing in seconds).

This module provides the same functions as `pycbc.frame.gwf`.
"""

import math

import h5py

_HDF5_SIGNATURE = b"\x89HDF\r\n\x1a\n"

# Times closer than this to a sample are treated as on the sample, as lalframe
# does
_FUZZ = 0.1 / 16384

# Channels are named '<IFO>:<NAME>', e.g. 'H1:STRAIN'
_DATASETS = {
    "STRAIN": "strain/Strain",
    "DQMASK": "quality/simple/DQmask",
    "INJMASK": "quality/injections/Injmask",
}


def is_format(path):
    """Test if a file is an HDF5 file by checking its signature."""
    try:
        with open(path, "rb") as f:
            return f.read(8) == _HDF5_SIGNATURE
    except OSError:
        return False


def span(path):
    """Return the GPS start time and duration of the data in a file."""
    with h5py.File(path, "r") as f:
        return float(f["meta/GPSstart"][()]), float(f["meta/Duration"][()])


def _dataset(h5file, channel):
    """
    Return the dataset for `channel`, which is '<IFO>:STRAIN',
    '<IFO>:DQMASK' or '<IFO>:INJMASK', where <IFO> must be the file's
    detector.
    """
    ifo = h5file["meta/Detector"].asstr()[()]
    prefix, _, name = channel.partition(":")
    if prefix != ifo or name not in _DATASETS:
        known = ", ".join(f"{ifo}:{name}" for name in _DATASETS)
        raise ValueError(
            f"Channel {channel} is not in {h5file.filename}; "
            f"GWOSC HDF5 files contain {known}"
        )
    return h5file[_DATASETS[name]]


def channel_info(path, channel):
    """Return the numpy dtype and sample rate of a channel in a file."""
    with h5py.File(path, "r") as f:
        dset = _dataset(f, channel)
        return dset.dtype, 1.0 / dset.attrs["Xspacing"]


def read(path, channel, start, end, check_integrity=False):
    """
    Read the samples of a channel between GPS times `start` and `end`,
    which must lie within the file. `check_integrity` is ignored, as GWOSC
    HDF5 files contain no checksums.

    Returns
    -------
    data : numpy.ndarray
    epoch : float
        GPS time of the first sample.
    delta_t : float
        Sample spacing in seconds.

    """
    with h5py.File(path, "r") as f:
        dset = _dataset(f, channel)
        epoch = float(dset.attrs["Xstart"])
        delta_t = float(dset.attrs["Xspacing"])
        # Choose samples as lalframe does (see LALFrStreamReadTS_source.c):
        # start at the first sample after `start`, or at a sample within
        # _FUZZ of it, and read int((end - start) / delta_t) samples.
        first = max(0, math.ceil((start - epoch - _FUZZ) / delta_t))
        count = int((end - start) / delta_t)
        return dset[first : first + count], epoch + first * delta_t, delta_t
