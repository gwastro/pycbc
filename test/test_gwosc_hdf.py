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
Unit tests for reading GWOSC HDF5 files through pycbc.frame
"""

import os
import tempfile
import unittest

import h5py
import numpy
from utils import parse_args_cpu_only, simple_exit

import pycbc.frame
from pycbc.io import get_file
from pycbc.types import TimeSeries

parse_args_cpu_only("GWOSC HDF5 I/O")

RATE = 256


def write_gwosc_hdf(path, start, duration, ifo="H1", seed=0):
    """
    Write a file with the layout of a GWOSC HDF5 strain file and return
    the strain and data quality arrays.
    """
    rng = numpy.random.default_rng(seed)
    strain = rng.normal(size=duration * RATE)
    dqmask = rng.integers(0, 128, size=duration).astype(numpy.uint32)
    with h5py.File(path, "w") as f:
        f["meta/Detector"] = ifo
        f["meta/GPSstart"] = start
        f["meta/Duration"] = duration
        for name, data, spacing in [
            ("strain/Strain", strain, 1.0 / RATE),
            ("quality/simple/DQmask", dqmask, 1.0),
            ("quality/injections/Injmask", dqmask, 1.0),
        ]:
            f[name] = data
            f[name].attrs["Xstart"] = start
            f[name].attrs["Xspacing"] = spacing
    return strain, dqmask


class TestGWOSCHDF(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.start = 1000000000
        self.paths = []
        self.strain = []
        self.dq = []
        for i in range(3):
            start = self.start + 16 * i
            path = os.path.join(self.tmp.name, f"H-H1_TEST-{start}-16.hdf5")
            strain, dq = write_gwosc_hdf(path, start, 16, seed=i)
            self.paths.append(path)
            self.strain.append(strain)
            self.dq.append(dq)
        self.all_strain = numpy.concatenate(self.strain)

    def tearDown(self):
        self.tmp.cleanup()

    def test_whole_file(self):
        ts = pycbc.frame.read_frame(self.paths[0], "H1:STRAIN")
        self.assertEqual(ts.start_time, self.start)
        self.assertEqual(ts.delta_t, 1.0 / RATE)
        numpy.testing.assert_array_equal(ts.numpy(), self.strain[0])

    def test_data_quality(self):
        ts = pycbc.frame.read_frame(self.paths[0], "H1:DQMASK")
        self.assertEqual(ts.delta_t, 1.0)
        numpy.testing.assert_array_equal(ts.numpy(), self.dq[0])
        ts = pycbc.frame.read_frame(self.paths[0], "H1:INJMASK")
        numpy.testing.assert_array_equal(ts.numpy(), self.dq[0])

    def test_slice_across_files(self):
        ts = pycbc.frame.read_frame(
            self.paths,
            "H1:STRAIN",
            start_time=self.start + 10,
            end_time=self.start + 40,
        )
        self.assertEqual(ts.start_time, self.start + 10)
        numpy.testing.assert_array_equal(
            ts.numpy(), self.all_strain[10 * RATE : 40 * RATE]
        )

    def test_start_between_samples(self):
        # As with lalframe, data starts at the first sample after start_time
        start = self.start + 1 + 0.3 / RATE
        ts = pycbc.frame.read_frame(
            self.paths[0], "H1:STRAIN", start_time=start, duration=1
        )
        self.assertEqual(ts.start_time, self.start + 1 + 1.0 / RATE)
        numpy.testing.assert_array_equal(
            ts.numpy(), self.strain[0][RATE + 1 : 2 * RATE + 1]
        )

    def test_glob_and_cache_file(self):
        ts = pycbc.frame.read_frame(os.path.join(self.tmp.name, "*.hdf5"), "H1:STRAIN")
        numpy.testing.assert_array_equal(ts.numpy(), self.all_strain)

        cache_path = os.path.join(self.tmp.name, "test.lcf")
        with open(cache_path, "w") as f:
            for i, path in enumerate(self.paths):
                start = self.start + 16 * i
                f.write(f"H H1_TEST {start} 16 file://localhost{path}\n")
        ts = pycbc.frame.read_frame(cache_path, "H1:STRAIN")
        numpy.testing.assert_array_equal(ts.numpy(), self.all_strain)

    def test_gap(self):
        files = [self.paths[0], self.paths[2]]
        with self.assertRaises(ValueError):
            pycbc.frame.read_frame(
                files, "H1:STRAIN", start_time=self.start, end_time=self.start + 48
            )
        # Without an end time, read up to the gap
        ts = pycbc.frame.read_frame(files, "H1:STRAIN")
        numpy.testing.assert_array_equal(ts.numpy(), self.strain[0])

    def test_unknown_channel(self):
        for channel in ["L1:STRAIN", "H1:GDS-CALIB_STRAIN", "STRAIN"]:
            with self.assertRaises(ValueError):
                pycbc.frame.read_frame(self.paths[0], channel)

    def test_nonstandard_file_name(self):
        path = os.path.join(self.tmp.name, "strain.hdf5")
        strain, _ = write_gwosc_hdf(path, self.start, 16, seed=0)
        ts = pycbc.frame.read_frame(path, "H1:STRAIN")
        self.assertEqual(ts.start_time, self.start)
        numpy.testing.assert_array_equal(ts.numpy(), strain)

    def test_mixed_with_gwf(self):
        gwf_path = os.path.join(self.tmp.name, f"H-H1_TEST-{self.start + 48}-16.gwf")
        gwf_strain = numpy.random.default_rng(5).normal(size=16 * RATE)
        pycbc.frame.write_frame(
            gwf_path,
            "H1:STRAIN",
            TimeSeries(gwf_strain, delta_t=1.0 / RATE, epoch=self.start + 48),
        )
        ts = pycbc.frame.read_frame([self.paths[2], gwf_path], "H1:STRAIN")
        numpy.testing.assert_array_equal(
            ts.numpy(), numpy.concatenate([self.strain[2], gwf_strain])
        )

    def test_data_buffer(self):
        buf = pycbc.frame.DataBuffer(
            self.paths,
            "H1:STRAIN",
            self.start + 4,
            max_buffer=8,
            force_update_cache=False,
        )
        ts = buf.advance(4)
        self.assertEqual(buf.raw_sample_rate, RATE)
        numpy.testing.assert_array_equal(
            ts.numpy(), self.strain[0][4 * RATE : 8 * RATE]
        )


class TestGWOSCRealData(unittest.TestCase):
    def test_same_as_gwf(self):
        """The HDF5 and GWF versions of a GWOSC file give the same data"""
        url = (
            "https://gwosc.org/eventapi/json/GWTC-1-confident/GW150914/v3/"
            "H-H1_GWOSC_4KHZ_R1-1126259447-32.{}"
        )
        hdf_path = get_file(url.format("hdf5"), cache=True)
        gwf_path = get_file(url.format("gwf"), cache=True)
        for start, end in [(None, None), (1126259450.1, 1126259453.4)]:
            hdf = pycbc.frame.read_frame(
                hdf_path, "H1:STRAIN", start_time=start, end_time=end
            )
            gwf = pycbc.frame.read_frame(
                gwf_path, "H1:GWOSC-4KHZ_R1_STRAIN", start_time=start, end_time=end
            )
            self.assertEqual(hdf.start_time, gwf.start_time)
            numpy.testing.assert_array_equal(hdf.numpy(), gwf.numpy())


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestGWOSCHDF))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestGWOSCRealData))

if __name__ == "__main__":
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
