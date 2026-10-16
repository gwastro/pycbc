# Copyright (C) 2012  Andrew Miller, Alex Nitz, Josh Willis
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

#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
"""
These are the unittests for the pycbc frame/cache reading and writing functions
for GWF and GWOSC HDF5 files.
"""

import os
import tempfile
import unittest

import numpy
from utils import parse_args_cpu_only, simple_exit

import pycbc
import pycbc.frame
from pycbc.io import get_file
from pycbc.types import TimeSeries

# Frame tests only need to happen on the CPU
parse_args_cpu_only("Frame I/O")

RATE = 256


class FrameTestBase(unittest.TestCase):
    __test__ = False

    def setUp(self):
        numpy.random.seed(1023)
        self.size = pow(2, 12)

        self.data1 = numpy.random.rand(self.size).astype(self.dtype)
        self.data2 = numpy.random.rand(self.size).astype(self.dtype)

        # If the dtype is complex, throw in some complex values as well
        if self.dtype in (numpy.complex64, numpy.complex128):
            self.data1 += numpy.random.rand(self.size) * 1j
            self.data2 += numpy.random.rand(self.size) * 1j

        self.delta_t = 0.5
        self.epoch = 123456.0
        self.expected_data1 = TimeSeries(
            self.data1, dtype=self.dtype, epoch=self.epoch, delta_t=self.delta_t
        )
        self.expected_data2 = TimeSeries(
            self.data2, dtype=self.dtype, epoch=self.epoch, delta_t=self.delta_t
        )

    def test_frame(self):
        # File in temp directory deleted when garbage collected
        url = "https://github.com/gwastro/pycbc-config/raw/master/"
        url += "test_data_files/frametest{}.gwf"
        filename = get_file(url.format(self.data1.dtype), cache=True)

        # Make sure we can run from one directory higher as well
        if not os.path.exists(filename):
            filename = "test/" + filename

        # Reading just one channel first
        ts1 = pycbc.frame.read_frame(filename, "channel1")
        # Checking all values
        self.assertEqual(ts1, self.expected_data1)
        # Now checking the start time
        self.assertEqual(ts1.start_time, self.epoch)
        # And the duration
        self.assertEqual(ts1.end_time - ts1.start_time, self.size * self.delta_t)

        # Now reading multiple channels
        ts2 = pycbc.frame.read_frame(filename, ["channel1", "channel2"])
        # We should get back a list
        self.assertTrue(type(ts2) is list)
        self.assertEqual(ts2[0], self.expected_data1)
        self.assertEqual(ts2[1], self.expected_data2)
        self.assertEqual(ts2[0].start_time, self.epoch)
        self.assertEqual(ts2[1].start_time, self.epoch)
        self.assertEqual(ts2[0].end_time - ts2[0].start_time, self.size * self.delta_t)
        self.assertEqual(ts2[1].end_time - ts2[1].start_time, self.size * self.delta_t)

        # These are the times and indices for the segment we will try to read
        start = self.epoch + 10
        end = self.epoch + 50
        startind = int(10 / self.delta_t)
        endind = int(50 / self.delta_t)

        # Now reading in a specific segment with an integer
        ts3 = pycbc.frame.read_frame(
            filename, "channel1", start_time=int(start), end_time=int(end)
        )

        # The same, but with a LIGOTimeGPS for the start and end times
        ts4 = pycbc.frame.read_frame(
            filename, "channel1", start_time=start, end_time=end
        )

        # Now we will check those two TimeSeries
        self.assertEqual(ts3, self.expected_data1[startind:endind])
        self.assertEqual(ts4, self.expected_data1[startind:endind])
        dur3 = float(ts3.end_time) - float(ts3.start_time)
        self.assertTrue(40 - dur3 < self.delta_t)
        self.assertEqual(ts3.start_time, start)
        dur4 = float(ts4.end_time) - float(ts4.start_time)
        self.assertTrue(40 - dur4 < self.delta_t)
        self.assertEqual(ts4.start_time, start)

        # And now some cases that should raise errors

        # There must be a span grater than 0
        self.assertRaises(
            ValueError,
            pycbc.frame.read_frame,
            filename,
            "channel1",
            start_time=self.epoch,
            end_time=self.epoch,
        )
        # The start must be before the end
        self.assertRaises(
            ValueError,
            pycbc.frame.read_frame,
            filename,
            "channel1",
            start_time=self.epoch + 1,
            end_time=self.epoch,
        )

    def test_write_and_read(self):
        """Test symmetric writing and reading of GWF and GWOSC HDF5 files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            gwf_path = os.path.join(tmpdir, "test.gwf")
            hdf_path = os.path.join(tmpdir, "test.hdf5")

            # Write single channel
            pycbc.frame.write_frame(gwf_path, "H1:STRAIN", self.expected_data1)
            pycbc.frame.write_frame(hdf_path, "H1:STRAIN", self.expected_data1)

            # Read back with read_frame
            ts_gwf = pycbc.frame.read_frame(gwf_path, "H1:STRAIN")
            ts_hdf = pycbc.frame.read_frame(hdf_path, "H1:STRAIN")

            self.assertEqual(ts_gwf, self.expected_data1)
            self.assertEqual(ts_hdf, self.expected_data1)
            # Cross-format verification: read what was written in both
            self.assertEqual(ts_gwf, ts_hdf)

            # Write multiple channels
            gwf_multi = os.path.join(tmpdir, "multi.gwf")
            hdf_multi = os.path.join(tmpdir, "multi.hdf5")
            pycbc.frame.write_frame(
                gwf_multi,
                ["H1:STRAIN", "H1:INJMASK"],
                [self.expected_data1, self.expected_data2],
            )
            pycbc.frame.write_frame(
                hdf_multi,
                ["H1:STRAIN", "H1:INJMASK"],
                [self.expected_data1, self.expected_data2],
            )

            ts_gwf_list = pycbc.frame.read_frame(gwf_multi, ["H1:STRAIN", "H1:INJMASK"])
            ts_hdf_list = pycbc.frame.read_frame(hdf_multi, ["H1:STRAIN", "H1:INJMASK"])

            self.assertEqual(ts_gwf_list[0], self.expected_data1)
            self.assertEqual(ts_gwf_list[1], self.expected_data2)
            self.assertEqual(ts_hdf_list[0], self.expected_data1)
            self.assertEqual(ts_hdf_list[1], self.expected_data2)
            self.assertEqual(ts_gwf_list[0], ts_hdf_list[0])
            self.assertEqual(ts_gwf_list[1], ts_hdf_list[1])

            # Writing without IFO prefix must raise ValueError
            no_ifo_path = os.path.join(tmpdir, "no_ifo.hdf5")
            with self.assertRaises(ValueError):
                pycbc.frame.write_frame(no_ifo_path, "STRAIN", self.expected_data1)

            # Writing channels with mismatched IFO must raise ValueError
            mismatched_path = os.path.join(tmpdir, "mismatched.hdf5")
            with self.assertRaises(ValueError):
                pycbc.frame.write_frame(
                    mismatched_path,
                    ["H1:STRAIN", "L1:STRAIN"],
                    [self.expected_data1, self.expected_data2],
                )


def write_gwosc_hdf(path, start, duration, ifo="H1", seed=0):
    """Write a GWOSC HDF5 strain file layout using write_frame."""
    rng = numpy.random.default_rng(seed)
    strain = rng.normal(size=duration * RATE)
    dqmask = rng.integers(0, 128, size=duration).astype(numpy.uint32)
    pycbc.frame.write_frame(
        path,
        [f"{ifo}:STRAIN", f"{ifo}:DQMASK", f"{ifo}:INJMASK"],
        [
            TimeSeries(strain, delta_t=1.0 / RATE, epoch=start),
            TimeSeries(dqmask, delta_t=1.0, epoch=start),
            TimeSeries(dqmask, delta_t=1.0, epoch=start),
        ],
    )
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


# We take a factory approach so we can test all possible dtypes we support
TestClasses = []
types = [numpy.float32, numpy.float64, numpy.complex64, numpy.complex128]

for ty in types:
    klass = type(f"{ty.__name__}_Test", (FrameTestBase,), {"dtype": ty})
    klass.__test__ = True
    vars()[klass.__name__] = klass
    TestClasses.append(klass)
    del klass

TestClasses.extend([TestGWOSCHDF, TestGWOSCRealData])

if __name__ == "__main__":
    suite = unittest.TestSuite()
    for klass in TestClasses:
        suite.addTest(unittest.TestLoader().loadTestsFromTestCase(klass))
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
