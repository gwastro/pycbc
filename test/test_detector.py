# Copyright (C) 2018 Alex Nitz
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
These are the unittests for the pycbc.waveform module
"""
from __future__ import print_function
import pycbc.detector as det
import unittest, numpy
from numpy.random import uniform, seed
seed(0)

# We require lal as a reference comparison
import lal

from utils import simple_exit

class TestDetector(unittest.TestCase):
    def setUp(self):
        self.d = [det.Detector(ifo)
                  for ifo, name in det.get_available_detectors()]

        # not distributed sanely, but should provide some good coverage
        N = 1000
        self.ra = uniform(0, numpy.pi * 2, size=N)
        self.dec = uniform(-numpy.pi, numpy.pi, size=N)
        self.pol = uniform(0, numpy.pi * 2, size=N)
        self.time = uniform(1126000000.0, 1336096017.0, size=N)

    def test_light_time(self):
        for d1 in self.d:
            for d2 in self.d:
                t1 = lal.LightTravelTime(d1.lal(), d2.lal()) * 1e-9
                t2 = d1.light_travel_time_to_detector(d2)
                self.assertAlmostEqual(t1, t2, 7)

    def test_antenna_pattern(self):
        vals = list(zip(self.ra, self.dec, self.pol, self.time))
        for ifo in self.d:
            fp = []
            fc = []
            for ra1, dec1, pol1, time1 in vals:
                gmst = lal.GreenwichMeanSiderealTime(time1)
                fp1, fc1 = tuple(lal.ComputeDetAMResponse(ifo.response, ra1, dec1, pol1, gmst))
                fp.append(fp1)
                fc.append(fc1)

            fp2, fc2 = ifo.antenna_pattern(self.ra, self.dec, self.pol, self.time)

            fp = numpy.array(fp)
            fc = numpy.array(fc)

            diff1 = fp - fp2
            diff2 = fc - fc2
            diff = abs(numpy.concatenate([diff1, diff2]))
            tolerance = 2e-4
            print("Max antenna diff:", ifo.name, diff.max())

            self.assertLess(diff.max(), tolerance)

    def test_delay_from_detector(self):
        ra, dec, time = self.ra[0:10], self.dec[0:10], self.time[0:10]
        for d1 in self.d:
            for d2 in self.d:
                time1 = []
                for ra1, dec1, tim1 in zip(ra, dec, time):
                    t1 = lal.ArrivalTimeDiff(d1.location, d2.location,
                                             ra1, dec1, tim1)
                    time1.append(t1)
                time1 = numpy.array(time1)
                time2 = d1.time_delay_from_detector(d2, ra, dec, time)
                self.assertLess(abs(time1 - time2).max(), 1e-3)

    def test_optimal_orientation(self):
        for d1 in self.d:
            ra, dec = d1.optimal_orientation(self.time[0])
            ra1 = d1.longitude + lal.GreenwichMeanSiderealTime(self.time[0]) % (numpy.pi *2)
            dec1 = d1.latitude

            self.assertAlmostEqual(ra, ra1, 3)
            self.assertAlmostEqual(dec, dec1, 7)


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestDetector))

if __name__ == '__main__':
    from astropy.utils import iers
    iers.conf.auto_download = False
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
