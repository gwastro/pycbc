# Copyright (C) 2026  Shichao Wu, Alex Nitz, Jacopo Tissino, Jan Harms
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
Regression tests for `pycbc.coordinates.reference_point`, PyCBC's SSB-hub
coordinate unification (arrival time, sky localization, polarization angle)
for an arbitrary fixed reference point ("Ref" frame).
"""
import numpy
import unittest
from scipy.optimize import fsolve

from astropy.constants import c

from pycbc.coordinates import space, reference_point as refpt
from utils import simple_exit


seed = 8202
numpy.random.seed(seed)

NUM_SAMPLES = 20


def _random_params(num=NUM_SAMPLES):
    t_ssb = numpy.random.uniform(1.0e8, 6.0e8, size=num)
    longitude = numpy.random.uniform(0.0, 2 * numpy.pi, size=num)
    latitude = numpy.random.uniform(-numpy.pi / 2, numpy.pi / 2, size=num)
    polarization = numpy.random.uniform(0.0, 2 * numpy.pi, size=num)
    return t_ssb, longitude, latitude, polarization


class TestReferencePointCoordinates(unittest.TestCase):
    def test_rotation_matrix_is_identity(self):
        r = refpt.rotation_matrix_ssb_to_ref()
        self.assertEqual(r.shape, (3, 3))
        self.assertTrue(numpy.array_equal(r, numpy.eye(3)))

    def test_t_ref_from_ssb_matches_fsolve(self):
        """The closed-form `t_ref_from_ssb` must agree with an
        independent, brute-force `scipy.optimize.fsolve` cross-check of
        the same light-travel-time relation -- a fixed point needs no
        root-finding, but the closed-form answer should still match a
        numerical solver applied to the same underlying equation.
        """
        ref_position = numpy.array([1.0e10, -2.0e10, 3.0e9])
        t_ssb, lon, lat, _ = _random_params(num=5)
        for i in range(5):
            k = space.localization_to_propagation_vector(
                lon[i], lat[i], use_astropy=False)

            def eqn(t_ref, t_ssb=t_ssb[i], k=k):
                return t_ref - (t_ssb + numpy.vdot(k, ref_position) / c.value)

            t_ref_fsolve = fsolve(eqn, t_ssb[i])[0]
            t_ref_closed = refpt.t_ref_from_ssb(
                t_ssb[i], lon[i], lat[i], ref_position)
            self.assertAlmostEqual(t_ref_closed, t_ref_fsolve, places=6)

    def test_t_ssb_from_t_ref_is_inverse(self):
        ref_position = numpy.array([5.0e9, 1.0e9, -4.0e9])
        t_ssb, lon, lat, _ = _random_params()
        t_ref = refpt.t_ref_from_ssb(t_ssb, lon, lat, ref_position)
        t_ssb_rt = refpt.t_ssb_from_t_ref(t_ref, lon, lat, ref_position)
        self.assertLess(numpy.max(numpy.abs(t_ssb_rt - t_ssb)), 1e-8)

    def test_zero_ref_position_degenerates_to_ssb(self):
        """`ref_position=(0,0,0)` must give t_ref == t_ssb exactly (no
        light-travel-time correction at all)."""
        ref_position = numpy.zeros(3)
        t_ssb, lon, lat, pol = _random_params()
        t_ref, lon_r, lat_r, pol_r = refpt.ssb_to_ref(
            t_ssb, lon, lat, pol, ref_position)
        self.assertTrue(numpy.array_equal(t_ref, t_ssb))
        # lon/lat/pol go through a vector<->angle round trip even under
        # the identity rotation, so ordinary floating-point trig error
        # applies (same tolerance as the identity-rotation design-lock
        # test above), not bit-exact equality.
        self.assertLess(numpy.max(numpy.abs(lon_r - lon)), 1e-10)
        self.assertLess(numpy.max(numpy.abs(lat_r - lat)), 1e-10)
        self.assertLess(numpy.max(numpy.abs(pol_r - pol)), 1e-10)

    def test_ssb_to_ref_round_trip(self):
        ref_position = numpy.array([2.0e10, 3.0e10, -1.0e10])
        t_ssb, lon, lat, pol = _random_params()
        t_ref, lon_r, lat_r, pol_r = refpt.ssb_to_ref(
            t_ssb, lon, lat, pol, ref_position)
        t_ssb_rt, lon_rt, lat_rt, pol_rt = refpt.ref_to_ssb(
            t_ref, lon_r, lat_r, pol_r, ref_position)
        self.assertLess(numpy.max(numpy.abs(t_ssb_rt - t_ssb)), 1e-6)
        self.assertLess(numpy.max(numpy.abs(lon_rt - lon)), 1e-8)
        self.assertLess(numpy.max(numpy.abs(lat_rt - lat)), 1e-8)
        self.assertLess(numpy.max(numpy.abs(pol_rt - pol)), 1e-8)

    def test_ssb_to_ref_design_lock_identity_rotation(self):
        """Because `rotation_matrix_ssb_to_ref` is the identity, sky
        localization/polarization must come out of `ssb_to_ref`
        numerically unchanged; only arrival time changes.
        """
        ref_position = numpy.array([1.0e10, 0.0, 0.0])
        t_ssb, lon, lat, pol = _random_params()
        t_ref, lon_r, lat_r, pol_r = refpt.ssb_to_ref(
            t_ssb, lon, lat, pol, ref_position)
        self.assertLess(numpy.max(numpy.abs(lon_r - lon)), 1e-10)
        self.assertLess(numpy.max(numpy.abs(lat_r - lat)), 1e-10)
        self.assertLess(numpy.max(numpy.abs(pol_r - pol)), 1e-10)
        self.assertFalse(numpy.array_equal(t_ref, t_ssb))

    def test_ref_to_geo_round_trip(self):
        ref_position = numpy.array([1.0e9, -2.0e9, 3.0e8])
        t_ssb, lon, lat, pol = _random_params(num=5)
        t_ref, lon_r, lat_r, pol_r = refpt.ssb_to_ref(
            t_ssb, lon, lat, pol, ref_position)
        t_geo, lon_g, lat_g, pol_g = refpt.ref_to_geo(
            t_ref, lon_r, lat_r, pol_r, ref_position)
        t_ref_rt, lon_r_rt, lat_r_rt, pol_r_rt = refpt.geo_to_ref(
            t_geo, lon_g, lat_g, pol_g, ref_position)
        self.assertLess(numpy.max(numpy.abs(t_ref_rt - t_ref)), 1e-2)
        self.assertLess(numpy.max(numpy.abs(lon_r_rt - lon_r)), 1e-6)
        self.assertLess(numpy.max(numpy.abs(lat_r_rt - lat_r)), 1e-6)
        self.assertLess(numpy.max(numpy.abs(pol_r_rt - pol_r)), 1e-6)

    def test_ref_to_lisa_round_trip(self):
        ref_position = numpy.array([1.0e9, -2.0e9, 3.0e8])
        t_ssb, lon, lat, pol = _random_params(num=5)
        t_ref, lon_r, lat_r, pol_r = refpt.ssb_to_ref(
            t_ssb, lon, lat, pol, ref_position)
        t_lisa, lon_l, lat_l, pol_l = refpt.ref_to_lisa(
            t_ref, lon_r, lat_r, pol_r, ref_position)
        t_ref_rt, lon_r_rt, lat_r_rt, pol_r_rt = refpt.lisa_to_ref(
            t_lisa, lon_l, lat_l, pol_l, ref_position)
        self.assertLess(numpy.max(numpy.abs(t_ref_rt - t_ref)), 1e-2)
        self.assertLess(numpy.max(numpy.abs(lon_r_rt - lon_r)), 1e-6)
        self.assertLess(numpy.max(numpy.abs(lat_r_rt - lat_r)), 1e-6)
        self.assertLess(numpy.max(numpy.abs(pol_r_rt - pol_r)), 1e-6)

    def test_ref_to_lisa_round_trip_with_taiji_orbit(self):
        """`orbit=` passed through ref_to_lisa/lisa_to_ref reaches
        `space.ssb_to_lisa`/`lisa_to_ssb` unchanged, generalizing to
        Taiji/TianQin/numeric orbits for free.
        """
        from pycbc.coordinates import space_orbit
        orbit = space_orbit.TaijiAnalyticOrbit()
        ref_position = numpy.array([1.0e9, -2.0e9, 3.0e8])
        t_ssb = 2.0e7
        lon, lat, pol = 1.3, -0.4, 0.9
        t_ref, lon_r, lat_r, pol_r = refpt.ssb_to_ref(
            t_ssb, lon, lat, pol, ref_position)
        t_taiji, lon_t, lat_t, pol_t = refpt.ref_to_lisa(
            t_ref, lon_r, lat_r, pol_r, ref_position, orbit=orbit)
        t_ref_rt, lon_r_rt, lat_r_rt, pol_r_rt = refpt.lisa_to_ref(
            t_taiji, lon_t, lat_t, pol_t, ref_position, orbit=orbit)
        self.assertAlmostEqual(t_ref_rt, t_ref, places=2)
        self.assertAlmostEqual(lon_r_rt, lon_r, places=6)
        self.assertAlmostEqual(lat_r_rt, lat_r, places=6)
        self.assertAlmostEqual(pol_r_rt, pol_r, places=6)


class TestOptionalLunarskyRef(unittest.TestCase):
    """Tests requiring the optional `lunarsky` package (Moon<->Ref
    round trip goes through `pycbc.coordinates.moon`, whose
    barycenter-only path doesn't need lunarsky, but keeping this
    consistent with `test_coordinates_moon.py`'s own skip pattern for
    the site-specific case).
    """
    def setUp(self):
        try:
            import lunarsky  # noqa: F401
        except ImportError:
            self.skipTest('lunarsky not installed; skipping Moon<->Ref '
                          'site-specific tests')

    def test_ref_to_moon_round_trip_with_site(self):
        ref_position = numpy.array([1.0e9, -2.0e9, 3.0e8])
        lon_site, lat_site = numpy.deg2rad(30.0), numpy.deg2rad(-89.0)
        t_ssb, lon, lat, pol = _random_params(num=5)
        t_ref, lon_r, lat_r, pol_r = refpt.ssb_to_ref(
            t_ssb, lon, lat, pol, ref_position)
        t_moon, lon_m, lat_m, pol_m = refpt.ref_to_moon(
            t_ref, lon_r, lat_r, pol_r, ref_position,
            longitude_site=lon_site, latitude_site=lat_site)
        t_ref_rt, lon_r_rt, lat_r_rt, pol_r_rt = refpt.moon_to_ref(
            t_moon, lon_m, lat_m, pol_m, ref_position,
            longitude_site=lon_site, latitude_site=lat_site)
        self.assertLess(numpy.max(numpy.abs(t_ref_rt - t_ref)), 1e-2)
        self.assertLess(numpy.max(numpy.abs(lon_r_rt - lon_r)), 1e-6)
        self.assertLess(numpy.max(numpy.abs(lat_r_rt - lat_r)), 1e-6)
        self.assertLess(numpy.max(numpy.abs(pol_r_rt - pol_r)), 1e-6)


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestReferencePointCoordinates))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestOptionalLunarskyRef))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
