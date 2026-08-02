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
Regression tests for `pycbc.coordinates.moon`, PyCBC's SSB-hub coordinate
unification (arrival time, sky localization, polarization angle) for a
lunar detector (e.g. LGWA). Does not test any antenna-pattern/response
function -- that layer does not exist yet (see `moon.py`'s module
docstring and `pycbc.detector.space._Generic_detector`).
"""
import numpy
import unittest

from pycbc.coordinates import space, moon, space_orbit
from utils import simple_exit


seed = 8202
numpy.random.seed(seed)

# the ~1 sidereal year real-ephemeris queries inside fsolve are not free, so
# these tests use fewer random samples than the older, all-analytic tests in
# test_coordinates_space.py, consistent with the scale used elsewhere in
# this session's own additions to test_coordinates_space_orbit.py.
NUM_SAMPLES = 20


def _random_params(num=NUM_SAMPLES):
    t_ssb = numpy.random.uniform(1.0e8, 6.0e8, size=num)
    longitude = numpy.random.uniform(0.0, 2 * numpy.pi, size=num)
    latitude = numpy.random.uniform(-numpy.pi / 2, numpy.pi / 2, size=num)
    polarization = numpy.random.uniform(0.0, 2 * numpy.pi, size=num)
    return t_ssb, longitude, latitude, polarization


class TestMoonCoordinates(unittest.TestCase):
    def test_rotation_matrix_is_identity(self):
        r = moon.rotation_matrix_ssb_to_moon()
        self.assertEqual(r.shape, (3, 3))
        self.assertTrue(numpy.array_equal(r, numpy.eye(3)))

    def test_ssb_to_moon_round_trip_center(self):
        """ssb_to_moon -> moon_to_ssb recovers the original
        (tc, longitude, latitude, polarization), for the Moon's
        barycenter (no lunarsky needed).
        """
        t_ssb, lon, lat, pol = _random_params()
        t_moon, lon_m, lat_m, pol_m = moon.ssb_to_moon(t_ssb, lon, lat, pol)
        t_ssb_rt, lon_rt, lat_rt, pol_rt = moon.moon_to_ssb(
            t_moon, lon_m, lat_m, pol_m)
        self.assertLess(numpy.max(numpy.abs(t_ssb_rt - t_ssb)), 1e-2)
        self.assertLess(numpy.max(numpy.abs(lon_rt - lon)), 1e-8)
        self.assertLess(numpy.max(numpy.abs(lat_rt - lat)), 1e-8)
        self.assertLess(numpy.max(numpy.abs(pol_rt - pol)), 1e-8)

    def test_ssb_to_moon_design_lock_identity_rotation(self):
        """Because `rotation_matrix_ssb_to_moon` is the identity (see
        module docstring, and Tissino et al. 2026 arXiv:2606.04918 Eq.
        7-9), the sky-localization/polarization values must come out of
        `ssb_to_moon` numerically unchanged -- only the arrival time
        changes. This locks that design choice down against an
        accidental future change (e.g. someone adding a non-trivial
        rotation without updating this test).
        """
        t_ssb, lon, lat, pol = _random_params()
        t_moon, lon_m, lat_m, pol_m = moon.ssb_to_moon(t_ssb, lon, lat, pol)
        # the identity rotation is applied exactly, but the vector <->
        # angle round trip through localization_to_propagation_vector/
        # propagation_vector_to_localization still incurs ordinary
        # floating-point trig error, hence a tight tolerance rather than
        # exact equality.
        self.assertLess(numpy.max(numpy.abs(lon_m - lon)), 1e-10)
        self.assertLess(numpy.max(numpy.abs(lat_m - lat)), 1e-10)
        self.assertLess(numpy.max(numpy.abs(pol_m - pol)), 1e-10)
        self.assertFalse(numpy.array_equal(t_moon, t_ssb))

    def test_moon_barycenter_distance_from_earth(self):
        """The Moon's barycenter, from `moon_site_position_ssb` (or
        equivalently `space_orbit._real_body_position_velocity(t,
        'moon')`), should be within the real Earth-Moon distance range
        (~356 000-407 000 km) of `earth_position_ssb`'s real Earth
        position -- both come from the same underlying JPL ephemeris, so
        this is a strong, real-data consistency check, not just a sanity
        bound.
        """
        times = numpy.random.uniform(1.0e8, 6.0e8, size=10)
        for t in times:
            p_moon = moon.moon_site_position_ssb(t).flatten()
            p_earth = space.earth_position_ssb(t)[0].flatten()
            dist = numpy.linalg.norm(p_moon - p_earth)
            self.assertGreater(dist, 3.5e8)
            self.assertLess(dist, 4.1e8)

    def test_moon_to_geo_round_trip(self):
        t_ssb, lon, lat, pol = _random_params()
        t_moon, lon_m, lat_m, pol_m = moon.ssb_to_moon(t_ssb, lon, lat, pol)
        t_geo, lon_g, lat_g, pol_g = moon.moon_to_geo(
            t_moon, lon_m, lat_m, pol_m)
        t_moon_rt, lon_m_rt, lat_m_rt, pol_m_rt = moon.geo_to_moon(
            t_geo, lon_g, lat_g, pol_g)
        self.assertLess(numpy.max(numpy.abs(t_moon_rt - t_moon)), 1e-2)
        self.assertLess(numpy.max(numpy.abs(lon_m_rt - lon_m)), 1e-6)
        self.assertLess(numpy.max(numpy.abs(lat_m_rt - lat_m)), 1e-6)
        self.assertLess(numpy.max(numpy.abs(pol_m_rt - pol_m)), 1e-6)

    def test_moon_to_lisa_round_trip(self):
        t_ssb, lon, lat, pol = _random_params()
        t_moon, lon_m, lat_m, pol_m = moon.ssb_to_moon(t_ssb, lon, lat, pol)
        t_lisa, lon_l, lat_l, pol_l = moon.moon_to_lisa(
            t_moon, lon_m, lat_m, pol_m)
        t_moon_rt, lon_m_rt, lat_m_rt, pol_m_rt = moon.lisa_to_moon(
            t_lisa, lon_l, lat_l, pol_l)
        self.assertLess(numpy.max(numpy.abs(t_moon_rt - t_moon)), 1e-2)
        self.assertLess(numpy.max(numpy.abs(lon_m_rt - lon_m)), 1e-6)
        self.assertLess(numpy.max(numpy.abs(lat_m_rt - lat_m)), 1e-6)
        self.assertLess(numpy.max(numpy.abs(pol_m_rt - pol_m)), 1e-6)

    def test_moon_to_lisa_round_trip_with_taiji_orbit(self):
        """`orbit=` passed through moon_to_lisa/lisa_to_moon reaches
        `space.ssb_to_lisa`/`lisa_to_ssb` unchanged, so it generalizes to
        Taiji/TianQin/numeric orbits for free -- no separate
        moon_to_taiji is needed.
        """
        orbit = space_orbit.TaijiAnalyticOrbit()
        t_ssb = 2.0e7
        lon, lat, pol = 1.3, -0.4, 0.9
        t_moon, lon_m, lat_m, pol_m = moon.ssb_to_moon(t_ssb, lon, lat, pol)
        t_taiji, lon_t, lat_t, pol_t = moon.moon_to_lisa(
            t_moon, lon_m, lat_m, pol_m, orbit=orbit)
        t_moon_rt, lon_m_rt, lat_m_rt, pol_m_rt = moon.lisa_to_moon(
            t_taiji, lon_t, lat_t, pol_t, orbit=orbit)
        self.assertAlmostEqual(t_moon_rt, t_moon, places=2)
        self.assertAlmostEqual(lon_m_rt, lon_m, places=6)
        self.assertAlmostEqual(lat_m_rt, lat_m, places=6)
        self.assertAlmostEqual(pol_m_rt, pol_m, places=6)

    def test_moon_to_geo_lal_convention_flag(self):
        """`lal_convention=False` should differ from the default
        (`lal_convention=True`, unchanged from before this flag existed)
        by exactly a +/-pi polarization flip, and should still round trip
        correctly through `geo_to_moon(lal_convention=False)`. This locks
        down the `lgwa_response`-facing convention (see
        `pycbc.detector.space._LGWA_detector`) against an accidental
        future change.
        """
        t_ssb, lon, lat, pol = _random_params(num=5)
        t_moon, lon_m, lat_m, pol_m = moon.ssb_to_moon(t_ssb, lon, lat, pol)
        t_g1, lon_g1, lat_g1, pol_g1 = moon.moon_to_geo(
            t_moon, lon_m, lat_m, pol_m, lal_convention=True)
        t_g2, lon_g2, lat_g2, pol_g2 = moon.moon_to_geo(
            t_moon, lon_m, lat_m, pol_m, lal_convention=False)
        self.assertTrue(numpy.array_equal(t_g1, t_g2))
        self.assertTrue(numpy.array_equal(lon_g1, lon_g2))
        self.assertTrue(numpy.array_equal(lat_g1, lat_g2))
        flip = numpy.mod(pol_g1 - pol_g2, 2 * numpy.pi)
        self.assertLess(numpy.max(numpy.abs(flip - numpy.pi)), 1e-10)

        t_m_rt, lon_m_rt, lat_m_rt, pol_m_rt = moon.geo_to_moon(
            t_g2, lon_g2, lat_g2, pol_g2, lal_convention=False)
        self.assertLess(numpy.max(numpy.abs(t_m_rt - t_moon)), 1e-2)
        self.assertLess(numpy.max(numpy.abs(lon_m_rt - lon_m)), 1e-6)
        self.assertLess(numpy.max(numpy.abs(lat_m_rt - lat_m)), 1e-6)
        self.assertLess(numpy.max(numpy.abs(pol_m_rt - pol_m)), 1e-6)

    def test_moon_site_position_requires_both_or_neither(self):
        with self.assertRaises(ValueError):
            moon.moon_site_position_ssb(2.0e7, longitude=0.1, latitude=None)
        with self.assertRaises(ValueError):
            moon.moon_site_position_ssb(2.0e7, longitude=None, latitude=0.2)


class TestOptionalLunarskySite(unittest.TestCase):
    """Tests requiring the optional `lunarsky` package (real lunar
    libration for a specific surface site). Skipped entirely if
    `lunarsky` is not installed.
    """
    def setUp(self):
        try:
            import lunarsky  # noqa: F401
        except ImportError:
            self.skipTest('lunarsky not installed; skipping site-specific '
                          'lunar position tests')

    def test_site_position_close_to_barycenter(self):
        """A specific surface site's position must be within the Moon's
        radius (~1737 km) of the barycenter position.
        """
        t = 2.0e7
        lon_site, lat_site = numpy.deg2rad(30.0), numpy.deg2rad(-89.0)
        p_site = moon.moon_site_position_ssb(
            t, longitude=lon_site, latitude=lat_site).flatten()
        p_center = moon.moon_site_position_ssb(t).flatten()
        dist = numpy.linalg.norm(p_site - p_center)
        self.assertLess(dist, 1.8e6)

    def test_site_position_is_deterministic(self):
        t = 2.0e7
        lon_site, lat_site = numpy.deg2rad(10.0), numpy.deg2rad(45.0)
        p1 = moon.moon_site_position_ssb(t, lon_site, lat_site)
        p2 = moon.moon_site_position_ssb(t, lon_site, lat_site)
        self.assertTrue(numpy.array_equal(p1, p2))

    def test_arrival_time_site_vs_center_within_light_travel_bound(self):
        """The arrival-time difference between a specific site and the
        Moon's barycenter must be bounded by (Moon's radius)/c
        (~5.8 ms), not e.g. accidentally computed in the wrong frame
        (which would give a much larger, unphysical difference -- as
        happened during development when the ICRS/ecliptic conversion
        was first checked).
        """
        t_ssb = 2.0e7
        lon, lat, pol = 1.1, -0.2, 0.7
        lon_site, lat_site = numpy.deg2rad(30.0), numpy.deg2rad(-89.0)
        t_moon_center, _, _, _ = moon.ssb_to_moon(t_ssb, lon, lat, pol)
        t_moon_site, _, _, _ = moon.ssb_to_moon(
            t_ssb, lon, lat, pol,
            longitude_site=lon_site, latitude_site=lat_site)
        self.assertLess(abs(t_moon_site - t_moon_center), 0.01)

    def test_ssb_to_moon_round_trip_with_site(self):
        t_ssb, lon, lat, pol = _random_params(num=5)
        lon_site, lat_site = numpy.deg2rad(30.0), numpy.deg2rad(-89.0)
        t_moon, lon_m, lat_m, pol_m = moon.ssb_to_moon(
            t_ssb, lon, lat, pol,
            longitude_site=lon_site, latitude_site=lat_site)
        t_ssb_rt, lon_rt, lat_rt, pol_rt = moon.moon_to_ssb(
            t_moon, lon_m, lat_m, pol_m,
            longitude_site=lon_site, latitude_site=lat_site)
        self.assertLess(numpy.max(numpy.abs(t_ssb_rt - t_ssb)), 1e-2)
        self.assertLess(numpy.max(numpy.abs(lon_rt - lon)), 1e-6)
        self.assertLess(numpy.max(numpy.abs(lat_rt - lat)), 1e-6)


class TestRealBodyPositionVelocityRefactor(unittest.TestCase):
    """`space_orbit._real_body_position_velocity` was generalized from
    the pre-existing `_real_earth_position_velocity` (adding a `body`
    argument); this checks that generalization did not change the
    existing Earth-specific behavior, and that the new 'moon' argument
    gives a sane result.
    """
    def test_earth_unchanged_by_generalization(self):
        t = numpy.random.uniform(1e6, 3.1e7, size=10)
        pos_new, vel_new = space_orbit._real_body_position_velocity(
            t, 'earth')
        pos_wrapper, vel_wrapper = space_orbit._real_earth_position_velocity(
            t)
        self.assertTrue(numpy.array_equal(pos_new, pos_wrapper))
        self.assertTrue(numpy.array_equal(vel_new, vel_wrapper))

    def test_moon_position_velocity_sane(self):
        t = numpy.random.uniform(1e6, 3.1e7, size=10)
        pos, vel = space_orbit._real_body_position_velocity(t, 'moon')
        dist_au = numpy.linalg.norm(pos, axis=-1) / 1.495978707e11
        # the Moon is ~1 AU from the SSB (it follows the Earth-Moon
        # barycenter's heliocentric orbit), same order as Earth's own
        # distance.
        self.assertTrue(numpy.all((dist_au > 0.9) & (dist_au < 1.1)))
        speed = numpy.linalg.norm(vel, axis=-1)
        # Earth's orbital speed is ~29.8 km/s; the Moon's is close to
        # this (dominated by the same heliocentric motion, +/- its own
        # ~1 km/s orbit around Earth).
        self.assertTrue(numpy.all((speed > 2.5e4) & (speed < 3.2e4)))


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestMoonCoordinates))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestOptionalLunarskySite))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestRealBodyPositionVelocityRefactor))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
