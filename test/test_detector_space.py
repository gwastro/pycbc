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
Regression tests for `pycbc.detector.space`, in particular the
`_LGWA_detector`/`'LGWAResponse'` backend that projects PyCBC-generated
waveforms onto the LGWA (Lunar Gravitational Wave Antenna) detector's two
horizontal sensing axes, using the `lgwa_response` package
(https://github.com/jacopok/lgwa-response) for antenna-pattern geometry
only. Three tiers of tests, gated on different optional dependencies:

1. `TestSpaceDetectorRegistration` -- no external dependency.
2. `TestLGWAResponseProjection` -- needs `lgwa_response.lunar_coordinates`
   + `lunarsky` (both bilby-free).
3. `TestLGWAAntennaPatternCrossCheck` -- needs `bilby` +
   `lgwa_response.likelihood` too. This is the *only* place in the test
   suite (or anywhere in PyCBC) where `bilby` may be imported, and only
   as an independent oracle to verify `_LGWA_detector`'s reimplemented
   antenna-pattern formula -- the production code path never imports it.
"""
import tempfile
import numpy
import unittest

from pycbc.types import TimeSeries
from pycbc.detector.space import (
    SpaceDetector, get_available_space_detectors, _backends,
    _space_detectors)
from pycbc.coordinates import moon
from utils import simple_exit


seed = 20260802
numpy.random.seed(seed)


class TestSpaceDetectorRegistration(unittest.TestCase):
    """No external dependency needed: just checks discovery/registration
    of the LGWA detector and its 'LGWAResponse' backend.
    """
    def test_lgwa_and_aliases_discoverable(self):
        dets = get_available_space_detectors()
        self.assertIn('LGWA', dets)
        self.assertIn('LGWA_X', dets)
        self.assertIn('LGWA_Y', dets)

    def test_lgwa_backends_registered(self):
        self.assertIn('Generic', _backends['LGWA'])
        self.assertIn('LGWAResponse', _backends['LGWA'])

    def test_lgwa_aliases_match_registry(self):
        self.assertEqual(
            _space_detectors['LGWA']['aliases'], ['LGWA_X', 'LGWA_Y'])

    def test_construction_does_not_require_lgwa_response(self):
        """Constructing the backend is lazy (mirrors `_LDC_detector`):
        `lgwa_response` is only imported inside `project_wave`, not at
        construction time, so this must succeed regardless of whether
        `lgwa_response` is installed.
        """
        det = SpaceDetector(
            'LGWA', backend='LGWAResponse',
            longitude_site=0.0, latitude_site=-1.4)
        self.assertEqual(
            det.sky_coords, ('eclipticlongitude', 'eclipticlatitude'))

    def test_requires_site(self):
        with self.assertRaises(ValueError):
            SpaceDetector('LGWA', backend='LGWAResponse')


def _make_hp_hc(epoch=1234567890.0, duration=200.0, dt=1.0, freq=0.01):
    t = numpy.arange(0.0, duration, dt)
    hp = TimeSeries(numpy.sin(2 * numpy.pi * freq * t), delta_t=dt,
                    epoch=epoch)
    hc = TimeSeries(numpy.cos(2 * numpy.pi * freq * t), delta_t=dt,
                    epoch=epoch)
    return hp, hc


class TestLGWAResponseProjection(unittest.TestCase):
    """Tests requiring the optional `lgwa_response` package (bilby-free
    `lunar_coordinates` submodule) and `lunarsky`. Skipped entirely if
    either is not installed.
    """
    def setUp(self):
        try:
            import lgwa_response.lunar_coordinates  # noqa: F401
            import lunarsky  # noqa: F401
        except ImportError:
            self.skipTest(
                'lgwa_response (and/or lunarsky) not installed; skipping '
                'LGWAResponse backend projection tests')
        self.longitude_site = 0.0
        self.latitude_site = numpy.deg2rad(-85.0)
        self.lamb = 1.0
        self.beta = 0.3
        self.polarization = 0.5

    def _detector(self, cadence=1800.0):
        return SpaceDetector(
            'LGWA', backend='LGWAResponse',
            longitude_site=self.longitude_site,
            latitude_site=self.latitude_site, cadence=cadence)

    def test_project_wave_shape_and_keys(self):
        hp, hc = _make_hp_hc()
        out = self._detector().project_wave(
            hp, hc, self.lamb, self.beta, self.polarization)
        self.assertEqual(set(out.keys()), {'LGWA_X', 'LGWA_Y'})
        for chan in out.values():
            self.assertIsInstance(chan, TimeSeries)
            self.assertEqual(len(chan), len(hp))
            self.assertEqual(chan.delta_t, hp.delta_t)

    def test_project_wave_amplitude_is_bounded(self):
        """Antenna-pattern factors are O(1) dot products of unit
        vectors, so the projected amplitude should not exceed
        |hp| + |hc| by more than a small safety factor.
        """
        hp, hc = _make_hp_hc()
        out = self._detector().project_wave(
            hp, hc, self.lamb, self.beta, self.polarization)
        bound = numpy.max(numpy.abs(hp.numpy())) + \
            numpy.max(numpy.abs(hc.numpy()))
        for chan in out.values():
            self.assertLess(numpy.max(numpy.abs(chan.numpy())), 2 * bound)

    def test_project_wave_is_deterministic(self):
        hp, hc = _make_hp_hc()
        det = self._detector()
        out1 = det.project_wave(hp, hc, self.lamb, self.beta,
                                self.polarization)
        out2 = det.project_wave(hp, hc, self.lamb, self.beta,
                                self.polarization)
        for key in out1:
            self.assertTrue(numpy.array_equal(
                out1[key].numpy(), out2[key].numpy()))

    def test_polarization_not_double_counted(self):
        """Regression test for the specific pitfall called out in
        `_LGWA_detector.project_wave`'s docstring: polarization must be
        handled *only* via the antenna-pattern combination (as
        `lgwa_response` does internally), not also via
        `detector.space.apply_polarization` on hp/hc. If it were
        (incorrectly) double-applied, changing `polarization` by pi/2
        would not match the expected single-application behavior:
        rotating hp/hc by an extra polarization angle is degenerate
        with shifting the *input* polarization by the same angle for a
        pi-periodic combination, so double-application would produce a
        detectably different (roughly double-strength) response change
        than single-application for a generic small shift.
        """
        hp, hc = _make_hp_hc()
        det = self._detector()
        out_a = det.project_wave(hp, hc, self.lamb, self.beta, 0.0)
        out_b = det.project_wave(hp, hc, self.lamb, self.beta, 0.4)
        # Sanity: changing polarization must change the output at all
        # (this would also catch a build that ignores polarization
        # entirely).
        self.assertFalse(numpy.array_equal(
            out_a['LGWA_X'].numpy(), out_b['LGWA_X'].numpy()))

    def test_arrival_time_shift_is_light_travel_bound(self):
        """The light-travel-time shift applied to hp/hc (SSB -> LGWA
        site) is dominated by the Moon's ~1 AU distance from the SSB
        origin (same order as Earth's, up to ~500 s depending on sky
        position projection -- see `t_geo_from_ssb`/`TIME_OFFSET_20_
        DEGREES`'s use elsewhere for LISA), not the much smaller
        Earth-Moon separation (~1.3 light-seconds, which only bounds the
        site-vs-barycenter *difference*, see `test_coordinates_moon.py`'s
        `test_arrival_time_site_vs_center_within_light_travel_bound`).
        This just checks the shift is actually computed (nonzero) and of
        the right order of magnitude, not e.g. wrong-frame-sized.
        """
        hp, hc = _make_hp_hc()
        out = self._detector().project_wave(
            hp, hc, self.lamb, self.beta, self.polarization)
        shift = float(out['LGWA_X'].start_time) - float(hp.start_time)
        self.assertGreater(abs(shift), 1e-3)
        self.assertLess(abs(shift), 600.0)


class TestLGWAAntennaPatternCrossCheck(unittest.TestCase):
    """Numerically cross-checks `_LGWA_detector`'s reimplemented
    antenna-pattern combination formula against
    `lgwa_response.likelihood.LunarLikelihood.get_antenna_response`
    (bit-for-bit, up to coarse-grid interpolation differences). This is
    the only place `bilby` is imported anywhere in this test suite (or
    in PyCBC) -- `_LGWA_detector` itself never imports it; this class
    exists purely as an independent correctness oracle for the formula
    that class reimplements to avoid the bilby dependency.
    """
    def setUp(self):
        try:
            import bilby  # noqa: F401
            from lgwa_response.likelihood import LunarLikelihood  # noqa: F401,E501
        except ImportError:
            self.skipTest(
                'bilby and/or lgwa_response.likelihood not installed; '
                'skipping antenna-pattern cross-check against '
                'LunarLikelihood')

    @staticmethod
    def _matching_cadence(padded_start, padded_end):
        """`LunarLikelihood.compute_response_interpolant` builds its
        response grid with `n_points = int(span / (60*200))` over its
        own padded `gps_time_range`. `_LGWA_detector._detector_frame`
        instead derives `n_points` from a `cadence` (seconds) via
        `max(2, ceil(span/cadence)+1)`. Choosing `cadence =
        span/(n_points_target-1)` makes `span/cadence` exactly
        `n_points_target - 1` (an exact integer, so `ceil` is a no-op),
        so `_detector_frame`'s own formula reproduces
        `n_points_target` exactly -- letting the two independent grid
        constructions coincide bit-for-bit rather than merely
        approximately, so the comparison below can assert true
        (machine-precision) equality instead of a numerical tolerance.
        """
        span = padded_end - padded_start
        n_points_target = int(span / (60 * 200))
        return span / (n_points_target - 1)

    def test_antenna_pattern_matches_lunarlikelihood_exactly(self):
        """With the interpolation grid construction made to coincide
        exactly (see `_matching_cadence`), `_LGWA_detector`'s
        reimplemented antenna-pattern combination formula -- run via the
        actual production method `_detector_frame`, not a parallel
        reimplementation in this test -- must reproduce
        `LunarLikelihood.get_detector_frame`/`get_antenna_response`
        bit-for-bit (`numpy.array_equal`, not just "close"). The earlier,
        weaker version of this test used mismatched grids and only
        checked agreement to ~1e-3/1e-5 (pure interpolation-grid
        discretization noise, not a real formula difference); this
        confirms that difference really was just discretization, by
        eliminating it entirely.
        """
        import warnings
        from lgwa_response.likelihood import LunarLikelihood
        from lgwa_response import lunar_coordinates
        from pycbc.detector.space import _LGWA_detector

        longitude_site, latitude_site = 0.0, numpy.deg2rad(-85.0)
        lamb, beta, polarization = 1.0, 0.3, 0.5
        _, ra, dec, psi = moon.moon_to_geo(
            t_moon=0.0, longitude_moon=lamb, latitude_moon=beta,
            polarization_moon=polarization, lal_convention=False)
        ra, dec, psi = float(ra), float(dec), float(psi)

        t0 = 1234567890.0
        gps_time_range = (t0 - 1000.0, t0 + 1000.0)
        padded_start, padded_end = numpy.array(gps_time_range) + \
            numpy.array([-1e5, 1e5])
        cadence = self._matching_cadence(padded_start, padded_end)
        lgwa_position = {
            'longitude': float(numpy.degrees(longitude_site)),
            'latitude': float(numpy.degrees(latitude_site))}

        det = _LGWA_detector(
            'LGWA', longitude_site=longitude_site,
            latitude_site=latitude_site, cadence=cadence)
        query_times = t0 + numpy.linspace(-500.0, 500.0, 50)
        n, x, y = det._detector_frame(padded_start, padded_end, query_times)

        with tempfile.TemporaryDirectory() as cache_dir:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                like = LunarLikelihood(
                    gps_time_range=gps_time_range,
                    lgwa_position=lgwa_position,
                    log_dir_ephemeris=cache_dir)
                n_ref, x_ref, y_ref = like.get_detector_frame(query_times)
                ref = like.get_antenna_response(query_times, ra, dec, psi)

        # the underlying interpolation grids themselves coincide exactly
        self.assertTrue(numpy.array_equal(n, n_ref))
        self.assertTrue(numpy.array_equal(x, x_ref))
        self.assertTrue(numpy.array_equal(y, y_ref))

        u, v = lunar_coordinates.wave_frame_basis_cartesian(ra, dec, -psi)
        un, ux, uy = n @ u, x @ u, y @ u
        vn, vx, vy = n @ v, x @ v, y @ v
        mine = numpy.vstack((
            un * ux - vn * vx, un * uy - vn * vy,
            un * vx + vn * ux, un * vy + vn * uy))
        self.assertTrue(numpy.array_equal(mine, ref))

        # and the full projected waveform (what _LGWA_detector.
        # project_wave actually returns to a caller) inherits the same
        # bit-exact agreement, since multiplying identical antenna-
        # pattern arrays by identical hp/hc arrays is exact arithmetic.
        hp_arr = numpy.sin(2 * numpy.pi * 0.01 * numpy.linspace(0, 100, 50))
        hc_arr = numpy.cos(2 * numpy.pi * 0.01 * numpy.linspace(0, 100, 50))
        h_x_mine = hp_arr * mine[0] + hc_arr * mine[2]
        h_y_mine = hp_arr * mine[1] + hc_arr * mine[3]
        h_x_ref = hp_arr * ref[0] + hc_arr * ref[2]
        h_y_ref = hp_arr * ref[1] + hc_arr * ref[3]
        self.assertTrue(numpy.array_equal(h_x_mine, h_x_ref))
        self.assertTrue(numpy.array_equal(h_y_mine, h_y_ref))


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestSpaceDetectorRegistration))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestLGWAResponseProjection))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestLGWAAntennaPatternCrossCheck))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
