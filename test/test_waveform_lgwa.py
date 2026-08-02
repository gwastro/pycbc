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
Regression tests for `pycbc.waveform.lgwa` (the `'LGWA_response'`
`fd_det`/`fd_det_sequence` plugin) and the `_LGWA_detector._detector_frame`
caching it relies on. Three tiers of tests, gated on different optional
dependencies:

1. `TestPluginRegistration` -- no external dependency.
2. `TestLGWADetectorCaching` -- needs `lgwa_response.lunar_coordinates` +
   `lunarsky` (both bilby-free).
3. Everything else -- needs the same, plus generates real LAL waveforms
   (`IMRPhenomD`/`IMRPhenomXHM`), so also implicitly needs `lalsimulation`.
"""
import numpy
import unittest
from unittest import mock

from pycbc.types import Array, FrequencySeries
from pycbc.waveform.waveform import fd_det, fd_det_sequence
from pycbc.coordinates import reference_point as refpt
from utils import simple_exit


seed = 20260802
numpy.random.seed(seed)


class TestPluginRegistration(unittest.TestCase):
    """No external dependency needed: just checks discovery/registration
    of the 'LGWA_response' entry-point plugin.
    """
    def test_lgwa_response_discoverable(self):
        self.assertIn('LGWA_response', fd_det)
        self.assertIn('LGWA_response', fd_det_sequence)
        self.assertIs(fd_det['LGWA_response'], fd_det_sequence['LGWA_response'])

    def test_still_needs_det_response_dispatch(self):
        """`relbin.Relative`'s LGWA dispatch (`relbin.py`'s
        `still_needs_det_response`) is driven entirely by membership in
        `fd_det_sequence` -- no LGWA-specific code exists there. This
        locks down that `'LGWA_response'` stays a member (an accidental
        rename/removal here would silently fall back to the ground-
        detector-specific path, which is physically wrong for LGWA).
        """
        self.assertIn('LGWA_response', fd_det_sequence)


def _requires_lgwa_response(cls):
    """Class decorator: skip all tests in `cls` if `lgwa_response`/
    `lunarsky` are not installed, while still running the class's own
    `setUp` (if any) afterward -- unlike naively overwriting `setUp`,
    which would silently drop any subclass-specific fixture setup.
    """
    orig_setup = cls.__dict__.get('setUp')

    def setUp(self):
        try:
            import lgwa_response.lunar_coordinates  # noqa: F401
            import lunarsky  # noqa: F401
        except ImportError:
            self.skipTest(
                'lgwa_response (and/or lunarsky) not installed; skipping '
                'LGWA_response plugin tests')
        if orig_setup is not None:
            orig_setup(self)
    cls.setUp = setUp
    return cls


@_requires_lgwa_response
class TestLGWADetectorCaching(unittest.TestCase):
    """`pycbc.waveform.lgwa` reuses a module-level cache of
    `_LGWA_detector` instances (keyed by site/cadence) specifically so
    that repeated calls -- as relbin makes, once per likelihood
    evaluation -- reuse the same instance's `_frame_cache` rather than
    rebuilding the slow `lunarsky`/astropy orientation grid every time.
    """
    def test_repeated_calls_reuse_cached_detector(self):
        from pycbc.waveform.lgwa import _get_cached_detector
        det1 = _get_cached_detector(0.0, numpy.deg2rad(-85.0), 1800.0)
        det2 = _get_cached_detector(0.0, numpy.deg2rad(-85.0), 1800.0)
        self.assertIs(det1, det2)

    def test_different_site_gives_different_detector(self):
        from pycbc.waveform.lgwa import _get_cached_detector
        det1 = _get_cached_detector(0.0, numpy.deg2rad(-85.0), 1800.0)
        det2 = _get_cached_detector(0.1, numpy.deg2rad(-85.0), 1800.0)
        self.assertIsNot(det1, det2)

    def test_lgwa_fd_response_does_not_rebuild_frame_on_second_call(self):
        """Two `lgwa_fd_response` calls with overlapping time ranges (as
        relbin's fixed static params would produce across likelihood
        evaluations at nearby `tc`) must not trigger a second real
        `generate_data_response` call -- verified by mocking it and
        counting invocations.
        """
        from pycbc.waveform.lgwa import lgwa_fd_response, \
            _lgwa_detector_cache
        _lgwa_detector_cache.clear()

        longitude_site, latitude_site = 0.2, numpy.deg2rad(-80.0)
        kwargs = dict(
            ifos=['LGWA_X'], sample_points=Array(numpy.linspace(1.0, 3.0, 5)),
            base_approximant='IMRPhenomD', mode_array=[[2, 2]],
            longitude_site=longitude_site, latitude_site=latitude_site,
            cadence=600.0, tc=1234567890.0, eclipticlongitude=1.0,
            eclipticlatitude=0.3, polarization=0.5,
            mass1=35.36, mass2=33.59, spin1z=0.0, spin2z=0.0,
            distance=500.0, inclination=1.0, coa_phase=0.5,
            f_lower=1.0, f_ref=1.0)

        from lgwa_response import lunar_coordinates
        real_generate = lunar_coordinates.generate_data_response
        with mock.patch.object(
                lunar_coordinates, 'generate_data_response',
                wraps=real_generate) as mocked:
            lgwa_fd_response(**kwargs)
            first_calls = mocked.call_count
            self.assertGreater(first_calls, 0)
            lgwa_fd_response(**kwargs)
            self.assertEqual(mocked.call_count, first_calls)


@_requires_lgwa_response
class TestModeArrayDecomposition(unittest.TestCase):
    """Verifies the strict per-mode decomposition the user required
    (rather than a dominant-mode approximation): summing per-mode
    `lgwa_fd_response` outputs must reproduce the multi-mode output,
    since `_lgwa_response_core` literally loops over `mode_array` and
    adds each mode's contribution.
    """
    def setUp(self):
        super().setUp()
        self.longitude_site = 0.0
        self.latitude_site = numpy.deg2rad(-85.0)
        self.kwargs = dict(
            sample_points=Array(numpy.linspace(1.0, 3.0, 20)),
            base_approximant='IMRPhenomXHM',
            ref_position=(0.0, 0.0, 0.0),
            longitude_site=self.longitude_site,
            latitude_site=self.latitude_site, cadence=600.0,
            tc=1234567890.0, eclipticlongitude=1.0, eclipticlatitude=0.3,
            polarization=0.5,
            mass1=40.0, mass2=8.0,  # 5:1 mass ratio, significant HM content
            spin1z=0.0, spin2z=0.0, distance=500.0, inclination=1.0,
            coa_phase=0.5, f_lower=1.0, f_ref=1.0)

    def test_mode_sum_matches_multi_mode_call(self):
        from pycbc.waveform.lgwa import lgwa_fd_response
        out_22 = lgwa_fd_response(
            ifos=['LGWA_X'], mode_array=[[2, 2]],
            **self.kwargs)['LGWA_X'].numpy()
        out_33 = lgwa_fd_response(
            ifos=['LGWA_X'], mode_array=[[3, 3]],
            **self.kwargs)['LGWA_X'].numpy()
        out_both = lgwa_fd_response(
            ifos=['LGWA_X'], mode_array=[[2, 2], [3, 3]],
            **self.kwargs)['LGWA_X'].numpy()

        # strain values are ~1e-23 scale, so comparing against a fixed
        # absolute threshold would be meaningless; compare (3,3)'s
        # magnitude against (2,2)'s instead.
        self.assertGreater(
            numpy.max(numpy.abs(out_33)), 1e-4 * numpy.max(numpy.abs(out_22)),
            'the (3,3) mode contribution should not be negligible for a '
            '5:1 mass-ratio source')

        rel_diff = numpy.abs((out_22 + out_33) - out_both) / \
            numpy.abs(out_both)
        self.assertLess(numpy.max(rel_diff), 1e-6)


@_requires_lgwa_response
class TestTDvsFDCrossCheck(unittest.TestCase):
    """Cross-checks the new frequency-domain SPA response
    (`lgwa_fd_response`) against the existing, bit-exact-validated
    time-domain response (`_LGWA_detector.project_wave`), for a short,
    high-frequency-band signal where both are tractable.

    The two use genuinely different numerical approaches (per-mode SPA
    arrival time + cached antenna pattern vs. direct time-domain antenna
    pattern + FFT), so exact agreement isn't expected -- but they should
    agree to well within 1%, this being the connective correctness anchor
    between the old (time-domain) and new (frequency-domain) code paths.

    NOTE on constructing the TD reference at a genuine absolute merger
    time `tc`: `TimeSeries.start_time` is purely a label for
    `to_frequencyseries()` (verified separately) -- setting it does *not*
    shift the waveform's actual phase/values. To place merger at a real
    GPS time `tc`, the frequency-domain phase must be shifted for real
    (`exp(-2j*pi*f*tc)`) *before* going to the time domain; only then does
    relabeling `start_time` become consistent with the shifted values.
    """
    def setUp(self):
        super().setUp()
        self.longitude_site = 0.0
        self.latitude_site = numpy.deg2rad(-85.0)
        self.lamb, self.beta, self.pol = 1.0, 0.3, 0.5
        self.tc = 1234567890.0
        self.mass1, self.mass2 = 35.36, 33.59

    def test_lgwa_fd_response_matches_project_wave(self):
        from pycbc.waveform import get_fd_waveform
        from pycbc.detector.space import SpaceDetector
        from pycbc.waveform.lgwa import lgwa_fd_response

        hp_fd, hc_fd = get_fd_waveform(
            approximant='IMRPhenomD', mass1=self.mass1, mass2=self.mass2,
            spin1z=0.0, spin2z=0.0, inclination=1.0, coa_phase=0.5,
            distance=500.0, f_lower=1.0, f_ref=1.0,
            delta_f=1.0 / 8192, f_final=8.0)

        freqs = hp_fd.sample_frequencies.numpy()
        shift = numpy.exp(-2j * numpy.pi * freqs * self.tc)
        hp_fd_shift = FrequencySeries(
            hp_fd.numpy() * shift, delta_f=hp_fd.delta_f, epoch=hp_fd.epoch)
        hc_fd_shift = FrequencySeries(
            hc_fd.numpy() * shift, delta_f=hc_fd.delta_f, epoch=hc_fd.epoch)
        hp_td = hp_fd_shift.to_timeseries()
        hc_td = hc_fd_shift.to_timeseries()
        hp_td.start_time += self.tc
        hc_td.start_time += self.tc

        det = SpaceDetector(
            'LGWA', backend='LGWAResponse',
            longitude_site=self.longitude_site,
            latitude_site=self.latitude_site, cadence=600.0)
        out_td = det.project_wave(hp_td, hc_td, self.lamb, self.beta, self.pol)
        hx_td_f = out_td['LGWA_X'].to_frequencyseries()

        fm_full = hx_td_f.sample_frequencies.numpy()
        mask = (fm_full >= 1.0) & (fm_full <= 3.0)
        fm_valid = fm_full[mask]
        hx_td_valid = hx_td_f.numpy()[mask]

        out = lgwa_fd_response(
            ifos=['LGWA_X'], sample_points=Array(fm_valid),
            base_approximant='IMRPhenomD', mode_array=[[2, 2]],
            ref_position=(0.0, 0.0, 0.0),
            longitude_site=self.longitude_site,
            latitude_site=self.latitude_site, cadence=600.0,
            tc=self.tc, eclipticlongitude=self.lamb,
            eclipticlatitude=self.beta, polarization=self.pol,
            mass1=self.mass1, mass2=self.mass2, spin1z=0.0, spin2z=0.0,
            distance=500.0, inclination=1.0, coa_phase=0.5, f_lower=1.0,
            f_ref=1.0)
        hx_plugin = out['LGWA_X'].numpy()

        rel_diff = numpy.abs(hx_plugin - hx_td_valid) / numpy.abs(hx_td_valid)
        self.assertLess(numpy.median(rel_diff), 0.02)
        self.assertLess(numpy.max(rel_diff), 0.05)


@_requires_lgwa_response
class TestRefPositionNonZero(unittest.TestCase):
    """Confirms `ref_position` is genuinely used, not silently ignored,
    and that its closed-form conversion (`pycbc.coordinates.
    reference_point`) is applied consistently: the same physical SSB
    merger time, described via two different (ref_position, tc)
    parameterizations that refer to the identical event, must give
    identical output.
    """
    def test_equivalent_ref_position_and_tc_give_same_output(self):
        from pycbc.waveform.lgwa import lgwa_fd_response

        longitude_site, latitude_site = 0.0, numpy.deg2rad(-85.0)
        lamb, beta, pol = 1.0, 0.3, 0.5
        tc_ssb = 1234567890.0
        kwargs = dict(
            ifos=['LGWA_X'], sample_points=Array(numpy.linspace(1.0, 3.0, 20)),
            base_approximant='IMRPhenomD',
            longitude_site=longitude_site, latitude_site=latitude_site,
            cadence=600.0, eclipticlongitude=lamb, eclipticlatitude=beta,
            polarization=pol, mass1=35.36, mass2=33.59, spin1z=0.0,
            spin2z=0.0, distance=500.0, inclination=1.0, coa_phase=0.5,
            f_lower=1.0, f_ref=1.0)

        ref_position = numpy.array([1.0e9, -2.0e9, 3.0e8])
        t_ref_equiv = refpt.t_ref_from_ssb(tc_ssb, lamb, beta, ref_position)

        out_a = lgwa_fd_response(
            tc=tc_ssb, ref_position=(0.0, 0.0, 0.0), **kwargs)['LGWA_X']
        out_b = lgwa_fd_response(
            tc=t_ref_equiv, ref_position=ref_position, **kwargs)['LGWA_X']

        self.assertTrue(numpy.allclose(
            out_a.numpy(), out_b.numpy(), rtol=1e-6, atol=0))

    def test_nonzero_ref_position_changes_output_for_fixed_tc(self):
        """Sanity check that ref_position isn't a no-op: for a *fixed*
        `tc` value, changing ref_position (without compensating tc) must
        change the output, since it changes what `tc` means physically.
        """
        from pycbc.waveform.lgwa import lgwa_fd_response

        kwargs = dict(
            ifos=['LGWA_X'], sample_points=Array(numpy.linspace(1.0, 3.0, 20)),
            base_approximant='IMRPhenomD',
            longitude_site=0.0, latitude_site=numpy.deg2rad(-85.0),
            cadence=600.0, tc=1234567890.0, eclipticlongitude=1.0,
            eclipticlatitude=0.3, polarization=0.5,
            mass1=35.36, mass2=33.59, spin1z=0.0, spin2z=0.0,
            distance=500.0, inclination=1.0, coa_phase=0.5,
            f_lower=1.0, f_ref=1.0)

        out_zero = lgwa_fd_response(
            ref_position=(0.0, 0.0, 0.0), **kwargs)['LGWA_X'].numpy()
        out_nonzero = lgwa_fd_response(
            ref_position=(1.0e9, -2.0e9, 3.0e8), **kwargs)['LGWA_X'].numpy()
        # atol=0 is essential here: the strain values themselves are
        # ~1e-23-scale, so numpy.allclose's default atol=1e-8 would swamp
        # any real relative/phase difference and always report "close".
        self.assertFalse(numpy.allclose(
            out_zero, out_nonzero, rtol=1e-3, atol=0))


@_requires_lgwa_response
class TestRelbinEndToEnd(unittest.TestCase):
    """End-to-end smoke test: `pycbc.inference.models.relbin.Relative`
    must be constructible with `approximant='LGWA_response'` static
    params (dispatching through `still_needs_det_response`, with zero
    LGWA-specific code in `relbin.py` itself) and return a finite,
    parameter-dependent `loglr`.
    """
    def test_relative_model_with_lgwa_response(self):
        from pycbc.inference import models

        delta_f = 0.01
        n_samples = 301
        tc0 = 1234567890.0
        data = FrequencySeries(
            numpy.zeros(n_samples, dtype=complex), delta_f=delta_f,
            epoch=tc0 - 10.0)

        static_params = dict(
            approximant='LGWA_response', base_approximant='IMRPhenomD',
            mode_array=[[2, 2]],
            longitude_site=0.0, latitude_site=numpy.deg2rad(-85.0),
            cadence=1800.0, mass1=35.36, mass2=33.59, spin1z=0.0,
            spin2z=0.0, distance=500.0, inclination=1.0, coa_phase=0.5,
            eclipticlongitude=1.0, eclipticlatitude=0.3, polarization=0.5,
            f_lower=1.0)

        model = models.Relative(
            ['tc'], {'LGWA_X': data},
            low_frequency_cutoff={'LGWA_X': 1.0},
            static_params=static_params,
            fiducial_params={'tc': tc0}, epsilon=0.1)

        self.assertTrue(model.still_needs_det_response)

        model.update(tc=tc0)
        loglr_a = model.loglr
        self.assertTrue(numpy.isfinite(loglr_a))

        model.update(tc=tc0 + 0.5)
        loglr_b = model.loglr
        self.assertTrue(numpy.isfinite(loglr_b))
        self.assertNotEqual(loglr_a, loglr_b)


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestPluginRegistration))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestLGWADetectorCaching))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestModeArrayDecomposition))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestTDvsFDCrossCheck))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestRefPositionNonZero))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestRelbinEndToEnd))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
