# Copyright (C) 2026  Alex Nitz
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
"""Tests moving the locked time region to wherever the peak has gone.

The region is locked around a reference waveform. As the parameters move
away from that reference the peak moves too, and once it leaves the region
the marginalization integrates over times the signal is not at. These drive
the reference away from the signal and check the likelihood survives it.
"""

import contextlib
import copy
import logging
import unittest

import numpy

from utils import simple_exit

from pycbc.detector import Detector
from pycbc.types import TimeSeries
from pycbc.distributions import JointDistribution, SinAngle, Uniform
from pycbc.inference import models
from pycbc.inference.models.tools import DistMarg
from pycbc.noise import noise_from_psd
from pycbc.psd import aLIGOZeroDetHighPower
from pycbc.waveform import get_td_waveform

TC = 1187008882.42
FLOW, SEGLEN, SRATE = 25., 64, 2048
SAMPLE_RATE = 4096
VARIABLE = ['distance', 'inclination', 'tc']
POINT = dict(distance=40., inclination=0.5)
# how far the reference is put from the signal. The last of these moves the
# peak by more than the width of a locked region.
OFFSETS = [0.0, 0.002, 0.005, 0.01]
SEARCH = dict(peak_lock_search_samples=246, peak_lock_search_decimate=8)


class TestPeakLockSearch(unittest.TestCase):

    MODEL = 'RelativeTimeDom'

    @classmethod
    def setUpClass(cls):
        hp, hc = get_td_waveform(approximant='TaylorF2', f_lower=FLOW,
                                 delta_t=1. / SRATE, mass1=1.4, mass2=1.35,
                                 distance=40., inclination=0.5, coa_phase=1.1)
        hp.start_time += TC
        hc.start_time += TC
        psd = aLIGOZeroDetHighPower(int(SRATE * SEGLEN / 2) + 1,
                                    1. / SEGLEN, FLOW)
        # one detector: a locked region can be narrower than the light
        # travel time between two, which is a separate problem from this one
        noise = noise_from_psd(int(SEGLEN * SRATE), 1. / SRATE, psd, seed=7)
        noise._epoch = TC - SEGLEN / 2
        signal = Detector('H1').project_wave(hp, hc, 1.7, -0.4, 0.3)
        cls.data = {'H1': noise.add_into(signal).to_frequencyseries()}
        cls.psds = {'H1': psd}
        cls.static = dict(mass1=1.4, mass2=1.35, f_lower=FLOW,
                          approximant='TaylorF2', ra=1.7, dec=-0.4,
                          polarization=0.3)
        cls.fiducial = dict(mass1=1.4, tc=TC, ra=1.7, dec=-0.4,
                            polarization=0.3)

    def model(self, mass_offset=0.0, data=None, **kwargs):
        numpy.random.seed(5)
        prior = JointDistribution(
            list(VARIABLE), SinAngle(inclination=None),
            Uniform(distance=(10, 200)), Uniform(tc=(TC - 0.1, TC + 0.1)))
        fiducial = dict(self.fiducial)
        fiducial['mass1'] = fiducial['mass1'] - mass_offset
        return getattr(models, self.MODEL)(
            list(VARIABLE), copy.deepcopy(self.data if data is None else data),
            low_frequency_cutoff={'H1': FLOW}, psds=self.psds,
            static_params=self.static, prior=prior,
            fiducial_params=fiducial, epsilon=0.1, marginalize_phase=True,
            marginalize_vector_params='tc', marginalize_vector_samples=2000,
            sample_rate=SAMPLE_RATE, **kwargs)

    def loglr(self, **kwargs):
        model = self.model(**kwargs)
        model.update(**POINT)
        return model.loglr

    @contextlib.contextmanager
    def logging_enabled(self):
        """assertLogs cannot see anything through a global logging.disable.

        test_fft_unthreaded silences WARNING when it is imported and never
        puts it back, so in a whole-suite run every module after it sees no
        records at all. Restore whatever was set rather than assuming.
        """
        previous = logging.root.manager.disable
        logging.disable(logging.NOTSET)
        try:
            yield
        finally:
            logging.disable(previous)

    def test_off_by_default(self):
        """A model that did not ask for it must not move its region."""
        model = self.model(mass_offset=0.01, peak_lock_snr=4.0)
        before = dict(model.tstart)
        model.update(**POINT)
        model.loglr
        for ifo in before:
            self.assertEqual(before[ifo], model.tstart[ifo])

    def test_the_region_follows_a_moving_peak(self):
        """The answer has to survive the reference drifting away.

        Without the search the region stops holding the peak and the
        likelihood collapses; with it the answer stays near the one an
        unrestricted region gives.
        """
        exact = self.loglr()
        for offset in OFFSETS:
            searched = self.loglr(mass_offset=offset, peak_lock_snr=4.0,
                                  **SEARCH)
            self.assertLess(
                abs(searched - exact), 50.,
                "with a reference %s off, searching gave %s against %s"
                % (offset, searched, exact))

    def test_it_matters(self):
        """The largest offset must actually break the unsearched case.

        Without this the test above could pass on a region that never
        needed moving.
        """
        exact = self.loglr()
        locked = self.loglr(mass_offset=OFFSETS[-1], peak_lock_snr=4.0)
        self.assertGreater(
            abs(locked - exact), 100.,
            "a reference %s off did not spill out of the region; the test "
            "proves nothing" % OFFSETS[-1])

    def test_the_region_depends_only_on_the_parameters(self):
        """Calling twice has to land the region in the same place.

        The offset is measured from the region as locked, not from where the
        last call left it, so repeated evaluation must not walk it along.
        """
        model = self.model(mass_offset=0.01, peak_lock_snr=4.0, **SEARCH)
        model.update(**POINT)
        model.loglr
        first = dict(model.tstart)
        for _ in range(3):
            model.update(**POINT)
            model.loglr
        for ifo in first:
            self.assertAlmostEqual(first[ifo], model.tstart[ifo], places=9)

    def test_a_matched_reference_barely_moves_the_region(self):
        """With nothing to follow, the region must stay where it was.

        The offset is measured from the reference peak, so a reference that
        is the signal has none to find and the region can move by at most
        the sample it is rounded to. Measuring from the middle of the region
        instead moves it by the coarse spacing, which the reference peak
        does not sit at the centre of.
        """
        model = self.model(peak_lock_snr=4.0, **SEARCH)
        before = dict(model.tstart)
        model.update(**POINT)
        model.loglr
        for ifo in before:
            self.assertLessEqual(
                abs(model.tstart[ifo] - before[ifo]) * SAMPLE_RATE, 1.0,
                "a matched reference moved %s by %s samples"
                % (ifo, (model.tstart[ifo] - before[ifo]) * SAMPLE_RATE))

    def test_the_region_lands_on_a_sample_of_its_own_grid(self):
        """The region may only move by whole samples.

        A fractional offset would leave the marginalization sampling the
        series at shifted phases, which costs more than placing the region
        more precisely gains. The reference is the signal here: that is
        where the coarse peak has a genuine fraction of a sample to it,
        rather than happening to land on one.
        """
        model = self.model(peak_lock_snr=4.0, **SEARCH)
        before = dict(model.tstart)
        model.update(**POINT)
        model.loglr
        for ifo in before:
            samples = (model.tstart[ifo] - before[ifo]) * SAMPLE_RATE
            self.assertAlmostEqual(
                samples, round(samples), places=6,
                msg="%s moved by %s samples" % (ifo, samples))

    def test_the_stride_is_shorter_than_the_region(self):
        """The search can only place the peak to within its own stride.

        A stride wider than the region can leave the peak outside it even
        though the search found it, which is the failure the search exists
        to prevent. The region here is 10 samples against a stride of 8;
        striding by 32 was measured to cost 157 nats.
        """
        model = self.model(peak_lock_snr=4.0, **SEARCH)
        for ifo in model.num_samples:
            self.assertLessEqual(
                SEARCH['peak_lock_search_decimate'], model.num_samples[ifo],
                "stride %s against a region of %s samples in %s"
                % (SEARCH['peak_lock_search_decimate'],
                   model.num_samples[ifo], ifo))

    def test_precalculated_points_hold_the_region_still(self):
        """Points drawn over the locked region cannot outlive it moving.

        With precalculate_marginalization_points the draw is a subset of
        points chosen once, over the region as locked. Moving the region
        after that leaves the points describing somewhere the model is no
        longer looking, and nothing about the result says so.
        """
        model = self.model(mass_offset=0.01, peak_lock_snr=4.0, **SEARCH)
        locked = dict(model.tstart)
        model.update(**POINT)
        model.loglr
        self.assertNotEqual(model.tstart['H1'], locked['H1'],
                            "the region never moved, so this proves nothing")

        # the waveforms the search would be given, and the region back where
        # it was locked, so the only thing left to stop it is the guard
        wfs = model.get_waveforms(model.current_params)
        model.tstart = dict(locked)
        model.premarg = {}
        model.follow_peak(wfs)
        self.assertEqual(model.tstart, locked,
                         "the region moved out from under precalculated "
                         "points")

    def test_regions_with_no_time_in_common_do_not_go_negative(self):
        """Locking narrowly can leave the detectors nothing to agree on.

        The regions are intersected pairwise, widened by the light travel
        time, so that a time in one has a sky position consistent with the
        others. Peaks further apart than that leave the intersection empty,
        and an empty region is a negative sample count that reaches the
        kernel as ValueError: negative dimensions are not allowed.
        """
        rate = 4096.0
        span = 0.1
        peaks = {'H1': 0.02, 'L1': 0.08}          # 60 ms apart, H1-L1 is 10
        snrs = {}
        for ifo, at in peaks.items():
            t = numpy.arange(int(span * rate)) / rate
            z = 30.0 * numpy.exp(-0.5 * ((t - at) / 2e-4) ** 2.0)
            snrs[ifo] = TimeSeries(z + 0j, delta_t=1.0 / rate, epoch=TC)

        class Two(DistMarg):
            data = {'H1': None, 'L1': None}
            marginalized_vector_priors = {
                'tc': Uniform(tc=(TC, TC + span))}

        two = Two()
        with self.logging_enabled(), \
                self.assertLogs(level='WARNING') as caught:
            two.setup_peak_lock(sample_rate=rate, snrs=snrs,
                                peak_lock_snr=5.0, peak_lock_ratio=20,
                                peak_lock_region=2)
        self.assertIn('no time in common', caught.records[0].getMessage())
        for ifo in peaks:
            self.assertGreater(two.num_samples[ifo], 0,
                               "%s region is %s samples"
                               % (ifo, two.num_samples[ifo]))

    def test_the_offset_comes_from_the_loudest_detector(self):
        """One detector supplies the offset, so it should be the best one.

        The peak of a louder detector is the less uncertain in time, and
        it is the one the likelihood is most sensitive to.
        """
        rate = 4096.0
        span = 0.1
        snrs, asked = {}, []
        for ifo, at, amp in (('H1', 0.04, 12.0), ('L1', 0.05, 30.0)):
            t = numpy.arange(int(span * rate)) / rate
            z = amp * numpy.exp(-0.5 * ((t - at) / 2e-4) ** 2.0)
            snrs[ifo] = TimeSeries(z + 0j, delta_t=1.0 / rate, epoch=TC)

        class Two(DistMarg):
            data = {'H1': None, 'L1': None}
            marginalized_vector_priors = {'tc': Uniform(tc=(TC, TC + span))}

            def coarse_series(self, ifo, wfs, tstart, delta_t, num_samples):
                asked.append(ifo)
                return numpy.zeros(num_samples)

        two = Two()
        two.setup_peak_lock(sample_rate=rate, snrs=snrs, peak_lock_snr=5.0,
                            peak_lock_ratio=1e4, **SEARCH)
        two.follow_peak({'H1': None, 'L1': None})
        self.assertEqual(asked, ['L1'],
                         "the offset was taken from %s, not the loudest"
                         % asked)

    def test_a_model_that_cannot_search_says_so(self):
        """Silently ignoring the option would be worse than not having it.

        The search needs the model to evaluate its predictor on a grid it
        chooses. A model whose reference series is fixed has no such method
        and no need of one -- SingleTemplate is one, and takes these same
        options through the same call -- but a run that asked for the
        search still has to be told it is not getting it.
        """
        self.assertFalse(hasattr(models.SingleTemplate, 'coarse_series'))

        class Bare(DistMarg):
            marginalized_vector_priors = {}

        bare = Bare()
        with self.logging_enabled(), \
                self.assertLogs(level='WARNING') as caught:
            bare.setup_peak_lock(sample_rate=SAMPLE_RATE, **SEARCH)
        self.assertIsNone(bare.peak_lock_search_samples)
        self.assertIn('peak_lock_search_samples',
                      caught.records[0].getMessage())


class TestPeakLockSearchRelativeTime(TestPeakLockSearch):
    """The same, through the class that keeps hp and hc apart.

    It has its own coarse_series, calling a different predictor, so the
    search is a separate path there and not covered by the one above.
    """

    MODEL = 'RelativeTime'


suite = unittest.TestSuite()
# not a loop variable: anything left bound at module scope here is a
# TestCase subclass, which pytest collects and runs a second time
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestPeakLockSearch))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestPeakLockSearchRelativeTime))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
