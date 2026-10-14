# Copyright (C) 2012  Alex Nitz, Josh Willis
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

import argparse
import unittest

import numpy
from utils import parse_args_all_schemes, simple_exit

from pycbc.types import Array, complex64, float32, zeros

_scheme, _context = parse_args_all_schemes("correlate")

from pycbc.catalog import Merger
from pycbc.filter import matched_filter_core, resample_to_delta_t, sigma
from pycbc.noise import noise_from_psd
from pycbc.psd import aLIGOZeroDetHighPower, interpolate, inverse_spectrum_truncation
from pycbc.vetoes import (
    SingleDetPowerChisq,
    chisq_accum_bin,
    power_chisq,
    power_chisq_at_points_from_precomputed,
)
from pycbc.vetoes.chisq import (
    power_chisq_bin_weights,
    power_chisq_bins_and_weights,
)
from pycbc.vetoes.chisq_cpu import chisq_accum_bin_numpy
from pycbc.waveform import get_fd_waveform

trusted_accum = chisq_accum_bin_numpy


class TestChisq(unittest.TestCase):
    def setUp(self, *args):
        self.context = _context
        self.scheme = _scheme
        self.tolerance = 1e-6
        xr = numpy.random.uniform(low=-1, high=1.0, size=2**20)
        xi = numpy.random.uniform(low=-1, high=1.0, size=2**20)
        self.x = Array(xr + xi * 1.0j, dtype=complex64)
        self.z = zeros(2**20, dtype=float32)
        for i in range(4):
            trusted_accum(self.z, self.x)

        m = Merger("GW170814")

        ifos = ["H1", "L1", "V1"]
        data = {}
        psd = {}
        for ifo in ifos:
            # Read in and condition the data and measure PSD
            ts = m.strain(ifo).highpass_fir(15, 512)
            data[ifo] = resample_to_delta_t(ts, 1.0 / 2048).crop(2, 2)
            p = data[ifo].psd(2)
            p = interpolate(p, data[ifo].delta_f)
            p = inverse_spectrum_truncation(
                p, int(2 * data[ifo].sample_rate), low_frequency_cutoff=15.0
            )
            psd[ifo] = p
        hp, _ = get_fd_waveform(
            approximant="IMRPhenomD",
            mass1=31.36,
            mass2=31.36,
            f_lower=20.0,
            delta_f=data[ifo].delta_f,
        )
        hp.resize(len(psd[ifo]))

        # For each ifo use this template to calculate the SNR time series
        snr = {}
        snr_unnorm = {}
        norm = {}
        corr = {}
        for ifo in ifos:
            snr_unnorm[ifo], corr[ifo], norm[ifo] = matched_filter_core(
                hp, data[ifo], psd=psd[ifo], low_frequency_cutoff=20
            )
            snr[ifo] = snr_unnorm[ifo] * norm[ifo]

        self.snr = snr
        self.snr_unnorm = snr_unnorm
        self.norm = norm
        self.corr = corr
        self.hp = hp
        self.data = data
        self.psd = psd
        self.ifos = ifos

    def test_accum(self):
        with self.context:
            z = zeros(2**20, dtype=float32)
            for i in range(4):
                chisq_accum_bin(z, self.x)
            self.assertTrue(self.z.almost_equal_elem(z, self.tolerance))

    def test_chisq(self):
        chisq_quick = {}
        chisq_full = {}
        chisq_ref = {}
        for ifo in self.ifos:
            nbins = 26
            dof = nbins * 2 - 2

            bins, weights = power_chisq_bins_and_weights(
                self.hp, nbins, self.psd[ifo], low_frequency_cutoff=20.0
            )
            dof = (len(bins) - 1) * 2 - 2
            chisq = power_chisq_at_points_from_precomputed(
                self.corr[ifo],
                self.snr_unnorm[ifo][27402:27492].data,
                self.norm[ifo],
                bins,
                indices=numpy.arange(27402, 27492),
                weights=weights,
            )
            chisq_quick[ifo] = chisq / dof

            chisq = power_chisq(
                self.hp, self.data[ifo], nbins, self.psd[ifo], low_frequency_cutoff=20.0
            )
            chisq_full[ifo] = chisq[27402:27492] / dof
            max_diff = max(abs(chisq_full[ifo] - chisq_quick[ifo]))
            self.assertTrue(max_diff < 1e-5)

    def test_power_chisq_return_bins_slice(self):
        ifo = self.ifos[0]
        nbins = 16
        sl = slice(27400, 27500)
        # 1. Full return_bins (user-side slicing)
        chisq_full, bins_full = power_chisq(
            self.hp,
            self.data[ifo],
            nbins,
            self.psd[ifo],
            low_frequency_cutoff=20.0,
            return_bins=True,
        )
        # 2. Optimized internal return_bins_slice
        chisq_sl, bins_sl = power_chisq(
            self.hp,
            self.data[ifo],
            nbins,
            self.psd[ifo],
            low_frequency_cutoff=20.0,
            return_bins_slice=sl,
        )

        self.assertEqual(len(bins_full), nbins)
        self.assertEqual(len(bins_sl), nbins)
        # Verify exact numerical equivalence across all bins
        for i in range(nbins):
            expected = bins_full[i][sl].numpy()
            actual = bins_sl[i]
            numpy.testing.assert_allclose(actual, expected, rtol=1e-5, atol=1e-7)


class TestPowerChisqOptions(unittest.TestCase):
    def parse(self, args):
        parser = argparse.ArgumentParser()
        SingleDetPowerChisq.insert_option_group(parser)
        return SingleDetPowerChisq.from_cli(parser.parse_args(args))

    def test_options_reach_the_estimator(self):
        est = self.parse(["--chisq-bins", "16", "--chisq-snr-threshold", "4.5"])
        self.assertTrue(est.do)
        self.assertEqual(est.num_bins, "16")
        self.assertEqual(est.snr_threshold, 4.5)

    def test_off_by_default(self):
        est = self.parse([])
        self.assertFalse(est.do)
        self.assertIsNone(est.snr_threshold)


class Params:
    # A bank template's params: SingleDetPowerChisq keys its cache of bins
    # on this object
    pass


class TestBinWeights(unittest.TestCase):
    def setUp(self):
        self.sr, duration = 2048, 16
        self.n = self.sr * duration
        self.df = 1.0 / duration
        self.psd = aLIGOZeroDetHighPower(self.n // 2 + 1, self.df, 15)
        hp, _ = get_fd_waveform(
            approximant="IMRPhenomD", mass1=30, mass2=30, f_lower=20, delta_f=self.df
        )
        hp.resize(self.n // 2 + 1)
        hp.params = Params()
        hp.f_lower = 20.0
        self.template = hp

    def trigger(self, num_bins, amp, seed):
        """
        Inject a signal that exactly matches the template, at SNR amp,
        into Gaussian noise, and compute the chisq at its time. Returns the
        chisq, its dof, the chisq with every bin weighted equally (the
        standard formula), and the bins, correlation, SNR normalization and
        index used.
        """
        data = noise_from_psd(self.n, 1.0 / self.sr, self.psd, seed=seed)
        data = data.to_frequencyseries()
        norm = amp / sigma(self.template, psd=self.psd, low_frequency_cutoff=20)
        data = data + (self.template * norm).cyclic_time_shift(8.0)
        snr, corr, snr_norm = matched_filter_core(
            self.template, data, psd=self.psd, low_frequency_cutoff=20
        )
        index = numpy.array([8 * self.sr])
        # The number of bins is given as text, as from --chisq-bins
        est = SingleDetPowerChisq(str(num_bins))
        chi, dof = est.values(
            corr, snr.numpy()[index], snr_norm, self.psd, index, self.template
        )
        bins = est.cached_chisq_bins(self.template, self.psd)
        equal = power_chisq_at_points_from_precomputed(
            corr, snr.numpy()[index], snr_norm, bins, index
        )
        return chi[0], dof[0], equal[0], bins, corr, snr_norm, index[0]

    def test_unequal_bin_formula(self):
        # Compute the weighted chisq directly, sum_j |z_j|^2 / f_j - |rho|^2,
        # where f_j is the fraction of the template's power in bin j and z_j
        # the SNR accumulated over bin j, and compare. With 256 bins over a
        # 16 s template each bin is only a few frequency samples wide, so the
        # f_j differ noticeably from 1 / 256.
        chi, dof, equal, bins, corr, snr_norm, i = self.trigger(256, 30, 0)
        # The template's power, |h|^2 / S, over the frequencies the bins cover
        power = numpy.zeros(len(self.psd))
        band = slice(bins[0], bins[-1])
        power[band] = abs(self.template.numpy()[band]) ** 2 / self.psd.numpy()[band]
        frac = numpy.array([power[a:b].sum() for a, b in zip(bins[:-1], bins[1:])])
        frac /= frac.sum()
        shift = numpy.exp(2j * numpy.pi * numpy.arange(len(corr)) * i / self.n)
        z = numpy.array(
            [
                numpy.sum((corr.numpy() * shift)[a:b])
                for a, b in zip(bins[:-1], bins[1:])
            ]
        )
        expect = (numpy.sum(abs(z) ** 2 / frac) - abs(z.sum()) ** 2) * snr_norm**2
        # The chisq is accumulated in single precision, so agreement to 1e-4
        # is expected
        self.assertLess(abs(chi / expect - 1), 1e-4)
        # The standard (equal-weight) chisq differs from it by 4e-3 here, so
        # this tolerance does distinguish the two
        self.assertGreater(abs(equal / expect - 1), 1e-3)
        self.assertEqual(dof, 2 * (len(bins) - 1) - 2)

    def test_loud_signal_is_central(self):
        # For a signal that exactly matches the template, the weighted chisq
        # divided by its dof (510) has mean 1 and standard deviation 0.063,
        # so the mean over 3 noise realizations lies within 0.15 (4 standard
        # deviations) of 1. The standard chisq instead grows with the square
        # of the SNR when the bins are unequal: here it averages 1.37 x dof
        # at SNR 200, which the second check confirms.
        vals = numpy.array([self.trigger(256, 200, s)[:3] for s in range(3)])
        chi, dof, equal = vals.T
        self.assertLess(abs(numpy.mean(chi / dof) - 1), 0.15)
        self.assertGreater(numpy.mean(equal / dof), 1.15)

    def test_empty_bins_are_merged(self):
        # Template power of 1 in each of samples 10-11 and 20-29 and zero
        # elsewhere, so that the bins 12-15, 15-20 and the last bin, 30-35,
        # contain no power
        series = numpy.cumsum(
            numpy.r_[
                numpy.zeros(10),
                numpy.ones(2),
                numpy.zeros(8),
                numpy.ones(10),
                numpy.zeros(5),
            ]
        )
        bins = numpy.array([10, 11, 12, 15, 20, 25, 30, 35])
        merged, weights = power_chisq_bin_weights(series, bins)
        # Each empty bin is merged into the bin that follows it; the empty
        # last bin is merged into the one before it
        numpy.testing.assert_array_equal(merged, [10, 11, 12, 25, 35])
        # The remaining bins contain 1, 1, 5 and 5 of the 12 units of power,
        # and each weight is the total power over the bin's power
        numpy.testing.assert_allclose(weights, 12.0 / numpy.array([1, 1, 5, 5]))


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestChisq))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestPowerChisqOptions))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestBinWeights))

if __name__ == "__main__":
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
