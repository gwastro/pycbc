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
import types
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
    power_chisq_from_precomputed,
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
        _, bins_full = power_chisq(
            self.hp,
            self.data[ifo],
            nbins,
            self.psd[ifo],
            low_frequency_cutoff=20.0,
            return_bins=True,
        )
        # 2. Optimized internal return_bins_slice
        _, bins_sl = power_chisq(
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


def bin_snrs(corr, bins, index):
    """
    Return z_j, the SNR accumulated over each chisq bin j at the time sample
    index, by summing the correlation over the bin's frequencies with the
    phase that shifts it to that time.
    """
    n = len(corr)
    shifted = corr.numpy() * numpy.exp(2j * numpy.pi * numpy.arange(n) * index / n)
    return numpy.array([shifted[a:b].sum() for a, b in zip(bins[:-1], bins[1:])])


def bin_power_fractions(template, psd, bins):
    """Return the fraction of the template's power, |h|^2 / S, in each bin."""
    band = slice(bins[0], bins[-1])
    power = numpy.zeros(len(psd))
    power[band] = abs(template.numpy()[band]) ** 2 / psd.numpy()[band]
    frac = numpy.array([power[a:b].sum() for a, b in zip(bins[:-1], bins[1:])])
    return frac / frac.sum()


class TestBinWeights(unittest.TestCase):
    """
    Each chisq bin j is weighted by 1 / f_j, where f_j is the fraction of
    the template's power it holds, so that

        chisq = (sum_j |z_j|^2 / f_j - |rho|^2) * snr_norm^2.

    With bins of equal power every 1 / f_j is the number of bins, which is
    the standard power chisq.
    """

    def setUp(self):
        self.sr, duration = 2048, 16
        self.n = self.sr * duration
        df = 1.0 / duration
        self.psd = aLIGOZeroDetHighPower(self.n // 2 + 1, df, 15)
        hp, _ = get_fd_waveform(
            approximant="IMRPhenomD", mass1=30, mass2=30, f_lower=20, delta_f=df
        )
        hp.resize(self.n // 2 + 1)
        self.template = hp
        self.index = 8 * self.sr

    def filter(self, snr, seed=0):
        """
        Return the SNR series, correlation and SNR normalization for Gaussian
        noise holding a signal that exactly matches the template, at the
        given SNR and at the time sample self.index.
        """
        data = noise_from_psd(self.n, 1.0 / self.sr, self.psd, seed=seed)
        data = data.to_frequencyseries()
        amp = snr / sigma(self.template, psd=self.psd, low_frequency_cutoff=20)
        data += (self.template * amp).cyclic_time_shift(self.index / self.sr)
        return matched_filter_core(
            self.template, data, psd=self.psd, low_frequency_cutoff=20
        )

    def test_empty_bins_are_merged(self):
        # Template power of 1 in each of samples 10-11 and 20-29, and zero
        # elsewhere, so bins 12-15, 15-20 and 30-35 hold no power
        power = numpy.zeros(35)
        power[10:12] = 1
        power[20:30] = 1
        bins = numpy.array([10, 11, 12, 15, 20, 25, 30, 35])
        merged, weights = power_chisq_bin_weights(numpy.cumsum(power), bins)
        # Empty bins join the bin after them, and the empty last bin joins
        # the one before it, leaving bins with 1, 1, 5 and 5 of the 12 units
        numpy.testing.assert_array_equal(merged, [10, 11, 12, 25, 35])
        numpy.testing.assert_allclose(weights, 12.0 / numpy.array([1, 1, 5, 5]))

    def test_weighted_chisq_formula(self):
        # Equal-power bins, with every other edge removed in the upper half so
        # that some bins hold twice the power of others
        bins, _ = power_chisq_bins_and_weights(
            self.template, 32, self.psd, low_frequency_cutoff=20
        )
        bins = numpy.concatenate([bins[:16], bins[16:-1:2], bins[-1:]])
        frac = bin_power_fractions(self.template, self.psd, bins)
        weights = 1.0 / frac

        snr, corr, snr_norm = self.filter(snr=30)
        z = bin_snrs(corr, bins, self.index)
        expected = (numpy.sum(abs(z) ** 2 / frac) - abs(z.sum()) ** 2) * snr_norm**2

        index = numpy.array([self.index])
        at_points = power_chisq_at_points_from_precomputed(
            corr, snr.numpy()[index], snr_norm, bins, index, weights
        )[0]
        by_fft = power_chisq_from_precomputed(
            corr, snr, snr_norm, bins, weights=weights
        )[self.index]
        # The chisq is accumulated in single precision
        self.assertAlmostEqual(at_points / expected, 1, delta=1e-4)
        self.assertAlmostEqual(by_fft / expected, 1, delta=1e-4)

        # Weighting every bin equally, as for bins of equal power, gives a
        # different answer for these bins
        equal = power_chisq_at_points_from_precomputed(
            corr, snr.numpy()[index], snr_norm, bins, index
        )[0]
        self.assertGreater(abs(equal / expected - 1), 0.1)

    def test_matched_signal_has_chisq_near_dof(self):
        # For a signal that exactly matches the template, chisq / dof has
        # mean 1 and, with 256 bins (dof 510), standard deviation 0.063, so
        # it lies within 0.25 (4 standard deviations) of 1. Weighting the
        # bins equally gives 1.7 to 2.0 here at SNR 300, which this rejects.
        # The template is given the attributes a bank template has, as
        # SingleDetPowerChisq reads them.
        self.template.params = types.SimpleNamespace()
        self.template.f_lower = 20.0
        snr, corr, snr_norm = self.filter(snr=300)
        index = numpy.array([self.index])
        # The number of bins is given as text, as from --chisq-bins
        chisq, dof = SingleDetPowerChisq("256").values(
            corr, snr.numpy()[index], snr_norm, self.psd, index, self.template
        )
        self.assertAlmostEqual(chisq[0] / dof[0], 1, delta=0.25)


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestChisq))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestPowerChisqOptions))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestBinWeights))

if __name__ == "__main__":
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
