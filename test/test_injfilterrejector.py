# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.

"""
Unit test for the optimal-SNR cut of InjFilterRejector.
"""

import unittest
from pycbc import DYN_RANGE_FAC
from pycbc.filter import sigma
from pycbc.inject.injfilterrejector import InjFilterRejector
from pycbc.psd import aLIGOZeroDetHighPower
from pycbc.waveform import get_td_waveform
from utils import parse_args_cpu_only, simple_exit

parse_args_cpu_only("InjFilterRejector")

DELTA_T = 1.0 / 4096
F_LOWER = 20.0


def bns_signal():
    hp, _ = get_td_waveform(approximant='TaylorT4', mass1=1.4, mass2=1.3,
                            delta_t=DELTA_T, f_lower=F_LOWER, distance=100)
    hp.start_time = 0
    return hp


class TestOptimalSNR(unittest.TestCase):
    def rejector(self, **kwargs):
        args = dict(chirp_time_window=None, match_threshold=None,
                    f_lower=F_LOWER, optimal_snr_threshold=4.0)
        args.update(kwargs)
        return InjFilterRejector('injections.hdf', **args)

    def test_matches_sigma(self):
        """The coarse in-job optimal SNR agrees with a direct sigma()."""
        signal = bns_signal()
        rej = self.rejector()
        rej.generate_short_inj_from_inj(signal.copy(), 0)
        tlen = 256 * int(1 / DELTA_T)
        delta_f = 1.0 / 256
        psd = aLIGOZeroDetHighPower(tlen // 2 + 1, delta_f, F_LOWER)
        padded = signal.copy()
        padded.resize(tlen)
        expected = sigma(padded, psd=psd, low_frequency_cutoff=F_LOWER)
        got = rej.optimal_snr(0, psd * DYN_RANGE_FAC ** 2)
        self.assertAlmostEqual(got / expected, 1.0, delta=0.01)

    def test_signal_untouched_without_match_test(self):
        """Only the optimal-SNR cut: the waveform handed in is not padded,
        since the caller adds it to the data afterwards."""
        signal = bns_signal()
        n = len(signal)
        self.rejector().generate_short_inj_from_inj(signal, 0)
        self.assertEqual(len(signal), n)

    def test_match_test_still_pads_in_place(self):
        """With the match test on, the existing in-place behaviour is kept."""
        signal = bns_signal()
        n = len(signal)
        rej = self.rejector(match_threshold=0.9)
        rej.generate_short_inj_from_inj(signal, 0)
        self.assertGreater(len(signal), n)
        self.assertIn(0, rej.short_injections)
        self.assertIn(0, rej.injection_power)

    def test_disabled_stores_nothing(self):
        rej = self.rejector(optimal_snr_threshold=None)
        self.assertFalse(rej.enabled)
        rej.generate_short_inj_from_inj(bns_signal(), 0)
        self.assertIsNone(rej.optimal_snr(0, None))


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestOptimalSNR))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
