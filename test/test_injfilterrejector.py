# Copyright (C) 2026 Alexander Harvey Nitz
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
Unit tests for all major features of InjFilterRejector:
1. Chirp time window filtering (keeping/rejecting based on tau0)
2. Coarse match threshold filtering (keeping/rejecting based on match)
3. Trigger window filtering (keeping/rejecting trigger intervals)
4. Optimal SNR threshold filtering (keeping/rejecting based on SNR)
5. Direct agreement of optimal_snr() with pycbc.filter.sigma()
"""

import unittest
from types import SimpleNamespace

import numpy as np
from utils import parse_args_cpu_only, simple_exit

from pycbc import DYN_RANGE_FAC
from pycbc.filter import sigma
from pycbc.inject.injfilterrejector import InjFilterRejector
from pycbc.psd import aLIGOZeroDetHighPower
from pycbc.types import FrequencySeries
from pycbc.waveform import get_td_waveform

DELTA_T = 1.0 / 4096
F_LOWER = 20.0


def bns_signal(tc=100.0, distance=100):
    hp, _ = get_td_waveform(approximant='TaylorT4', mass1=1.4, mass2=1.3,
                            delta_t=DELTA_T, f_lower=F_LOWER,
                            distance=distance)
    hp.start_time = tc
    return hp


class TestInjFilterRejector(unittest.TestCase):
    """Tests verifying that InjFilterRejector keeps and rejects as expected."""

    def test_chirp_time_rejection(self):
        """Verify chirp-time window keeps near templates and rejects."""
        rej = InjFilterRejector('dummy.hdf', chirp_time_window=2.0,
                                match_threshold=None, f_lower=F_LOWER,
                                seg_buffer=0)
        inj_table = [SimpleNamespace(mass1=1.4, mass2=1.3,
                                     geocent_end_time=100.0,
                                     geocent_end_time_ns=0,
                                     simulation_id=0)]
        rej.injection_params = SimpleNamespace(table=inj_table)
        rej.injection_ids = [0]

        class DummyBank:
            def __init__(self, m1, m2):
                self.table = [{'mass1': m1, 'mass2': m2}]

        seg = SimpleNamespace(start_time=90.0, end_time=110.0, psd=None)

        # Same masses: tau0 difference is 0 <= 2.0s -> KEPT
        keep = rej.template_segment_checker(DummyBank(1.4, 1.3), 0, seg)
        self.assertTrue(keep)

        # Distant masses (30, 30): tau0 diff >> 2.0s -> REJECTED
        reject = rej.template_segment_checker(DummyBank(30.0, 30.0), 0, seg)
        self.assertFalse(reject)

    def test_match_threshold_rejection(self):
        """Verify coarse match threshold keeps matching and rejects."""
        rej = InjFilterRejector('dummy.hdf', chirp_time_window=None,
                                match_threshold=0.8, f_lower=F_LOWER,
                                seg_buffer=0)
        sig = bns_signal(tc=100.0)
        rej.generate_short_inj_from_inj(sig.copy(), 0)
        inj_table = [SimpleNamespace(mass1=1.4, mass2=1.3,
                                     geocent_end_time=100.0,
                                     geocent_end_time_ns=0,
                                     simulation_id=0)]
        rej.injection_params = SimpleNamespace(table=inj_table)
        rej.injection_ids = [0]

        class DummyBankMatching:
            table = [SimpleNamespace(f_lower=F_LOWER)]

            def generate_with_delta_f_and_max_freq(self, t_num, fmax, deltaf,
                                                   low_frequency_cutoff,
                                                   cached_mem=None):
                return rej.short_injections[0].copy()

        class DummyBankMismatch:
            table = [SimpleNamespace(f_lower=F_LOWER)]

            def generate_with_delta_f_and_max_freq(self, t_num, fmax, deltaf,
                                                   low_frequency_cutoff,
                                                   cached_mem=None):
                data = np.zeros(len(rej.short_injections[0]),
                                dtype=np.complex64)
                data[len(data) // 2] = 1.0
                return FrequencySeries(data, delta_f=deltaf)

        tlen = 256 * int(1 / DELTA_T)
        delta_f = 1.0 / 256
        psd = aLIGOZeroDetHighPower(tlen // 2 + 1, delta_f,
                                    F_LOWER) * DYN_RANGE_FAC ** 2
        seg = SimpleNamespace(start_time=90.0, end_time=110.0, psd=psd)

        # Matching waveform: match ~ 1.0 > 0.8 -> KEPT
        keep = rej.template_segment_checker(DummyBankMatching(), 0, seg)
        self.assertTrue(keep)

        # Mismatching waveform: match < 0.8 -> REJECTED
        rej._short_template_id = None
        reject = rej.template_segment_checker(DummyBankMismatch(), 0, seg)
        self.assertFalse(reject)

    def test_trigger_window_rejection(self):
        """Verify trigger window keeps triggers near injections and rejects."""
        rej = InjFilterRejector('dummy.hdf', chirp_time_window=None,
                                match_threshold=None, f_lower=F_LOWER,
                                inj_trigger_window=1.0)
        rej.injection_params = SimpleNamespace(
            end_times=lambda: [100.0, 200.0]
        )
        trigs = np.array([50.0, 99.5, 100.2, 150.0, 200.5, 250.0])
        mask = rej.find_indices_in_injection_intervals(trigs)

        # Triggers at 99.5, 100.2, 200.5 are within 1.0s of 100.0 or 200.0
        expected_mask = np.array([False, True, True, False, True, False])
        self.assertTrue(np.array_equal(mask, expected_mask))
        self.assertEqual(list(trigs[mask]), [99.5, 100.2, 200.5])

    def test_optimal_snr_threshold_rejection(self):
        """Verify optimal SNR threshold keeps high SNR and rejects low SNR."""
        rej = InjFilterRejector('dummy.hdf', chirp_time_window=None,
                                match_threshold=None, f_lower=F_LOWER,
                                optimal_snr_threshold=10.0, seg_buffer=0)
        sig = bns_signal(tc=100.0, distance=100)  # SNR ~ 32
        rej.generate_short_inj_from_inj(sig.copy(), 0)
        inj_table = [SimpleNamespace(mass1=1.4, mass2=1.3,
                                     geocent_end_time=100.0,
                                     geocent_end_time_ns=0,
                                     simulation_id=0)]
        rej.injection_params = SimpleNamespace(table=inj_table)
        rej.injection_ids = [0]

        class DummyBank:
            table = [{'mass1': 1.4, 'mass2': 1.3}]

        tlen = 256 * int(1 / DELTA_T)
        delta_f = 1.0 / 256
        psd = aLIGOZeroDetHighPower(tlen // 2 + 1, delta_f,
                                    F_LOWER) * DYN_RANGE_FAC ** 2
        seg = SimpleNamespace(start_time=90.0, end_time=110.0, psd=psd)

        # Threshold 10.0: SNR ~ 32 >= 10.0 -> KEPT
        keep = rej.template_segment_checker(DummyBank(), 0, seg)
        self.assertTrue(keep)

        # Threshold 50.0: SNR ~ 32 < 50.0 -> REJECTED
        rej.optimal_snr_threshold = 50.0
        reject = rej.template_segment_checker(DummyBank(), 0, seg)
        self.assertFalse(reject)

    def test_optimal_snr_calculation(self):
        """Verify optimal_snr() calculates coarse SNR directly."""
        sig = bns_signal(tc=0.0)
        rej = InjFilterRejector('dummy.hdf', chirp_time_window=None,
                                match_threshold=None, f_lower=F_LOWER,
                                optimal_snr_threshold=4.0)
        rej.generate_short_inj_from_inj(sig.copy(), 0)
        tlen = 256 * int(1 / DELTA_T)
        delta_f = 1.0 / 256
        psd = aLIGOZeroDetHighPower(tlen // 2 + 1, delta_f, F_LOWER)
        padded = sig.copy()
        padded.resize(tlen)
        expected = sigma(padded, psd=psd, low_frequency_cutoff=F_LOWER)
        got = rej.optimal_snr(0, psd * DYN_RANGE_FAC ** 2)

        # Coarse SNR agrees with full sigma within ~5% band truncation
        self.assertAlmostEqual(got / expected, 1.0, delta=0.05)

    def test_disabled_rejector(self):
        """Verify disabled rejector allows all templates and triggers."""
        rej = InjFilterRejector(None, chirp_time_window=None,
                                match_threshold=None, f_lower=F_LOWER)
        self.assertFalse(rej.enabled)
        # Allows all templates
        self.assertTrue(rej.template_segment_checker(None, 0, None))
        # Allows all triggers
        trigs = np.array([10.0, 20.0, 30.0])
        self.assertEqual(rej.find_indices_in_injection_intervals(trigs),
                         slice(None))
        # Returns None for optimal_snr
        self.assertIsNone(rej.optimal_snr(0, None))


suite = unittest.TestSuite()
loader = unittest.TestLoader()
suite.addTest(loader.loadTestsFromTestCase(TestInjFilterRejector))

if __name__ == '__main__':
    parse_args_cpu_only("InjFilterRejector")
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
