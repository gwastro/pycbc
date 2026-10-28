# Copyright (C) 2026 Alex Nitz
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
"""Unit tests for hierarchical ratio/FIR matched filtering and follow-up."""

import os
import tempfile
import unittest

import h5py
import numpy as np

import pycbc
import pycbc.psd
from pycbc.filter import sigma
from pycbc.filter.matched_filter_ratio import MatchedFilterRatioControl
from pycbc.types import FrequencySeries, TimeSeries
from pycbc.vetoes.autochisq import SingleDetAutoChisq
from pycbc.vetoes.chain import SignalConsistencyChain
from pycbc.vetoes.chisq import SingleDetPowerChisq
from pycbc.vetoes.sgchisq import SingleDetSGChisq
from pycbc.waveform import compress, get_fd_waveform
from pycbc.waveform.bank import RatioFilterBank

try:
    import matchedfilter  # noqa: F401

    HAS_MATCHEDFILTER = True
except ImportError:
    HAS_MATCHEDFILTER = False


def _make_synthetic_fir_bank(path, num_levels=2, sample_rate=2048.0, seg_len=64.0):
    """Write a minimal synthetic N-level RatioFilterBank HDF5 file and return
    ``(htilde, delta_f, flen, f_low)``.
    """
    n_samples = round(sample_rate * seg_len)
    flen = n_samples // 2 + 1
    delta_f = 1.0 / seg_len
    f_low = 30.0

    hp, _ = get_fd_waveform(
        approximant="TaylorF2",
        mass1=1.4,
        mass2=1.4,
        delta_f=delta_f,
        f_lower=f_low,
        f_final=sample_rate / 2.0,
    )
    htilde = FrequencySeries(np.zeros(flen, dtype=np.complex64), delta_f=delta_f)
    m = min(len(hp), flen)
    htilde[:m] = (hp[:m] * pycbc.DYN_RANGE_FAC).astype(np.complex64)

    with h5py.File(path, "w") as f:
        f.create_dataset("mass1", data=np.array([1.4], dtype=np.float64))
        f.create_dataset("mass2", data=np.array([1.4], dtype=np.float64))
        f.create_dataset("spin1z", data=np.array([0.0], dtype=np.float64))
        f.create_dataset("spin2z", data=np.array([0.0], dtype=np.float64))
        f.create_dataset("f_lower", data=np.array([f_low], dtype=np.float64))
        f.create_dataset(
            "f_final", data=np.array([sample_rate / 2.0], dtype=np.float64)
        )
        f.create_dataset("template_hash", data=np.array([1001], dtype=np.int64))
        f.create_dataset("approximant", data=np.array([b"TaylorF2"], dtype="S16"))

        fir_grp = f.create_group("fir_data")
        fir_grp.attrs["num_levels"] = int(num_levels)
        fir_grp.attrs["sample_rate"] = float(sample_rate)

        # Level 1: leaf impulse filter (scaled by 2 * sample_rate)
        lvl1 = fir_grp.create_group("level_1")
        g1 = lvl1.create_group("0")
        taps1 = np.zeros((1, 33), dtype=np.complex64)
        taps1[0, 16] = 2.0 * sample_rate + 0.0j
        g1.create_dataset("taps", data=taps1)
        g1.create_dataset("actual_tap_count", data=np.array([33], dtype=np.int32))
        g1.create_dataset("fine_bank_index", data=np.array([0], dtype=np.int64))
        g1.create_dataset("sample_shift", data=np.array([2], dtype=np.int64))
        g1.create_dataset("sigmas", data=np.array([[1.0, 1.5]], dtype=np.float64))

        # Higher levels: unit impulse filters
        for lvl in range(2, num_levels + 1):
            g = fir_grp.create_group(f"level_{lvl}").create_group("0")
            taps_u = np.zeros((1, 17), dtype=np.float32)
            taps_u[0, 8] = 1.0
            g.create_dataset("taps", data=taps_u)
            g.create_dataset("actual_tap_count", data=np.array([17], dtype=np.int32))
            g.create_dataset("fine_bank_index", data=np.array([0], dtype=np.int64))
            g.create_dataset("sample_shift", data=np.array([3], dtype=np.int64))
            g.create_dataset("sigmas", data=np.array([[1.0, 2.0]], dtype=np.float64))

        # Compress and store the top-level reference waveform
        top_lvl = fir_grp[f"level_{num_levels}"]
        kmax = int(np.nonzero(abs(htilde.numpy()))[0][-1])
        hc = compress.compress_waveform(
            htilde,
            np.array([f_low, kmax * delta_f], float),
            1e-3,
            "inline_quadratic",
            "single",
            psd=None,
        )
        hc.write_to_hdf(top_lvl, "0")

    return htilde, delta_f, flen, f_low


class TestMatchedFilterRatio(unittest.TestCase):
    @unittest.skipUnless(HAS_MATCHEDFILTER, "matchedfilter is not installed")
    def test_hierarchical_search_and_followup(self):
        sample_rate = 2048.0
        seg_len = 64.0
        n_samples = round(sample_rate * seg_len)
        t_inj = 32.0

        with tempfile.TemporaryDirectory() as tmpdir:
            bank_path = os.path.join(tmpdir, "bank.hdf")
            htilde, delta_f, flen, f_low = _make_synthetic_fir_bank(
                bank_path, num_levels=2, sample_rate=sample_rate, seg_len=seg_len
            )
            bank = RatioFilterBank(
                bank_path,
                filter_length=flen,
                delta_f=delta_f,
                dtype=np.complex64,
                low_frequency_cutoff=f_low,
            )

            # Verify multi-level shift and sigma accumulation on RatioFilterBank
            np.testing.assert_allclose(bank.sigma_rescale([0]), [3.0], rtol=1e-6)
            _, _, _, shifts = bank.get_firs(0, level=1, return_shifts=True)
            np.testing.assert_array_equal(shifts, [5])

            # Build colored noise + injected signal (SNR 15 in H1, SNR 4.5 in L1)
            psd = (
                pycbc.psd.aLIGOZeroDetHighPower(flen, delta_f, f_low)
                * (pycbc.DYN_RANGE_FAC**2)
            ).astype(np.float32)
            phase = np.exp(-2j * np.pi * np.arange(flen) * delta_f * t_inj).astype(
                np.complex64
            )
            sig_unit = htilde / float(
                sigma(htilde, psd=psd, low_frequency_cutoff=f_low)
            )

            rng = np.random.default_rng(42)
            segments = {}
            for ifo, snr_inj in (("H1", 15.0), ("L1", 4.5)):
                noise_fd = TimeSeries(
                    rng.standard_normal(n_samples).astype(np.float32),
                    delta_t=1.0 / sample_rate,
                ).to_frequencyseries()
                w_fd = noise_fd.numpy() / np.sqrt(n_samples * sample_rate / 2.0)
                stilde = FrequencySeries(
                    (
                        w_fd * np.sqrt(np.maximum(psd.numpy(), 0.0))
                        + snr_inj * sig_unit.numpy() * phase
                    ).astype(np.complex64),
                    delta_f=delta_f,
                    epoch=1000.0,
                )
                stilde[0] = 0.0
                stilde[len(stilde) - 1] = 0.0
                stilde.psd = psd
                stilde.analyze = slice(4096, n_samples - 4096)
                stilde.cumulative_index = 4096
                segments[ifo] = [stilde]

            # 1. Primary + follow-up search with power, auto, and sine-Gaussian chisq
            veto_chain = SignalConsistencyChain(
                SingleDetPowerChisq(num_bins=16),
                SingleDetAutoChisq(stride=2, num_points=10),
                SingleDetSGChisq(bank, "16", None, ["mtotal>1:20-15,20-30"]),
                steps=["power", "auto", "sg"],
                window=int(sample_rate),
            )
            control = MatchedFilterRatioControl(
                bank,
                segments,
                snr_threshold=6.0,
                high_frequency_cutoff=sample_rate * 0.45,
                cluster_window=1.0,
                followup_snr_threshold=3.0,
                followup_bin_width=0.1,
                followup_num_bins=10,
            )
            triggers, has_followup = control.run(veto_callback=veto_chain.values)

            # Verify primary trigger in H1
            h1 = triggers["H1"]
            peak = int(np.argmax(np.abs(h1["snr"])))
            self.assertGreater(abs(h1["snr"][peak]), 12.0)
            self.assertAlmostEqual(
                h1["time_index"][peak] / sample_rate, t_inj, delta=0.01
            )
            self.assertLess(h1["chisq"][peak] / h1["chisq_dof"][peak], 2.5)

            # Verify secondary follow-up trigger in L1 uses L1's own correlation
            self.assertGreater(has_followup["L1"], 0)
            l1 = triggers["L1"]
            self.assertTrue(np.any(l1["is_followup"]))
            l1_peak = int(np.argmax(np.abs(l1["snr"])))
            self.assertGreater(abs(l1["snr"][l1_peak]), 3.0)
            self.assertLess(l1["chisq"][l1_peak] / l1["chisq_dof"][l1_peak], 2.5)
            self.assertNotAlmostEqual(
                float(l1["chisq"][l1_peak]), float(h1["chisq"][peak]), places=3
            )

            # 2. Omitting follow-up options runs only the primary search
            primary_only = MatchedFilterRatioControl(
                bank, segments, snr_threshold=6.0, cluster_window=1.0
            )
            _, no_followup = primary_only.run()
            self.assertEqual(no_followup["L1"], 0)
            bank.close_filehandler()


if __name__ == "__main__":
    unittest.main()
