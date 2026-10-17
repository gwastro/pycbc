"""Known-answer checks for the CPU pointwise chi-square working precision."""
import unittest

import numpy as np

from pycbc.scheme import CPUScheme
from pycbc.types import Array
from pycbc.vetoes.chisq import power_chisq_at_points_from_precomputed
from pycbc.vetoes.chisq_cpu import shift_sum, point_chisq_code


class TestPointChisqArithmetic(unittest.TestCase):
    def test_long_rotation_preserves_single_coefficient_power(self):
        # One nonzero coefficient has unit power at every time, regardless
        # of the phase constant. Zeros force the actual recurrence to run.
        corr = np.zeros(2**20, dtype=np.complex64)
        corr[2**18 - 1] = 1
        bins = [0, 2**18, 2**19]
        for dtype in (np.complex64, np.complex128):
            with self.subTest(dtype=dtype), CPUScheme():
                value = Array(corr.astype(dtype), copy=False)
                power = shift_sum(value, [100000], bins)
                expected_dtype = np.empty((), dtype=dtype).real.dtype
                self.assertEqual(power.dtype, expected_dtype)
                np.testing.assert_allclose(power, 1, rtol=2e-7, atol=0)
                chisq = power_chisq_at_points_from_precomputed(
                    value, np.ones(1, dtype=dtype), 1., bins, [100000])
                np.testing.assert_allclose(chisq, 1, rtol=4e-7, atol=0)

    def test_small_terms_survive_within_a_bin(self):
        # At t=0 every phase is exactly one: this isolates accumulation.
        n = 2**16
        corr = np.full(n, 2.**-26, dtype=np.complex64)
        corr[0] = 1
        expected = (1 + (n - 1) * 2.**-26)**2
        with CPUScheme():
            actual = shift_sum(Array(corr, copy=False), [0], [0, n])
        np.testing.assert_allclose(actual, expected, rtol=1e-7, atol=0)

    def test_small_powers_survive_across_bins(self):
        # One coefficient per bin removes both recurrence and within-bin sums.
        n = 257
        corr = np.full(n, 2.**-13, dtype=np.complex64)
        corr[0] = 1
        expected = 1 + (n - 1) * 2.**-26
        with CPUScheme():
            actual = shift_sum(Array(corr, copy=False), [0], np.arange(n + 1))
        np.testing.assert_allclose(actual, expected, rtol=1e-7, atol=0)

    def test_direct_call_preserves_initial_output(self):
        output = np.array([7, 11], dtype=np.float64)
        corr = np.array([1, 2], dtype=np.complex128)
        shifts = np.array([0], dtype=np.float64)
        bins = np.array([0, 2], dtype=np.uint32)
        point_chisq_code(output, corr, 1, 2, shifts, bins, 1)
        np.testing.assert_array_equal(output, [16, 11])

    def test_empty_points(self):
        corr = Array(np.ones(2, dtype=np.complex64), copy=False)
        with CPUScheme():
            output = shift_sum(corr, [], [0, 2])
        self.assertEqual(output.dtype, np.float32)
        self.assertEqual(len(output), 0)

    def test_short_bins_agree_with_direct_sum(self):
        rng = np.random.default_rng(182)
        data = (rng.normal(size=128) + 1j*rng.normal(size=128)) / 10
        bins = [3, 19, 54, 64]
        points = np.array([0, 1.5, 31])
        for dtype in (np.complex64, np.complex128):
            corr = data.astype(dtype)
            expected = np.zeros(len(points))
            for lo, hi in zip(bins[:-1], bins[1:]):
                phase = np.exp(
                    2j*np.pi*np.outer(points, np.arange(lo, hi))/128)
                expected += abs(phase @ corr[lo:hi].astype(np.complex128))**2
            with self.subTest(dtype=dtype), CPUScheme():
                actual = shift_sum(Array(corr, copy=False), points, bins)
                # The independent truncated-pi defect remains on this branch.
                np.testing.assert_allclose(actual, expected,
                                           rtol=3e-7, atol=1e-8)


if __name__ == '__main__':
    unittest.main()
