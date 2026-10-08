import unittest
import numpy as np
from pycbc.types import TimeSeries, FrequencySeries
from pycbc.waveform import get_fd_waveform
from pycbc.vetoes import power_chisq

class TestChisqSlicing(unittest.TestCase):
    def test_power_chisq_return_bins_slice(self):
        delta_t = 1.0 / 2048
        N = 2048 * 8
        data = TimeSeries(np.random.normal(0, 1, N).astype(np.float64), delta_t=delta_t)
        psd = FrequencySeries(np.ones(N // 2 + 1, dtype=np.float64), delta_f=data.delta_f)
        hp, _ = get_fd_waveform(approximant="IMRPhenomD", mass1=30, mass2=30,
                                 f_lower=20.0, delta_f=data.delta_f)
        hp.resize(len(psd))

        nbins = 16
        sl = slice(100, 200)

        # Full return_bins
        chisq_full, bins_full = power_chisq(
            hp, data, nbins, psd, low_frequency_cutoff=20.0, return_bins=True)
        # Optimized internal return_bins_slice
        chisq_sl, bins_sl = power_chisq(
            hp, data, nbins, psd, low_frequency_cutoff=20.0, return_bins_slice=sl)

        self.assertEqual(len(bins_full), nbins)
        self.assertEqual(len(bins_sl), nbins)
        for i in range(nbins):
            expected = bins_full[i][sl].numpy()
            actual = bins_sl[i]
            np.testing.assert_allclose(actual, expected, rtol=1e-5, atol=1e-7)

if __name__ == '__main__':
    unittest.main()
