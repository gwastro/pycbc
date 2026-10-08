import unittest
import numpy as np
from pycbc.waveform import get_fd_waveform, compress
from pycbc.filter import match


class TestWaveformCompress(unittest.TestCase):
    def test_compress_waveform_accuracy(self):
        hp, _ = get_fd_waveform(approximant='TaylorF2', mass1=1.4, mass2=1.4,
                                f_lower=30.0, delta_f=1.0/16)
        fmin = 30.0
        sample_points = np.array([fmin, 400.0], dtype=float)
        cwave = compress.compress_waveform(
            hp, sample_points, tolerance=1e-3,
            interpolation='inline_linear', precision='single')
        self.assertGreater(len(cwave.sample_points), 2)

        decomp = cwave.decompress(df=hp.delta_f)
        decomp.resize(len(hp))
        m, _ = match(hp, decomp, low_frequency_cutoff=30.0)
        self.assertGreater(m, 0.98)

    def test_sample_point_deduplication_allows_adjacent(self):
        hp, _ = get_fd_waveform(approximant='TaylorF2', mass1=1.4, mass2=1.4,
                                f_lower=30.0, delta_f=1.0/16)
        fmin = 30.0
        sample_points = np.array([fmin, 200.0, 400.0], dtype=float)
        cwave = compress.compress_waveform(
            hp, sample_points, tolerance=5e-4,
            interpolation='inline_linear', precision='single')
        diffs = np.diff(cwave.sample_points)
        # Verify all points are strictly monotonically increasing (no duplicate points)
        self.assertTrue(np.all(diffs > 0))


if __name__ == '__main__':
    unittest.main()
