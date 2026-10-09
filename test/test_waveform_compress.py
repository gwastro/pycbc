import unittest
import numpy as np
from pycbc.waveform import get_fd_waveform, compress
from pycbc.filter import match


class TestWaveformCompress(unittest.TestCase):
    def test_compress_waveform_accuracy(self):
        # Test representative cases (BNS TaylorF2 and BBH IMRPhenomD)
        # Verifies sample point counts and matches against expected physical bounds
        cases = [
            # (approximant, mass1, mass2, f_lower, delta_f, f_max, min_pts, max_pts, min_match)
            ('TaylorF2', 1.4, 1.4, 30.0, 1.0/16, 400.0, 180, 260, 0.985),
            ('TaylorF2', 1.4, 1.4, 30.0, 1.0/8, 300.0, 200, 280, 0.975),
            ('IMRPhenomD', 30.0, 30.0, 20.0, 1.0/8, 300.0, 18, 35, 0.995),
        ]
        for approx, m1, m2, flow, df, fmax, min_pts, max_pts, min_match in cases:
            hp, _ = get_fd_waveform(approximant=approx, mass1=m1, mass2=m2,
                                    f_lower=flow, delta_f=df)
            sample_points = np.array([flow, fmax], dtype=float)
            cwave = compress.compress_waveform(
                hp, sample_points, tolerance=1e-3,
                interpolation='inline_linear', precision='single')
            self.assertGreaterEqual(
                len(cwave.sample_points), min_pts,
                f"Sample points below expected floor for {approx}")
            self.assertLessEqual(
                len(cwave.sample_points), max_pts,
                f"Sample points above expected ceiling for {approx}")

            decomp = cwave.decompress(df=hp.delta_f)
            decomp.resize(len(hp))
            m, _ = match(hp, decomp, low_frequency_cutoff=flow)
            self.assertGreater(
                m, min_match,
                f"Match {m} below threshold {min_match} for {approx}")

    def test_compression_interpolation_orders(self):
        hp, _ = get_fd_waveform(approximant='TaylorF2', mass1=1.4, mass2=1.4,
                                f_lower=30.0, delta_f=1.0/8)
        interp_schemes = ['inline_linear', 'inline_quadratic', 'inline_cubic',
                          'linear', 'cubic']
        for interp in interp_schemes:
            sample_points = np.linspace(30.0, 300.0, 60, dtype=float)
            cwave = compress.compress_waveform(
                hp, sample_points, tolerance=1e-3,
                interpolation=interp, precision='single')
            decomp = cwave.decompress(df=hp.delta_f)
            decomp.resize(len(hp))
            m, _ = match(hp, decomp, low_frequency_cutoff=30.0)
            self.assertGreater(m, 0.975, f"Failed for interpolation={interp}")

    def test_sample_point_deduplication_allows_adjacent(self):
        hp, _ = get_fd_waveform(approximant='TaylorF2', mass1=1.4, mass2=1.4,
                                f_lower=30.0, delta_f=1.0/16)
        # Test deduplication and ordering even if caller passes unsorted/duplicate points
        sample_points = np.array([400.0, 30.0, 200.0, 200.0], dtype=float)
        cwave = compress.compress_waveform(
            hp, sample_points, tolerance=5e-4,
            interpolation='inline_linear', precision='single')
        diffs = np.diff(cwave.sample_points)
        self.assertTrue(np.all(diffs > 0))

    def test_narrow_slice_and_adjacent_points_no_singularity(self):
        hp, _ = get_fd_waveform(approximant='TaylorF2', mass1=1.4, mass2=1.4,
                                f_lower=30.0, delta_f=1.0/16)
        df = hp.delta_f
        # Test _vecdiff directly on adjacent frequencies (<= 1 sample apart)
        diff = compress._vecdiff(hp, hp, 30.0, 30.0 + df)
        self.assertAlmostEqual(diff, 0.0, places=10)

        # Test compression with adjacent sample points in the initial set
        sample_points = np.array([30.0, 30.0 + df, 30.0 + 2*df, 100.0], dtype=float)
        cwave = compress.compress_waveform(
            hp, sample_points, tolerance=1e-3,
            interpolation='inline_linear', precision='single')
        diffs = np.diff(cwave.sample_points)
        self.assertTrue(np.all(diffs > 0))


if __name__ == '__main__':
    unittest.main()
