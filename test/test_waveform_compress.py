import unittest
import numpy as np
from pycbc.waveform import get_fd_waveform, compress
from pycbc.filter import match


class TestWaveformCompress(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Ensure delta_f is chosen so that the full physical signal length
        # fits within 1 / delta_f to avoid time-domain aliasing.
        # BNS duration from 30 Hz is ~59s, so delta_f = 1/64 (T = 64s > 59s).
        cls.bns_flow = 30.0
        cls.bns_df = 1.0 / 64
        cls.hp_bns, _ = get_fd_waveform(
            approximant='TaylorF2', mass1=1.4, mass2=1.4,
            f_lower=cls.bns_flow, delta_f=cls.bns_df)

        # BBH duration from 20 Hz is ~1.1s, so delta_f = 1/4 (T = 4s > 1.1s).
        cls.bbh_flow = 20.0
        cls.bbh_df = 1.0 / 4
        cls.hp_bbh, _ = get_fd_waveform(
            approximant='IMRPhenomD', mass1=30.0, mass2=30.0,
            f_lower=cls.bbh_flow, delta_f=cls.bbh_df)

    def test_compress_waveform_accuracy(self):
        cases = [
            (self.hp_bns, self.bns_flow, 400.0, 180, 260, 0.985),
            (self.hp_bbh, self.bbh_flow, 300.0, 15, 35, 0.995),
        ]
        for hp, flow, fmax, min_pts, max_pts, min_match in cases:
            sample_points = np.array([flow, fmax], dtype=float)
            cwave = compress.compress_waveform(
                hp, sample_points, tolerance=1e-3,
                interpolation='inline_linear', precision='single')
            self.assertGreaterEqual(
                len(cwave.sample_points), min_pts,
                f"Sample points below expected floor for flow={flow}")
            self.assertLessEqual(
                len(cwave.sample_points), max_pts,
                f"Sample points above expected ceiling for flow={flow}")

            decomp = cwave.decompress(df=hp.delta_f)
            decomp.resize(len(hp))
            m, _ = match(hp, decomp, low_frequency_cutoff=flow)
            self.assertGreater(
                m, min_match,
                f"Match {m} below threshold {min_match} for flow={flow}")

    def test_compression_interpolation_orders(self):
        interp_schemes = ['inline_linear', 'inline_quadratic', 'inline_cubic',
                          'linear', 'cubic']
        for interp in interp_schemes:
            sample_points = np.linspace(self.bns_flow, 300.0, 60, dtype=float)
            cwave = compress.compress_waveform(
                self.hp_bns, sample_points, tolerance=1e-3,
                interpolation=interp, precision='single')
            decomp = cwave.decompress(df=self.hp_bns.delta_f)
            decomp.resize(len(self.hp_bns))
            m, _ = match(self.hp_bns, decomp, low_frequency_cutoff=self.bns_flow)
            self.assertGreater(m, 0.975, f"Failed for interpolation={interp}")

    def test_sample_point_deduplication_allows_adjacent(self):
        # Verify deduplication and sorting even if caller passes unsorted/duplicate points
        sample_points = np.array([400.0, self.bns_flow, 200.0, 200.0], dtype=float)
        cwave = compress.compress_waveform(
            self.hp_bns, sample_points, tolerance=5e-4,
            interpolation='inline_linear', precision='single')
        diffs = np.diff(cwave.sample_points)
        self.assertTrue(np.all(diffs > 0))

    def test_narrow_slice_and_adjacent_points_no_singularity(self):
        df = self.hp_bns.delta_f
        # Test _vecdiff directly on adjacent frequencies (<= 1 sample apart)
        diff = compress._vecdiff(self.hp_bns, self.hp_bns, self.bns_flow, self.bns_flow + df)
        self.assertAlmostEqual(diff, 0.0, places=10)

        # Test compression with adjacent sample points in the initial set
        sample_points = np.array([self.bns_flow, self.bns_flow + df, self.bns_flow + 2*df, 100.0], dtype=float)
        cwave = compress.compress_waveform(
            self.hp_bns, sample_points, tolerance=1e-3,
            interpolation='inline_linear', precision='single')
        diffs = np.diff(cwave.sample_points)
        self.assertTrue(np.all(diffs > 0))


if __name__ == '__main__':
    unittest.main()
