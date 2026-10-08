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

    def test_compression_representative_cases(self):
        # Case 1: BNS system (TaylorF2, 1.4 + 1.4 Msun)
        hp_bns, _ = get_fd_waveform(approximant='TaylorF2', mass1=1.4, mass2=1.4,
                                    f_lower=30.0, delta_f=1.0/8)
        cwave_bns = compress.compress_waveform(
            hp_bns, np.array([30.0, 300.0], dtype=float), tolerance=1e-3,
            interpolation='inline_linear', precision='single')
        decomp_bns = cwave_bns.decompress(df=hp_bns.delta_f)
        decomp_bns.resize(len(hp_bns))
        m_bns, _ = match(hp_bns, decomp_bns, low_frequency_cutoff=30.0)
        self.assertGreater(m_bns, 0.97)

        # Case 2: BBH system (IMRPhenomD, 30.0 + 30.0 Msun)
        hp_bbh, _ = get_fd_waveform(approximant='IMRPhenomD', mass1=30.0, mass2=30.0,
                                    f_lower=20.0, delta_f=1.0/8)
        cwave_bbh = compress.compress_waveform(
            hp_bbh, np.array([20.0, 300.0], dtype=float), tolerance=1e-3,
            interpolation='inline_linear', precision='single')
        decomp_bbh = cwave_bbh.decompress(df=hp_bbh.delta_f)
        decomp_bbh.resize(len(hp_bbh))
        m_bbh, _ = match(hp_bbh, decomp_bbh, low_frequency_cutoff=20.0)
        self.assertGreater(m_bbh, 0.98)

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
            self.assertGreater(m, 0.97, f"Failed for interpolation={interp}")

    def test_sample_point_deduplication_allows_adjacent(self):
        hp, _ = get_fd_waveform(approximant='TaylorF2', mass1=1.4, mass2=1.4,
                                f_lower=30.0, delta_f=1.0/16)
        fmin = 30.0
        sample_points = np.array([fmin, 200.0, 400.0], dtype=float)
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
        self.assertEqual(diff, 0.0)

        # Test compression with adjacent sample points in the initial set
        sample_points = np.array([30.0, 30.0 + df, 30.0 + 2*df, 100.0], dtype=float)
        cwave = compress.compress_waveform(
            hp, sample_points, tolerance=1e-3,
            interpolation='inline_linear', precision='single')
        diffs = np.diff(cwave.sample_points)
        self.assertTrue(np.all(diffs > 0))


if __name__ == '__main__':
    unittest.main()
