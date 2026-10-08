import unittest
import tempfile
import os
import h5py
import numpy as np
from pycbc.events import eventmgr


class TestEventMgr(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.h5_path = os.path.join(self.temp_dir.name, "test_events.hdf")

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_h5filesyntsugar_append_and_context(self):
        # 1. Write initial dataset using group and default mode='a'
        with eventmgr.H5FileSyntSugar(self.h5_path, group='H1') as f:
            f['snr'] = np.array([5.5, 6.2, 7.8], dtype=np.float32)

        # 2. Append dataset for second detector using legacy prefix alias
        with eventmgr.H5FileSyntSugar(self.h5_path, prefix='L1') as f:
            f['snr'] = np.array([4.2, 8.1], dtype=np.float32)

        # 3. Overwrite dataset in first detector (del + recreate) using positional group
        with eventmgr.H5FileSyntSugar(self.h5_path, 'H1') as f:
            f['snr'] = np.array([10.0], dtype=np.float32)

        # Verify content
        with h5py.File(self.h5_path, 'r') as hf:
            self.assertIn('H1/snr', hf)
            self.assertIn('L1/snr', hf)
            np.testing.assert_allclose(hf['H1/snr'][:], [10.0], rtol=1e-5)
            np.testing.assert_allclose(hf['L1/snr'][:], [4.2, 8.1], rtol=1e-5)


if __name__ == '__main__':
    unittest.main()
