"""Unit tests for pycbc.events.eventmgr module.
"""

import os
import tempfile
import unittest
import h5py
import numpy as np

from pycbc.events import eventmgr


class MockTemplate(object):
    """Mock template object providing attributes required by EventManager."""
    def __init__(self, template_hash=12345, template_duration=1.0):
        self.template_hash = template_hash
        self.template_duration = template_duration


class MockOptions(object):
    """Mock CLI options namespace providing attributes for EventManager."""
    def __init__(self, channel_name='H1:STRAIN'):
        self.channel_name = channel_name
        self.sample_rate = 2048
        self.gps_start_time = 1000000000
        self.gps_end_time = 1000000064
        self.segment_start_pad = 0
        self.segment_end_pad = 0
        self.trig_start_time = 1000000000
        self.trig_end_time = 1000000064
        self.autochi_number_points = 10
        self.autochi_onesided = None
        self.autochi_two_phase = False
        self.autochi_max_valued_dof = None
        self.psdvar_segment = None
        self.chisq_threshold = None
        self.chisq_bins = None
        self.chisq_delta = 0
        self.newsnr_threshold = None
        self.keep_loudest_interval = None
        self.cluster_window = 0


class TestEventManager(unittest.TestCase):
    """Test EventManager multi-detector HDF5 output and lifecycle."""

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.h5_path = os.path.join(self.temp_dir.name, "events.hdf")
        self.columns = [
            'snr', 'chisq', 'bank_chisq', 'bank_chisq_dof',
            'cont_chisq', 'time_index', 'sigmasq'
        ]
        self.column_types = [
            np.float32, np.float32, np.float32, np.float32,
            np.float32, np.int64, np.float32
        ]

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_multi_ifo_eventmgr_end_to_end(self):
        # 1. Instantiate EventManager for H1 with ifo='H1'
        opt_h1 = MockOptions(channel_name='H1:STRAIN')
        mgr_h1 = eventmgr.EventManager(
            opt_h1, self.columns, self.column_types, ifo='H1'
        )
        self.assertEqual(mgr_h1.ifo, 'H1')

        mgr_h1.new_template(tmplt=MockTemplate(template_hash=101))
        h1_events = [
            np.array([8.5, 12.0], dtype=np.float32),  # snr
            np.array([1.1, 0.95], dtype=np.float32),  # chisq
            np.array([0.0, 0.0], dtype=np.float32),   # bank_chisq
            np.array([0.0, 0.0], dtype=np.float32),   # bank_chisq_dof
            np.array([0.0, 0.0], dtype=np.float32),   # cont_chisq
            np.array([2048, 4096], dtype=np.int64),   # time_index
            np.array([50.0, 60.0], dtype=np.float32),  # sigmasq
        ]
        mgr_h1.add_template_events(self.columns, h1_events)
        mgr_h1.finalize_template_events()
        mgr_h1.write_events(self.h5_path)

        # 2. Instantiate EventManager for L1 and append to the same HDF file
        opt_l1 = MockOptions(channel_name='L1:STRAIN')
        mgr_l1 = eventmgr.EventManager(
            opt_l1, self.columns, self.column_types, ifo='L1'
        )
        self.assertEqual(mgr_l1.ifo, 'L1')

        mgr_l1.new_template(tmplt=MockTemplate(template_hash=202))
        l1_events = [
            np.array([6.5], dtype=np.float32),        # snr
            np.array([1.05], dtype=np.float32),       # chisq
            np.array([0.0], dtype=np.float32),        # bank_chisq
            np.array([0.0], dtype=np.float32),        # bank_chisq_dof
            np.array([0.0], dtype=np.float32),        # cont_chisq
            np.array([3000], dtype=np.int64),         # time_index
            np.array([45.0], dtype=np.float32),       # sigmasq
        ]
        mgr_l1.add_template_events(self.columns, l1_events)
        mgr_l1.finalize_template_events()
        mgr_l1.write_events(self.h5_path)

        # 3. Instantiate via from_multi_ifo_interface for V1
        opt_multi = MockOptions()
        opt_multi.channel_name = {'V1': 'V1:STRAIN'}
        mgr_v1 = eventmgr.EventManager.from_multi_ifo_interface(
            opt_multi, 'V1', self.columns, self.column_types
        )
        self.assertEqual(mgr_v1.ifo, 'V1')
        mgr_v1.new_template(tmplt=MockTemplate(template_hash=303))
        mgr_v1.finalize_template_events()
        mgr_v1.write_events(self.h5_path)

        # 4. Verify the multi-detector HDF5 structure
        with h5py.File(self.h5_path, 'r') as hf:
            # Check groups exist
            self.assertIn('H1', hf)
            self.assertIn('L1', hf)
            self.assertIn('V1', hf)

            # Check H1 datasets
            self.assertIn('H1/snr', hf)
            np.testing.assert_allclose(
                hf['H1/snr'][:], [8.5, 12.0], rtol=1e-5
            )
            np.testing.assert_allclose(
                hf['H1/template_hash'][:], [101, 101]
            )
            self.assertIn('H1/search/start_time', hf)

            # Check L1 datasets
            self.assertIn('L1/snr', hf)
            np.testing.assert_allclose(
                hf['L1/snr'][:], [6.5], rtol=1e-5
            )
            np.testing.assert_allclose(
                hf['L1/template_hash'][:], [202]
            )
            self.assertIn('L1/search/start_time', hf)

            # Check V1 metadata
            self.assertIn('V1/search/start_time', hf)

    def test_eventmgr_single_ifo_fallback(self):
        # Verify fallback to opt.channel_name[0:2] when ifo is not passed
        opt = MockOptions(channel_name='H1:STRAIN')
        mgr = eventmgr.EventManager(opt, self.columns, self.column_types)
        self.assertIsNone(mgr.ifo)

        mgr.new_template(tmplt=MockTemplate(template_hash=999))
        mgr.finalize_template_events()
        mgr.write_events(self.h5_path)

        with h5py.File(self.h5_path, 'r') as hf:
            self.assertIn('H1/search/start_time', hf)


if __name__ == '__main__':
    unittest.main()
