"""Unit tests for bounded timeslide coincidence and multi-detector trigger reading."""

import os
import tempfile
import unittest

import h5py
import numpy as np

from pycbc.events import coinc, significance
from pycbc.io.hdf import ReadByTemplate
from pycbc.results.layout import two_column_layout
from pycbc.workflow.core import File as WorkflowFile


class TestBoundedCoincidence(unittest.TestCase):
    """Test bounded timeslide coincidence and multi-detector trigger reading."""

    def test_time_coincidence_max_slide(self):
        """Verify time_coincidence restricts slides to |slide| <= max_slide."""
        t1 = np.array([100.0, 200.0])
        t2 = np.array([100.0, 100.5, 102.0])

        # Unbounded: slides at -20, -5, 0 (from t1=100) and 980, 995, 1000 (from t1=200)
        _, _, slide_all = coinc.time_coincidence(t1, t2, window=0.01, slide_step=0.1)
        self.assertEqual(set(slide_all.tolist()), {-20, -5, 0, 980, 995, 1000})

        # Bounded to max_slide = 10 -> slides outside [-10, 10] are excluded
        _, _, slide = coinc.time_coincidence(
            t1, t2, window=0.01, slide_step=0.1, max_slide=10
        )
        self.assertEqual(set(slide.tolist()), {-5, 0})
        self.assertTrue(np.all(np.abs(slide) <= 10))

        # Foreground only (max_slide = 0)
        idx1_z, idx2_z, slide_z = coinc.time_coincidence(
            t1, t2, window=0.01, slide_step=0.1, max_slide=0
        )
        self.assertEqual(slide_z.tolist(), [0])
        self.assertEqual(idx1_z.tolist(), [0])
        self.assertEqual(idx2_z.tolist(), [0])

    def test_time_multi_coincidence_max_slide(self):
        """Verify time_multi_coincidence forwards max_slide to pairwise coincidence."""
        times = {
            "H1": np.array([100.0, 200.0]),
            "L1": np.array([100.0, 100.5, 102.0]),
            "V1": np.array([100.0, 105.0, 110.0]),
        }
        ids, slide = coinc.time_multi_coincidence(
            times, slide_step=0.1, slop=0.01, pivot="L1", fixed="H1", max_slide=10
        )
        self.assertTrue(np.all(np.abs(slide) <= 10))
        self.assertEqual(set(slide.tolist()), {0, 5})
        self.assertEqual(len(ids), 3)

    def test_read_by_template_multi_ifo_file(self):
        """Verify ReadByTemplate selects the requested IFO from a combined trigger file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            trig_path = os.path.join(tmpdir, "combined.hdf")
            with h5py.File(trig_path, "w") as tf:
                for ifo, snr_val in (("H1", 6.5), ("L1", 7.5)):
                    grp = tf.create_group(ifo)
                    dset = grp.create_dataset(
                        "snr", data=np.array([snr_val], dtype=np.float32)
                    )
                    ref_dt = h5py.special_dtype(ref=h5py.RegionReference)
                    ref_ds = grp.create_dataset("snr_template", (1,), dtype=ref_dt)
                    ref_ds[0] = dset.regionref[0:1]
                    sgrp = grp.create_group("search")
                    sgrp["start_time"] = np.array([90], dtype=np.int64)
                    sgrp["end_time"] = np.array([110], dtype=np.int64)

            rbt_l1 = ReadByTemplate(trig_path, ifo="L1")
            self.assertEqual(rbt_l1.ifo, "L1")
            np.testing.assert_allclose(rbt_l1.get_data("snr", 0), [7.5])

    def test_bounded_background_time_and_decimation_scaling(self):
        """Verify bounded timeslide duration and non-divisible decimation weights."""
        # Two disjoint 500 s coinc segments -> 1000 s foreground time
        c_start = np.array([1000.0, 2000.0])
        c_end = np.array([1500.0, 2500.0])
        max_slide = 12
        ts_interval = 0.1
        offsets = np.arange(1, max_slide + 1) * ts_interval
        bg_time = 2.0 * float(
            coinc.timeslide_durations(c_start, c_start, c_end, c_end, offsets).sum()
        )
        # Each slide s loses s * 0.1 s per disjoint segment (2 segments)
        expected_bg = 2.0 * sum(1000.0 - 2.0 * s * ts_interval for s in range(1, 13))
        self.assertAlmostEqual(bg_time, expected_bg, places=6)

        # Non-divisible decimation factor F=5 on max_slide=12 keeps 2*(12//5)=4 slides
        # out of 24; exact weight 12 / (12 // 5) = 6.0 recovers the full-slide rate
        factor = 5
        eff_weight = max_slide / (max_slide // factor)
        self.assertAlmostEqual(eff_weight, 6.0)
        _, fg_far, _ = significance.get_far(
            np.array([8.0, 9.0, 10.0, 11.0]),
            np.array([8.5]),
            np.full(4, eff_weight),
            bg_time,
        )
        self.assertAlmostEqual(float(fg_far[0]), (3 * eff_weight + 1.0) / bg_time)

    def test_two_column_layout_custom_filename_preserves_dir(self):
        """Verify two_column_layout embeds target directory when filename is given."""
        with tempfile.TemporaryDirectory() as tmpdir:
            target_dir = os.path.join(tmpdir, "subsection")
            os.makedirs(target_dir)
            out_well = os.path.join(tmpdir, "well_custom.html")
            wf_file = WorkflowFile("H1", "TEST_PLOT", (1000, 100), "test.png")
            two_column_layout(target_dir, [(wf_file,)], filename=out_well)
            with open(out_well, encoding="utf-8") as fh:
                content = fh.read()
            self.assertIn(os.path.join(target_dir, "test.png"), content)


if __name__ == "__main__":
    unittest.main()
