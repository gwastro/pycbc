"""Unit tests for high-frequency cutoff handling in PSD and range plots."""

import os
import subprocess
import sys
import tempfile
import unittest

import numpy as np

import pycbc
import pycbc.filter
import pycbc.types
import pycbc.waveform
from pycbc.io.hdf import HFile

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


class TestRangeAndPSDCutoff(unittest.TestCase):
    """Test high-frequency cutoff in sigma and PSD/range plotting scripts."""

    def test_stopband_cutoff_in_sigma_and_plots(self):
        """Verify --high-frequency-cutoff excludes anti-aliasing stopbands."""
        df = 1.0
        n_pts = 513  # 0 to 512 Hz
        flow = 30.0
        fhigh = 400.0

        # In-band scaled PSD with steep stopband rolloff above 460 Hz
        psd_data = np.full(n_pts, 1e8, dtype=np.float32)
        psd_data[int(460 / df) :] = 1e-4
        psd = pycbc.types.FrequencySeries(psd_data, delta_f=df)

        out = pycbc.types.zeros(n_pts, dtype=np.complex64)
        htilde = pycbc.waveform.get_waveform_filter(
            out,
            mass1=1.4,
            mass2=1.4,
            approximant="SPAtmplt",
            f_lower=flow,
            delta_f=df,
            delta_t=1.0 / (2 * 512),
            distance=1.0 / pycbc.DYN_RANGE_FAC,
        ).astype(np.complex64)

        sig_unbounded = pycbc.filter.sigma(htilde, psd=psd, low_frequency_cutoff=flow)
        sig_bounded = pycbc.filter.sigma(
            htilde, psd=psd, low_frequency_cutoff=flow, high_frequency_cutoff=fhigh
        )
        self.assertGreater(sig_unbounded, 10.0 * sig_bounded)
        self.assertGreater(sig_bounded, 0.0)

        # Write minimal merged-PSD HDF5 file and exercise both plotting scripts
        with tempfile.TemporaryDirectory() as tmpdir:
            psd_path = os.path.join(tmpdir, "H1-PSD.hdf")
            with HFile(psd_path, "w") as h5f:
                h5f.attrs["low_frequency_cutoff"] = flow
                h5f.attrs["dynamic_range_factor"] = pycbc.DYN_RANGE_FAC
                ds = h5f.create_dataset("H1/psds/0", data=psd_data)
                ds.attrs["delta_f"] = df
                ds.attrs["epoch"] = 1000000000
                h5f["H1/start_time"] = np.array([1000000000], dtype=np.uint32)
                h5f["H1/end_time"] = np.array([1000000512], dtype=np.uint32)

            spec_png = os.path.join(tmpdir, "spectrum.png")
            subprocess.run(
                [
                    sys.executable,
                    os.path.join(REPO_ROOT, "bin", "plotting", "pycbc_plot_psd_file"),
                    "--psd-files",
                    psd_path,
                    "--high-frequency-cutoff",
                    str(fhigh),
                    "--output-file",
                    spec_png,
                ],
                check=True,
            )
            self.assertGreater(os.path.getsize(spec_png), 0)

            range_png = os.path.join(tmpdir, "range.png")
            subprocess.run(
                [
                    sys.executable,
                    os.path.join(REPO_ROOT, "bin", "plotting", "pycbc_plot_range"),
                    "--psd-files",
                    psd_path,
                    "--mass1",
                    "1.4",
                    "--mass2",
                    "1.4",
                    "--approximant",
                    "SPAtmplt",
                    "--high-frequency-cutoff",
                    str(fhigh),
                    "--output-file",
                    range_png,
                ],
                check=True,
            )
            self.assertGreater(os.path.getsize(range_png), 0)


if __name__ == "__main__":
    unittest.main()
