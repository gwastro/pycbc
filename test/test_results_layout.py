"""Tests for the well pages written by pycbc.results.layout"""

import os
import tempfile
import unittest
from types import SimpleNamespace

from utils import parse_args_cpu_only, simple_exit

from pycbc.results import layout

parse_args_cpu_only("results.layout")


class TestTwoColumnLayout(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.addCleanup(os.chdir, os.getcwd())
        os.chdir(tmp.name)
        self.results_dir = os.path.join(os.getcwd(), "results")
        self.cols = [(SimpleNamespace(name="H1-PLOT-1-2.png"), None)]

    def expected_reference(self):
        return f"path_exists('{self.results_dir}/H1-PLOT-1-2.png')"

    def test_default_location(self):
        layout.two_column_layout(self.results_dir, self.cols)
        with open(os.path.join(self.results_dir, "well.html")) as well:
            self.assertIn(self.expected_reference(), well.read())

    def test_filename_in_working_directory(self):
        # A workflow job is given a bare file name and writes to its working
        # directory, but the plots are not there: the well is staged out to
        # sit next to them in the results directory.
        filename = "well_H1_TEST-1-2.html"
        layout.two_column_layout(self.results_dir, self.cols, filename=filename)
        with open(filename) as well:
            self.assertIn(self.expected_reference(), well.read())


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestTwoColumnLayout))

if __name__ == "__main__":
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
