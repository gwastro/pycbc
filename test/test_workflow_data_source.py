"""Unit tests for checking a workflow job's source of data."""

import sys
import unittest

from pycbc.workflow.configuration import WorkflowConfigParser
from pycbc.workflow.core import Executable


def make_node(section_options=None):
    cp = WorkflowConfigParser()
    cp.add_section("workflow")
    cp.add_section("executables")
    cp.set("executables", "inspiral", sys.executable)
    cp.add_section("pegasus_profile")
    cp.add_section("inspiral")
    for opt, value in (section_options or {}).items():
        cp.set("inspiral", opt, value)
    exe = Executable(cp, "inspiral", ifos="H1", out_dir=".")
    return exe.create_node()


class TestOneDataSource(unittest.TestCase):
    def assert_rejected(self, node, *names):
        with self.assertRaises(ValueError) as caught:
            node._finalize()
        for name in names:
            self.assertIn(name, str(caught.exception))

    def test_single_source(self):
        node = make_node({"frame-type": "H1_HOFT_C00"})
        node._finalize()

    def test_no_source(self):
        node = make_node()
        node.add_opt("--sample-rate", 1024)
        node._finalize()

    def test_unrelated_similar_option(self):
        node = make_node({"fake-strain": "zeroNoise", "fake-strain-seed": "1"})
        node._finalize()

    def test_two_sources(self):
        node = make_node({"frame-type": "H1_HOFT_C00"})
        node.add_list_opt("--frame-files", ["a.gwf", "b.gwf"])
        self.assert_rejected(node, "--frame-type", "--frame-files")

    def test_repeated_source_through_add_arg(self):
        node = make_node()
        node.add_list_opt("--frame-files", ["a.gwf"])
        node.add_arg("--frame-files b.gwf")
        self.assert_rejected(node, "--frame-files")

    def test_multi_ifo_single_source(self):
        node = make_node()
        node.add_arg("--frame-files H1:a.gwf H1:b.gwf L1:c.gwf")
        node._finalize()

    def test_multi_ifo_different_ifos_may_differ(self):
        node = make_node()
        node.add_arg("--frame-files H1:a.gwf")
        node.add_list_opt("--fake-strain", ["L1:zeroNoise"])
        node._finalize()

    def test_multi_ifo_same_ifo_two_sources(self):
        node = make_node()
        node.add_arg("--frame-files H1:a.gwf L1:b.gwf")
        node.add_list_opt("--fake-strain", ["L1:zeroNoise"])
        self.assert_rejected(node, "for L1", "--frame-files", "--fake-strain")

    def test_multi_ifo_raw_options(self):
        node = make_node()
        node.add_raw_arg("--frame-files")
        node.add_raw_arg(" ")
        node.add_raw_arg("H1:a.gwf")
        node.add_raw_arg(" ")
        node.add_list_opt("--frame-type", ["H1:H1_HOFT_C00"])
        self.assert_rejected(node, "for H1")


if __name__ == "__main__":
    unittest.main()
