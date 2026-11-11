# Copyright (C) 2026 AerosaneDev
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

"""Unit tests for minifollowup workflow job setup."""

import tempfile
import unittest

import igwn_segments as segments

from pycbc.workflow.configuration import WorkflowConfigParser
from pycbc.workflow.minifollowups import make_audio_file


class _DummyWorkflow:
    """Lightweight workflow container for testing job setup without Pegasus."""

    def __init__(self, cp, ifos, analysis_time):
        self.cp = cp
        self.ifos = list(ifos)
        self.analysis_time = analysis_time
        self.nodes = []

    def __iadd__(self, node):
        self.nodes.append(node)
        return self


def _make_config(ifo="H1", pad_data=8):
    cp = WorkflowConfigParser()
    cp.add_section("workflow")
    cp.set("workflow", "file-retention-level", "all_files")
    cp.add_section("workflow-ifos")
    cp.set("workflow-ifos", ifo.lower(), "")
    cp.add_section("executables")
    cp.set("executables", "page_audio", "/bin/true")
    cp.add_section("pegasus_profile")
    cp.add_section("page_audio")
    cp.set("page_audio", "time-window", "2,1")
    cp.set("page_audio", "pad-data", str(pad_data))
    sec = f"page_audio-{ifo.lower()}"
    cp.add_section(sec)
    cp.set(sec, "channel-name", f"{ifo}:GDS-CALIB_STRAIN")
    return cp


def _opt(node, name):
    opts = node._options
    return opts[opts.index(name) + 1]


class TestMakeAudioFile(unittest.TestCase):
    def setUp(self):
        self.analysis_time = segments.segment(1000000000, 1000001000)
        self.data_segs = segments.segmentlist([self.analysis_time])
        self.out_dir = tempfile.mkdtemp()
        self.workflow = _DummyWorkflow(_make_config(), ["H1"], self.analysis_time)

    def test_outputs_and_options(self):
        files = make_audio_file(
            self.workflow,
            "H1",
            1000000500.25,
            self.out_dir,
            data_segments=self.data_segs,
            tags=["FG", "1"],
        )
        # only the HTML snippet goes to the page layout
        self.assertEqual(len(files), 1)
        self.assertTrue(files[0].name.endswith(".html"))
        self.assertEqual(len(self.workflow.nodes), 1)
        node = self.workflow.nodes[0]
        # but the wav file must still be an output of the node
        exts = sorted(f.name.rsplit(".", 1)[1] for f in node.output_files)
        self.assertEqual(exts, ["html", "wav"])
        self.assertIn("FG_1", files[0].name)
        self.assertEqual(_opt(node, "--center-time"), "1000000500.25")
        self.assertEqual(_opt(node, "--gps-start-time"), "1000000400")
        self.assertEqual(_opt(node, "--gps-end-time"), "1000000600")
        self.assertEqual(_opt(node, "--time-window"), "2,1")
        self.assertEqual(_opt(node, "--channel-name"), "H1:GDS-CALIB_STRAIN")

    def test_clamps_to_data_segment_with_padding(self):
        make_audio_file(
            self.workflow,
            "H1",
            1000000010.0,
            self.out_dir,
            data_segments=self.data_segs,
        )
        node = self.workflow.nodes[0]
        self.assertEqual(_opt(node, "--gps-start-time"), "1000000008")
        self.assertEqual(_opt(node, "--gps-end-time"), "1000000110")

    def test_no_trigger_sentinel_adds_nothing(self):
        files = make_audio_file(
            self.workflow, "H1", -1.0, self.out_dir, data_segments=self.data_segs
        )
        self.assertEqual(len(files), 0)
        self.assertEqual(len(self.workflow.nodes), 0)

    def test_time_outside_segments_raises(self):
        with self.assertRaises(ValueError):
            make_audio_file(
                self.workflow,
                "H1",
                2000000000.0,
                self.out_dir,
                data_segments=self.data_segs,
            )


if __name__ == "__main__":
    unittest.main()
