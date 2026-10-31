# Copyright (C) 2026 Alex Nitz
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

"""Unit tests for multi-detector matched-filter workflow job setup."""

import os
import tempfile
import unittest

import igwn_segments as segments
import lal

from pycbc.workflow.configuration import WorkflowConfigParser
from pycbc.workflow.core import File, FileList
from pycbc.workflow.matched_filter import setup_matchedfltr_workflow


class _DummyWorkflow:
    """Lightweight workflow container for testing job setup without Pegasus."""

    def __init__(self, cp, ifos, analysis_time):
        self.cp = cp
        self.ifos = list(ifos)
        self.analysis_time = analysis_time
        self.nodes = []

    def add_node(self, node):
        self.nodes.append(node)


def _make_config(
    ifos=("H1", "L1", "V1"),
    method="WORKFLOW_MULTIPLE_IFOS",
    segment_length=512,
    start_pad=144,
    end_pad=16,
    pad_data=8,
    min_analysis_segs=2,
    max_analysis_segs=10,
    min_analysis_length=1040,
    allow_zero_padding=False,
):
    """Build a WorkflowConfigParser for matched-filter job setup tests."""
    cp = WorkflowConfigParser()
    cp.add_section("workflow")
    cp.set("workflow", "file-retention-level", "all_triggers")
    cp.set("workflow", "start-time", "10000")
    cp.set("workflow", "end-time", "30000")

    cp.add_section("workflow-ifos")
    for ifo in ifos:
        cp.set("workflow-ifos", ifo.lower(), "")
        cp.set("workflow", f"{ifo.lower()}-channel-name", f"{ifo}:GDS-CALIB_STRAIN")

    cp.add_section("executables")
    inspiral_exe = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "..", "bin", "pycbc_inspiral")
    )
    cp.set("executables", "inspiral", inspiral_exe)

    cp.add_section("workflow-matchedfilter")
    cp.set("workflow-matchedfilter", "matchedfilter-method", method)
    if min_analysis_segs is not None:
        cp.set(
            "workflow-matchedfilter",
            "min-analysis-segments",
            str(min_analysis_segs),
        )
    if max_analysis_segs is not None:
        cp.set(
            "workflow-matchedfilter",
            "max-analysis-segments",
            str(max_analysis_segs),
        )
    if min_analysis_length is not None:
        cp.set(
            "workflow-matchedfilter",
            "min-analysis-length",
            str(min_analysis_length),
        )

    cp.add_section("inspiral")
    cp.set("inspiral", "segment-length", str(segment_length))
    cp.set("inspiral", "segment-start-pad", str(start_pad))
    cp.set("inspiral", "segment-end-pad", str(end_pad))
    cp.set("inspiral", "pad-data", str(pad_data))
    if allow_zero_padding:
        cp.set("inspiral", "allow-zero-padding", "")

    for ifo in ifos:
        sec = f"inspiral-{ifo.lower()}"
        cp.add_section(sec)
        cp.set(sec, "channel-name", f"{ifo}:GDS-CALIB_STRAIN")

    return cp


def _make_banks(ifos, analysis_time, num_banks=1):
    banks = FileList([])
    for idx in range(num_banks):
        tags = [f"BANK{idx}"] if num_banks > 1 else []
        bank = File(
            list(ifos),
            "TMPLTBANK",
            analysis_time,
            file_url=f"file:///tmp/bank_{idx}.hdf",
            tags=tags,
        )
        banks.append(bank)
    return banks


def _make_datafind_outs(ifos, analysis_time):
    df_outs = FileList([])
    for ifo in ifos:
        df = File(
            ifo,
            "DATAFIND",
            analysis_time,
            file_url=f"file:///tmp/{ifo}_data.gwf",
        )
        df_outs.append(df)
    return df_outs


class TestMultiDetectorWorkflowJobSetup(unittest.TestCase):
    """Verify multi-detector workflow job preparation loses no analyzable time."""

    def setUp(self):
        self.tmp_dir = tempfile.TemporaryDirectory()
        self.analysis_time = segments.segment(10000, 30000)

    def tearDown(self):
        self.tmp_dir.cleanup()

    def _run_both_modes(
        self,
        science_segs,
        num_banks=1,
        with_datafind=True,
        **config_kwargs,
    ):
        ifos = sorted(science_segs.keys())
        banks = _make_banks(ifos, self.analysis_time, num_banks=num_banks)
        df_outs = (
            _make_datafind_outs(ifos, self.analysis_time)
            if with_datafind
            else FileList([])
        )

        cp_indep = _make_config(
            ifos=ifos, method="WORKFLOW_INDEPENDENT_IFOS", **config_kwargs
        )
        wf_indep = _DummyWorkflow(cp_indep, ifos, self.analysis_time)
        out_indep = setup_matchedfltr_workflow(
            wf_indep,
            science_segs,
            df_outs,
            banks,
            output_dir=os.path.join(self.tmp_dir.name, "indep"),
        )

        cp_multi = _make_config(
            ifos=ifos, method="WORKFLOW_MULTIPLE_IFOS", **config_kwargs
        )
        wf_multi = _DummyWorkflow(cp_multi, ifos, self.analysis_time)
        out_multi = setup_matchedfltr_workflow(
            wf_multi,
            science_segs,
            df_outs,
            banks,
            output_dir=os.path.join(self.tmp_dir.name, "multi"),
        )
        return wf_indep, out_indep, wf_multi, out_multi

    def _assert_zero_time_lost_and_valid_segments(
        self,
        science_segs,
        out_indep,
        out_multi,
        min_data_length=1040,
        left_pad=152,
        right_pad=24,
        num_banks=1,
    ):
        ifos = sorted(science_segs.keys())
        bank_tags = [f"BANK{i}" for i in range(num_banks)] if num_banks > 1 else [None]

        for btag in bank_tags:
            for ifo in ifos:
                indep_segs = [
                    f.segment
                    for f in out_indep
                    if ifo in f.ifo_list and (btag is None or btag in f.tags)
                ]
                multi_segs = [
                    f.segment
                    for f in out_multi
                    if ifo in f.ifo_list and (btag is None or btag in f.tags)
                ]

                # 1. Check no overlapping valid segments in multi-IFO mode
                total_multi_dur = sum(abs(s) for s in multi_segs)
                coalesced_multi = segments.segmentlist(multi_segs).coalesce()
                self.assertEqual(
                    total_multi_dur,
                    abs(coalesced_multi),
                    f"Overlapping trigger segments found for {ifo} (bank={btag}): "
                    f"{multi_segs}",
                )

                # 2. Check zero analyzable time is lost (exact segmentlist equality
                #    and exact total duration equality against independent-IFO mode)
                coalesced_indep = segments.segmentlist(indep_segs).coalesce()
                self.assertEqual(
                    abs(coalesced_multi),
                    abs(coalesced_indep),
                    f"Analyzable time lost or changed for {ifo}: "
                    f"multi={abs(coalesced_multi)}s vs indep={abs(coalesced_indep)}s",
                )
                self.assertEqual(
                    coalesced_multi,
                    coalesced_indep,
                    f"Valid trigger segments for {ifo} differ between "
                    f"WORKFLOW_MULTIPLE_IFOS ({coalesced_multi}) and "
                    f"WORKFLOW_INDEPENDENT_IFOS ({coalesced_indep})",
                )

        # 3. Verify every job's data segment is valid for all active IFOs
        for f in out_multi:
            data_seg = f.metadata["data_seg"]
            valid_seg = f.segment
            self.assertGreaterEqual(
                abs(data_seg),
                min_data_length,
                f"Job data_seg {data_seg} is shorter than min_data_length "
                f"({min_data_length})",
            )
            needed_seg = segments.segment(
                valid_seg[0] - left_pad, valid_seg[1] + right_pad
            )
            self.assertIn(
                needed_seg,
                data_seg,
                f"Required padded segment {needed_seg} not contained in "
                f"job data_seg {data_seg}",
            )
            for ifo in f.ifo_list:
                self.assertTrue(
                    any(data_seg in sci_seg for sci_seg in science_segs[ifo]),
                    f"Job data_seg {data_seg} for {f.ifo_list} lies outside "
                    f"{ifo} science segments {science_segs[ifo]}",
                )

    def test_identical_science_segments(self):
        """All detectors share identical science segments and are grouped together."""
        common = segments.segmentlist(
            [
                segments.segment(10000, 15000),
                segments.segment(16000, 20000),
            ]
        )
        science_segs = {
            "H1": segments.segmentlist(common),
            "L1": segments.segmentlist(common),
            "V1": segments.segmentlist(common),
        }
        _, out_indep, _, out_multi = self._run_both_modes(science_segs)
        self._assert_zero_time_lost_and_valid_segments(
            science_segs, out_indep, out_multi
        )
        # Every multi-IFO job should analyze all 3 detectors together
        for f in out_multi:
            self.assertEqual(f.ifo_list, ["H1", "L1", "V1"])
        # max_valid_length is 3520s (10 PSD segments): [10152, 14976] (4824s)
        # and [16152, 19976] (3824s) each require ceil(len / 3520) = 2 jobs.
        self.assertEqual(len(out_multi), 4)

    def test_staggered_detector_turn_on_and_off(self):
        """Staggered start/end times across H1, L1, V1 lose zero valid time."""
        # L1 starts 500s after H1 and ends 500s before H1.
        # V1 starts only 180s after H1 (< min_data_length=1040).
        science_segs = {
            "H1": segments.segmentlist([segments.segment(10000, 20000)]),
            "L1": segments.segmentlist([segments.segment(10500, 19500)]),
            "V1": segments.segmentlist([segments.segment(10180, 20000)]),
        }
        _, out_indep, _, out_multi = self._run_both_modes(science_segs)
        self._assert_zero_time_lost_and_valid_segments(
            science_segs, out_indep, out_multi
        )
        # Each transition slice shorter than max_valid_length (3520s) is covered
        # by a single job rather than tiled into multiple jobs:
        # [10152, 10332] (H1, 1 job), [10332, 10652] (H1+V1, 1 job),
        # [10652, 19476] (H1+L1+V1, 3 jobs), [19476, 19976] (H1+V1, 1 job).
        self.assertEqual(len(out_multi), 6)

    def test_detector_dropout_and_short_unanalyzable_segment(self):
        """A detector dropping out and returning mid-run causes no time loss in other IFOs."""
        science_segs = {
            "H1": segments.segmentlist([segments.segment(10000, 22000)]),
            "L1": segments.segmentlist([segments.segment(10000, 22000)]),
            "V1": segments.segmentlist(
                [
                    segments.segment(10000, 14000),
                    # Short segment (< min_data_length=1040) that cannot be analyzed
                    segments.segment(14100, 14500),
                    segments.segment(15000, 22000),
                ]
            ),
        }
        _, out_indep, _, out_multi = self._run_both_modes(science_segs)
        self._assert_zero_time_lost_and_valid_segments(
            science_segs, out_indep, out_multi
        )

    def test_short_science_segment_handoff(self):
        """Handoff where two science segments overlap by less than min_data_length."""
        # H1 [10000, 12000] and L1 [11400, 13400] overlap on [11400, 12000] (600s),
        # which is shorter than min_data_length (1040s), but both H1 and L1 have
        # valid triggers on [11552, 11976]. They must be split into separate jobs
        # on that slice so neither loses valid trigger time.
        science_segs = {
            "H1": segments.segmentlist([segments.segment(10000, 12000)]),
            "L1": segments.segmentlist([segments.segment(11400, 13400)]),
        }
        _, out_indep, _, out_multi = self._run_both_modes(science_segs)
        self._assert_zero_time_lost_and_valid_segments(
            science_segs, out_indep, out_multi
        )

    def test_split_banks_and_cli_arguments(self):
        """Multiple split template banks and per-IFO frame/channel arguments."""
        science_segs = {
            "H1": segments.segmentlist([segments.segment(10000, 15000)]),
            "L1": segments.segmentlist([segments.segment(10500, 15000)]),
        }
        _, out_indep, wf_multi, out_multi = self._run_both_modes(
            science_segs, num_banks=2, with_datafind=True
        )
        self._assert_zero_time_lost_and_valid_segments(
            science_segs, out_indep, out_multi, num_banks=2
        )

        # Verify CLI arguments on multi-IFO nodes
        multi_ifo_nodes = [
            n for n in wf_multi.nodes if n.executable.ifo_list == ["H1", "L1"]
        ]
        self.assertGreater(len(multi_ifo_nodes), 0)
        for node in multi_ifo_nodes:
            self.assertIn("--channel-name", node._options)
            ch_idx = node._options.index("--channel-name")
            self.assertEqual(
                node._options[ch_idx + 1],
                "H1:GDS-CALIB_STRAIN L1:GDS-CALIB_STRAIN",
            )
            raw_str = "".join(str(x) for x in node._raw_options)
            self.assertIn("--frame-files", raw_str)
            self.assertIn("H1:", raw_str)
            self.assertIn("L1:", raw_str)

    def test_allow_zero_padding(self):
        """Zero-padding mode preserves all valid time and extends data segments."""
        science_segs = {
            "H1": segments.segmentlist([segments.segment(10000, 15000)]),
            "L1": segments.segmentlist([segments.segment(10300, 14800)]),
        }
        _, out_indep, _, out_multi = self._run_both_modes(
            science_segs, allow_zero_padding=True
        )
        # With allow_zero_padding=True, left_pad = right_pad = pad_data = 8
        self._assert_zero_time_lost_and_valid_segments(
            science_segs,
            out_indep,
            out_multi,
            min_data_length=720,
            left_pad=8,
            right_pad=8,
        )

    def test_ligo_time_gps_multi_job_slicing(self):
        """LIGOTimeGPS segment endpoints split cleanly across multiple jobs."""
        gps_seg = segments.segment(lal.LIGOTimeGPS(10000), lal.LIGOTimeGPS(20000))
        science_segs = {
            "H1": segments.segmentlist([gps_seg]),
            "L1": segments.segmentlist([gps_seg]),
        }
        _, out_indep, _, out_multi = self._run_both_modes(science_segs)
        self.assertGreater(len(out_multi), 1)
        self._assert_zero_time_lost_and_valid_segments(
            science_segs, out_indep, out_multi
        )


if __name__ == "__main__":
    unittest.main()
