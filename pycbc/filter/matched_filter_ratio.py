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
"""Hierarchical ratio/FIR matched filtering and multi-detector follow-up."""

import collections
import logging
import multiprocessing as mp

import numpy as np

from pycbc.filter.matchedfilter import get_cutoff_indices
from pycbc.libutils import import_optional
from pycbc.types import FrequencySeries, complex64, float32

logger = logging.getLogger(__name__)

mf_search = import_optional("matchedfilter.search")
mf_td = import_optional("matchedfilter.time_domain")

# Module-level reference and worker function used by fork-based multiprocessing
# to dispatch individual top-level template hierarchies across worker processes.
_WORKER_CONTROL = None


def _run_top_worker(top_idx):
    return _WORKER_CONTROL.run_top(top_idx, veto_callback=_WORKER_CONTROL.veto_callback)


class MatchedFilterRatioControl:
    """
    Control class for hierarchical ratio/FIR matched filtering and
    multi-detector follow-up using ``matchedfilter``.

    Parameters
    ----------
    bank : RatioFilterBank
        The hierarchical ratio/FIR template bank.
    segments : dict of list of FrequencySeries
        Mapping from detector name to list of frequency-domain data segments
        (each with ``.psd`` and ``.analyze`` attached).
    snr_threshold : float
        SNR threshold for primary trigger generation.
    high_frequency_cutoff : float, optional
        High frequency cutoff in Hz (defaults to Nyquist).
    cluster_window : float, optional
        Length of time-clustering window in seconds (default: 0.0).
    followup_snr_threshold : float, optional
        Secondary detector SNR threshold for follow-up triggers. Multi-detector
        follow-up is enabled only when ``followup_snr_threshold``,
        ``followup_bin_width``, and ``followup_num_bins`` are all provided.
    followup_bin_width : float, optional
        Width in seconds of each bin in the follow-up window. Required together
        with ``followup_snr_threshold`` and ``followup_num_bins`` to enable
        follow-up.
    followup_num_bins : int, optional
        Number of bins in the follow-up window. Required together with
        ``followup_snr_threshold`` and ``followup_bin_width`` to enable
        follow-up.
    false_dismissal : float, optional
        Tolerated probability of missing a signal at ``snr_threshold`` in the
        hierarchical leaf filter (default: 1e-3).
    cores : int, optional
        Number of CPU cores to use for multiprocessing over top-level templates
        (default: 1).
    """

    def __init__(
        self,
        bank,
        segments,
        snr_threshold,
        high_frequency_cutoff=None,
        cluster_window=0.0,
        followup_snr_threshold=None,
        followup_bin_width=None,
        followup_num_bins=None,
        false_dismissal=1e-3,
        cores=1,
    ):
        self.bank = bank
        self.segments = segments
        self.instruments = list(segments.keys())
        self.seg_ifos = [ifo for ifo in self.instruments if len(segments[ifo]) > 0]
        self.threshold = float(snr_threshold)
        self.delta_f = float(bank.delta_f)
        self.sample_rate = (
            float(round(1.0 / segments[self.seg_ifos[0]][0].delta_t))
            if self.seg_ifos
            else float(bank.sample_rate)
        )
        self.f_high = high_frequency_cutoff
        self.followup_threshold = followup_snr_threshold
        self.followup_bin_width = followup_bin_width
        self.followup_num_bins = followup_num_bins
        self.do_followup = (
            followup_snr_threshold is not None
            and followup_bin_width is not None
            and followup_num_bins is not None
        )
        self.cores = max(1, int(cores or 1))
        self.veto_callback = None

        # Verify that the bank's FIR tap sample rate is an integer multiple of data rate
        exact_dec = float(bank.sample_rate) / self.sample_rate
        self.decimation = round(exact_dec)
        if abs(exact_dec - self.decimation) > 1e-5 or self.decimation < 1:
            raise ValueError(
                f"Bank sample rate ({bank.sample_rate} Hz) must be an integer "
                f"multiple of data sample rate ({self.sample_rate} Hz)."
            )

        # Configure leaf TimeDomainFilterBank options and pre-overwhiten segments
        self.leaf_options = {
            "tap_sample_rate": int(bank.sample_rate),
            "data_sample_rate": int(self.sample_rate),
            "false_dismissal": float(false_dismissal),
            "binsize": round(float(cluster_window) * int(self.sample_rate))
            if cluster_window
            else None,
        }
        self.overwhitened = {
            ifo: [st / st.psd for st in segments[ifo]] for ifo in self.seg_ifos
        }

    @staticmethod
    def insert_option_group(parser):
        """Add Ratio/FIR matched-filter and multi-detector follow-up options."""
        group = parser.add_argument_group("Ratio Filter Options")
        group.add_argument(
            "--snr-threshold",
            type=float,
            default=5.5,
            help="SNR threshold for primary trigger generation (default: 5.5).",
        )
        group.add_argument(
            "--high-frequency-cutoff",
            type=float,
            help="Optional high frequency cutoff in Hz (defaults to Nyquist).",
        )
        group.add_argument(
            "--cluster-window",
            type=float,
            default=0.0,
            help="Optional length of clustering window in seconds (default: 0.0).",
        )
        group.add_argument(
            "--ratio-filter-false-dismissal",
            type=float,
            default=1e-3,
            metavar="PROB",
            help="Optional tolerated probability of missing a signal at the "
            "SNR threshold (default: 1e-3).",
        )
        group.add_argument(
            "--followup-snr-threshold",
            type=float,
            help="Optional secondary detector SNR threshold for follow-up "
            "triggers. Multi-detector follow-up runs only when "
            "--followup-snr-threshold, --followup-bin-width, and "
            "--followup-num-bins are all provided.",
        )
        group.add_argument(
            "--followup-bin-width",
            type=float,
            help="Optional width in seconds of each bin in the follow-up "
            "window (required with --followup-snr-threshold and "
            "--followup-num-bins to enable follow-up).",
        )
        group.add_argument(
            "--followup-num-bins",
            type=int,
            help="Optional number of bins in the follow-up window (required "
            "with --followup-snr-threshold and --followup-bin-width to enable "
            "follow-up).",
        )
        group.add_argument(
            "--cores",
            type=int,
            default=1,
            help="Optional number of CPU cores to use for multiprocessing "
            "over top templates (default: 1).",
        )
        return group

    @classmethod
    def from_cli(cls, args, bank, segments):
        """Construct ``MatchedFilterRatioControl`` from parsed CLI options."""
        return cls(
            bank,
            segments,
            snr_threshold=args.snr_threshold,
            high_frequency_cutoff=args.high_frequency_cutoff,
            cluster_window=args.cluster_window,
            followup_snr_threshold=args.followup_snr_threshold,
            followup_bin_width=args.followup_bin_width,
            followup_num_bins=args.followup_num_bins,
            false_dismissal=args.ratio_filter_false_dismissal,
            cores=args.cores,
        )

    def run(self, veto_callback=None, top_indices=None):
        """Run the hierarchical FIR search across ``top_indices`` (defaults to
        ``self.bank.top_indices``) and return ``(triggers, has_followup)``.
        """
        triggers = {ifo: collections.defaultdict(list) for ifo in self.instruments}
        has_followup = dict.fromkeys(self.instruments, 0)
        if not self.seg_ifos:
            return triggers, has_followup

        if top_indices is None:
            top_indices = [int(x) for x in self.bank.top_indices]
        else:
            top_indices = [int(x) for x in top_indices]

        def _collect(results_iter):
            for worker_trigs, worker_followup in results_iter:
                for ifo in self.instruments:
                    for col, arr in worker_trigs[ifo].items():
                        if len(arr) > 0:
                            triggers[ifo][col].append(arr)
                    has_followup[ifo] += worker_followup[ifo]

        n_workers = min(self.cores, max(1, len(top_indices)))
        if n_workers > 1:
            global _WORKER_CONTROL
            logger.info(
                "Filtering %d top templates across %d cores...",
                len(top_indices),
                n_workers,
            )
            self.veto_callback = veto_callback
            _WORKER_CONTROL = self
            try:
                self.bank.close_filehandler()
                with mp.get_context("fork").Pool(processes=n_workers) as pool:
                    _collect(pool.imap_unordered(_run_top_worker, top_indices))
            finally:
                _WORKER_CONTROL = None
                self.veto_callback = None
        else:
            logger.info("Filtering %d top templates...", len(top_indices))
            _collect(
                self.run_top(idx, veto_callback=veto_callback) for idx in top_indices
            )

        collapsed = {
            ifo: {
                col: np.concatenate(arrs) if arrs else np.empty(0)
                for col, arrs in triggers[ifo].items()
            }
            for ifo in self.instruments
        }
        return collapsed, has_followup

    def run_top(self, top_idx, veto_callback=None):
        """Filter a single top-level reference template hierarchy across all
        segments and detectors, returning ``(triggers, followup_counts)``.
        """
        self.veto_callback = veto_callback
        self._triggers = {
            ifo: collections.defaultdict(list) for ifo in self.instruments
        }
        self._followup = dict.fromkeys(self.instruments, 0)

        # Build the Inner/Leaf filter tree below this top-level reference template
        children, self._groups = self.bank.build_children(int(top_idx))
        if not self._groups:
            return {ifo: {} for ifo in self.instruments}, self._followup

        # Decompress and band-limit the top-level reference waveform
        self.top_template = self.bank.get_coarse_template(int(top_idx))
        first_seg = self.segments[self.seg_ifos[0]][0]
        n_series = (len(first_seg) - 1) * 2
        kmin, kmax = get_cutoff_indices(
            self.top_template.f_lower,
            self.f_high,
            self.delta_f,
            n_series,
        )
        self.top_template[:kmin] = 0
        self.top_template[kmax:] = 0

        # Calibrate the root power profile and construct the HierarchicalSearch
        profile = mf_td.reference_profile(
            self.top_template,
            first_seg.psd,
            float(self.top_template.f_lower or 0.0),
            float(self.f_high or 0.0),
            self.delta_f,
        )
        self.root = mf_search.Root(
            self.top_template,
            (kmin, kmax),
            children=children,
            profile=(profile, self.delta_f),
            node_id=int(top_idx),
        )
        self.search = mf_search.HierarchicalSearch(
            [self.root],
            threshold=self.threshold,
            leaf_options=self.leaf_options,
        )

        self._fine_correlations = {}
        self._ref_blocks = {}
        self._batch = None

        # Normalize each segment's root C2R IFFT output to unit-variance SNR
        inputs = {
            ifo: mf_search.BlockInput(
                self.overwhitened[ifo],
                [(st.analyze.start, st.analyze.stop) for st in segs],
                lambda _root, chunk, segs=segs: (
                    2.0
                    * segs[chunk].delta_f
                    * segs[chunk].delta_t
                    / (
                        np.sqrt(self.top_template.sigmasq(segs[chunk].psd))
                        * self.decimation
                    )
                ),
                tags=list(segs),
            )
            for ifo in self.seg_ifos
            for segs in (self.segments[ifo],)
        }
        self.search.run(inputs, on_peaks=self._on_peaks)

        # Concatenate per-batch trigger lists into flat arrays per detector
        collapsed = {
            ifo: {
                col: np.concatenate(arrs) if arrs else np.empty(0)
                for col, arrs in self._triggers[ifo].items()
            }
            for ifo in self.instruments
        }
        return collapsed, self._followup

    def fine_correlation(self, j, reach=0):
        """Return ``(corr, index)`` for trigger ``j`` in the current batch over
        one leaf FFT block around the trigger, re-centering if ``reach`` extends
        into the FIR wrap-around region.
        """
        ifo, leaf_bank, cid, local_idxs, t_idxs, tstarts, ref_snr = self._batch
        i = int(local_idxs[j])
        t = int(t_idxs[j])
        tstart = int(tstarts[j])
        n = leaf_bank.get_block_length(i)
        half = int(leaf_bank.effective_data_counts[i]) // 2

        # Re-center the FFT block if the requested reach enters the FIR wrap region
        if t - reach - tstart < half or t + reach - tstart >= n - half:
            tstart = min(max(t - n // 2, 0), len(ref_snr) - n)
        if t - reach - tstart < half or t + reach - tstart >= n - half:
            return None, None

        # Multiply the leaf FIR spectrum by the FFT of the parent SNR block
        key = (ifo, cid, i, tstart)
        if key not in self._fine_correlations:
            block_key = (ifo, cid, tstart, n)
            if block_key not in self._ref_blocks:
                self._ref_blocks[block_key] = np.fft.fft(
                    ref_snr[tstart : tstart + n]
                ).astype(np.complex64)
            corr = FrequencySeries(
                leaf_bank.get_filter_f(i) * self._ref_blocks[block_key] / n,
                delta_f=self.sample_rate / n,
            )
            self._fine_correlations[key] = corr
        return self._fine_correlations[key], t - tstart

    def _on_peaks(self, ctx):
        """Record primary triggers for each detector and follow up surviving
        primaries in secondary detectors.
        """
        self._fine_correlations.clear()
        self._ref_blocks.clear()
        primaries = {
            ifo: self._record_triggers(ctx, ifo, res) for ifo, res in ctx.peaks.items()
        }
        if self.do_followup and len(ctx.channels) > 1:
            self._follow_up_secondaries(ctx, primaries)

    def _record_triggers(self, ctx, ifo, res, is_followup=False):
        """Apply ``self.veto_callback`` (if set) to ``res`` and append surviving
        triggers for ``ifo``.
        """
        ids_local = res.template_indices
        t_idxs = res.sample_indices
        if len(ids_local) == 0:
            return ids_local, t_idxs

        # Map leaf-local template indices and sample lags to fine-bank indices and times
        cid = ctx.leaf.node_id
        fine_indices, shifts, _ = self._groups[cid]
        stilde = ctx.tags[ifo]
        snrs = res.snr
        tstarts = res.block_starts
        phys_t = t_idxs - (shifts[ids_local] if shifts is not None else 0)
        gids = fine_indices[ids_local]

        # Evaluate signal-consistency vetoes on the trigger batch and filter survivors
        vals = {}
        if self.veto_callback is not None:
            leaf_bank = self.search.leaf_bank(self.root, ctx.leaf)
            self._batch = (
                ifo,
                leaf_bank,
                cid,
                ids_local,
                t_idxs,
                tstarts,
                ctx.series(ifo),
            )

            def _template(j):
                tmplt = self.bank.reconstruct_fine_template(
                    self.top_template,
                    cid,
                    ids_local[j],
                    leaf_bank,
                )
                step = round(tmplt.delta_f / stilde.delta_f)
                tmplt.psd = stilde.psd[::step][: len(tmplt)]
                return tmplt

            kept, vals = self.veto_callback(
                snrs,
                phys_t,
                gids,
                _template,
                self.fine_correlation,
                stilde=stilde,
            )
            if not np.all(kept):
                ids_local = ids_local[kept]
                t_idxs = t_idxs[kept]
                snrs = snrs[kept]
                phys_t = phys_t[kept]
                gids = gids[kept]
                vals = {k: v[kept] for k, v in vals.items()}
                if len(ids_local) == 0:
                    return ids_local, t_idxs

        # Convert segment-relative sample indices to cumulative indices and record
        offset = (
            stilde.cumulative_index - stilde.analyze.start
            if getattr(stilde, "analyze", None) is not None
            else 0
        )
        sigmasq = self.top_template.sigmasq(stilde.psd) * (
            self.bank.sigma_rescale(gids) ** 2
        )
        batch = {
            "template_id": gids,
            "time_index": (phys_t + offset).astype(np.int64),
            "snr": snrs.astype(complex64),
            "sigmasq": sigmasq.astype(float32),
            "is_followup": np.full(len(ids_local), is_followup, dtype=bool),
            **vals,
        }
        for col, arr in batch.items():
            self._triggers[ifo][col].append(arr)
        return ids_local, t_idxs

    def _follow_up_secondaries(self, ctx, primaries):
        """Refilter surviving primary triggers from each detector in all other
        active detectors at ``self.followup_threshold``.
        """
        bs = max(1, round(self.followup_bin_width * self.sample_rate))
        wlen = self.followup_num_bins * bs
        half = (self.followup_num_bins // 2) * bs + bs // 2

        # Follow up surviving primary triggers from each detector in all other detectors
        for ifo_a in ctx.channels:
            ids, t_idxs = primaries[ifo_a]
            if len(ids) == 0:
                continue
            st_a = ctx.tags[ifo_a]
            t0_a, dt_a = float(st_a.start_time), float(st_a.delta_t)
            for ifo_b in ctx.channels:
                if ifo_b == ifo_a:
                    continue
                st_b = ctx.tags[ifo_b]
                t0_b, dt_b = float(st_b.start_time), float(st_b.delta_t)
                v0, v1 = st_b.analyze.start, st_b.analyze.stop
                for local_id, t_a in zip(ids, t_idxs, strict=False):
                    # Shift the follow-up window to stay inside ifo_b's valid region
                    t_center = round((t0_a + int(t_a) * dt_a - t0_b) / dt_b)
                    start = min(max(t_center - half, v0), max(v0, v1 - wlen))
                    stop = min(start + wlen, v1)
                    if start >= stop:
                        continue
                    self._followup[ifo_b] += 1
                    res = ctx.refilter(
                        ifo_b,
                        windows=slice(start, stop),
                        binsize=bs,
                        threshold=self.followup_threshold,
                        template_index=int(local_id),
                    )
                    if res.snr.size:
                        self._record_triggers(ctx, ifo_b, res, is_followup=True)
