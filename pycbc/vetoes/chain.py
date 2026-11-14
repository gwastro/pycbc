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
"""Run signal-consistency tests and thresholds on trigger batches."""

import logging
import time
import types
from collections import defaultdict

import numpy as np

from pycbc.events import cuts, ranking
from pycbc.io.hdf import chisq_choices, get_chisq_from_file_choice
from pycbc.vetoes.autochisq import SingleDetAutoChisq
from pycbc.vetoes.bank_chisq import SingleDetBankVeto
from pycbc.vetoes.chisq import SingleDetPowerChisq
from pycbc.vetoes.sgchisq import SingleDetSGChisq

logger = logging.getLogger("pycbc.vetoes.chain")


class SignalConsistencyChain:
    """
    Run single-detector signal-consistency tests and intermediate thresholds
    in a configurable sequence.

    Each step is either a test name (``'power'``, ``'bank'``, ``'auto'``,
    ``'sg'``) or a threshold ``'PARAMETER:VALUE:LIMIT:ACTION'`` where
    ``ACTION`` is ``'remove'`` (drop failing primary triggers from output;
    treated as ``'skip'`` on secondary follow-up triggers) or ``'skip'``
    (keep failing triggers with default values for later tests).

    Parameters
    ----------
    *estimators : object
        Active signal-consistency test instances.
    steps : list of str, optional
        Ordered sequence of stage names and threshold specifications.
    window : int, optional
        Sample window within which triggers of the same template share a
        chi-square evaluation (default: 0).
    max_iterations : int, optional
        Maximum chi-square evaluations per stage per batch (default: None).
    input_columns : iterable of str, optional
        Pre-existing trigger columns available to threshold steps.
    defaults : dict, optional
        Override default column fill values for skipped triggers.
    **named_estimators : object
        Additional signal-consistency test instances passed by keyword.
    """

    TESTS = (
        SingleDetPowerChisq,
        SingleDetBankVeto,
        SingleDetAutoChisq,
        SingleDetSGChisq,
    )
    columns = {col: dtype for cls in TESTS for col, dtype in cls.columns.items()}

    def __init__(
        self,
        *estimators,
        steps=None,
        window=0,
        max_iterations=None,
        input_columns=(),
        defaults=None,
        **named_estimators,
    ):
        self.stages = {
            est.name: est
            for est in (*estimators, *named_estimators.values())
            if est is not None and est.do
        }

        steps = steps or []
        self.description = " ".join(steps)

        self.window = window
        self.max_iterations = max_iterations
        self.defaults = defaults or {}
        self.stage_times = {cls.name: 0.0 for cls in self.TESTS}

        available = {"snr"} | set(input_columns)
        self.steps = []
        for step in steps:
            if ":" not in step:
                if step not in self.stages:
                    raise ValueError(
                        f"Stage '{step}' is in the chain but its options are not set"
                    )
                stage = self.stages[step]
                self.steps.append(stage)
                available |= set(stage.columns)
            else:
                param, func, value, action = self._parse_threshold(step)
                missing = set(self._required(param)) - available
                if missing:
                    raise ValueError(
                        f"'{step}' needs {sorted(missing)}, "
                        "which no earlier step provides"
                    )
                self.steps.append((param, func, value, action, step))
        logger.info("Signal consistency chain: %s", self.description)

    @staticmethod
    def _parse_threshold(step):
        fields = step.split(":")
        if len(fields) != 4 or fields[3].lower() not in ("remove", "skip"):
            raise ValueError(
                "Expected PARAMETER:VALUE:LIMIT:ACTION with ACTION "
                f"'remove' or 'skip'; got '{step}'"
            )
        action = fields.pop().lower()
        choices = (
            list(ranking.sngls_ranking_function_dict)
            + list(chisq_choices)
            + ["snr", "psd_var_val"]
        )
        (param, func), value = next(
            iter(cuts.convert_inputstr(":".join(fields), choices).items())
        )
        return param, func, value, action

    @staticmethod
    def _required(param):
        if param in ranking.sngls_ranking_function_dict:
            return ranking.reqd_datasets.get(param, [])
        if param in chisq_choices:
            probe = defaultdict(lambda: np.full(1, 2.0))
            get_chisq_from_file_choice(probe, param)
            return list(probe)
        return [param]

    def values(
        self,
        snr,
        time_index,
        template_id,
        template,
        correlation,
        stilde=None,
        columns=None,
        is_followup=False,
    ):
        """
        Evaluate the chain on a batch of triggers from one data segment.

        Parameters
        ----------
        snr : array_like
            Complex SNR values for each trigger in the batch.
        time_index : array_like
            Sample indices for each trigger in the batch.
        template_id : array_like
            Template indices for each trigger in the batch.
        template : callable
            Callback ``template(j)`` returning the frequency-domain template.
        correlation : callable
            Callback ``correlation(j, reach)`` returning ``(corr, index)``.
        stilde : FrequencySeries, optional
            Segment frequency series with attached ``.psd``.
        columns : dict, optional
            Mapping of pre-existing trigger column arrays.
        is_followup : bool, optional
            If True, ``'remove'`` threshold steps skip remaining stages
            without dropping the sub-threshold follow-up triggers (default:
            False).

        Returns
        -------
        keep : numpy.ndarray
            Boolean mask of retained triggers.
        vals : dict of str to numpy.ndarray
            Evaluated signal-consistency output column arrays.
        """
        batch = types.SimpleNamespace(
            snr=np.asarray(snr),
            time_index=np.asarray(time_index),
            template_id=np.asarray(template_id),
            template=template,
            correlation=correlation,
            stilde=stilde,
            columns=columns or {},
        )
        n = len(batch.snr)
        keep = np.ones(n, dtype=bool)
        alive = np.ones(n, dtype=bool)
        vals = {
            col: np.full(
                n,
                self.defaults.get(col, src.default_values[col]),
                dtype=dtype,
            )
            for src in (*self.TESTS, *self.stages.values())
            for col, dtype in src.columns.items()
        }
        for step in self.steps:
            if isinstance(step, tuple):
                idx = np.flatnonzero(alive)
                failed = idx[~self._threshold(step, batch, vals, idx)]
                alive[failed] = False
                remove = step[3] == "remove" and not is_followup
                if remove:
                    keep[failed] = False
                logger.debug(
                    "%s: %d of %d triggers %s",
                    step[4],
                    len(failed),
                    len(idx),
                    "removed" if remove else "skipped",
                )
            else:
                self._run_stage(step, batch, vals, alive)
        return keep, vals

    def _threshold(self, step, batch, vals, idx):
        param, func, value = step[:3]
        trigs = {"snr": np.abs(batch.snr[idx])}
        for col, arr in vals.items():
            trigs[col] = arr[idx]
        for col, arr in batch.columns.items():
            trigs[col] = np.asarray(arr)[idx]
        if "chisq_dof" in trigs:
            trigs["chisq_dof"] = trigs["chisq_dof"] / 2 + 1

        if param in ranking.sngls_ranking_function_dict:
            stat = ranking.get_sngls_ranking_from_trigs(trigs, param)
        elif param in chisq_choices:
            stat = get_chisq_from_file_choice(trigs, param)
        else:
            stat = trigs[param]
        return func(stat, value)

    def _run_stage(self, stage, batch, vals, alive):
        if not np.any(alive):
            return
        start = time.perf_counter()
        priority = np.abs(batch.snr)
        times = batch.time_index
        tids = batch.template_id
        assigned = ~alive
        computed = 0
        for j in np.argsort(priority)[::-1]:
            if not alive[j]:
                continue
            if self.max_iterations is not None and computed >= self.max_iterations:
                break
            if assigned[j]:
                continue
            out = stage.evaluate_batch(batch, j)
            if out is None:
                continue
            share = (
                (np.abs(times - times[j]) <= self.window)
                & (tids == tids[j])
                & ~assigned
            )
            for col, val in out.items():
                vals[col][share] = val
            assigned |= share
            computed += 1
        alive &= assigned
        elapsed = time.perf_counter() - start
        self.stage_times[stage.name] = self.stage_times.get(stage.name, 0.0) + elapsed
        logger.debug(
            "%s chisq: %d evaluations in %.3f s",
            stage.name,
            computed,
            elapsed,
        )

    @classmethod
    def insert_option_group(cls, parser):
        names = ", ".join(f"'{t.name}'" for t in cls.TESTS)
        group = parser.add_argument_group("Signal consistency tests")
        group.add_argument(
            "--chisq-window",
            type=float,
            default=1.0,
            help="Triggers of the same template within this window (s) share "
            "a computed chi-square value.",
        )
        group.add_argument(
            "--chisq-iterations",
            type=int,
            default=10000,
            help="Maximum evaluations per test per batch, loudest trigger first.",
        )
        group.add_argument(
            "--signal-consistency-chain",
            nargs="+",
            metavar="STEP",
            help=f"Optional ordered steps of the chain: test names ({names}) "
            "or thresholds PARAMETER:VALUE:LIMIT:ACTION with ACTION 'remove' "
            "(drop failing primary triggers; treated as 'skip' on secondary "
            "follow-up triggers) or 'skip' (e.g. "
            "'snr:5.25:lower:skip power newsnr:5:lower:remove sg').",
        )
        for test_cls in cls.TESTS:
            test_cls.insert_option_group(parser)

    @classmethod
    def from_cli(cls, args, bank, sample_rate, defaults=None):
        tests = [test_cls.from_cli(args, bank=bank) for test_cls in cls.TESTS]
        extra = ["psd_var_val"] if getattr(args, "psdvar_segment", None) else []
        return cls(
            *tests,
            steps=args.signal_consistency_chain,
            window=int(args.chisq_window * sample_rate),
            max_iterations=args.chisq_iterations,
            input_columns=extra,
            defaults=defaults,
        )
