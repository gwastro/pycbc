"""Mock simulation to easily test and profile PyCBC Live's coincidence code."""

import unittest
from types import SimpleNamespace
import numpy as np
import copy
import logging
from astropy.utils.data import download_file
from pycbc import gps_now
from pycbc.events.coinc import LiveCoincTimeslideBackgroundEstimator as Coincer
from utils import simple_exit
import validation_code.old_coinc as old_coinc

OriginalCoincer = old_coinc.LiveCoincTimeslideBackgroundEstimator

np.random.seed(0)

class SingleDetTrigSimulator:
    """An object that simulates single-detector triggers in the same format
    as produced by the matched-filtering processes of PyCBC Live.
    """
    def __init__(self, num_templates, analysis_chunk, detectors, num_trigs_per_block):
        self.num_templates = num_templates
        self.detectors = detectors
        self.analysis_chunk = analysis_chunk
        self.start_time = int(1448294502) #gps_now()
        self.num_trigs = num_trigs_per_block

    def get_trigs(self):
        trigs = {}
        for det in self.detectors:
            rand_end = np.random.randint(
                self.start_time*4096,
                (self.start_time + self.analysis_chunk)*4096,
                size=self.num_trigs
            )
            rand_end = (rand_end / 4096.).astype(np.float64)

            trigs[det] = {
                "snr": np.random.uniform(4.5, 10, size=self.num_trigs).astype(np.float32),
                "end_time": rand_end,
                "chisq": np.random.uniform(0.5, 1.5, size=self.num_trigs).astype(np.float32),
                "chisq_dof": np.ones(self.num_trigs, dtype=np.int32) * 10,
                "coa_phase": np.random.uniform(0, 2*np.pi, size=self.num_trigs).astype(np.float32),
                "sigmasq": np.ones(self.num_trigs, dtype=np.float32),  # FIXME (maybe)
                "template_id": np.random.uniform(
                    0,
                    self.num_templates,
                    size=self.num_trigs
                ).astype(np.int32),
                "mass1": np.random.uniform(2.0, 100.0, size=self.num_trigs).astype(np.float32),
                "mass2": np.random.uniform(2.0, 100.0, size=self.num_trigs).astype(np.float32)
            }
        self.start_time += self.analysis_chunk
        return trigs


class TestPyCBCLiveCoinc(unittest.TestCase):
    def setUp(self, *args):
        # Uncomment for more verbosity
        # logging.basicConfig(format="%(asctime)s %(message)s",
        #                     level=logging.INFO)

        # simulate the `args` object we normally get from the command line arguments

        url = 'https://github.com/gwastro/pycbc-config/raw/master/'
        url += 'test_data_files/{}-PTA_HISTOGRAM.hdf'
        stat_file_paths = [
            download_file(url.format("H1L1"), cache=True),
        ]
        args = SimpleNamespace(
            sngl_ranking="snr",
            ranking_statistic="phasetd",
            statistic_files=[stat_file_paths],
            statistic_keywords=None,
            statistic_features=None,
            timeslide_interval=0.1,
            background_ifar_limit=100,
            store_background=True,
            coinc_window_pad=0.002,
            statistic_refresh_rate=None,
        )

        # save args and some parameters for possible replay tracing
        self._test_args = args
        self._test_analysis_chunk = None
        self._test_detectors = None

        # number of templates in the bank
        self.num_templates = 10

        # duration of analysis segment
        analysis_chunk = 2000

        # combination of two detectors to analyze
        detectors = ["H1", "L1"]

        # store for replay mapping
        self._test_analysis_chunk = analysis_chunk
        self._test_detectors = detectors

        # number of single-detector triggers per detector per chunk
        num_single_trigs = 400

        self.num_iterations = 15

        # create the single-detector trigger simulator
        single_det_trig_sim = SingleDetTrigSimulator(
            self.num_templates, analysis_chunk, detectors, num_single_trigs
        )

        self.new_trigs = [single_det_trig_sim.get_trigs()
                          for _ in range(self.num_iterations)]

        # create the current "coincer" object
        self.new_coincer = Coincer.from_cli(args, self.num_templates,
                                            analysis_chunk, detectors)

        # create the validation "coincer" object
        self.old_coincer = OriginalCoincer.from_cli(args, self.num_templates,
                                                    analysis_chunk, detectors)


    def test_coincer_runs(self):
        # the following loop simulates the "infinite" analysis loop
        # (though we only do a few iterations here)

        def assess_same_output(idx, newout, oldout, new_coincer, old_coincer, trigs):
            checkkeys = [
                'background/time',
                'background/count',
                'background/stat',
                'foreground/ifar',
                'foreground/stat',
                'foreground/type'
            ]

            for ifo in ['H1', 'L1']:
                checkkeys += [
                    f'foreground/{ifo}/snr',
                    f'foreground/{ifo}/end_time',
                    f'foreground/{ifo}/chisq',
                    f'foreground/{ifo}/chisq_dof',
                    f'foreground/{ifo}/coa_phase',
                    f'foreground/{ifo}/sigmasq',
                    f'foreground/{ifo}/template_id',
                    f'foreground/{ifo}/stat'
                ]

            for key in checkkeys:
                if key not in newout:
                    self.assertTrue(key not in oldout)
                else:
                    self.assertTrue(key in oldout)
                    print(i, key)
                    if type(newout[key]) is np.ndarray:
                        print(len(newout[key]), len(oldout[key]))
                        self.assertTrue(len(newout[key]) == len(oldout[key]))
                        diff = newout[key] - oldout[key]
                        print(diff.mean(), diff.std())
                        not_close = np.logical_not(np.isclose(newout[key], oldout[key]))
                        print(sum(not_close), np.flatnonzero(not_close))
                        print(newout[key][not_close], oldout[key][not_close])

                        # Dump deeper internal state to help debug where the
                        # divergence originates. This will compare single-det
                        # statistic outputs and per-template buffers.
                        def debug_state():
                            print('\n==== DEBUG STATE DUMP ====', flush=True)
                            # Compare single-det statistic outputs for each ifo
                            for ifo in ['H1', 'L1']:
                                print(f'-- IFO {ifo} single-stat comparison --')
                                try:
                                    trigsc = copy.copy(trigs[ifo])
                                    trigsc['ifo'] = ifo
                                    trigsc['chisq'] = trigs[ifo]['chisq'] * trigs[ifo]['chisq_dof']
                                    trigsc['chisq_dof'] = (trigs[ifo]['chisq_dof'] + 2) / 2
                                except Exception as e:
                                    print('  failed preparing trigsc:', e)
                                    continue

                                try:
                                    new_single = new_coincer.stat_calculator.single(trigsc)
                                except Exception as e:
                                    print('  new stat_calculator.single() failed:', e)
                                    new_single = None
                                try:
                                    old_single = old_coincer.stat_calculator.single(trigsc)
                                except Exception as e:
                                    print('  old stat_calculator.single() failed:', e)
                                    old_single = None

                                if new_single is None or old_single is None:
                                    print('  one of the single outputs is None')
                                else:
                                    try:
                                        eq = np.isclose(new_single, old_single)
                                        print('  single shapes:', new_single.shape, old_single.shape)
                                        if eq.all():
                                            print('  single arrays match (all close)')
                                        else:
                                            nc = np.flatnonzero(np.logical_not(eq))
                                            print('  single arrays differ at', len(nc), 'positions; first few indices:', nc[:10])
                                            print('  new_single[diffs]:', new_single[nc[:10]])
                                            print('  old_single[diffs]:', old_single[nc[:10]])
                                    except Exception as e:
                                        print('  error comparing single arrays:', e)

                            # Compare per-template buffers for each ifo
                            for ifo in ['H1', 'L1']:
                                print(f'-- IFO {ifo} per-template buffers --')
                                if ifo not in new_coincer.singles or ifo not in old_coincer.singles:
                                    print('  singles buffer missing for ifo in one of the coincers')
                                    continue
                                for temp in range(self.num_templates):
                                    try:
                                        newbuf = new_coincer.singles[ifo].data(temp)
                                    except Exception as e:
                                        newbuf = None
                                    try:
                                        oldbuf = old_coincer.singles[ifo].data(temp)
                                    except Exception as e:
                                        oldbuf = None

                                    if newbuf is None or oldbuf is None:
                                        if newbuf is None and oldbuf is None:
                                            continue
                                        print(f'  template {temp}: buffer missing in one side (new {newbuf is None}, old {oldbuf is None})')
                                        continue

                                    # Print lengths and compare key fields if present
                                    print(f'  template {temp}: new len {len(newbuf)}, old len {len(oldbuf)}')
                                    for field in getattr(newbuf, 'dtype', {}).names or []:
                                        if field in getattr(oldbuf, 'dtype', {}).names:
                                            try:
                                                nn = newbuf[field]
                                                oo = oldbuf[field]
                                                if len(nn) == 0 and len(oo) == 0:
                                                    continue
                                                close = np.isclose(nn, oo)
                                                if not close.all():
                                                    nc = np.flatnonzero(np.logical_not(close))
                                                    print(f'    field {field} differs at {len(nc)} positions; first indices: {nc[:5]}')
                                                    print('    new:', nn[nc[:5]])
                                                    print('    old:', oo[nc[:5]])
                                                    # Stop after first differing template for brevity
                                                    raise StopIteration
                                            except StopIteration:
                                                break
                                            except Exception as e:
                                                print(f'    comparing field {field} raised:', e)

                            # Compare coinc buffer summary
                            try:
                                print('-- coinc buffer summary --')
                                print(' new coincs len:', len(new_coincer.coincs.data))
                                print(' old coincs len:', len(old_coincer.coincs.data))
                                print(' new coincs first 10:', new_coincer.coincs.data[:10])
                                print(' old coincs first 10:', old_coincer.coincs.data[:10])
                            except Exception as e:
                                print('  failed reading coinc buffers:', e)

                            # If we found differing entries, retrigger the
                            # new coincer with a trace flag set so the
                            # statistic code will print intermediate
                            # binning/lookups for the first differing
                            # index. This is a one-off debug run.
                            try:
                                if 'not_close' in locals() and np.any(not_close):
                                    first_idx = int(np.flatnonzero(not_close)[0])
                                    print('Triggering phasetd trace for index', first_idx)
                                    try:
                                        new_coincer.stat_calculator.trace_idx = first_idx
                                        # Re-run add_singles to exercise the
                                        # statistic code with the trace flag
                                        new_coincer.add_singles(trigs)
                                    except Exception as e:
                                        print('  failed to retrigger new_coincer for trace:', e)
                                    # Also try the non-invasive debug helper that
                                    # returns per-candidate internals for a
                                    # requested global index. We pass a
                                    # zeroed shift vector and zero to_shift
                                    # multipliers which mirrors the no-timeslide
                                    # behaviour used in this test.
                                    try:
                                        shift = np.zeros(len(new_coincer.coincs.data), dtype=np.float64)
                                        to_shift = [0 for _ in new_coincer.stat_calculator.ifos]
                                        print('Calling debug_trace_candidate for NEW stat...')
                                        ddn = new_coincer.stat_calculator.debug_trace_candidate(trigs, shift, to_shift, first_idx)
                                        print(' debug_trace_candidate (new) keys:', list(ddn.keys()))
                                    except Exception as e:
                                        print('  debug_trace_candidate (new) failed:', e)
                                    try:
                                        print('Calling debug_trace_candidate for OLD stat...')
                                        ddo = old_coincer.stat_calculator.debug_trace_candidate(trigs, shift, to_shift, first_idx)
                                        print(' debug_trace_candidate (old) keys:', list(ddo.keys()))
                                    except Exception as e:
                                        print('  debug_trace_candidate (old) failed:', e)
                                    
                                    # If the trace didn't map the global index, try a deterministic
                                    # replay to map the global coinc buffer index to a per-call
                                    # rtype element. This creates fresh coincers, replays the
                                    # same sequence of single-det triggers and scans the
                                    # saved _last_logsignalrate_data entries produced by the
                                    # statistic to locate which per-call rtype element
                                    # produced the global mismatched value.
                                    try:
                                        # Find the first global mismatch index between
                                        # the two coinc buffers (if present)
                                        gdiff = np.flatnonzero(np.logical_not(np.isclose(
                                            new_coincer.coincs.data, old_coincer.coincs.data)))
                                        if len(gdiff):
                                            global_idx = int(gdiff[0])
                                            target_val = float(new_coincer.coincs.data[global_idx])
                                            print('Attempting replay mapping for global index', global_idx, 'target_val', target_val)

                                            # Create fresh replay coincers (deterministic RNG)
                                            replay_new = Coincer.from_cli(self._test_args, self.num_templates,
                                                                          self._test_analysis_chunk, self._test_detectors)
                                            replay_old = OriginalCoincer.from_cli(self._test_args, self.num_templates,
                                                                                 self._test_analysis_chunk, self._test_detectors)

                                            found = False
                                            # iterate same deterministic single-det trig sequence
                                            for j, rtrigs in enumerate(self.new_trigs):
                                                replay_new.add_singles(rtrigs)
                                                replay_old.add_singles(rtrigs)

                                                # Examine saved internals from the new replay stat
                                                sc_new = replay_new.stat_calculator
                                                if not hasattr(sc_new, '_last_logsignalrate_data'):
                                                    continue
                                                for ent in sc_new._last_logsignalrate_data:
                                                    re = ent.get('ref_entry', {})
                                                    rtype = re.get('rate_indices', [])
                                                    if len(rtype) == 0:
                                                        continue
                                                    # For each rtype element, call debug_trace_candidate
                                                    for ridx in rtype:
                                                        try:
                                                            ddn_r = sc_new.debug_trace_candidate(rtrigs, np.zeros(len(replay_new.coincs.data), dtype=np.float64), [0 for _ in sc_new.ifos], int(ridx))
                                                        except Exception:
                                                            continue
                                                        rvn = float(ddn_r.get('rate_value', float('nan')))
                                                        if np.isclose(rvn, target_val):
                                                            print('Mapped global index', global_idx, 'to replay iteration', j, 'ridx', int(ridx))
                                                            print('Replayed NEW rate_value:', rvn)
                                                            # Now get old stat trace for same ridx
                                                            sc_old = replay_old.stat_calculator
                                                            try:
                                                                ddo_r = sc_old.debug_trace_candidate(rtrigs, np.zeros(len(replay_old.coincs.data), dtype=np.float64), [0 for _ in sc_old.ifos], int(ridx))
                                                            except Exception as e:
                                                                print('  old debug_trace_candidate failed during replay mapping:', e)
                                                                ddo_r = None
                                                            print('Replayed OLD rate_value:', None if ddo_r is None else float(ddo_r.get('rate_value', float('nan'))))
                                                            found = True
                                                            break
                                                    if found:
                                                        break
                                                if found:
                                                    break
                                            if not found:
                                                    # Quick scan the currently-saved internals from the original run
                                                    try:
                                                        sc_cur = new_coincer.stat_calculator
                                                        if hasattr(sc_cur, '_last_logsignalrate_data'):
                                                            for ent in sc_cur._last_logsignalrate_data:
                                                                re = ent.get('ref_entry', {})
                                                                vals = re.get('rate_values', None)
                                                                rinds = re.get('rate_indices', None)
                                                                if vals is None or rinds is None:
                                                                    continue
                                                                for pos, rv in enumerate(vals):
                                                                    if np.isclose(rv, target_val, rtol=1e-06, atol=1e-08):
                                                                        mapped = int(rinds[pos])
                                                                        print('Quick-mapped global index', global_idx, 'to saved ref_entry ridx', mapped)
                                                                        try:
                                                                            ddn_q = sc_cur.debug_trace_candidate(trigs, np.zeros(len(new_coincer.coincs.data), dtype=np.float64), [0 for _ in sc_cur.ifos], mapped)
                                                                            print(' quick NEW rate_value:', ddn_q.get('rate_value'))
                                                                        except Exception as e:
                                                                            print(' quick NEW debug_trace_candidate failed:', e)
                                                                        try:
                                                                            sc_old_cur = old_coincer.stat_calculator
                                                                            ddo_q = sc_old_cur.debug_trace_candidate(trigs, np.zeros(len(old_coincer.coincs.data), dtype=np.float64), [0 for _ in sc_old_cur.ifos], mapped)
                                                                            print(' quick OLD rate_value:', ddo_q.get('rate_value'))
                                                                        except Exception as e:
                                                                            print(' quick OLD debug_trace_candidate failed:', e)
                                                                        found = True
                                                                        break
                                                                if found:
                                                                    break
                                                    except Exception as e:
                                                        print(' quick scan of saved internals failed:', e)
                                                    if not found:
                                                        print('Replay mapping did not find a matching rtype for target value')
                                    except Exception as e:
                                        print('  replay mapping failed:', e)
                            except Exception as e:
                                print('  trace trigger check failed:', e)
                            # Also print any saved internals from the last
                            # logsignalrate call so we can inspect the
                            # binned arrays and lookup indices post-mortem.
                            try:
                                if hasattr(new_coincer.stat_calculator, '_last_logsignalrate_data'):
                                    print('Saved logsignalrate internals (new):', len(new_coincer.stat_calculator._last_logsignalrate_data))
                                    for ent in new_coincer.stat_calculator._last_logsignalrate_data[-1:]:
                                        re = ent.get('ref_entry', {})
                                        print('  ref_ifo:', re.get('ref_ifo'))
                                        r = re.get('rtype', [])
                                        print('  rtype length:', len(r))
                                        if len(r):
                                            print('  rtype min/max:', (int(r.min()), int(r.max())))
                                        b = re.get('binned', [])
                                        if b:
                                            tb, pb, sb = b[0]
                                            print('  sample tbin[:5]:', tb[:5])
                                            print('  sample pbin[:5]:', pb[:5])
                                            print('  sample sbin[:5]:', sb[:5])
                                        print('  loc (sample):', None if re.get('loc') is None else re.get('loc')[:5])
                                if hasattr(old_coincer.stat_calculator, '_last_logsignalrate_data'):
                                    print('Saved logsignalrate internals (old):', len(old_coincer.stat_calculator._last_logsignalrate_data))
                                    # Try a deeper per-position comparison across saved internals
                                    try:
                                        print('Attempting per-position rate comparison across saved internals...')
                                        shift = np.zeros(len(new_coincer.coincs.data), dtype=np.float64)
                                        to_shift = [0 for _ in new_coincer.stat_calculator.ifos]
                                        found = False
                                        for ent in new_coincer.stat_calculator._last_logsignalrate_data:
                                            re = ent.get('ref_entry', {})
                                            rtype = re.get('rtype', [])
                                            if len(rtype) == 0:
                                                continue
                                            for ridx in rtype:
                                                try:
                                                    # ridx is a per-call rtype index; try tracing it on both stats
                                                    dn = new_coincer.stat_calculator.debug_trace_candidate(trigs, shift, to_shift, int(ridx))
                                                except Exception as e:
                                                    # skip if cannot trace
                                                    continue
                                                try:
                                                    do = old_coincer.stat_calculator.debug_trace_candidate(trigs, shift, to_shift, int(ridx))
                                                except Exception as e:
                                                    continue
                                                # Compare rate values if present
                                                if 'rate_value' in dn and 'rate_value' in do:
                                                    rvn = float(dn['rate_value'])
                                                    rvo = float(do['rate_value'])
                                                    if not np.isclose(rvn, rvo):
                                                        print('Per-position mismatch found for rtype element', int(ridx))
                                                        print('  new rate_value:', rvn)
                                                        print('  old rate_value:', rvo)
                                                        print('  ref_ifo (new):', dn.get('ref_ifo'))
                                                        print('  pos_in_rtype (new):', dn.get('pos_in_rtype'))
                                                        print('  nbinned_pos (new):', dn.get('nbinned_pos'))
                                                        print('  loc_pos (new):', dn.get('loc_pos'))
                                                        found = True
                                                        break
                                            if found:
                                                break
                                    except Exception as e:
                                        print('  per-position deep compare failed:', e)
                            except Exception as e:
                                print('  printing saved internals failed:', e)

                            print('==== END DEBUG DUMP ====', flush=True)

                        # Run the debug dump
                        try:
                            debug_state()
                        except Exception as e:
                            print('debug_state failed with exception:', e)

                        self.assertTrue(
                            np.isclose(newout[key], oldout[key]).all()
                        )
                    else:
                        self.assertTrue(newout[key] == oldout[key])

        for i in range(self.num_iterations):
            logging.info("Iteration %d", i)
            single_det_trigs = self.new_trigs[i]
            cres = self.new_coincer.add_singles(single_det_trigs)
            ocres = self.old_coincer.add_singles(single_det_trigs)
            assess_same_output(i, cres, ocres, self.new_coincer, self.old_coincer, single_det_trigs)

        # Are they the same coincs now?
        new_coincer = self.new_coincer
        old_coincer = self.old_coincer
        self.assertTrue(len(new_coincer.coincs.data) == len(old_coincer.coincs.data))
        self.assertTrue(np.isclose(new_coincer.coincs.data, old_coincer.coincs.data, rtol=1e-06).all())

        for ifo in new_coincer.singles:
            lgc = True
            for temp in range(self.num_templates):
                # Check that all singles, for all templates, are identical
                lgc = lgc & (new_coincer.singles[ifo].data(temp) == old_coincer.singles[ifo].data(temp)).all()
            self.assertTrue(lgc)

suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestPyCBCLiveCoinc))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
