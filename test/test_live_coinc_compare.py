"""Mock simulation to easily test and profile PyCBC Live's coincidence code."""

import unittest
import os
import copy
import tempfile
from types import SimpleNamespace
import numpy as np
import h5py
import logging
from pycbc.io import get_file
from pycbc.events.coinc import LiveCoincTimeslideBackgroundEstimator as Coincer
from utils import simple_exit
import validation_code.old_coinc as old_coinc

OriginalCoincer = old_coinc.LiveCoincTimeslideBackgroundEstimator

# This seed is chosen because the impelentations agree here.
# They should only differ due to different numerical precission
SEED = int(os.environ.get('PYCBC_LIVE_COINC_SEED', 0))
START_TIME = 1187008882

class SingleDetTrigSimulator:
    """An object that simulates single-detector triggers in the same format
    as produced by the matched-filtering processes of PyCBC Live.
    """
    def __init__(self, num_templates, analysis_chunk, detectors, num_trigs_per_block):
        self.num_templates = num_templates
        self.detectors = detectors
        self.analysis_chunk = analysis_chunk
        self.start_time = START_TIME
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


def add_loud_trigger_pair(trigs, template_id, end_time, snr=50.0):
    """Append one obviously loud, exactly-coincident trigger to each
    detector's trigger dict (in place), guaranteeing a genuine (zerolag)
    H1-L1 coincidence with a very high ranking statistic value.

    Parameters
    ----------
    trigs: dict of dict
        Per-ifo trigger dicts, as produced by SingleDetTrigSimulator.get_trigs.
    template_id: int
        Template index shared by both detectors' extra trigger, so that
        they are eligible to form a coincidence.
    end_time: float
        GPS end time shared by both detectors' extra trigger (zero time
        difference guarantees it lands within the coincidence window).
    snr: float
        SNR to give the extra trigger in both detectors. Kept well above
        the normal simulated range (4.5-10) so the resulting coinc is
        unambiguously the loudest thing around.
    """
    extra = {
        'snr': snr,
        'end_time': end_time,
        'chisq': 0.5,
        'chisq_dof': 10,
        'coa_phase': 0.0,
        'sigmasq': 1.0,
        'template_id': template_id,
        'mass1': 30.0,
        'mass2': 30.0,
    }
    for det in trigs:
        for key, value in trigs[det].items():
            trigs[det][key] = np.append(
                value, np.array([extra[key]], dtype=value.dtype)
            )


class TestPyCBCLiveCoinc(unittest.TestCase):
    def setUp(self, *args):
        np.random.seed(SEED)

        # Uncomment for more verbosity
        # logging.basicConfig(format="%(asctime)s %(message)s",
        #                     level=logging.INFO)

        # simulate the `args` object we normally get from the command line arguments

        url = 'https://github.com/gwastro/pycbc-config/raw/master/'
        url += 'test_data_files/{}-PTA_HISTOGRAM.hdf'
        stat_file_paths = [
            get_file(url.format("H1L1"), cache=True),
        ]
        # kept on self so other tests can build coincers with variations
        # (e.g. a different ifar_remove_threshold)
        self.args = args = SimpleNamespace(
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
            ifar_remove_threshold=None,
        )

        # number of templates in the bank
        self.num_templates = 10

        # duration of analysis segment
        analysis_chunk = 2000
        self.analysis_chunk = analysis_chunk

        # combination of two detectors to analyze
        detectors = ["H1", "L1"]
        self.detectors = detectors

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

        def assess_same_output(newout, oldout):
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

                    a = newout[key]
                    b = oldout[key]

                    if key == 'foreground/stat':
                        self.assertIsInstance(a, np.ndarray)
                        self.assertEqual(a.ndim, 1)
                        self.assertEqual(len(a), 1)
                        self.assert_foreground_stat_hdf_readable(a)
                        self.assertEqual(len(a), len(np.atleast_1d(b)))
                        self.assertTrue(np.isclose(a, np.atleast_1d(b)).all())
                        continue

                    if isinstance(a, np.ndarray):
                        # compare shapes and values
                        self.assertEqual(len(a), len(b))

                        a_comp = a
                        b_comp = b

                        # For background/stat, order by time as the sort is not stable
                        if key == 'background/stat' and len(a) > 1:
                            tnew = newout.get('background/time', None)
                            told = oldout.get('background/time', None)
                            idx_new = np.argsort(tnew, kind='stable')
                            idx_old = np.argsort(told, kind='stable')
                            a_comp = a[idx_new]
                            b_comp = b[idx_old]

                        self.assertTrue(np.isclose(a_comp, b_comp).all())
                    else:
                        self.assertEqual(a,b)

        for i in range(self.num_iterations):
            logging.info("Iteration %d", i)
            single_det_trigs = self.new_trigs[i]
            cres = self.new_coincer.add_singles(single_det_trigs)
            ocres = self.old_coincer.add_singles(single_det_trigs)
            assess_same_output(cres, ocres)

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

    def test_ifar_remove_threshold(self):
        """With a nonzero ifar_remove_threshold, a loud zerolag coincidence
        should mark its analysis chunk as loud: the chunk is then excluded
        from the background time and future background coincidences, while
        the zerolag candidate itself is still found and reported exactly
        as it would be without the threshold set.

        The `old_coinc` validation code does not implement this feature at
        all, so this test compares two instances of the *new* coincer
        (with and without the threshold) against each other, rather than
        against `old_coinc` as in test_coincer_runs.
        """
        threshold = 1.0  # years

        args_thresh = copy.copy(self.args)
        args_thresh.ifar_remove_threshold = threshold
        args_nothresh = copy.copy(self.args)
        args_nothresh.ifar_remove_threshold = None

        coincer_thresh = Coincer.from_cli(
            args_thresh, self.num_templates, self.analysis_chunk,
            self.detectors
        )
        coincer_nothresh = Coincer.from_cli(
            args_nothresh, self.num_templates, self.analysis_chunk,
            self.detectors
        )

        # A few chunks of ordinary noise triggers to establish a
        # background, fed identically to both coincers.
        num_warmup = 4
        for i in range(num_warmup):
            trigs = self.new_trigs[i]
            coincer_thresh.add_singles(copy.deepcopy(trigs))
            coincer_nothresh.add_singles(copy.deepcopy(trigs))

        # Engineer an unambiguous, very loud zerolag coincidence: identical,
        # very high SNR triggers for H1 and L1 in the same template at
        # exactly the same time.
        loud_trigs = copy.deepcopy(self.new_trigs[num_warmup])
        loud_time = loud_trigs['H1']['end_time'][0]
        add_loud_trigger_pair(loud_trigs, template_id=0, end_time=loud_time)

        res_thresh = coincer_thresh.add_singles(copy.deepcopy(loud_trigs))
        res_nothresh = coincer_nothresh.add_singles(copy.deepcopy(loud_trigs))

        # The candidate itself is found and reported the same way whether
        # or not removal is enabled.
        self.assertIn('foreground/ifar', res_thresh)
        self.assertIn('foreground/ifar', res_nothresh)
        self.assertGreater(res_thresh['foreground/ifar'], threshold)

        # Only the thresholded coincer marks its chunk loud.
        expected_chunk = int(loud_time // self.analysis_chunk)
        self.assertEqual(len(coincer_nothresh.loud_chunks), 0)
        self.assertIn(expected_chunk, coincer_thresh.loud_chunks)

        # A loud chunk is excluded from both ifos' contribution to the
        # background time, so it must be strictly smaller than in the
        # unfiltered run, even though both saw identical triggers.
        self.assertLess(coincer_thresh.background_time,
                        coincer_nothresh.background_time)

        # On the very update that creates the loud chunk, the background
        # coincs formed from it are stripped out of the surviving
        # (post-clustering) winners before either coincer's buffer is
        # updated, so the thresholded coincer must not have gained more
        # background coincs than the unfiltered one on this update.
        # (This direct comparison only holds right at the point the chunk
        # is first marked loud: on later updates, excluding loud-chunk
        # coincs *before* clustering can shift which coincs win each
        # cluster, so background counts are no longer simply ordered.)
        self.assertLessEqual(
            len(coincer_thresh.coincs.data), len(coincer_nothresh.coincs.data)
        )

        # Continue for a few more chunks, comfortably inside the lookback
        # window so the loud chunk isn't pruned yet: the exclusion should
        # keep reducing the thresholded coincer's background time.
        last = min(num_warmup + 4, self.num_iterations)
        for i in range(num_warmup + 1, last):
            trigs = self.new_trigs[i]
            coincer_thresh.add_singles(copy.deepcopy(trigs))
            coincer_nothresh.add_singles(copy.deepcopy(trigs))

        self.assertIn(expected_chunk, coincer_thresh.loud_chunks)
        self.assertLess(coincer_thresh.background_time,
                        coincer_nothresh.background_time)

    def test_foreground_stat_hdf_contract(self):
        self.assert_foreground_stat_hdf_readable(np.array([12.5]))

    def assert_foreground_stat_hdf_readable(self, stat):
        """Check live HDF output keeps foreground/stat slice-readable."""
        fd, path = tempfile.mkstemp(suffix='.hdf')
        os.close(fd)
        try:
            with h5py.File(path, 'w') as fp:
                fp['foreground/stat'] = stat
            with h5py.File(path, 'r') as fp:
                saved = fp['foreground/stat'][:]
        finally:
            os.remove(path)
        self.assertEqual(saved.shape, (1,))
        self.assertTrue(np.isclose(saved, stat).all())

suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestPyCBCLiveCoinc))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
