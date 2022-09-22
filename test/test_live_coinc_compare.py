"""Mock simulation to easily test and profile PyCBC Live's coincidence code."""

import unittest
from types import SimpleNamespace
import numpy as np
# Some duplicate imports, but I want to copy code without changing it!
import numpy, logging, pycbc.pnutils, pycbc.conversions, copy, lal
import cProfile
from astropy.utils.data import download_file
from pycbc import gps_now
from pycbc.events.coinc import LiveCoincTimeslideBackgroundEstimator as Coincer
import pycbc.events.coinc
from utils import simple_exit
import validation_code.old_coinc as old_coinc

OriginalCoincer = old_coinc.LiveCoincTimeslideBackgroundEstimator

class SingleDetTrigSimulator:
    """An object that simulates single-detector triggers in the same format
    as produced by the matched-filtering processes of PyCBC Live.
    """
    def __init__(self, num_templates, analysis_chunk, detectors, num_trigs_per_block):
        self.num_templates = num_templates
        self.detectors = detectors
        self.analysis_chunk = analysis_chunk
        self.start_time = gps_now()
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
                ).astype(np.int32)
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
            timeslide_interval=0.1,
            background_ifar_limit=100,
            store_background=True
        )

        # number of templates in the bank
        self.num_templates = 10

        # duration of analysis segment
        analysis_chunk = 2000

        # combination of two detectors to analyze
        detectors = ["H1", "L1"]

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
                    if type(newout[key]) is np.ndarray:
                        self.assertTrue(len(newout[key]) == len(oldout[key]))
                        self.assertTrue(
                            numpy.isclose(newout[key], oldout[key]).all()
                        )
                    else:
                        self.assertTrue(newout[key] == oldout[key])

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
        self.assertTrue(numpy.isclose(new_coincer.coincs.data, old_coincer.coincs.data, rtol=1e-06).all())

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
