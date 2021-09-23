"""Unit test for coincident ranking statistic implementations."""

import unittest
import numpy as np
from utils import parse_args_cpu_only, simple_exit
from pycbc.events.stat import statistic_dict


# this test only needs to happen on the CPU
parse_args_cpu_only('coinc stats')

class CoincStatTest(unittest.TestCase):
    def setUp(self):
        # one could loop over all single rankings
        # and detector combinations for a complete test
        self.single_rank = 'snr'
        self.detectors = ['H1', 'L1']

        # simulate some single-detector triggers from Gaussian noise
        self.num_trigs = n = 100
        self.single_trigs = {}
        for d in self.detectors:
            self.single_trigs[d] = {
                'snr': np.random.chisquare(2, size=n) ** 0.5,
                'coa_phase': np.random.uniform(0, 2 * np.pi, size=n),
                'end_time': 1295441120 + np.random.uniform(-0.01, 0.01, size=n),
                'sigmasq': np.random.uniform(1, 10, size=n)
            }

        self.time_slide_shift = 0.1
        self.time_slide_vector = np.zeros(n, dtype=int)

        # here we should also prepare some files needed by some stats,
        # like the PTA histograms or fake trigger fits
        self.aux_files = []


# dynamically insert a test case method for each available statistic
for stat_name in statistic_dict:
    # FIXME until we can fake the required files,
    # do not test stats that require them
    if 'phasetd' in stat_name or 'exp_fit' in stat_name \
            or 'ExpFit' in stat_name:
        continue

    def stat_test_method(self, sn=stat_name):
        # instantiate an object for this stat
        stat = statistic_dict[sn](self.single_rank,
                                  files=self.aux_files,
                                  ifos=self.detectors)

        # get single-detector statistic at each detector from the fake triggers
        single_ranks = [(d, stat.single(self.single_trigs[d]))
                        for d in self.detectors]

        # pretend the fake triggers are coincident, and rank the coincidences
        coinc_ranks = stat.rank_stat_coinc(
            single_ranks,
            self.time_slide_vector,
            self.time_slide_shift,
            (0, 0) # does this make sense?
        )

        # basic sanity check on the produced ranks
        self.assertEqual(len(coinc_ranks), self.num_trigs)
        self.assertFalse(np.isnan(coinc_ranks).any())

        # here one could add a statistical test based on known analytic
        # behavior of a particular statistic in Gaussian noise

    setattr(CoincStatTest, 'test_' + stat_name, stat_test_method)

# create and populate unittest's test suite
suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(CoincStatTest))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
