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
        self.sngl_rank = 'snr'
        self.detectors = ['H1', 'L1']

        # fake single-detector triggers
        self.sngl_trigs = [(d, {'snglstat': np.random.uniform(0.1, 100, size=10)})
                           for d in self.detectors]

        # here we should also prepare some files needed by some stats,
        # like the PTA histograms or fake trigger fits

    def test_stats(self):
        # FIXME it seems this is only necessary on Python 2
        if not hasattr(self, 'subTest'):
            class DummyCM:
                def __init__(self, msg=''):
                    print('Testing ' + msg)
                def __enter__(self):
                    return self
                def __exit__(self, ex_type, ex_value, ex_tb):
                    return False
            self.subTest = DummyCM

        for stat in statistic_dict:
            # FIXME until we can fake the required files,
            # do not test stats that require them
            if 'phasetd' in stat or 'exp_fit' in stat or 'ExpFit' in stat:
                continue

            with self.subTest(msg=stat):
                # instantiate an object for this stat
                instance = statistic_dict[stat](self.sngl_rank,
                                                ifos=self.detectors)

                # try to rank the fake single triggers
                ranks = instance.rank_stat_coinc(self.sngl_trigs,
                                                 None, None, None)

                # basic sanity check on the produced ranks
                self.assertFalse(np.isnan(ranks).any())


suite = unittest.TestSuite()

# it would be nicer to add each stat as a separate test here
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(CoincStatTest))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
