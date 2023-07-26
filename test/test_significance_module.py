"""Unit test for converting sets of statistic values into significances."""

import unittest
import argparse
import itertools
import copy
import numpy as np
from utils import simple_exit
from pycbc.events import significance

def parse_args(args):
    # Helper function to convert a list of flags/options into
    # an arguments structure
    parser = argparse.ArgumentParser()
    significance.insert_significance_option_group(parser)
    return parser, parser.parse_args(args)

ifos = ['H1','L1','V1']
# What combinations of the ifos can we make?
combos = []
for l in np.arange(len(ifos)) + 1:
            combos += [''.join(c)
                       for c in itertools.combinations(ifos, l)]

class SignificanceParserTest(unittest.TestCase):
    def setUp(self):
        # Set up some things we will want to use in the tests:
        self.maxDiff = None
        self.ifos = copy.copy(ifos)
        self.combos = copy.copy(combos)

# Tuples of inputs and the errors they should create
tests_which_sysexit = []

# Try to use a calculation method which doesn't exist
tests_which_sysexit.append((['--far-calculation-method',
                             'H1L1:nonexistent_method'],
                            'method_doesnt_exist'))

# Try to set a fit threshold when using n_louder method
tests_which_sysexit.append((['--fit-threshold',
                             'H1L1:6'],
                            'n_louder_with_threshold'))

# Try to set a fit function when using n_louder method
tests_which_sysexit.append((['--fit-function',
                             'H1L1:exponential'],
                            'function_with_n_louder'))

# Try to set a fit threshold which is not a number
tests_which_sysexit.append((['--fit-threshold',
                             'H1L1:not_a_number'],
                            'threshold_not_a_number'))


# Try to set a fit function which isn't expected
tests_which_sysexit.append((['--far-calculation-method',
                             'H1L1:trigger_fit',
                             '--fit-function',
                             'H1L1:spanish_inquisition'],
                            'function_doesnt_exist'))

# Dynamically add sysexit tests into the class
for test_sysexit in tests_which_sysexit:
    parser, args = parse_args(test_sysexit[0])
    def check_sysexit_test(self, a=args, p=parser):
        with self.assertRaises(SystemExit):
            significance.check_significance_options(
                a, p)
    setattr(SignificanceParserTest,
        'test_parser_' + test_sysexit[1],
        check_sysexit_test)


# Set up the default values of the output dictionary, we will edit this
# for each test
default_dict = {}
# Default Values
for combo in combos:
    default_dict[combo] = copy.deepcopy(significance._default_opt_dict)

tests_which_pass = []

# Does passing no options return the default dictionary?
tests_which_pass.append(([], default_dict, 'default_vals'))

# Try to add a detector combination which does not exist
# - should return dictionary including the nonexistent combination
# as we want to be able to give combos to scripts where they aren't valid
extra_combo_dict = copy.deepcopy(default_dict)
extra_combo_dict['H1G1'] = {}
extra_combo_dict['H1G1']['method'] = 'trigger_fit'
extra_combo_dict['H1G1']['fit_function'] = None
extra_combo_dict['H1G1']['fit_threshold'] = 6.
tests_which_pass.append((['--far-calculation-method',
                          'H1G1:trigger_fit',
                          '--fit-threshold', 'H1G1:6'],
                         extra_combo_dict,
                         'extra_combo'))

# Supply different methods for the different combinations, and check that
# they are taken in properly
test_dict = copy.deepcopy(default_dict)
test_dict['H1L1']['method'] = 'trigger_fit'
test_dict['H1']['method'] = 'trigger_fit'
test_dict['L1']['method'] = 'trigger_fit'
test_dict['H1L1']['fit_function'] = 'power'
test_dict['H1']['fit_function'] = 'rayleigh'
test_dict['L1']['fit_function'] = 'exponential'
test_dict['H1L1']['fit_threshold'] = 6
test_dict['H1']['fit_threshold'] = 5.5
test_dict['L1']['fit_threshold'] = 5

calc_methods = ['H1L1:trigger_fit', 'H1:trigger_fit', 'L1:trigger_fit']
functions = ['H1L1:power', 'H1:rayleigh', 'L1:exponential']
thresholds = ['H1L1:6', 'H1:5.5', 'L1:5']
tests_which_pass.append((['--far-calculation-method'] + calc_methods +
                         ['--fit-function'] + functions + 
                         ['--fit-threshold'] + thresholds,
                         test_dict,
                         'different_combos'))

# Dynamically add value tests for the parser
for test_values in tests_which_pass:
    parser, args = parse_args(test_values[0])
    def digest_values_test(self, a=args, tv=test_values[1]):
        method_dict = significance.digest_significance_options(
            self.combos, a)
        self.assertEqual(method_dict, tv)


    setattr(SignificanceParserTest,
            'test_parser_values_' + test_values[2],
            digest_values_test)

class SignificanceMethodTest(unittest.TestCase):
    def setUp(self):
        self.test_fg_stat = np.random.normal(loc=5,scale=2,size=50)
        self.test_bg_stat = np.random.normal(loc=5,scale=2,size=500)
        self.dec_facs = np.ones_like(self.test_bg_stat)


method_functions = {
    'n_louder': [None],
    'trigger_fit': ['exponential','rayleigh']
}

# Dynamically add method tests into the class
for method in significance._significance_meth_dict:
    for function in method_functions[method]:
        method_dict = {}
        method_dict['method'] = method
        method_dict['fit_function'] = function
        method_dict['fit_threshold'] = None if not function else 0

        def meth_test(self, md=method_dict):
            back_cnum, fnlouder = significance.get_n_louder(
                self.test_bg_stat,
                self.test_fg_stat,
                self.dec_facs,
                **method_dict)

            back_stat_sort = np.argsort(self.test_bg_stat)
            back_far_sort = np.argsort(back_cnum)

            fore_stat_sort = np.argsort(self.test_fg_stat)
            fore_far_sort = np.argsort(fnlouder)

            # Basic sanity check - there should be one n_louder value
            # per stat value
            self.assertEqual(len(back_cnum), len(self.test_bg_stat))
            self.assertEqual(len(fnlouder), len(self.test_fg_stat))

            # None of the output should be NaN or infinite
            self.assertTrue(np.isfinite(back_cnum).all())
            self.assertTrue(np.isfinite(fnlouder).all())

            # The background stat value order should be the reverse of the
            # n_louder order
            back_stat_sort = np.argsort(self.test_bg_stat)
            back_far_sort = np.argsort(back_cnum)
            self.assertTrue(np.array_equal(back_stat_sort,
                                           back_far_sort[::-1]))

            fore_stat_sort = np.argsort(self.test_fg_stat)
            fore_far_sort = np.argsort(fnlouder)
            # As fg events could have an equal number of louder bg events,
            # argsort be the opposite way round for the far sort and stat
            # sort. So we need to use the recovered n_louder as the equal
            # equality test array
            self.assertTrue(np.array_equal(fnlouder[fore_stat_sort],
                                           fnlouder[fore_far_sort][::-1]))

        setattr(SignificanceMethodTest,
                'test_%s_%s' % (method, function),
                meth_test)

# create and populate unittest's test suite
suite = unittest.TestSuite()
test_loader = unittest.TestLoader()
suite.addTest(test_loader.loadTestsFromTestCase(SignificanceMethodTest))
suite.addTest(test_loader.loadTestsFromTestCase(SignificanceParserTest))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
