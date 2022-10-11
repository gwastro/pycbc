"""
Unit tests for cuts being applied to trigger and templates
"""

import unittest
import argparse
import copy
import numpy as np
from utils import simple_exit
from pycbc.events import cuts, ranking
from pycbc.tmpltbank import bank_conversions


def parse_args(args):
    # Helper function to convert a list of flags/options into
    # an arguments structure
    parser = argparse.ArgumentParser()

    cuts.insert_cuts_option_group(parser)
    return parser.parse_args(args)


class CutsErrorsTest(unittest.TestCase):
    def setUp(self):
        # Set up some things we will want to use in the tests:
        self.maxDiff = None

# Tuples of inputs for tests which should fail
test_parser_raises_error = []

test_parser_raises_error.append((['--template-cuts',
                                 'this_is_not_right'],
                                 ValueError,
                                 'incorrect_format_template_cut'))

test_parser_raises_error.append((['--trigger-cuts',
                                 'this_is_not_right'],
                                 ValueError,
                                 'incorrect_format_trigger_cut'))

test_parser_raises_error.append((['--trigger-cuts',
                                 'chi_eff:0.9:upper'],
                                 NotImplementedError,
                                 'template_param_given_as_trigger_cut'))

test_parser_raises_error.append((['--template-cuts',
                                 'newsnr:5.5:lower'],
                                 NotImplementedError,
                                 'trigger_param_given_as_template_cut'))

test_parser_raises_error.append((['--trigger-cuts',
                                 'newsnr:5.5:nonsense'],
                                 NotImplementedError,
                                 'nonsense_limit'))

test_parser_raises_error.append((['--trigger-cuts',
                                 'newsnr:notafloat:upper'],
                                 ValueError,
                                 'threshold_not_a_float'))

# Dynamically add error tests into the class
for test_parser_error in test_parser_raises_error:
    args = parse_args(test_parser_error[0])
    def check_sysexit_test(self, a=args, te=test_parser_error[1]):
        with self.assertRaises(te):
            cuts.ingest_cuts_option_group(a)
    setattr(CutsErrorsTest,
            'test_parser_error_' + test_parser_error[2],
            check_sysexit_test)

class CutsParserTest(unittest.TestCase):
    def setUp(self):
        # Set up some things we will want to use in the tests:
        self.maxDiff = None

# Tuples of inputs for tests where we are checking the output
tests_parser_output = []

# No cuts: empty dicts as output
tests_parser_output.append(([], ({}, {}), "no_cuts_given"))

# Multiple cuts on the same parameter / cut type,
# should only use the strictest cut
tests_parser_output.append((['--trigger-cuts',
                             'snr:4:lower', 'snr:5:lower'],
                            ({('snr', np.greater): 5}, {}),
                            "multiple_similar_cuts"))

tests_parser_output.append((['--trigger-cuts',
                             'snr:5:lower', 'snr:4:lower'],
                            ({('snr', np.greater): 5}, {}),
                            "multiple_similar_cuts_2"))

# Multiple cuts exactly the same, should give the warning but still complete
tests_parser_output.append((['--trigger-cuts',
                             'snr:5:lower', 'snr:5:lower'],
                            ({('snr', np.greater): 5}, {}),
                            "multiple_same_cuts"))

# Dynamically add value tests for the parser
for test_values in tests_parser_output:
    args = parse_args(test_values[0])
    def ingest_values_test(self, a=args, tv=test_values[1]):
        cuts_dicts = cuts.ingest_cuts_option_group(a)
        self.assertEqual(cuts_dicts, tv)


    setattr(CutsParserTest,
            'test_parser_values_' + test_values[2],
            ingest_values_test)

# Set up some random datasets to test the cuts:

# Make this reproducible
np.random.seed(1865)

# Set values for triggers:
n_triggers = 10000
trigger_dset = {}
trigger_dset['snr'] = np.random.random(n_triggers) * 10
trigger_dset['chisq'] = np.random.random(n_triggers) * 100
trigger_dset['chisq_dof'] = np.ceil(np.random.random(n_triggers) * 100) + 2
trigger_dset['sg_chisq'] = np.random.random(n_triggers) * 2

# Set values for templates:
template_dset = {}
n_templates = 10000

# For one of the tests, we want to use an exact value of mass1,
# so set it manually
mass1 = list(np.random.random(n_templates - 1) * 100) + [59]
mass2 = list(np.random.random(n_templates - 1) * 100) + [10]

# Assert mass1 > mass2
template_dset['mass1'] = np.maximum(mass1, mass2)
template_dset['mass2'] = np.minimum(mass1, mass2)

# Spins must be in the range -1, 1, (but not exactly 1):
def make_random_spins(n_templates):
    spins = np.random.normal(size=n_templates, scale=0.4)
    spins = np.maximum(spins, -0.998)
    spins = np.minimum(spins, 0.998)
    return spins

template_dset['spin1z'] = make_random_spins(n_templates)
template_dset['spin2z'] = make_random_spins(n_templates)

template_dset['template_duration'] = np.random.random(n_templates) * 100


class CutsTest(unittest.TestCase):
    def setUp(self):
        # Set up some things we will want to use in the tests:
        self.maxDiff = None

# Set up values to be tested:
test_cut_output = []


# Lower limits on a trigger parameter we can read directly
def snr_lower_test(trigs, trig_idx, **kwargs):
    return all(trigs['snr'][trig_idx] > 4)


test_cut_output.append((['--trigger-cuts', 'snr:4:lower'],
                        snr_lower_test,
                        'snr_cut_lower'))


# Upper limits on a derived trigger parameter
def newsnr_sgveto_upper_test(trigs, trig_idx, **kwargs):
    rwsnr = ranking.get_newsnr_sgveto(trigs)
    return all(rwsnr[trig_idx] < 10)

test_cut_output.append((['--trigger-cuts', 'newsnr_sgveto:10:upper'],
                        newsnr_sgveto_upper_test,
                        'nsnr_sgveto_cut_upper'))

# Lower limit on directly-read template bank parameter
def template_duration_test(temps, temp_idx, **kwargs):
    return all(temps['template_duration'][temp_idx] > 10)

test_cut_output.append((['--template-cuts',
                         'template_duration:10:lower'],
                        template_duration_test,
                        'template_duration_lower'))

# make sure the "or equal to" is working properly
def mass1_ge_test(temps, temp_idx, **kwargs):
    return any(temps['mass1'][temp_idx] == 59)

test_cut_output.append((['--template-cuts',
                         'mass1:59:lower_inc'],
                        mass1_ge_test,
                        'mass1_orequal'))

# Upper and lower limits on a derived parameter from the bank
def chi_eff_upper_lower_test(temps, temp_idx, **kwargs):
    chi_eff = bank_conversions.get_bank_property('chi_eff', temps, temp_idx)
    return all(np.logical_and(chi_eff > -0.5,
                              chi_eff <= 0.8))


test_cut_output.append((['--template-cuts', 'chi_eff:0.8:upper',
                         'chi_eff:-0.5:lower_inc'],
                         chi_eff_upper_lower_test,
                         'chi_eff_cut_upper_lower'))


# FIXME: Once we can fake the statistic files, create a test
# for template fits cuts

# Dynamically add value tests for the parser
for test_values in test_cut_output:
    args = parse_args(test_values[0])
    def cut_values_test(self, a=args, test_func=test_values[1]):
        # Copy the global variables, as we need to use local in the tests
        triggers = copy.deepcopy(trigger_dset)
        bank = copy.deepcopy(template_dset)

        # Take in the arguments which define the cuts
        trigger_cut_dict, template_cut_dict = cuts.ingest_cuts_option_group(a)

        # Get out the indices which meet the cut criteria
        templates_idx = cuts.apply_template_cuts(bank, template_cut_dict)
        triggers_idx = cuts.apply_trigger_cuts(triggers, trigger_cut_dict)

        # Using the checks from the definition tuple, is the cut working?
        self.assertTrue(test_func(temps=bank,
                                  temp_idx=templates_idx,
                                  trigs=trigger_dset,
                                  trig_idx=triggers_idx))

    setattr(CutsTest,
            'test_cuts_correct_' + test_values[2],
            cut_values_test)

# create and populate unittest's test suite
suite = unittest.TestSuite()
test_loader = unittest.TestLoader()
suite.addTest(test_loader.loadTestsFromTestCase(CutsErrorsTest))
suite.addTest(test_loader.loadTestsFromTestCase(CutsParserTest))
suite.addTest(test_loader.loadTestsFromTestCase(CutsTest))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
