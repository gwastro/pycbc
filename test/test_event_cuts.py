import numpy as np
from pycbc.events import cuts
from pycbc.bank import bank_conversions as bank_conv
from utils import simple_exit
import copy
import unittest
import argparse

# Use same seed to make the test reproducible
np.random.seed(1865)

def parse_args(args):
    # Helper function to convert a list of flags/options into
    # an arguments structure
    parser = argparse.ArgumentParser()
    cuts.insert_cuts_option_group(parser)
    return parser, parser.parse_args(args)

bank_size = 100
min_mass = 1
max_mass = 100
min_spin = -0.998
max_spin = 0.998

thresh_mass1 = 60
thresh_mass2 = 40
thresh_spin1 = 0.1
thresh_spin2 = 0.05

mass_range = (max_mass - min_mass)
mass1 = min_mass + np.random.rand(bank_size) * mass_range
mass2 = min_mass + np.random.rand(bank_size) * mass_range
# add a template with an exact number to test the greater than or
# equal comparisons
mass1 = np.concatenate((mass1, [thresh_mass1]))
mass2 = np.concatenate((mass2, [thresh_mass2]))

primary_mass = np.maximum(mass1, mass2)
secondary_mass = np.minimum(mass2, mass2)

spin_range = max_spin - min_spin
spin1 = min_spin + np.random.rand(bank_size) * spin_range
spin2 = min_spin + np.random.rand(bank_size) * spin_range

# add a template with an exact number to test the greater than or
# equal comparisons
spin1 = np.concatenate((spin1, [thresh_spin1]))
spin2 = np.concatenate((spin2, [thresh_spin2]))

bank_size += 1

bank = {
    'mass1': primary_mass,
    'mass2': secondary_mass,
    'spin1z': spin1,
    'spin2z': spin2,
    'f_lower': 15. * np.ones_like(spin1),
}

thresh_bank = {
    'mass1': np.array([thresh_mass1]),
    'mass2': np.array([thresh_mass2]),
    'spin1z': np.array([thresh_spin1]),
    'spin2z': np.array([thresh_spin2]),
    'f_lower': np.array([15.]),
}

class CutsTest(unittest.TestCase):
    def setUp(self):
        self.bank = bank
        self.bank_size = self.bank['mass1'].size

# Tuples of inputs and the name of the test
tests = []
# List which parameters to test
template_params_to_test = bank_conv.conversion_options

for param_name in template_params_to_test:
    # Get the value of the threshold
    threshold = bank_conv.get_bank_property(param_name, thresh_bank, [0])[-1]
    # Add a test that the upper limits work
    tests.append((["--template-cuts", "{}:{}:upper".format(param_name, threshold)],
                 param_name + '_upper'))
    # Add a test that the lower limits work
    tests.append((["--template-cuts", "{}:{}:lower".format(param_name, threshold)],
                 param_name + '_lower'))

for test in tests:
    parser, args = parse_args(test[0])
    def template_cuts_test(self, a=args):
        # Take in the arguments
        _, template_cut_dict = cuts.ingest_cuts_option_group(a)
        # Apply the cuts
        idx_keep = cuts.apply_template_cuts(self.bank, template_cut_dict)

        # For each cut, test that the correct templates pass
        for key_value_cut in a.template_cuts:
            key, value, cut = key_value_cut.split(':')
            value = float(value)
            param = bank_conv.get_bank_property(key, self.bank,
                                                np.arange(self.bank_size))
            if cut == 'upper':
                self.assertTrue(all(param[idx_keep] < value))
            elif cut == 'lower':
                self.assertTrue(all(param[idx_keep] > value))

        # Edit the tests to be inclusive of the cut threshold
        b = copy.deepcopy(a)
        template_cuts_new = []
        for key_value_cut in b.template_cuts:
            # Adding '_inc' means that the threshold is inclusive
            template_cuts_new.append(key_value_cut + '_inc')
        # Reassign the list of cuts in the arguments
        b.template_cuts = template_cuts_new
        # Take in the arguments
        trigger_cut_dict_b, template_cut_dict_b = cuts.ingest_cuts_option_group(b)
        # Apply the cuts
        idx_keep_b = cuts.apply_template_cuts(self.bank, template_cut_dict_b)

        # For each cut, test that the correct templates pass
        for key_value_cut in a.template_cuts:
            key, value, cut = key_value_cut.split(':')
            value = float(value)
            param = bank_conv.get_bank_property(key, self.bank,
                                                np.arange(self.bank_size))
            if cut == 'upper_inc':
                self.assertTrue(all(param[idx_keep_b] <= value))
            elif cut == 'lower_inc':
                self.assertTrue(all(param[idx_keep_b] >= value))

        # The cut threshold is designed to be equal to one template,
        # so the inclusive cut should keep one more than the exclusive
        self.assertTrue(len(idx_keep) == len(idx_keep_b) - 1)

    setattr(CutsTest, "test_cuts_" + test[1], template_cuts_test)

# TODO:
# Add trigger cuts tests (need to work out how to make the triggers dataset/file)
# Add template fit cuts tests    

# create and populate unittest's test suite
suite = unittest.TestSuite()
test_loader = unittest.TestLoader()
suite.addTest(test_loader.loadTestsFromTestCase(CutsTest))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
