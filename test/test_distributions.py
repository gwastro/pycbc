# Copyright (C) 2017 Christopher M. Biwer
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

"""
These are the unittests for distributions in the pycbc.distribtions subpackage.
"""
import unittest
import numpy
import os
from pycbc import distributions
from pycbc.inference import option_utils
from utils import parse_args_cpu_only
from utils import simple_exit

# some distributions are not checked in this unit test
EXCLUDE_DIST_NAMES = ["fromfile", "arbitrary"]

# tests only need to happen on the CPU
parse_args_cpu_only("Distributions")

class TestDistributions(unittest.TestCase):

    def setUp(self, *args):

        # set random seed
        numpy.random.seed(1024)

        # path to example configuration file for testing
        config_path = "/".join([os.path.dirname(os.path.realpath(__file__)),
                                "../examples/distributions/example.ini"])

        # get a set of simulated command line options for
        # configuration file reading
        class Arguments(object):
            config_overrides = []
            config_files = [config_path]
        self.opts = Arguments()

        # read configuration files
        self.cp = option_utils.config_parser_from_cli(self.opts)
        args = option_utils.read_args_from_config(self.cp)
        self.variable_args, self.static_args, self.contraints = args

        # read distributions
        self.dists = distributions.read_distributions_from_config(self.cp)

        # check that all distriubtions will be tested
        for dname, dclass in distributions.distribs.iteritems():
            if (not numpy.any([isinstance(dist, dclass)
                               for dist in self.dists])
                and dname not in EXCLUDE_DIST_NAMES):
                raise ValueError("There is no test for {}".format(dname))

    def test_pdf_rvs(self):

        # 
        for dist in self.dists:
            pass

suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestDistributions))

if __name__ == "__main__":
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
