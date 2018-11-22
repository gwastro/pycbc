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

import numpy
import unittest
from pycbc import transforms
from utils import parse_args_cpu_only
from utils import simple_exit

# list of transforms without an inverse function and to ignore
IGNORE = [t.name for t in transforms.common_cbc_transforms
          if t.inverse is None]

# ranges to draw random numbers for each parameter
RANGES = {
    "mass1" : (1.0, 100.0),
    "mass2" : (1.0, 100.0),
    "mchirp" : (1.0, 20.0),
    "q" : (1.0, 10.0),
    "spin1_a" : (0.0, 1.0),
    "spin1_polar" : (0, numpy.pi),
    "spin1_azimuthal" : (0.0, 2 * numpy.pi),
    "spin2_a" : (0.0, 1.0),
    "spin2_polar" : (0, numpy.pi),
    "spin2_azimuthal" : (0.0, 2 * numpy.pi),
    "chi_eff" : (-1.0, 1.0),
    "chi_a" : (0.0, 1.0),
    "chi_p" : (0.0, 1.0),
    "phi_s" : (0.0, 2 * numpy.pi),
    "phi_a" : (0.0, 2 * numpy.pi),
    "xi1" : (0.0, 1.0),
    "xi2" : (0.0, 1.0),
    "chirp_distance" : (2.0, 10.0),
}

# tests only need to happen on the CPU
parse_args_cpu_only("Transforms")

class TestTransforms(unittest.TestCase):

    def setUp(self):

        # set random seed
        numpy.random.seed(1024)

    def test_inverse(self):

        # set threshold how similar values must be
        threshold = 0.001

        # loop over forward CBC transforms
        for trans in transforms.common_cbc_forward_transforms:

            # check if inverse exists
            if trans.name in IGNORE:
                continue
            inv = trans.inverse()

            # generate some random points
            in_map = {p : numpy.random.uniform(*RANGES[p])
                      for p in trans.inputs}

            # transforms to and back from inverse transform
            intermediate_map = trans.transform(in_map)
            out_map = inv.transform(intermediate_map)

            # check that input equals outputs to some threshold
            in_arr = numpy.array([in_map[p] for p in trans.inputs])
            out_arr = numpy.array([out_map[p] for p in trans.inputs])
            if not numpy.all(1.0 - in_arr / out_arr < threshold):
                raise ValueError(
                "Transform {} does not map back to itself.".format(trans.name))

suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestTransforms))

if __name__ == "__main__":
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)

