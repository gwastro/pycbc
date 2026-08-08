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
    "tc" : (1126259462.43, 1526259462.43),
    "ra" : (0.0, 2 * numpy.pi),
    "dec" : (-numpy.pi / 2, numpy.pi / 2),
    "eclipticlongitude" : (0.0, 2 * numpy.pi),
    "eclipticlatitude" : (-numpy.pi / 2, numpy.pi / 2),
    "polarization" : (0.0, 2 * numpy.pi),
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
            if trans.name == 'spherical_to_cartesian':
                # spherical to cartesian requires the cartesian and spherical
                # parameter names to be specified, which we can get from
                # the inputs and outputs
                inv = trans.inverse(*trans._outputs+trans._inputs)
            else:
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

class TestSpaceDetectorTransformsSetOutputs(unittest.TestCase):
    """`SSBToLISA`, `GEOToLISA` and `SSBToGEO` are "inverse-direction"
    subclasses (of `LISAToSSB`, `LISAToGEO`, `GEOToSSB` respectively) that
    fully redefine `__init__` rather than extending the parent's, purely for
    code sharing of `transform`/`inverse_transform` via class-attribute
    assignment. Calling `super().__init__()` from any of these three would
    silently invoke the *parent* class's `__init__` (next in MRO) instead of
    `BaseTransform.__init__`, re-running the parent's default-parameter
    substitution logic and clobbering any custom parameter names just set.
    This checks that custom parameter names survive construction, and that
    `.inputs`/`.outputs` (needed by `from_config`) still get set correctly.
    """
    def test_ssb_to_lisa_custom_names(self):
        t = transforms.SSBToLISA(
            tc_lisa_param='tc_lisa', longitude_lisa_param='lon_lisa',
            latitude_lisa_param='lat_lisa', polarization_lisa_param='pol_lisa')
        self.assertEqual(t.tc_lisa_param, 'tc_lisa')
        self.assertEqual(
            t.outputs, {'tc_lisa', 'lon_lisa', 'lat_lisa', 'pol_lisa'})
        self.assertEqual(
            t.inputs, {'tc', 'eclipticlongitude', 'eclipticlatitude',
                      'polarization'})

    def test_geo_to_lisa_custom_names(self):
        t = transforms.GEOToLISA(tc_lisa_param='tc_lisa')
        self.assertEqual(t.tc_lisa_param, 'tc_lisa')
        self.assertIn('tc_lisa', t.outputs)

    def test_ssb_to_geo_custom_names(self):
        t = transforms.SSBToGEO(tc_geo_param='tc_geo')
        self.assertEqual(t.tc_geo_param, 'tc_geo')
        self.assertIn('tc_geo', t.outputs)

    def test_ssb_to_moon_custom_names(self):
        t = transforms.SSBToMoon(
            tc_moon_param='tc_moon', longitude_moon_param='lon_moon',
            latitude_moon_param='lat_moon', polarization_moon_param='pol_moon')
        self.assertEqual(t.tc_moon_param, 'tc_moon')
        self.assertEqual(
            t.outputs, {'tc_moon', 'lon_moon', 'lat_moon', 'pol_moon'})
        self.assertEqual(
            t.inputs, {'tc', 'eclipticlongitude', 'eclipticlatitude',
                      'polarization'})

    def test_geo_to_moon_custom_names(self):
        t = transforms.GEOToMoon(tc_moon_param='tc_moon')
        self.assertEqual(t.tc_moon_param, 'tc_moon')
        self.assertIn('tc_moon', t.outputs)

    def test_lisa_to_moon_custom_names(self):
        t = transforms.LISAToMoon(tc_moon_param='tc_moon')
        self.assertEqual(t.tc_moon_param, 'tc_moon')
        self.assertIn('tc_moon', t.outputs)


class TestMoonTransformsRoundTrip(unittest.TestCase):
    """Functional round-trip checks for the Moon<->SSB/GEO/LISA
    transforms, mirroring `TestTransforms.test_inverse`'s pattern for the
    pre-existing space-detector transforms.
    """
    def setUp(self):
        numpy.random.seed(1024)
        self.in_map = {
            'tc': numpy.random.uniform(*RANGES['tc']),
            'eclipticlongitude': numpy.random.uniform(
                *RANGES['eclipticlongitude']),
            'eclipticlatitude': numpy.random.uniform(
                *RANGES['eclipticlatitude']),
            'polarization': numpy.random.uniform(*RANGES['polarization']),
        }

    def _check_round_trip(self, forward_cls, inverse_cls):
        forward = forward_cls()
        inverse = inverse_cls()
        intermediate = forward.transform(self.in_map)
        out_map = inverse.transform(intermediate)
        for key in self.in_map:
            self.assertAlmostEqual(
                out_map[key], self.in_map[key], places=2,
                msg=f'{forward_cls.name} -> {inverse_cls.name} does not '
                    f'recover {key}')

    def test_ssb_moon_round_trip(self):
        self._check_round_trip(transforms.SSBToMoon, transforms.MoonToSSB)

    def test_ssb_moon_round_trip_other_direction(self):
        self._check_round_trip(transforms.MoonToSSB, transforms.SSBToMoon)


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestTransforms))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestSpaceDetectorTransformsSetOutputs))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestMoonTransformsRoundTrip))

if __name__ == "__main__":
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
