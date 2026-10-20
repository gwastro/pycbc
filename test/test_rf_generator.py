"""
Unit tests for the detector-frame waveform generators when the radiation
frame (RF) is used, alone or alongside other detectors.
"""

import unittest

import numpy
from utils import parse_args_cpu_only, simple_exit

from pycbc.detector import Detector
from pycbc.inference.models.gated_gaussian_noise import BaseGatedGaussian
from pycbc.waveform import generator

parse_args_cpu_only("Detector-frame generators")

# fiducial values to test RF detector below
LOCATION = {"tc": 3.1, "ra": 1.37, "dec": -1.26, "polarization": 2.76}
STATIC = {
    "mass1": 38.6,
    "mass2": 29.3,
    "inclination": 0.4,
    "coa_phase": 0.3,
    "distance": 400.0,
    "f_lower": 20.0,
    "delta_f": 1.0 / 8,
    "approximant": "IMRPhenomD",
}


class TestRFDetFrameGenerator(unittest.TestCase):
    """Check that the det frame generators can handle a radiation-frame (RF)
    detector, and properly does/does not require orientation params.
    """

    def _generator(self, cls, detectors, location=LOCATION, **kwargs):
        static = STATIC.copy()
        static.update(kwargs)
        static.update(location)
        return cls(generator.FDomainCBCGenerator, 0.0, detectors=detectors, **static)

    def _assert_close(self, a, b):
        numpy.testing.assert_allclose(
            a.numpy(), b.numpy(), rtol=1e-12, atol=1e-12 * abs(b).max()
        )

    def test_rf_only_needs_no_location(self):
        """An RF-only generator does not require sky location parameters"""
        for dets in [None, ["RF"]]:
            gen = self._generator(generator.FDomainDetFrameGenerator, dets, location={})
            self.assertEqual(gen.detector_names, ["RF"])
            self.assertFalse(gen.has_response)
            self.assertIn("RF", gen.generate())

    def test_mixed_needs_location(self):
        """Real detectors alongside RF still require location parameters"""
        with self.assertRaises(ValueError):
            self._generator(
                generator.FDomainDetFrameGenerator, ["H1", "RF"], location={"tc": 3.1}
            )

    def test_mixed_one_pol(self):
        """Check that DetFrame correctly skips location params only on an
        RF detector instance; the behaviour should be unchanged for regular
        detectors."""
        # combo of known detectors and rf
        mixed = self._generator(
            generator.FDomainDetFrameGenerator, ["H1", "L1", "RF"]
        ).generate()
        # just known detectors
        dets = self._generator(
            generator.FDomainDetFrameGenerator, ["H1", "L1"]
        ).generate()
        # just the rf detector
        rf = self._generator(generator.FDomainDetFrameGenerator, ["RF"]).generate()
        # check wf keys are read in properly
        self.assertEqual(sorted(mixed), ["H1", "L1", "RF"])
        # no RF waveform is returned unless RF is requested
        self.assertEqual(sorted(dets), ["H1", "L1"])
        # waveforms for known detectors should be unchanged with RF...
        for det in ["H1", "L1"]:
            self._assert_close(mixed[det], dets[det])
        # ...and vice versa; RF should not depend on other detectors
        self._assert_close(mixed["RF"], rf["RF"])
        # the RF waveform is the plus polarization, which differs from the
        # projected waveforms with non-trivial location params
        diff = abs(mixed["RF"] - mixed["H1"]).max()
        self.assertGreater(diff, 0.1 * abs(mixed["RF"]).max())

    def test_rf_tc_ref_frame(self):
        """If tc is given in a detector reference frame, check that specifying
        RF correctly translates it back to the geocentric frame.
        """
        # known det plus RF
        mixed = self._generator(
            generator.FDomainDetFrameGenerator, ["H1", "RF"], tc_ref_frame="H1"
        ).generate()
        # just one known det
        h1 = self._generator(
            generator.FDomainDetFrameGenerator, ["H1"], tc_ref_frame="H1"
        ).generate()
        # convert tc (in the H1 frame) to geocenter (the RF frame)
        geotc = LOCATION["tc"] - Detector("H1").time_delay_from_earth_center(
            LOCATION["ra"], LOCATION["dec"], LOCATION["tc"]
        )
        # just the RF det, with tc converted to geocenter
        rf = self._generator(
            generator.FDomainDetFrameGenerator, ["RF"], location={"tc": geotc}
        ).generate()
        # check that the waveforms are close between models, i.e. the shift
        # was applied properly
        self._assert_close(mixed["H1"], h1["H1"])
        self._assert_close(mixed["RF"], rf["RF"])
        # check that the shift actually changed the waveform
        unshifted = self._generator(
            generator.FDomainDetFrameGenerator, ["RF"]
        ).generate()
        diff = abs(mixed["RF"] - unshifted["RF"]).max()
        self.assertGreater(diff, 0.1 * abs(unshifted["RF"]).max())

    def test_mixed_two_pol(self):
        """Check that the TwoPol generator generates the same wfs regardless
        of if RF is included.
        """
        # this generator does not use the polarization
        loc = {k: LOCATION[k] for k in ["tc", "ra", "dec"]}

        # known det plus RF
        mixed = self._generator(
            generator.FDomainDetFrameTwoPolGenerator, ["H1", "RF"], location=loc
        ).generate()
        # just one known det
        dets = self._generator(
            generator.FDomainDetFrameTwoPolGenerator, ["H1"], location=loc
        ).generate()
        # just the RF det
        rf = self._generator(
            generator.FDomainDetFrameTwoPolGenerator, ["RF"], location=loc
        ).generate()
        # check that polarizations match with and without extra dets
        for ii in range(2):
            self._assert_close(mixed["H1"][ii], dets["H1"][ii])
            self._assert_close(mixed["RF"][ii], rf["RF"][ii])

    def test_mixed_modes(self):
        """Check that the Modes generator produces the same waveforms whether
        or not RF is specified.
        """

        def gen(dets):
            return generator.FDomainDetFrameModesGenerator(
                generator.FDomainCBCModesGenerator,
                0.0,
                detectors=dets,
                **dict(
                    STATIC,
                    approximant="IMRPhenomXHM",
                    **{k: LOCATION[k] for k in ["tc", "ra", "dec"]},
                ),
            ).generate()

        # known det plus RF
        mixed = gen(["H1", "RF"])
        # just one known det
        dets = gen(["H1"])
        # just the RF det
        rf = gen(["RF"])
        # check that the modes match across models
        for mode in dets["H1"]:
            for ii in range(2):
                self._assert_close(mixed["H1"][mode][ii], dets["H1"][mode][ii])
                self._assert_close(mixed["RF"][mode][ii], rf["RF"][mode][ii])


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestRFDetFrameGenerator))

if __name__ == "__main__":
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
