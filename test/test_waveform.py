# Copyright (C) 2012  Alex Nitz, Josh Willis
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

#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
"""
These are the unittests for the pycbc.waveform module
"""
import unittest
import numpy
from numpy import sqrt, cos, sin
from pycbc.scheme import CPUScheme
from pycbc.waveform import (get_td_waveform, get_fd_waveform,
                            get_fd_waveform_sequence)
from pycbc.waveform import generator
from utils import parse_args_all_schemes, simple_exit
from pycbc.types import Array

_scheme, _context = parse_args_all_schemes("Waveform")

# We only check a few as some require auxiliary files
good_waveforms = ['IMRPhenomD', 'TaylorF2', 'SpinTaylorT5',
                  'IMRPhenomPv2', 'IMRPhenomPv3HM',
                  'IMRPhenomPv3']

class TestWaveform(unittest.TestCase):
    def setUp(self,*args):
        self.context = _context
        self.scheme = _scheme

    def test_generation(self):
        with self.context:
            for waveform in good_waveforms:
                print(waveform)
                hc,hp = get_td_waveform(approximant=waveform,mass1=20,mass2=20,delta_t=1.0/4096,f_lower=40)
                self.assertTrue(len(hc)> 0)
            for waveform in good_waveforms:
                print(waveform)
                htilde, g = get_fd_waveform(approximant=waveform,mass1=20,mass2=20,delta_f=1.0/256,f_lower=40)
                self.assertTrue(len(htilde)> 0)

    def test_frequency_sequence(self):
        sample_points = numpy.geomspace(10, 400, 50)
        hp, hc = get_fd_waveform_sequence(approximant="IMRPhenomXAS", mass1=20, mass2=20, sample_points=sample_points)
        hp_ref, hc_ref = get_fd_waveform_sequence(approximant="IMRPhenomXAS", mass1=20, mass2=20, sample_points=Array(sample_points))
        self.assertEqual(hp, hp_ref)
        self.assertEqual(hc, hc_ref)

    def test_spintaylorf2GPU(self):

        print(type(self.context))
        if isinstance(self.context, CPUScheme):
            return

        fl = 25
        delta_f = 1.0 / 256

        for m1 in [3, 5, 15]:
               for m2 in [1., 2., 3.]:
                   for s1 in [0.001, 1.0, 10]:
                       for s1Ctheta in [-1.,0.,0.5,1.]:
                           for s1phi in [0,2.09,4.18]:
                               for inclination in [0.2,1.2]:
                                   s1x = s1 * sqrt(1-s1Ctheta**2) * cos(s1phi)
                                   s1y = s1 * sqrt(1-s1Ctheta**2) * sin(s1phi)
                                   s1z = s1 * s1Ctheta
                                   # Generate SpinTaylorF2 from lalsimulation
                                   hpLAL,hcLAL = get_fd_waveform( mass1=m1, mass2=m2, spin1x=s1x, spin1y=s1y,spin1z=s1z, delta_f=delta_f, f_lower=fl,approximant="SpinTaylorF2", amplitude_order=0, phase_order=7, inclination=inclination )

                                   #Generate SpinTaylorF2 from SpinTaylorF2.py
                                   with self.context:
                                        hp,hc = get_fd_waveform( mass1=m1, mass2=m2, spin1x=s1x, spin1y=s1y,spin1z=s1z, delta_f=delta_f, f_lower=fl,approximant="SpinTaylorF2", amplitude_order=0, phase_order=7, inclination=inclination )

                                   o =  overlap(hpLAL, hp)
                                   self.assertAlmostEqual(1.0, o, places=4)
                                   o =  overlap(hcLAL, hc)
                                   self.assertAlmostEqual(1.0, o, places=4)

                                   ampPLAL=numpy.abs(hpLAL.data)
                                   ampP=numpy.abs(hp.data)
                                   phasePLAL=numpy.unwrap(numpy.angle(hpLAL.data))
                                   phaseP=numpy.unwrap(numpy.angle(hp.data))
                                   ampCLAL=numpy.abs(hcLAL.data)
                                   ampC=numpy.abs(hc.data)
                                   phaseCLAL=numpy.unwrap(numpy.angle(hcLAL.data))
                                   phaseC=numpy.unwrap(numpy.angle(hc.data))
                                   indexampP=numpy.where( ampPLAL!= 0)
                                   indexphaseP=numpy.where( phasePLAL!= 0)
                                   indexampC=numpy.where( ampCLAL!= 0)
                                   indexphaseC=numpy.where( phaseCLAL!= 0)
                                   AmpDiffP = max(abs ( (ampP[indexampP]-ampPLAL[indexampP]) / ampPLAL[indexampP] ) )
                                   PhaseDiffP = max(abs ( (phaseP[indexphaseP] - phasePLAL[indexphaseP]) / phasePLAL[indexphaseP] ) )
                                   AmpDiffC = max(abs ( (ampC[indexampP]-ampCLAL[indexampP]) / ampCLAL[indexampP] ) )
                                   PhaseDiffC = max(abs ( (phaseC[indexphaseP] - phaseCLAL[indexphaseP]) / phaseCLAL[indexphaseP] ) )
                                   self.assertTrue(AmpDiffP < 0.00001)
                                   self.assertTrue(PhaseDiffP < 0.00001)
                                   self.assertTrue(AmpDiffC < 0.00001)
                                   self.assertTrue(PhaseDiffC < 0.00001)
                                   print("..checked m1: %s m2:: %s s1x: %s s1y: %s s1z: %s Inclination: %s" % (m1, m2, s1x, s1y, s1z, inclination))

    def test_errors(self):
        func = get_fd_waveform
        self.assertRaises(ValueError,func,approximant="BLAH")
        self.assertRaises(ValueError,func,approximant="SpinTaylorF2",mass1=3)
        self.assertRaises(ValueError,func,approximant="SpinTaylorF2",mass1=3,mass2=3)
        self.assertRaises(ValueError,func,approximant="SpinTaylorF2",mass1=3,mass2=3,phase_order=7)
        self.assertRaises(ValueError,func,approximant="SpinTaylorF2",mass1=3,mass2=3,phase_order=7)
        self.assertRaises(ValueError,func,approximant="SpinTaylorF2",mass1=3)

        func = get_fd_waveform
        self.assertRaises(ValueError,func,approximant="BLAH")
        self.assertRaises(ValueError,func,approximant="TaylorF2",mass1=3)
        self.assertRaises(ValueError,func,approximant="TaylorF2",mass1=3,mass2=3)
        self.assertRaises(ValueError,func,approximant="TaylorF2",mass1=3,mass2=3,phase_order=7)
        self.assertRaises(ValueError,func,approximant="TaylorF2",mass1=3,mass2=3,phase_order=7)
        self.assertRaises(ValueError,func,approximant="TaylorF2",mass1=3)

        for func in [get_fd_waveform,get_td_waveform]:
            self.assertRaises(ValueError,func,approximant="BLAH")
            self.assertRaises(ValueError,func,approximant="IMRPhenomB",mass1=3)
            self.assertRaises(ValueError,func,approximant="IMRPhenomB",mass1=3,mass2=3)
            self.assertRaises(ValueError,func,approximant="IMRPhenomB",mass1=3,mass2=3,phase_order=7)
            self.assertRaises(ValueError,func,approximant="IMRPhenomB",mass1=3,mass2=3,phase_order=7)
            self.assertRaises(ValueError,func,approximant="IMRPhenomB",mass1=3)


# fiducial values to test RF detector below
LOCATION = {'tc': 3.1, 'ra': 1.37, 'dec': -1.26, 'polarization': 2.76}
STATIC = {'mass1': 38.6, 'mass2': 29.3, 'inclination': 0.4,
          'coa_phase': 0.3, 'distance': 400., 'f_lower': 20.,
          'delta_f': 1./8, 'approximant': 'IMRPhenomD'}


class TestRFDetFrameGenerator(unittest.TestCase):
    """Check that the det frame generators can handle a radiation-frame (RF)
    detector, and properly does/does not require orientation params.
    """
    def _generator(self, cls, detectors, location=LOCATION, **kwargs):
        static = STATIC.copy()
        static.update(kwargs)
        static.update(location)
        return cls(generator.FDomainCBCGenerator, 0., detectors=detectors,
                   **static)

    def _assert_close(self, a, b):
        numpy.testing.assert_allclose(a.numpy(), b.numpy(), rtol=1e-12,
                                      atol=1e-12 * abs(b).max())

    def test_rf_only_needs_no_location(self):
        """An RF-only generator does not require sky location parameters"""
        for dets in [None, ['RF']]:
            gen = self._generator(generator.FDomainDetFrameGenerator, dets,
                                  location={})
            self.assertEqual(gen.detector_names, ['RF'])
            self.assertFalse(gen.has_response)
            self.assertIn('RF', gen.generate())

    def test_mixed_needs_location(self):
        """Real detectors alongside RF still require location parameters"""
        with self.assertRaises(ValueError):
            self._generator(generator.FDomainDetFrameGenerator,
                            ['H1', 'RF'], location={'tc': 3.1})

    def test_mixed_one_pol(self):
        """Check that DetFrame correctly skips location params only on an
        RF detector instance; the behaviour should be unchanged for regular
        detectors."""
        # combo of known detectors and rf
        mixed = self._generator(generator.FDomainDetFrameGenerator,
                                ['H1', 'L1', 'RF']).generate()
        # just known detectors
        dets = self._generator(generator.FDomainDetFrameGenerator,
                               ['H1', 'L1']).generate()
        # just the rf detector
        rf = self._generator(generator.FDomainDetFrameGenerator,
                             ['RF']).generate()
        # check wf keys are read in properly
        self.assertEqual(sorted(mixed), ['H1', 'L1', 'RF'])
        # no RF waveform is returned unless RF is requested
        self.assertEqual(sorted(dets), ['H1', 'L1'])
        # waveforms for known detectors should be unchanged with RF...
        for det in ['H1', 'L1']:
            self._assert_close(mixed[det], dets[det])
        # ...and vice versa; RF should not depend on other detectors
        self._assert_close(mixed['RF'], rf['RF'])
        # the RF waveform is the plus polarization, which differs from the
        # projected waveforms with non-trivial location params
        diff = abs(mixed['RF'] - mixed['H1']).max()
        self.assertGreater(diff, 0.1 * abs(mixed['RF']).max())

    def test_rf_tc_ref_frame(self):
        """If tc is given in a detector reference frame, check that specifying
        RF correctly translates it back to the geocentric frame.
        """
        from pycbc.detector import Detector
        # known det plus RF
        mixed = self._generator(generator.FDomainDetFrameGenerator,
                                ['H1', 'RF'], tc_ref_frame='H1').generate()
        # just one known det
        h1 = self._generator(generator.FDomainDetFrameGenerator,
                             ['H1'], tc_ref_frame='H1').generate()
        # convert tc (in the H1 frame) to geocenter (the RF frame)
        geotc = LOCATION['tc'] - Detector('H1').time_delay_from_earth_center(
            LOCATION['ra'], LOCATION['dec'], LOCATION['tc'])
        # just the RF det, with tc converted to geocenter
        rf = self._generator(generator.FDomainDetFrameGenerator, ['RF'],
                             location={'tc': geotc}).generate()
        # check that the waveforms are close between models, i.e. the shift
        # was applied properly
        self._assert_close(mixed['H1'], h1['H1'])
        self._assert_close(mixed['RF'], rf['RF'])
        # check that the shift actually changed the waveform
        unshifted = self._generator(generator.FDomainDetFrameGenerator,
                                    ['RF']).generate()
        diff = abs(mixed['RF'] - unshifted['RF']).max()
        self.assertGreater(diff, 0.1 * abs(unshifted['RF']).max())

    def test_mixed_two_pol(self):
        """Check that the TwoPol generator generates the same wfs regardless
        of if RF is included.
        """
        # this generator does not use the polarization
        loc = {k: LOCATION[k] for k in ['tc', 'ra', 'dec']}

        # known det plus RF
        mixed = self._generator(generator.FDomainDetFrameTwoPolGenerator,
                                ['H1', 'RF'], location=loc).generate()
        # just one known det
        dets = self._generator(generator.FDomainDetFrameTwoPolGenerator,
                               ['H1'], location=loc).generate()
        # just the RF det
        rf = self._generator(generator.FDomainDetFrameTwoPolGenerator,
                             ['RF'], location=loc).generate()
        # check that polarizations match with and without extra dets
        for ii in range(2):
            self._assert_close(mixed['H1'][ii], dets['H1'][ii])
            self._assert_close(mixed['RF'][ii], rf['RF'][ii])

    def test_mixed_modes(self):
        """Check that the Modes generator produces the same waveforms whether
        or not RF is specified.
        """
        def gen(dets):
            return generator.FDomainDetFrameModesGenerator(
                generator.FDomainCBCModesGenerator, 0., detectors=dets,
                **dict(STATIC, approximant='IMRPhenomXHM',
                       **{k: LOCATION[k] for k in ['tc', 'ra', 'dec']})
                ).generate()
        # known det plus RF
        mixed = gen(['H1', 'RF'])
        # just one known det
        dets = gen(['H1'])
        # just the RF det
        rf = gen(['RF'])
        # check that the modes match across models
        for mode in dets['H1']:
            for ii in range(2):
                self._assert_close(mixed['H1'][mode][ii],
                                   dets['H1'][mode][ii])
                self._assert_close(mixed['RF'][mode][ii], rf['RF'][mode][ii])


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestWaveform))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestRFDetFrameGenerator))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
