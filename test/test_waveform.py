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
import sys
import pycbc
import unittest
from pycbc.types import *
from pycbc.scheme import *
from pycbc.filter import *
from pycbc.waveform import *
import pycbc.fft
import numpy
from numpy import sqrt, cos, sin
from utils import parse_args_all_schemes, simple_exit

_scheme, _context = parse_args_all_schemes("Waveform")

class TestWaveform(unittest.TestCase):
    def setUp(self,*args):
        self.context = _context
        self.scheme = _scheme

    def test_generation(self):
        with self.context:
            for waveform in td_approximants():
                print waveform
                hc,hp = get_td_waveform(approximant=waveform,mass1=20,mass2=20,delta_t=1.0/4096,f_lower=40)
                self.assertTrue(len(hc)> 0)
            for waveform in fd_approximants():
                print waveform
                htilde, g = get_fd_waveform(approximant=waveform,mass1=20,mass2=20,delta_f=1.0/256,f_lower=40)
                self.assertTrue(len(htilde)> 0)


    def test_spintaylorf2GPU(self):
    
        print type(self.context)
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
                                   print "..checked m1: %s m2:: %s s1x: %s s1y: %s s1z: %s Inclination: %s" % (m1, m2, s1x, s1y, s1z, inclination)

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


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestWaveform))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
