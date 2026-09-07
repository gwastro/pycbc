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
from pycbc.types import zeros, complex64
from pycbc.filter import overlap
from pycbc.waveform import get_fd_waveform, get_waveform_filter
from pycbc.waveform.spa_tmplt import findchirp_chirptime
from pycbc import pnutils
from utils import parse_args_all_schemes, simple_exit

import lal
import lalsimulation

_scheme, _context = parse_args_all_schemes("Waveform")

class TestSPAtmplt(unittest.TestCase):
    def setUp(self,*args):
        self.context = _context
        self.scheme = _scheme

    def test_spatmplt(self):
        fl = 25
        delta_f = 1.0 / 256

        for m1 in [1, 1.4, 20]:
            for m2 in [1.4, 20]:
                for s1 in  [-2, -1, -0.5, 0, 0.5, 1, 2]:
                    for s2 in [-2, -1, -0.5, 0, 0.5, 1, 2]:
                        # Generate TaylorF2 from lalsimulation, restricting to the capabilities of spatmplt
                        hpr,_ = get_fd_waveform( mass1=m1, mass2=m2, spin1z=s1, spin2z=s2,
                                                 delta_f=delta_f, f_lower=fl,
                                                 approximant="TaylorF2", amplitude_order=0,
                                                 spin_order=-1, phase_order=-1)
                        hpr=hpr.astype(complex64)

                        with self.context:
                            # Generate the spatmplt waveform
                            out = zeros(len(hpr), dtype=complex64)
                            hp = get_waveform_filter(out, mass1=m1, mass2=m2, spin1z=s1, spin2z=s2,
                                                     delta_f=delta_f, f_lower=fl, approximant="SPAtmplt",
                                                     amplitude_order=0, spin_order=-1, phase_order=-1)

                            # Check the diff is sane
                            mag = abs(hpr).sum()
                            diff = abs(hp - hpr).sum() / mag
                            self.assertTrue(diff < 0.01)

                            # Point to point overlap (no phase or time maximization)
                            o =  overlap(hp, hpr)
                            self.assertAlmostEqual(1.0, o, places=4)

                            print("checked m1: %s m2:: %s s1z: %s s2z: %s] overlap = %s, diff = %s" % (m1, m2, s1, s2, o, diff))


class TestChirpTime(unittest.TestCase):
    """Tests for ``findchirp_chirptime``, which derives its PN coefficients
    from LAL's TaylorF2 aligned-spin phasing.
    """

    def test_vs_taylorf2_phase(self):
        """The chirp time is, by construction, t_coalescence - t(f) with
        t(f) = (1 / 2 pi) dPsi/df in the stationary phase approximation.
        Cross-check it against the derivative of the actual LAL TaylorF2
        frequency-domain phase: the time between two in-band frequencies
        t(f1) - t(f2) must match the difference of the chirp times.
        """
        delta_f = 1.0 / 256
        f_lower = 15.0
        f1, f2 = 30.0, 55.0
        cases = [
            (1.2, 1.3, 0.0, 0.0),
            (1.4, 1.4, 0.7, 0.3),
            (1.4, 1.4, -0.7, -0.3),
            (2.5, 1.4, 0.4, 0.0),
            (10.0, 10.0, 0.0, 0.0),
            (10.0, 1.4, 0.6, 0.0),
            (20.0, 20.0, 0.5, 0.5),
            (30.0, 20.0, -0.6, 0.2),
            (5.0, 5.0, -0.8, -0.8),
        ]
        for m1, m2, s1z, s2z in cases:
            hp, _ = get_fd_waveform(
                approximant="TaylorF2", mass1=m1, mass2=m2,
                spin1z=s1z, spin2z=s2z, delta_f=delta_f, f_lower=f_lower,
                f_final=512.0, phase_order=-1, spin_order=-1)
            freq = hp.sample_frequencies.numpy()
            phase = numpy.unwrap(numpy.angle(hp.numpy()))

            def t_of_f(f0, half=1.0):
                sel = (freq > f0 - half) & (freq < f0 + half)
                # local cubic fit removes the (strong) phase curvature bias
                coeffs = numpy.polyfit(freq[sel] - f0, phase[sel], 3)
                return coeffs[2] / (2 * numpy.pi)

            dt_phase = -(t_of_f(f2) - t_of_f(f1))
            dt_formula = (findchirp_chirptime(m1, m2, f1, 7, s1z=s1z, s2z=s2z)
                          - findchirp_chirptime(m1, m2, f2, 7, s1z=s1z, s2z=s2z))
            self.assertAlmostEqual(dt_formula / dt_phase, 1.0, places=3,
                                   msg=f"m1={m1} m2={m2} s1z={s1z} s2z={s2z}")

    def test_roughly_matches_reduced_spin(self):
        """For light-to-moderate mass, moderately spinning systems the
        aligned-spin chirp time should roughly agree with LAL's
        reduced-spin TaylorF2 chirp time.
        """
        f_lower = 20.0
        for m1 in [1.2, 1.4, 3.0, 6.0, 10.0]:
            for m2 in [1.3, 1.4, 5.0, 10.0]:
                if m1 + m2 > 25.0:
                    continue
                for s1z in [-0.7, -0.3, 0.0, 0.3, 0.7]:
                    for s2z in [-0.7, 0.0, 0.7]:
                        chi = lalsimulation.SimInspiralTaylorF2ReducedSpinComputeChi(
                            m1, m2, s1z, s2z)
                        red = lalsimulation.SimInspiralTaylorF2ReducedSpinChirpTime(
                            f_lower, m1 * lal.MSUN_SI, m2 * lal.MSUN_SI, chi, 7)
                        mine = findchirp_chirptime(
                            m1, m2, f_lower, 7, s1z=s1z, s2z=s2z)
                        self.assertAlmostEqual(
                            mine / red, 1.0, delta=0.02,
                            msg=f"m1={m1} m2={m2} s1z={s1z} s2z={s2z}")

    def test_nonspinning_unchanged(self):
        """Spin defaults to zero, and the non-spinning result is the plain
        TaylorF2 chirp time (regression guard against accidental spin leakage).
        """
        for m1, m2 in [(1.4, 1.4), (10.0, 1.4), (25.0, 25.0)]:
            no_spin = findchirp_chirptime(m1, m2, 20.0, 7)
            explicit_zero = findchirp_chirptime(m1, m2, 20.0, 7, s1z=0.0, s2z=0.0)
            self.assertEqual(no_spin, explicit_zero)
            # aligned spin lengthens, anti-aligned shortens the inspiral
            self.assertGreater(
                findchirp_chirptime(m1, m2, 20.0, 7, s1z=0.6, s2z=0.6), no_spin)
            self.assertLess(
                findchirp_chirptime(m1, m2, 20.0, 7, s1z=-0.6, s2z=-0.6), no_spin)

    def test_get_inspiral_tf_uses_spin(self):
        """pnutils.get_inspiral_tf passes spins through to the TaylorF2
        time-frequency track.
        """
        kwargs = dict(tc=0.0, mass1=15.0, mass2=1.4, f_low=20.0,
                      approximant="TaylorF2")
        t_aligned, _ = pnutils.get_inspiral_tf(spin1=0.8, spin2=0.0, **kwargs)
        t_anti, _ = pnutils.get_inspiral_tf(spin1=-0.8, spin2=0.0, **kwargs)
        t_nospin, _ = pnutils.get_inspiral_tf(spin1=0.0, spin2=0.0, **kwargs)
        # at the lowest tracked frequency the aligned-spin system has a longer
        # inspiral, so it starts further before the coalescence time
        self.assertLess(t_aligned[0], t_nospin[0])
        self.assertGreater(t_anti[0], t_nospin[0])


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestSPAtmplt))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestChirpTime))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
