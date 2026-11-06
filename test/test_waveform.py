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
from numpy import cos, sin, sqrt
from utils import parse_args_all_schemes, simple_exit

from pycbc.filter import overlap
from pycbc.scheme import CPUScheme
from pycbc.types import Array, complex64, zeros
from pycbc.waveform import (
    get_fd_waveform,
    get_fd_waveform_sequence,
    get_td_waveform,
    get_waveform_filter,
)

_scheme, _context = parse_args_all_schemes("Waveform")

# We only check a few as some require auxiliary files
good_waveforms = [
    "IMRPhenomD",
    "TaylorF2",
    "SpinTaylorT5",
    "IMRPhenomPv2",
    "IMRPhenomPv3HM",
    "IMRPhenomPv3",
]


class TestWaveform(unittest.TestCase):
    def setUp(self, *args):
        self.context = _context
        self.scheme = _scheme

    def test_generation(self):
        with self.context:
            for waveform in good_waveforms:
                print(waveform)
                hp, _ = get_td_waveform(
                    approximant=waveform,
                    mass1=20,
                    mass2=20,
                    delta_t=1.0 / 4096,
                    f_lower=40,
                )
                self.assertTrue(len(hp) > 0)
            for waveform in good_waveforms:
                print(waveform)
                htilde, _ = get_fd_waveform(
                    approximant=waveform,
                    mass1=20,
                    mass2=20,
                    delta_f=1.0 / 256,
                    f_lower=40,
                )
                self.assertTrue(len(htilde) > 0)

    def test_frequency_sequence(self):
        sample_points = numpy.geomspace(10, 400, 50)
        hp, hc = get_fd_waveform_sequence(
            approximant="IMRPhenomXAS", mass1=20, mass2=20, sample_points=sample_points
        )
        hp_ref, hc_ref = get_fd_waveform_sequence(
            approximant="IMRPhenomXAS",
            mass1=20,
            mass2=20,
            sample_points=Array(sample_points),
        )
        self.assertEqual(hp, hp_ref)
        self.assertEqual(hc, hc_ref)

    def test_waveform_filter_taper(self):
        """Time-domain templates are tapered as given by their taper fields"""
        if self.scheme != "cpu":
            self.skipTest("constant tapering needs the CPU")

        def htilde(**taper):
            out = zeros(16 * 4096 // 2 + 1, dtype=complex64)
            return get_waveform_filter(
                out,
                approximant="TaylorT4",
                mass1=15,
                mass2=15,
                f_lower=20,
                delta_f=1.0 / 16,
                delta_t=1.0 / 4096,
                **taper,
            ).numpy()

        untapered = htilde()
        for taper in [
            {"taper": "start"},
            {"taper": "start", "taper_method": "constant", "taper_window": 0.5},
        ]:
            diff = abs(htilde(**taper) - untapered).max()
            self.assertGreater(diff / abs(untapered).max(), 1e-3, taper)

    def test_spintaylorf2_gpu(self):

        print(type(self.context))
        if isinstance(self.context, CPUScheme):
            return

        fl = 25
        delta_f = 1.0 / 256

        for m1 in [3, 5, 15]:
            for m2 in [1.0, 2.0, 3.0]:
                for s1 in [0.001, 1.0, 10]:
                    for s1_cos_theta in [-1.0, 0.0, 0.5, 1.0]:
                        for s1phi in [0, 2.09, 4.18]:
                            for inclination in [0.2, 1.2]:
                                s1x = s1 * sqrt(1 - s1_cos_theta**2) * cos(s1phi)
                                s1y = s1 * sqrt(1 - s1_cos_theta**2) * sin(s1phi)
                                s1z = s1 * s1_cos_theta
                                # Generate SpinTaylorF2 from lalsimulation
                                hp_lal, hc_lal = get_fd_waveform(
                                    mass1=m1,
                                    mass2=m2,
                                    spin1x=s1x,
                                    spin1y=s1y,
                                    spin1z=s1z,
                                    delta_f=delta_f,
                                    f_lower=fl,
                                    approximant="SpinTaylorF2",
                                    amplitude_order=0,
                                    phase_order=7,
                                    inclination=inclination,
                                )

                                # Generate SpinTaylorF2 from spintaylorf2.py
                                with self.context:
                                    hp, hc = get_fd_waveform(
                                        mass1=m1,
                                        mass2=m2,
                                        spin1x=s1x,
                                        spin1y=s1y,
                                        spin1z=s1z,
                                        delta_f=delta_f,
                                        f_lower=fl,
                                        approximant="SpinTaylorF2",
                                        amplitude_order=0,
                                        phase_order=7,
                                        inclination=inclination,
                                    )

                                o = overlap(hp_lal, hp)
                                self.assertAlmostEqual(1.0, o, places=4)
                                o = overlap(hc_lal, hc)
                                self.assertAlmostEqual(1.0, o, places=4)

                                amp_p_lal = numpy.abs(hp_lal.data)
                                amp_p = numpy.abs(hp.data)
                                phase_p_lal = numpy.unwrap(numpy.angle(hp_lal.data))
                                phase_p = numpy.unwrap(numpy.angle(hp.data))
                                amp_c_lal = numpy.abs(hc_lal.data)
                                amp_c = numpy.abs(hc.data)
                                phase_c_lal = numpy.unwrap(numpy.angle(hc_lal.data))
                                phase_c = numpy.unwrap(numpy.angle(hc.data))
                                index_amp_p = numpy.where(amp_p_lal != 0)
                                index_phase_p = numpy.where(phase_p_lal != 0)
                                index_amp_c = numpy.where(amp_c_lal != 0)
                                index_phase_c = numpy.where(phase_c_lal != 0)
                                amp_diff_p = max(
                                    abs(
                                        (amp_p[index_amp_p] - amp_p_lal[index_amp_p])
                                        / amp_p_lal[index_amp_p]
                                    )
                                )
                                phase_diff_p = max(
                                    abs(
                                        (
                                            phase_p[index_phase_p]
                                            - phase_p_lal[index_phase_p]
                                        )
                                        / phase_p_lal[index_phase_p]
                                    )
                                )
                                amp_diff_c = max(
                                    abs(
                                        (amp_c[index_amp_p] - amp_c_lal[index_amp_p])
                                        / amp_c_lal[index_amp_p]
                                    )
                                )
                                phase_diff_c = max(
                                    abs(
                                        (
                                            phase_c[index_phase_p]
                                            - phase_c_lal[index_phase_p]
                                        )
                                        / phase_c_lal[index_phase_p]
                                    )
                                )
                                self.assertTrue(amp_diff_p < 0.00001)
                                self.assertTrue(phase_diff_p < 0.00001)
                                self.assertTrue(amp_diff_c < 0.00001)
                                self.assertTrue(phase_diff_c < 0.00001)
                                print(
                                    f"..checked m1: {m1} m2:: {m2} s1x: {s1x} s1y: "
                                    f"{s1y} s1z: {s1z} Inclination: {inclination}"
                                )

    def test_errors(self):
        func = get_fd_waveform
        self.assertRaises(ValueError, func, approximant="BLAH")
        self.assertRaises(ValueError, func, approximant="SpinTaylorF2", mass1=3)
        self.assertRaises(
            ValueError, func, approximant="SpinTaylorF2", mass1=3, mass2=3
        )
        self.assertRaises(
            ValueError,
            func,
            approximant="SpinTaylorF2",
            mass1=3,
            mass2=3,
            phase_order=7,
        )
        self.assertRaises(
            ValueError,
            func,
            approximant="SpinTaylorF2",
            mass1=3,
            mass2=3,
            phase_order=7,
        )
        self.assertRaises(ValueError, func, approximant="SpinTaylorF2", mass1=3)

        func = get_fd_waveform
        self.assertRaises(ValueError, func, approximant="BLAH")
        self.assertRaises(ValueError, func, approximant="TaylorF2", mass1=3)
        self.assertRaises(ValueError, func, approximant="TaylorF2", mass1=3, mass2=3)
        self.assertRaises(
            ValueError, func, approximant="TaylorF2", mass1=3, mass2=3, phase_order=7
        )
        self.assertRaises(
            ValueError, func, approximant="TaylorF2", mass1=3, mass2=3, phase_order=7
        )
        self.assertRaises(ValueError, func, approximant="TaylorF2", mass1=3)

        for func in [get_fd_waveform, get_td_waveform]:
            self.assertRaises(ValueError, func, approximant="BLAH")
            self.assertRaises(ValueError, func, approximant="IMRPhenomB", mass1=3)
            self.assertRaises(
                ValueError, func, approximant="IMRPhenomB", mass1=3, mass2=3
            )
            self.assertRaises(
                ValueError,
                func,
                approximant="IMRPhenomB",
                mass1=3,
                mass2=3,
                phase_order=7,
            )
            self.assertRaises(
                ValueError,
                func,
                approximant="IMRPhenomB",
                mass1=3,
                mass2=3,
                phase_order=7,
            )
            self.assertRaises(ValueError, func, approximant="IMRPhenomB", mass1=3)


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestWaveform))

if __name__ == "__main__":
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
