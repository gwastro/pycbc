# Copyright (C) 2026  Shichao Wu, Alex Nitz
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
Regression tests for `pycbc.coordinates.space_orbit`.

Fed the circular LISA orbit that `pycbc.coordinates.space` assumes, the
generic orbit-provider machinery must reproduce those hard-coded functions.
The reference orbit is the LISA Data Challenge manual's "Equal arm analytic
orbit" (LISA-LCST-SGS-MAN-001, Sec. 8.1.1, Eq. 48-52), written out below
rather than imported from `lisaorbits`, so these tests have no optional
dependency and are not comparing the code against itself.
"""
import numpy
import unittest
from astropy.constants import au

from pycbc.coordinates import space
from pycbc.coordinates import space_orbit
from utils import simple_exit


seed = 8202
numpy.random.seed(seed)

# Reference constants for the LDC "Equal arm analytic orbit" (LISA flavor).
# The orbit *formula* below is deliberately re-derived rather than imported,
# so these tests compare two independent implementations; the physical
# constants it is evaluated with are shared, not re-typed.
OMEGA_0 = space.EARTH_ORBIT_ANGULAR_FREQUENCY
ARMLENGTH = 2.5e9  # matches pycbc.detector.space._space_detectors['LISA']
SEMI_MAJOR_AXIS = au.value
ECCENTRICITY = ARMLENGTH / (2 * SEMI_MAJOR_AXIS * numpy.sqrt(3))
T0 = space.TIME_OFFSET_20_DEGREES


def _random_sky_position(with_polarization=False):
    """Shared random (lam, beta[, pol]) sky position/polarization, drawn
    fresh each call from this module's seeded RNG. Used throughout this
    file instead of hard-coded literals, so a passing test isn't
    accidentally relying on some property of that one specific value.
    """
    lam = numpy.random.uniform(0.0, 2 * numpy.pi)
    beta = numpy.random.uniform(-numpy.pi / 2, numpy.pi / 2)
    if with_polarization:
        pol = numpy.random.uniform(0.0, 2 * numpy.pi)
        return lam, beta, pol
    return lam, beta


class AnalyticEqualArmOrbit:
    """The LDC manual's "Equal arm analytic orbit" (LISA-LCST-SGS-MAN-001,
    Sec. 8.1.1, Eq. 48-52; Rubbo, Cornish & Poujade 2004, Phys. Rev. D 69,
    082003), LISA flavour, written out independently of `pycbc`.

    Parameters
    ----------
    eccentricity : float, optional
        Constellation eccentricity. Default `ECCENTRICITY`, from LISA's own
        arm length.
    """
    def __init__(self, eccentricity=None):
        self.eccentricity = (
            ECCENTRICITY if eccentricity is None else eccentricity)

    def _phase(self, t):
        return OMEGA_0 * (t + T0)

    def compute_position(self, t, sc=(1, 2, 3)):
        t = numpy.atleast_1d(numpy.asarray(t, dtype=float))
        sc = numpy.atleast_1d(sc)
        alpha = self._phase(t)
        out = numpy.empty((len(t), len(sc), 3))
        for k, n in enumerate(sc):
            beta_n = (n - 1) * 2 * numpy.pi / 3.0
            out[:, k, 0] = SEMI_MAJOR_AXIS * numpy.cos(alpha) \
                + SEMI_MAJOR_AXIS * self.eccentricity * (
                    numpy.sin(alpha) * numpy.cos(alpha) * numpy.sin(beta_n)
                    - (1 + numpy.sin(alpha) ** 2) * numpy.cos(beta_n))
            out[:, k, 1] = SEMI_MAJOR_AXIS * numpy.sin(alpha) \
                + SEMI_MAJOR_AXIS * self.eccentricity * (
                    numpy.sin(alpha) * numpy.cos(alpha) * numpy.cos(beta_n)
                    - (1 + numpy.cos(alpha) ** 2) * numpy.sin(beta_n))
            out[:, k, 2] = -numpy.sqrt(3) * SEMI_MAJOR_AXIS \
                * self.eccentricity * numpy.cos(alpha - beta_n)
        return out


def _arm_lengths(orbit, t):
    pos = orbit.compute_position(t)
    r1, r2, r3 = pos[:, 0], pos[:, 1], pos[:, 2]
    return (numpy.linalg.norm(r1 - r2, axis=-1),
            numpy.linalg.norm(r2 - r3, axis=-1),
            numpy.linalg.norm(r1 - r3, axis=-1))


class TestConstellationFrame(unittest.TestCase):
    """Some tests here pin the new machinery against the existing
    `lisa_position_ssb`, `rotation_matrix_ssb_to_lisa`, `t_lisa_from_ssb`
    and `t_ssb_from_t_lisa`, showing it reduces to current behaviour for a
    circular orbit. The rest are self-contained property checks (arm
    length, orthonormality, round trip, Earth-trailing angle).
    """
    def setUp(self):
        self.orbit = AnalyticEqualArmOrbit()
        self.precision = 1e-6  # relative to ~1 AU length scale / O(1) rotation
        self.times = numpy.random.uniform(0.0, 3.15e7, size=50)

    def test_centroid_matches_lisa_position_ssb(self):
        """The constellation centroid derived from the 3 spacecraft
        positions must match the existing (eccentricity-independent)
        guiding-center formula `lisa_position_ssb`, since the average of the
        Equal Arm orbit's eccentricity terms over the 3 spacecraft is
        exactly zero at every order.
        """
        centroid, _ = space_orbit.constellation_frame(self.times, self.orbit)
        expected = numpy.array([
            space.lisa_position_ssb(t, T0)[0].flatten().astype(float)
            for t in self.times
        ])
        maxdiff = numpy.max(numpy.abs(centroid - expected))
        self.assertLess(maxdiff, self.precision * SEMI_MAJOR_AXIS,
                        f'constellation centroid does not match '
                        f'lisa_position_ssb; max abs diff {maxdiff}')

    def test_rotation_matches_rotation_matrix_ssb_to_lisa(self):
        """The rotation matrix derived from the 3 spacecraft positions must
        match the existing analytic `rotation_matrix_ssb_to_lisa`.
        """
        _, rotation = space_orbit.constellation_frame(self.times, self.orbit)
        for i, t in enumerate(self.times):
            alpha = OMEGA_0 * (t + T0)
            expected = space.rotation_matrix_ssb_to_lisa(alpha)
            maxdiff = numpy.max(numpy.abs(rotation[i] - expected))
            self.assertLess(maxdiff, self.precision,
                            f'rotation matrix at t={t} does not match '
                            f'rotation_matrix_ssb_to_lisa; max abs diff '
                            f'{maxdiff}')

    def test_reference_arrival_time_functions_round_trip(self):
        """The two tests below pin the new functions against
        `t_lisa_from_ssb`/`t_ssb_from_t_lisa`, so that pair has to be
        self-consistent for the comparison to mean anything. On master it
        has no round-trip test of its own -- only the indirect one through
        `test_round_robin` in test_coordinates_space.py, which goes through
        the four-parameter transforms -- so check it directly here first.
        """
        lam, beta = _random_sky_position()
        for t_ssb in self.times[:10]:
            t_lisa = space.t_lisa_from_ssb(t_ssb, lam, beta, T0)
            recovered = space.t_ssb_from_t_lisa(t_lisa, lam, beta, T0)
            self.assertLess(
                abs(recovered - t_ssb), 1e-10,
                f't_lisa_from_ssb/t_ssb_from_t_lisa do not round trip at '
                f't_ssb={t_ssb}')

    def test_t_detector_from_ssb_matches_t_lisa_from_ssb(self):
        lam, beta = _random_sky_position()
        k_ssb = space.localization_to_propagation_vector(
            lam, beta, use_astropy=False).flatten()
        for t_ssb in self.times[:10]:
            expected = space.t_lisa_from_ssb(t_ssb, lam, beta, T0)
            derived = space_orbit.t_detector_from_ssb(
                t_ssb, k_ssb, self.orbit)
            self.assertLess(abs(derived - expected), 1e-10,
                msg=f't_detector_from_ssb does not match t_lisa_from_ssb '
                    f'at t_ssb={t_ssb}')

    def test_t_ssb_from_t_detector_matches_t_ssb_from_t_lisa(self):
        lam, beta = _random_sky_position()
        k_ssb = space.localization_to_propagation_vector(
            lam, beta, use_astropy=False).flatten()
        for t_lisa in self.times[:10]:
            expected = space.t_ssb_from_t_lisa(t_lisa, lam, beta, T0)
            derived = space_orbit.t_ssb_from_t_detector(
                t_lisa, k_ssb, self.orbit)
            self.assertLess(abs(derived - expected), 1e-10,
                msg=f't_ssb_from_t_detector does not match '
                    f't_ssb_from_t_lisa at t_lisa={t_lisa}')

    def test_arm_length_matches_design_value(self):
        for d in _arm_lengths(self.orbit, self.times):
            self.assertLess(
                numpy.max(numpy.abs(d - ARMLENGTH)), 1.0,
                'LISA arm length deviates from the 2.5e6 km design value '
                'by more than 1 m')

    def test_trails_real_earth_by_documented_lag(self):
        """The guiding center must trail the real Earth (not a zero-phase
        circular proxy) by the 19-23 degree range documented for
        `TIME_OFFSET_20_DEGREES`, at any time -- not just near the
        particular epoch that constant happens to have been tuned at.
        """
        centroid, _ = space_orbit.constellation_frame(self.times, self.orbit)
        for i, t in enumerate(self.times):
            earth_pos = space.earth_position_ssb(t)[0].flatten().astype(
                float)
            cos_angle = numpy.dot(centroid[i], earth_pos) / (
                numpy.linalg.norm(centroid[i]) * numpy.linalg.norm(
                    earth_pos))
            angle_deg = numpy.rad2deg(numpy.arccos(
                numpy.clip(cos_angle, -1, 1)))
            self.assertTrue(
                17.0 < angle_deg < 24.0,
                f'LISA-to-real-Earth angle {angle_deg} deg at t={t} is '
                f'outside the documented ~19-23 degree range')

    def test_rotation_orthonormal(self):
        _, rotation = space_orbit.constellation_frame(self.times, self.orbit)
        for r in rotation:
            self.assertLess(
                numpy.max(numpy.abs(r @ r.T - numpy.eye(3))), 1e-8)

    def test_round_trip_time_delay(self):
        lam, beta = _random_sky_position()
        k_ssb = space.localization_to_propagation_vector(
            lam, beta, use_astropy=False).flatten()
        for t_ssb in self.times[:10]:
            t_det = space_orbit.t_detector_from_ssb(t_ssb, k_ssb, self.orbit)
            t_ssb_roundtrip = space_orbit.t_ssb_from_t_detector(
                t_det, k_ssb, self.orbit)
            self.assertLess(abs(t_ssb_roundtrip - t_ssb), 1e-10)


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestConstellationFrame))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
