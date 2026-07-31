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

The core claim being tested is backward compatibility: the generic,
orbit-provider-based machinery in `space_orbit` must reproduce the existing,
hard-coded circular-orbit LISA functions in `pycbc.coordinates.space` exactly
(to floating point precision) when fed the same analytic circular orbit that
those functions assume. The reference orbit used here is the "Equal arm
analytic orbit" defined in the LISA Data Challenge manual
(LISA-LCST-SGS-MAN-001, Sec. 8.1.1, Eq. 48-52); this is the same closed-form
per-spacecraft trajectory implemented by `lisaorbits.EqualArmlengthOrbits`,
reimplemented here directly (not imported from `lisaorbits`) so this test has
no optional dependency.

Beyond the LISA special case, this file also checks that `constellation_frame`
generalizes correctly to constellations with genuinely different geometry:
Taiji (heliocentric, like LISA but leading the Earth by 20 degrees instead of
trailing it, with a different arm length/eccentricity) and TianQin (a fast,
rigidly-rotating triangle around a geocentric guiding center, in a plane
fixed in inertial space rather than precessing over a year). The Taiji
fixture uses the same first-order-in-eccentricity Keplerian expansion as the
LDC manual's Eq. 48-52 above (Rubbo, Cornish & Poujade 2004, Phys. Rev. D 69,
082003), just with Taiji's own arm length/eccentricity and a +20 degree lead
angle instead of LISA's trailing offset. The TianQin fixture follows Hu et al
2018 (Class. Quantum Grav. 35, 095008), simplified to a pure circular
heliocentric guiding center coincident with the Earth (consistent with how
the LISA/Taiji guiding center above is also treated as the zero-eccentricity
limit of the same family of orbits).
"""
import os
import tempfile
import h5py
import numpy
import unittest
from astropy.constants import au

from pycbc import transforms
from pycbc.coordinates import space
from pycbc.coordinates import space_orbit
from utils import simple_exit


seed = 8202
numpy.random.seed(seed)

OMEGA_0 = 1.99098659277e-7  # LISA/Taiji-like orbital angular frequency [rad/s]
ARMLENGTH = 2.5e9  # matches pycbc.detector.space._space_detectors['LISA']
SEMI_MAJOR_AXIS = au.value
ECCENTRICITY = ARMLENGTH / (2 * SEMI_MAJOR_AXIS * numpy.sqrt(3))
T0 = space.TIME_OFFSET_20_DEGREES

# Real Earth's ecliptic longitude at GPS t=0, taken directly from pycbc's own
# (astropy-based) `earth_position_ssb` -- used below as the phase-zero
# reference for the Taiji/TianQin fixtures' guiding center, so that "leads/
# coincides with the Earth" checks can be made against the real Earth
# (`earth_position_ssb` itself) rather than an arbitrarily-phased toy
# circular reference. No external ephemeris/orbital-element constants are
# used here; this is the one already-existing real-ephemeris function in
# pycbc.coordinates.space, evaluated once.
PHI0_REAL = space.earth_position_ssb(0.0)[1]


class AnalyticEqualArmOrbit:
    """Reference implementation of the LDC manual's 'Equal arm analytic
    orbit' (LISA-LCST-SGS-MAN-001, Sec. 8.1.1, Eq. 48-52). Used only to
    validate `space_orbit.constellation_frame` and related functions against
    the existing analytic functions in `pycbc.coordinates.space`, which
    assume this same orbit but only track the (eccentricity-independent)
    guiding center.
    """
    def compute_position(self, t, sc=(1, 2, 3)):
        t = numpy.atleast_1d(numpy.asarray(t, dtype=float))
        sc = numpy.atleast_1d(sc)
        alpha = OMEGA_0 * (t + T0)
        out = numpy.empty((len(t), len(sc), 3))
        for k, n in enumerate(sc):
            beta_n = (n - 1) * 2 * numpy.pi / 3.0
            out[:, k, 0] = SEMI_MAJOR_AXIS * numpy.cos(alpha) \
                + SEMI_MAJOR_AXIS * ECCENTRICITY * (
                    numpy.sin(alpha) * numpy.cos(alpha) * numpy.sin(beta_n)
                    - (1 + numpy.sin(alpha) ** 2) * numpy.cos(beta_n))
            out[:, k, 1] = SEMI_MAJOR_AXIS * numpy.sin(alpha) \
                + SEMI_MAJOR_AXIS * ECCENTRICITY * (
                    numpy.sin(alpha) * numpy.cos(alpha) * numpy.cos(beta_n)
                    - (1 + numpy.cos(alpha) ** 2) * numpy.sin(beta_n))
            out[:, k, 2] = -numpy.sqrt(3) * SEMI_MAJOR_AXIS * ECCENTRICITY \
                * numpy.cos(alpha - beta_n)
        return out


# Taiji shares the LISA fixture's functional form (the same
# first-order-in-eccentricity Keplerian expansion as the LDC manual's
# Eq. 48-52 above, per Rubbo, Cornish & Poujade 2004, Phys. Rev. D 69,
# 082003), differing only in arm length/eccentricity and in leading (rather
# than trailing) the Earth-like guiding center by 20 degrees. Its phase is
# anchored to PHI0_REAL (real Earth's longitude at t=0) so that "leads the
# Earth by 20 degrees" can be checked against the real Earth position, not
# just an arbitrarily-phased circular proxy.
TAIJI_ARMLENGTH = 3.0e9  # matches pycbc.detector.space._space_detectors
TAIJI_ECCENTRICITY = TAIJI_ARMLENGTH / (2 * SEMI_MAJOR_AXIS * numpy.sqrt(3))
TAIJI_LEAD_ANGLE = numpy.deg2rad(20.0)


class AnalyticTaijiOrbit:
    """Reference implementation of the Taiji heliocentric orbit (first order
    in eccentricity, per Rubbo, Cornish & Poujade 2004, Phys. Rev. D 69,
    082003). Used only to check that `constellation_frame`/`NumericOrbits`
    handle a constellation genuinely different from the LISA fixture above
    (different arm length, eccentricity, and heliocentric lead angle).
    """
    def compute_position(self, t, sc=(1, 2, 3)):
        t = numpy.atleast_1d(numpy.asarray(t, dtype=float))
        sc = numpy.atleast_1d(sc)
        alpha = OMEGA_0 * t + PHI0_REAL + TAIJI_LEAD_ANGLE
        out = numpy.empty((len(t), len(sc), 3))
        for k, n in enumerate(sc):
            beta_n = (n - 1) * 2 * numpy.pi / 3.0
            out[:, k, 0] = SEMI_MAJOR_AXIS * numpy.cos(alpha) \
                + SEMI_MAJOR_AXIS * TAIJI_ECCENTRICITY * (
                    numpy.sin(alpha) * numpy.cos(alpha) * numpy.sin(beta_n)
                    - (1 + numpy.sin(alpha) ** 2) * numpy.cos(beta_n))
            out[:, k, 1] = SEMI_MAJOR_AXIS * numpy.sin(alpha) \
                + SEMI_MAJOR_AXIS * TAIJI_ECCENTRICITY * (
                    numpy.sin(alpha) * numpy.cos(alpha) * numpy.cos(beta_n)
                    - (1 + numpy.cos(alpha) ** 2) * numpy.sin(beta_n))
            out[:, k, 2] = -numpy.sqrt(3) * SEMI_MAJOR_AXIS \
                * TAIJI_ECCENTRICITY * numpy.cos(alpha - beta_n)
        return out


# TianQin: a fast-rotating rigid triangle (per Hu et al 2018, Class.
# Quantum Grav. 35, 095008) whose plane is fixed in inertial space (pointing
# towards the calibration source RX J0806.3+1527), around a guiding center
# that coincides with the Earth. The guiding center itself is simplified to
# a pure circular heliocentric orbit (Earth's real orbital eccentricity is
# irrelevant to what most tests below check), but its phase is anchored to
# PHI0_REAL (real Earth's longitude at t=0) so that "coincides with the
# Earth" can also be checked against the real Earth position, not just an
# arbitrarily-phased circular proxy.
TIANQIN_ARMLENGTH = 1.7e8  # matches pycbc.detector.space._space_detectors
TIANQIN_LAMBDA_S = numpy.deg2rad(120.5)  # direction to RX J0806.3+1527
TIANQIN_BETA_S = numpy.deg2rad(-4.7)
TIANQIN_F_SC = 1.0 / (3.65 * 86400.0)  # ~3.65 day rotation period


class AnalyticTianQinOrbit:
    """Reference implementation of the TianQin geocentric orbit (per Hu et
    al 2018, Class. Quantum Grav. 35, 095008). Used to check that
    `constellation_frame` correctly derives a *non-precessing* constellation
    plane, in contrast to LISA/Taiji above.
    """
    def compute_position(self, t, sc=(1, 2, 3), initial_orbit_phase=0.0):
        t = numpy.atleast_1d(numpy.asarray(t, dtype=float))
        sc = numpy.atleast_1d(sc)
        alpha_earth = OMEGA_0 * t + PHI0_REAL
        earth = SEMI_MAJOR_AXIS * numpy.stack([
            numpy.cos(alpha_earth), numpy.sin(alpha_earth),
            numpy.zeros_like(t)], axis=-1)
        out = numpy.empty((len(t), len(sc), 3))
        for k, n in enumerate(sc):
            kappa_n = 2 * numpy.pi / 3.0 * (n - 1) + initial_orbit_phase
            phase = 2 * numpy.pi * TIANQIN_F_SC * t + kappa_n
            xn = (TIANQIN_ARMLENGTH / numpy.sqrt(3)) * (
                numpy.sin(TIANQIN_BETA_S) * numpy.cos(TIANQIN_LAMBDA_S)
                * numpy.sin(phase)
                + numpy.sin(TIANQIN_LAMBDA_S) * numpy.cos(phase))
            yn = (TIANQIN_ARMLENGTH / numpy.sqrt(3)) * (
                numpy.sin(TIANQIN_BETA_S) * numpy.sin(TIANQIN_LAMBDA_S)
                * numpy.sin(phase)
                - numpy.cos(TIANQIN_LAMBDA_S) * numpy.cos(phase))
            zn = -(TIANQIN_ARMLENGTH / numpy.sqrt(3)) \
                * numpy.cos(TIANQIN_BETA_S) * numpy.sin(phase)
            out[:, k, 0] = earth[:, 0] + xn
            out[:, k, 1] = earth[:, 1] + yn
            out[:, k, 2] = earth[:, 2] + zn
        return out


def _arm_lengths(orbit, t):
    pos = orbit.compute_position(t)
    r1, r2, r3 = pos[:, 0], pos[:, 1], pos[:, 2]
    return (numpy.linalg.norm(r1 - r2, axis=-1),
            numpy.linalg.norm(r2 - r3, axis=-1),
            numpy.linalg.norm(r1 - r3, axis=-1))


class TestConstellationFrame(unittest.TestCase):
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

    def test_t_detector_from_ssb_matches_t_lisa_from_ssb(self):
        lam, beta = 1.234, 0.456
        k_ssb = space.localization_to_propagation_vector(
            lam, beta, use_astropy=False).flatten()
        for t_ssb in self.times[:10]:
            expected = space.t_lisa_from_ssb(t_ssb, lam, beta, T0)
            derived = space_orbit.t_detector_from_ssb(
                t_ssb, k_ssb, self.orbit)
            self.assertAlmostEqual(derived, expected, places=3,
                msg=f't_detector_from_ssb does not match t_lisa_from_ssb '
                    f'at t_ssb={t_ssb}')

    def test_t_ssb_from_t_detector_matches_t_ssb_from_t_lisa(self):
        lam, beta = 1.234, 0.456
        k_ssb = space.localization_to_propagation_vector(
            lam, beta, use_astropy=False).flatten()
        for t_lisa in self.times[:10]:
            expected = space.t_ssb_from_t_lisa(t_lisa, lam, beta, T0)
            derived = space_orbit.t_ssb_from_t_detector(
                t_lisa, k_ssb, self.orbit)
            self.assertAlmostEqual(derived, expected, places=3,
                msg=f't_ssb_from_t_detector does not match '
                    f't_ssb_from_t_lisa at t_lisa={t_lisa}')


class TestLinkVector(unittest.TestCase):
    """`link_vector` is the per-arm building block (unit vector + arm
    length) that single-link response formulas need, on top of the
    whole-constellation quantities `constellation_frame` already provides.
    Checked here against all three constellation fixtures (LISA, Taiji,
    TianQin), since a per-arm quantity should behave identically regardless
    of the specific orbit's geometry.
    """
    def setUp(self):
        self.times = numpy.random.uniform(0.0, 3.15e7, size=20)

    def _check_arm_length_and_normalization(self, orbit, expected_armlength):
        for i, j in [(1, 2), (2, 3), (1, 3)]:
            unit_vector, arm_length = space_orbit.link_vector(
                self.times, orbit, i, j)
            self.assertEqual(unit_vector.shape, (len(self.times), 1, 3))
            self.assertEqual(arm_length.shape, (len(self.times), 1))
            self.assertLess(
                numpy.max(numpy.abs(arm_length - expected_armlength)), 1.0,
                f'link ({i},{j}) length deviates from the design arm '
                f'length by more than 1 m')
            norms = numpy.linalg.norm(unit_vector, axis=-1)
            self.assertLess(numpy.max(numpy.abs(norms - 1.0)), 1e-10,
                            'link unit vector is not normalized')

    def test_lisa_arm_length_and_normalization(self):
        self._check_arm_length_and_normalization(
            AnalyticEqualArmOrbit(), ARMLENGTH)

    def test_taiji_arm_length_and_normalization(self):
        self._check_arm_length_and_normalization(
            AnalyticTaijiOrbit(), TAIJI_ARMLENGTH)

    def test_tianqin_arm_length_and_normalization(self):
        self._check_arm_length_and_normalization(
            AnalyticTianQinOrbit(), TIANQIN_ARMLENGTH)

    def test_antiparallel_and_symmetric_length(self):
        orbit = AnalyticTaijiOrbit()
        for i, j in [(1, 2), (2, 3), (1, 3)]:
            uv_ij, len_ij = space_orbit.link_vector(self.times, orbit, i, j)
            uv_ji, len_ji = space_orbit.link_vector(self.times, orbit, j, i)
            self.assertLess(
                numpy.max(numpy.abs(uv_ij + uv_ji)), 1e-10,
                f'link ({i},{j}) and ({j},{i}) unit vectors are not '
                f'antiparallel')
            self.assertLess(numpy.max(numpy.abs(len_ij - len_ji)), 1e-6,
                            f'link ({i},{j}) and ({j},{i}) lengths differ')

    def test_vectorized_multi_link_call_matches_individual_calls(self):
        orbit = AnalyticTianQinOrbit()
        sc_emitter = [1, 2, 3]
        sc_receiver = [2, 3, 1]
        uv_batch, len_batch = space_orbit.link_vector(
            self.times, orbit, sc_emitter, sc_receiver)
        self.assertEqual(uv_batch.shape, (len(self.times), 3, 3))
        for k, (i, j) in enumerate(zip(sc_emitter, sc_receiver)):
            uv_single, len_single = space_orbit.link_vector(
                self.times, orbit, i, j)
            self.assertLess(
                numpy.max(numpy.abs(uv_batch[:, k] - uv_single[:, 0])), 1e-10)
            self.assertLess(
                numpy.max(numpy.abs(len_batch[:, k] - len_single[:, 0])),
                1e-6)


class TestNumericOrbits(unittest.TestCase):
    def setUp(self):
        self.orbit = AnalyticEqualArmOrbit()
        # coarse grid, similar density to a typical lisaorbits orbit file
        self.t_grid = numpy.linspace(0.0, 3.15e7, 400)
        self.positions = self.orbit.compute_position(self.t_grid)
        self.numeric = space_orbit.NumericOrbits(self.t_grid, self.positions)
        self.t_query = numpy.random.uniform(1e5, 3.1e7, size=20)

    def test_position_interpolation_accuracy(self):
        """Interpolated positions must reproduce the analytic orbit to a
        small fraction of the ~1 AU length scale, at off-grid query times.
        """
        interpolated = self.numeric.compute_position(self.t_query)
        exact = self.orbit.compute_position(self.t_query)
        maxdiff = numpy.max(numpy.abs(interpolated - exact))
        self.assertLess(maxdiff, 1e-6 * SEMI_MAJOR_AXIS,
                        f'NumericOrbits position interpolation error too '
                        f'large: {maxdiff} m')

    def test_constellation_frame_consistent_with_numeric_orbits(self):
        """Feeding `constellation_frame` a `NumericOrbits` instance built
        from a dense sampling of the analytic orbit must reproduce the
        result of feeding it the analytic orbit provider directly.
        """
        centroid_num, rotation_num = space_orbit.constellation_frame(
            self.t_query, self.numeric)
        centroid_ana, rotation_ana = space_orbit.constellation_frame(
            self.t_query, self.orbit)
        self.assertLess(
            numpy.max(numpy.abs(centroid_num - centroid_ana)),
            1e-6 * SEMI_MAJOR_AXIS)
        self.assertLess(
            numpy.max(numpy.abs(rotation_num - rotation_ana)), 1e-6)

    def test_rejects_mismatched_shapes(self):
        with self.assertRaises(ValueError):
            space_orbit.NumericOrbits(self.t_grid, self.positions[:-1])
        with self.assertRaises(ValueError):
            space_orbit.NumericOrbits(
                self.t_grid, self.positions,
                velocities=self.positions[:, :, :2])

    def test_velocity_interpolation_accuracy(self):
        """`compute_velocity`, auto-derived from the position spline (no
        velocities given), must reproduce the analytic orbit's own exact
        velocity at off-grid query times.
        """
        exact_orbit = space_orbit.LisaAnalyticOrbit()
        positions = exact_orbit.compute_position(self.t_grid)
        numeric = space_orbit.NumericOrbits(self.t_grid, positions)
        interpolated = numeric.compute_velocity(self.t_query)
        exact = exact_orbit.compute_velocity(self.t_query)
        maxdiff = numpy.max(numpy.abs(interpolated - exact))
        self.assertLess(maxdiff, 1e-3,
                        f'NumericOrbits velocity interpolation error too '
                        f'large: {maxdiff} m/s')

    def test_acceleration_interpolation_accuracy(self):
        """`compute_acceleration`, derived as the second analytic
        derivative of the position spline (no velocities or accelerations
        given), must reproduce the analytic orbit's own exact acceleration
        at off-grid query times.
        """
        exact_orbit = space_orbit.LisaAnalyticOrbit()
        positions = exact_orbit.compute_position(self.t_grid)
        numeric = space_orbit.NumericOrbits(self.t_grid, positions)
        interpolated = numeric.compute_acceleration(self.t_query)
        exact = exact_orbit.compute_acceleration(self.t_query)
        maxdiff = numpy.max(numpy.abs(interpolated - exact))
        self.assertLess(maxdiff, 1e-6,
                        f'NumericOrbits acceleration interpolation error '
                        f'too large: {maxdiff} m/s^2')


class TestLisaOrbit(unittest.TestCase):
    """Standalone sanity checks for the LISA fixture, mirroring
    `TestTaijiOrbit`/`TestTianQinOrbit` below. `TestConstellationFrame`
    above already gives LISA a stronger check (byte-level agreement with
    the pre-existing hardcoded `pycbc.coordinates.space` formulas), but
    that is a different kind of test; these are the same direct,
    self-contained checks the other two missions get, for consistency.

    Unlike Taiji's fixture (whose 20 degree lead is a bare angular offset
    added directly to its own phase, checked against an equally-invented
    zero-phase circular proxy for the Earth), this fixture's `T0` is a time
    shift, not an angle, and it does not correspond to a clean angle offset
    from an arbitrary zero-phase circular reference -- comparing against
    one gives a constant but physically meaningless ~84 degrees. Comparing
    instead against the real Earth position (`space.earth_position_ssb`,
    real ephemeris) does give the intended ~20 degree lag, matching the
    documented 19-23 degree range for `TIME_OFFSET_20_DEGREES`.
    """
    def setUp(self):
        self.orbit = AnalyticEqualArmOrbit()
        self.times = numpy.random.uniform(0.0, 3.15e7, size=50)

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
        lam, beta = 1.1, -0.4
        k_ssb = space.localization_to_propagation_vector(
            lam, beta, use_astropy=False).flatten()
        for t_ssb in self.times[:10]:
            t_det = space_orbit.t_detector_from_ssb(t_ssb, k_ssb, self.orbit)
            t_ssb_roundtrip = space_orbit.t_ssb_from_t_detector(
                t_det, k_ssb, self.orbit)
            self.assertAlmostEqual(t_ssb_roundtrip, t_ssb, places=3)


class TestTaijiOrbit(unittest.TestCase):
    def setUp(self):
        self.orbit = AnalyticTaijiOrbit()
        self.times = numpy.random.uniform(0.0, 3.15e7, size=50)

    def test_arm_length_matches_design_value(self):
        for d in _arm_lengths(self.orbit, self.times):
            self.assertLess(
                numpy.max(numpy.abs(d - TAIJI_ARMLENGTH)), 1.0,
                'Taiji arm length deviates from the 3e6 km design value '
                'by more than 1 m')

    def test_leads_earth_like_guiding_center_by_20_degrees(self):
        """Taiji's guiding center must lead a coincident, zero-eccentricity
        circular orbit anchored to the same real-Earth phase reference
        (PHI0_REAL) as the Taiji fixture's own zeroth-order term, by exactly
        20 degrees, at every time (alpha'' = alpha - beta + 20 deg).
        """
        centroid, _ = space_orbit.constellation_frame(self.times, self.orbit)
        alpha_earth = OMEGA_0 * self.times + PHI0_REAL
        earth_like = SEMI_MAJOR_AXIS * numpy.stack([
            numpy.cos(alpha_earth), numpy.sin(alpha_earth),
            numpy.zeros_like(self.times)], axis=-1)
        cos_angle = numpy.sum(centroid * earth_like, axis=-1) / (
            numpy.linalg.norm(centroid, axis=-1)
            * numpy.linalg.norm(earth_like, axis=-1))
        angle_deg = numpy.rad2deg(numpy.arccos(numpy.clip(cos_angle, -1, 1)))
        self.assertLess(numpy.max(numpy.abs(angle_deg - 20.0)), 1e-6)

    def test_leads_real_earth_by_approximately_20_degrees(self):
        """Because the guiding center's phase is anchored to real Earth's
        longitude at t=0 (PHI0_REAL), Taiji must also lead the *real* Earth
        position (`space.earth_position_ssb`, real ephemeris, not the
        idealized circular proxy above) by approximately 20 degrees -- not
        exactly, since this fixture's guiding center is still an idealized
        circle while the real Earth's orbit is eccentric, but within a few
        degrees, and stably so over time (checked here over one year;
        verified separately to hold over multiple decades).
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
                15.0 < angle_deg < 25.0,
                f'Taiji-to-real-Earth angle {angle_deg} deg at t={t} is '
                f'outside the expected ~20 +/- 5 degree range')

    def test_rotation_orthonormal(self):
        _, rotation = space_orbit.constellation_frame(self.times, self.orbit)
        for r in rotation:
            self.assertLess(
                numpy.max(numpy.abs(r @ r.T - numpy.eye(3))), 1e-8)

    def test_round_trip_time_delay(self):
        lam, beta = 0.7, -0.3
        k_ssb = space.localization_to_propagation_vector(
            lam, beta, use_astropy=False).flatten()
        for t_ssb in self.times[:10]:
            t_det = space_orbit.t_detector_from_ssb(t_ssb, k_ssb, self.orbit)
            t_ssb_roundtrip = space_orbit.t_ssb_from_t_detector(
                t_det, k_ssb, self.orbit)
            self.assertAlmostEqual(t_ssb_roundtrip, t_ssb, places=3)


class TestTianQinOrbit(unittest.TestCase):
    def setUp(self):
        self.orbit = AnalyticTianQinOrbit()
        self.times = numpy.random.uniform(0.0, 3.15e7, size=50)

    def test_arm_length_matches_design_value(self):
        for d in _arm_lengths(self.orbit, self.times):
            self.assertLess(
                numpy.max(numpy.abs(d - TIANQIN_ARMLENGTH)), 1.0,
                'TianQin arm length deviates from the 1.7e5 km design '
                'value by more than 1 m')

    def test_centroid_matches_earthlike_guiding_center(self):
        """TianQin's 3 spacecraft are equally spaced in phase around their
        guiding center, so their average must cancel the fast-rotation term
        exactly, leaving just the (here, circular, real-Earth-phase-anchored
        via PHI0_REAL) Earth-like guiding center.
        """
        centroid, _ = space_orbit.constellation_frame(self.times, self.orbit)
        alpha_earth = OMEGA_0 * self.times + PHI0_REAL
        earth_like = SEMI_MAJOR_AXIS * numpy.stack([
            numpy.cos(alpha_earth), numpy.sin(alpha_earth),
            numpy.zeros_like(self.times)], axis=-1)
        self.assertLess(numpy.max(numpy.abs(centroid - earth_like)), 1e-3)

    def test_centroid_close_to_real_earth_position(self):
        """Because the guiding center's phase is anchored to real Earth's
        longitude at t=0 (PHI0_REAL), TianQin's centroid must also stay
        close to the *real* Earth position (`space.earth_position_ssb`),
        not just the idealized circular proxy above -- within a few percent
        of 1 AU, the expected residual from treating the guiding center as
        a perfect circle instead of Earth's real, slightly eccentric orbit.
        """
        centroid, _ = space_orbit.constellation_frame(self.times, self.orbit)
        for i, t in enumerate(self.times):
            earth_pos = space.earth_position_ssb(t)[0].flatten().astype(
                float)
            diff = numpy.linalg.norm(centroid[i] - earth_pos)
            self.assertLess(
                diff, 0.06 * SEMI_MAJOR_AXIS,
                f'TianQin centroid at t={t} is more than 6% of 1 AU away '
                f'from the real Earth position (diff={diff:.3e} m)')

    def test_constellation_plane_normal_is_time_independent(self):
        """Unlike LISA/Taiji (whose constellation plane precesses over a
        year), TianQin's plane is fixed in inertial space (pointing towards
        RX J0806.3+1527): the normal derived by `constellation_frame` must
        be the same at every sampled time, even though the in-plane axes
        rotate rapidly (~3.65 day period).
        """
        year_times = numpy.linspace(0.0, 3.15e7, 50)
        _, rotation = space_orbit.constellation_frame(year_times, self.orbit)
        normal = rotation[:, 2, :]
        x_axis = rotation[:, 0, :]
        self.assertLess(
            numpy.max(numpy.abs(normal - normal[0])), 1e-6,
            'TianQin constellation plane normal is not time-independent')
        # the in-plane axis, in contrast, must rotate substantially over a
        # year (~100 fast rotations), i.e. this is not a degenerate check
        self.assertGreater(numpy.max(numpy.abs(x_axis - x_axis[0])), 0.5)

    def test_rotation_orthonormal(self):
        _, rotation = space_orbit.constellation_frame(self.times, self.orbit)
        for r in rotation:
            self.assertLess(
                numpy.max(numpy.abs(r @ r.T - numpy.eye(3))), 1e-8)

    def test_round_trip_time_delay(self):
        lam, beta = 2.1, 0.5
        k_ssb = space.localization_to_propagation_vector(
            lam, beta, use_astropy=False).flatten()
        for t_ssb in self.times[:10]:
            t_det = space_orbit.t_detector_from_ssb(t_ssb, k_ssb, self.orbit)
            t_ssb_roundtrip = space_orbit.t_ssb_from_t_detector(
                t_det, k_ssb, self.orbit)
            self.assertAlmostEqual(t_ssb_roundtrip, t_ssb, places=3)


class TestProductionAnalyticOrbits(unittest.TestCase):
    """`space_orbit.LisaAnalyticOrbit`/`TaijiAnalyticOrbit`/
    `TianQinAnalyticOrbit` are the production (not test-only) analytic
    reference orbits: unlike the hand-written `AnalyticEqualArmOrbit`/
    `AnalyticTaijiOrbit`/`AnalyticTianQinOrbit` fixtures above (used purely
    to give `constellation_frame` etc. an independent implementation to
    check against), these are importable, user-facing classes intended
    for prototyping single-link response/TDI work before an official
    numeric orbit product exists.

    The key check here is that each *separately implemented* production
    class agrees with the hand-written test fixture bit for bit when given
    the same reference phase -- i.e. that promoting the fixtures' math to
    a real, documented, parameterized class did not introduce any
    transcription error. `LisaAnalyticOrbit` additionally must reproduce
    `pycbc.coordinates.space`'s pre-existing hardcoded functions exactly
    with its default (no-argument) `t0`, since that default is meant to be
    a drop-in stand-in for the existing `orbit=None` code path -- unlike
    Taiji/TianQin, it does not default to a real-Earth-anchored epoch.
    """
    def test_lisa_matches_hardcoded_space_functions(self):
        orbit = space_orbit.LisaAnalyticOrbit()
        times = numpy.random.uniform(0.0, 3.15e7, size=20)
        centroid, rotation = space_orbit.constellation_frame(times, orbit)
        expected_centroid = numpy.array([
            space.lisa_position_ssb(t, orbit.t0)[0].flatten().astype(float)
            for t in times
        ])
        self.assertLess(
            numpy.max(numpy.abs(centroid - expected_centroid)),
            1e-6 * SEMI_MAJOR_AXIS)
        for i, t in enumerate(times):
            alpha = OMEGA_0 * (t + orbit.t0)
            expected_rotation = space.rotation_matrix_ssb_to_lisa(alpha)
            self.assertLess(
                numpy.max(numpy.abs(rotation[i] - expected_rotation)), 1e-6)

    def test_lisa_matches_independent_fixture_implementation(self):
        production = space_orbit.LisaAnalyticOrbit()
        fixture = AnalyticEqualArmOrbit()
        times = numpy.random.uniform(0.0, 3.15e7, size=20)
        pos_production = production.compute_position(times)
        pos_fixture = fixture.compute_position(times)
        self.assertLess(
            numpy.max(numpy.abs(pos_production - pos_fixture)), 1e-3,
            'LisaAnalyticOrbit does not match the independent '
            'AnalyticEqualArmOrbit test fixture')

    def test_lisa_default_t0_matches_time_offset_20_degrees(self):
        orbit = space_orbit.LisaAnalyticOrbit()
        self.assertEqual(orbit.t0, space.TIME_OFFSET_20_DEGREES)

    def test_lisa_custom_t0_overrides_default(self):
        orbit = space_orbit.LisaAnalyticOrbit(t0=0.0)
        self.assertEqual(orbit.t0, 0.0)
        default_orbit = space_orbit.LisaAnalyticOrbit()
        self.assertGreater(
            numpy.max(numpy.abs(
                orbit.compute_position([1e7])
                - default_orbit.compute_position([1e7]))),
            1e6)

    def test_taiji_matches_independent_fixture_implementation(self):
        production = space_orbit.TaijiAnalyticOrbit(kappa0=PHI0_REAL)
        fixture = AnalyticTaijiOrbit()
        times = numpy.random.uniform(0.0, 3.15e7, size=20)
        pos_production = production.compute_position(times)
        pos_fixture = fixture.compute_position(times)
        self.assertLess(
            numpy.max(numpy.abs(pos_production - pos_fixture)), 1e-3,
            'TaijiAnalyticOrbit does not match the independent '
            'AnalyticTaijiOrbit test fixture')

    def test_tianqin_matches_independent_fixture_implementation(self):
        production = space_orbit.TianQinAnalyticOrbit(kappa0=PHI0_REAL)
        fixture = AnalyticTianQinOrbit()
        times = numpy.random.uniform(0.0, 3.15e7, size=20)
        pos_production = production.compute_position(times)
        pos_fixture = fixture.compute_position(times)
        self.assertLess(
            numpy.max(numpy.abs(pos_production - pos_fixture)), 1e-3,
            'TianQinAnalyticOrbit does not match the independent '
            'AnalyticTianQinOrbit test fixture')

    def test_default_kappa0_anchors_to_real_earth(self):
        for cls in (space_orbit.TaijiAnalyticOrbit,
                    space_orbit.TianQinAnalyticOrbit):
            orbit = cls()
            self.assertAlmostEqual(orbit.kappa0, PHI0_REAL, places=9)

    def test_custom_kappa0_overrides_default(self):
        for cls in (space_orbit.TaijiAnalyticOrbit,
                    space_orbit.TianQinAnalyticOrbit):
            orbit = cls(kappa0=0.0)
            self.assertEqual(orbit.kappa0, 0.0)
            default_orbit = cls()
            # different kappa0 must give a different position at a fixed
            # time (sanity check that the parameter actually does
            # something, not a silently-ignored constructor argument)
            self.assertGreater(
                numpy.max(numpy.abs(
                    orbit.compute_position([1e7])
                    - default_orbit.compute_position([1e7]))),
                1e6)

    def test_usable_as_orbit_provider(self):
        """Both classes must work as a drop-in `orbit=` argument to the
        already-generalized `pycbc.coordinates.space` functions, exactly
        like any other orbit provider (NumericOrbits, a test fixture, or a
        real lisaorbits.Orbits instance).
        """
        lam, beta, pol = 1.3, -0.2, 0.8
        for cls in (space_orbit.LisaAnalyticOrbit,
                    space_orbit.TaijiAnalyticOrbit,
                    space_orbit.TianQinAnalyticOrbit):
            orbit = cls()
            t_ssb = 2.0e7
            t_det, lam_det, beta_det, pol_det = space.ssb_to_lisa(
                t_ssb, lam, beta, pol, orbit=orbit)
            t_rt, lam_rt, beta_rt, pol_rt = space.lisa_to_ssb(
                t_det, lam_det, beta_det, pol_det, orbit=orbit)
            self.assertAlmostEqual(t_rt, t_ssb, places=3)
            self.assertAlmostEqual(lam_rt, lam, places=6)
            self.assertAlmostEqual(beta_rt, beta, places=6)
            self.assertAlmostEqual(pol_rt, pol, places=6)

    def test_arm_lengths_match_design_values(self):
        cases = [
            (space_orbit.LisaAnalyticOrbit(), ARMLENGTH),
            (space_orbit.TaijiAnalyticOrbit(), TAIJI_ARMLENGTH),
            (space_orbit.TianQinAnalyticOrbit(), TIANQIN_ARMLENGTH),
        ]
        times = numpy.random.uniform(0.0, 3.15e7, size=20)
        for orbit, expected_armlength in cases:
            for d in _arm_lengths(orbit, times):
                self.assertLess(
                    numpy.max(numpy.abs(d - expected_armlength)), 1.0)

    def test_velocity_matches_finite_difference_of_position(self):
        """`compute_velocity` is an exact analytic derivative, not a
        finite-difference approximation -- but a central finite difference
        of `compute_position`, at a small enough step, must still agree
        with it to high precision. This is independent of `lisaorbits`.
        """
        dt = 1.0  # s
        for orbit in (space_orbit.LisaAnalyticOrbit(),
                      space_orbit.TaijiAnalyticOrbit(),
                      space_orbit.TianQinAnalyticOrbit()):
            t = numpy.random.uniform(1e5, 3.1e7, size=10)
            pos_plus = orbit.compute_position(t + dt)
            pos_minus = orbit.compute_position(t - dt)
            finite_diff_vel = (pos_plus - pos_minus) / (2 * dt)
            analytic_vel = orbit.compute_velocity(t)
            self.assertLess(
                numpy.max(numpy.abs(analytic_vel - finite_diff_vel)), 1e-3)

    def test_acceleration_matches_finite_difference_of_velocity(self):
        """Same rationale as `test_velocity_matches_finite_difference_of_
        position`, one derivative order up.
        """
        dt = 1.0  # s
        for orbit in (space_orbit.LisaAnalyticOrbit(),
                      space_orbit.TaijiAnalyticOrbit(),
                      space_orbit.TianQinAnalyticOrbit()):
            t = numpy.random.uniform(1e5, 3.1e7, size=10)
            vel_plus = orbit.compute_velocity(t + dt)
            vel_minus = orbit.compute_velocity(t - dt)
            finite_diff_acc = (vel_plus - vel_minus) / (2 * dt)
            analytic_acc = orbit.compute_acceleration(t)
            self.assertLess(
                numpy.max(numpy.abs(analytic_acc - finite_diff_acc)), 1e-6)


class TestKeplerianOrbits(unittest.TestCase):
    """`LisaKeplerianOrbit`/`TaijiKeplerianOrbit` implement the second-
    order-in-eccentricity tilted-Kepler-ellipse equal-arm constellation
    (three spacecraft on independent, common-inclination, common-
    eccentricity Kepler ellipses, 120 degrees apart in both mean anomaly
    and ascending node) -- a standard "tilted formation" construction
    (the same physics `lisaorbits.KeplerianOrbits` implements),
    independently implemented here (own code, not ported from any other
    implementation) and validated against it.
    """
    def test_velocity_matches_finite_difference_of_position(self):
        dt = 1.0  # s
        for orbit in (space_orbit.LisaKeplerianOrbit(),
                      space_orbit.TaijiKeplerianOrbit()):
            t = numpy.random.uniform(1e5, 3.1e7, size=10)
            pos_plus = orbit.compute_position(t + dt)
            pos_minus = orbit.compute_position(t - dt)
            finite_diff_vel = (pos_plus - pos_minus) / (2 * dt)
            analytic_vel = orbit.compute_velocity(t)
            self.assertLess(
                numpy.max(numpy.abs(analytic_vel - finite_diff_vel)), 1e-3)

    def test_acceleration_matches_finite_difference_of_velocity(self):
        dt = 1.0  # s
        for orbit in (space_orbit.LisaKeplerianOrbit(),
                      space_orbit.TaijiKeplerianOrbit()):
            t = numpy.random.uniform(1e5, 3.1e7, size=10)
            vel_plus = orbit.compute_velocity(t + dt)
            vel_minus = orbit.compute_velocity(t - dt)
            finite_diff_acc = (vel_plus - vel_minus) / (2 * dt)
            analytic_acc = orbit.compute_acceleration(t)
            self.assertLess(
                numpy.max(numpy.abs(analytic_acc - finite_diff_acc)), 1e-6)

    def test_mean_arm_length_close_to_design_value(self):
        """Unlike `LisaAnalyticOrbit` (whose specific choice of orbital
        phases happens to keep the arm length essentially exactly
        constant at the design value for LISA's small eccentricity --
        confirmed here against real lisaorbits.EqualArmlengthOrbits data,
        spread ~1e-5 m), the true-Kepler-ellipse construction here trades
        that near-perfect constancy for matching a real Keplerian orbit's
        physics: real lisaorbits.KeplerianOrbits itself has a ~0.2%
        arm-length variation over a year for LISA's parameters, and a mean
        arm length that undershoots the nominal design value by a similar
        amount (confirmed against real lisaorbits.KeplerianOrbits data).
        This just checks the mean stays within that same, expected,
        percent-level range of the design value -- not that this
        construction is "more accurate" in an absolute sense, which it is
        not for this quantity.
        """
        t = numpy.linspace(0.0, 3.15e7, 200)
        orbit = space_orbit.LisaKeplerianOrbit(t0=0.0)
        _, length = space_orbit.link_vector(t, orbit, 1, 2)
        self.assertLess(abs(length.mean() - ARMLENGTH), 0.01 * ARMLENGTH)

    def test_usable_as_orbit_provider(self):
        lam, beta, pol = 1.1, -0.2, 0.7
        for orbit in (space_orbit.LisaKeplerianOrbit(),
                      space_orbit.TaijiKeplerianOrbit()):
            t_ssb = 2.0e7
            t_det, lam_det, beta_det, pol_det = space.ssb_to_lisa(
                t_ssb, lam, beta, pol, orbit=orbit)
            t_rt, lam_rt, beta_rt, pol_rt = space.lisa_to_ssb(
                t_det, lam_det, beta_det, pol_det, orbit=orbit)
            self.assertAlmostEqual(t_rt, t_ssb, places=3)
            self.assertAlmostEqual(lam_rt, lam, places=6)
            self.assertAlmostEqual(beta_rt, beta, places=6)
            self.assertAlmostEqual(pol_rt, pol, places=6)

    def test_matches_lisaorbits_keplerian_orbits_exactly(self):
        """Direct, machine-precision cross-check of position/velocity/
        acceleration against a real `lisaorbits.KeplerianOrbits` instance,
        using its own angular frequency (as in the analogous first-order
        test) to isolate the comparison from the unrelated
        EARTH_ORBIT_ANGULAR_FREQUENCY-vs-sqrt(GM_SUN/a**3) difference.
        """
        try:
            import lisaorbits
        except ImportError:
            self.skipTest('lisaorbits not installed; skipping exact '
                          'Keplerian-orbit cross-check')
        ref = lisaorbits.KeplerianOrbits(L=ARMLENGTH)
        ref_ecliptic = space_orbit.ICRSOrbitAdapter(ref)

        t = numpy.random.uniform(1e5, 3.1e7, size=20)
        sc = (1, 2, 3)
        alpha = ref.n * t

        self.assertAlmostEqual(
            space_orbit._kepler_orbit_elements(ARMLENGTH, ref.a)[0], ref.e,
            places=12)

        pos_mine = space_orbit._kepler_orbit_position(
            alpha, ARMLENGTH, ref.a, sc)
        vel_mine = space_orbit._kepler_orbit_velocity(
            alpha, ref.n, ARMLENGTH, ref.a, sc)
        acc_mine = space_orbit._kepler_orbit_acceleration(
            alpha, ref.n, ARMLENGTH, ref.a, sc)

        pos_ref = ref_ecliptic.compute_position(t, sc)
        vel_ref = ref_ecliptic.compute_velocity(t, sc)
        acc_ref = ref_ecliptic.compute_acceleration(t, sc)

        self.assertLess(numpy.max(numpy.abs(pos_mine - pos_ref)), 1e-2)
        self.assertLess(numpy.max(numpy.abs(vel_mine - vel_ref)), 1e-8)
        self.assertLess(numpy.max(numpy.abs(acc_mine - acc_ref)), 1e-14)


class TestNumericOrbitsGeneralizesBeyondLisa(unittest.TestCase):
    """`NumericOrbits` must interpolate Taiji- and TianQin-shaped orbits
    accurately too, not just the LISA special case in `TestNumericOrbits`
    above. TianQin's much shorter (~3.65 day) rotation period requires a
    denser sampling grid than LISA/Taiji's yearly one to interpolate well.
    """
    def test_taiji_interpolation_accuracy(self):
        orbit = AnalyticTaijiOrbit()
        t_grid = numpy.linspace(0.0, 3.15e7, 400)
        numeric = space_orbit.NumericOrbits(
            t_grid, orbit.compute_position(t_grid))
        t_query = numpy.random.uniform(1e5, 3.1e7, size=20)
        interpolated = numeric.compute_position(t_query)
        exact = orbit.compute_position(t_query)
        self.assertLess(numpy.max(numpy.abs(interpolated - exact)),
                        1e-6 * SEMI_MAJOR_AXIS)

    def test_tianqin_interpolation_accuracy(self):
        orbit = AnalyticTianQinOrbit()
        # dense enough to resolve the ~3.65 day rotation period
        t_grid = numpy.linspace(0.0, 3.15e7, 20000)
        numeric = space_orbit.NumericOrbits(
            t_grid, orbit.compute_position(t_grid))
        t_query = numpy.random.uniform(1e5, 3.1e7, size=20)
        interpolated = numeric.compute_position(t_query)
        exact = orbit.compute_position(t_query)
        self.assertLess(numpy.max(numpy.abs(interpolated - exact)),
                        1e-6 * TIANQIN_ARMLENGTH)


class TestSpaceAcceptsOrbitProvider(unittest.TestCase):
    """`pycbc.coordinates.space`'s public, PE-facing functions
    (`t_lisa_from_ssb`, `t_ssb_from_t_lisa`, `ssb_to_lisa`, `lisa_to_ssb`)
    accept an optional `orbit` argument that swaps out the hard-coded
    circular LISA orbit for any orbit provider from `space_orbit`. This
    checks that wiring end to end, both against `space_orbit`'s own
    functions directly and (for the LISA-equivalent fixture) against the
    default, `orbit=None` code path these functions already had.
    """
    def setUp(self):
        self.lisa_orbit = AnalyticEqualArmOrbit()
        self.taiji_orbit = AnalyticTaijiOrbit()
        self.tianqin_orbit = AnalyticTianQinOrbit()
        self.times = numpy.random.uniform(1e5, 3.1e7, size=10)

    def test_t_lisa_from_ssb_delegates_to_space_orbit(self):
        lam, beta = 1.1, -0.2
        k_ssb = space.localization_to_propagation_vector(
            lam, beta, use_astropy=False).flatten()
        for orbit in (self.lisa_orbit, self.taiji_orbit, self.tianqin_orbit):
            for t_ssb in self.times[:5]:
                via_space = space.t_lisa_from_ssb(
                    t_ssb, lam, beta, orbit=orbit)
                via_space_orbit = space_orbit.t_detector_from_ssb(
                    t_ssb, k_ssb, orbit)
                self.assertAlmostEqual(via_space, via_space_orbit, places=6)

    def test_orbit_path_matches_default_for_lisa_equivalent_fixture(self):
        """Feeding `ssb_to_lisa`/`lisa_to_ssb` the analytic orbit that
        `pycbc.coordinates.space`'s hard-coded circular LISA path already
        assumes must reproduce that default (`orbit=None`) path's output,
        not just be self-consistent with itself.
        """
        lam, beta, pol = 0.9, -0.3, 2.4
        for t_ssb in self.times[:5]:
            default = space.ssb_to_lisa(t_ssb, lam, beta, pol, t0=T0)
            via_orbit = space.ssb_to_lisa(
                t_ssb, lam, beta, pol, t0=T0, orbit=self.lisa_orbit)
            for d, o in zip(default, via_orbit):
                self.assertAlmostEqual(d, o, places=3)

    def test_ssb_to_lisa_round_trip_with_taiji_and_tianqin(self):
        """The full ssb -> detector -> ssb round trip must recover the
        original parameters, for constellations genuinely different from
        the LISA special case checked above.
        """
        lam, beta, pol = 2.5, 0.1, 1.0
        for orbit in (self.taiji_orbit, self.tianqin_orbit):
            for t_ssb in self.times[:5]:
                t_det, lam_det, beta_det, pol_det = space.ssb_to_lisa(
                    t_ssb, lam, beta, pol, orbit=orbit)
                t_rt, lam_rt, beta_rt, pol_rt = space.lisa_to_ssb(
                    t_det, lam_det, beta_det, pol_det, orbit=orbit)
                self.assertAlmostEqual(t_rt, t_ssb, places=3)
                self.assertAlmostEqual(lam_rt, lam, places=6)
                self.assertAlmostEqual(beta_rt, beta, places=6)
                self.assertAlmostEqual(pol_rt, pol, places=6)

    def test_lisa_to_geo_and_back_passes_orbit_through(self):
        """`lisa_to_geo`/`geo_to_lisa` compose `lisa_to_ssb`/`ssb_to_lisa`
        with `ssb_to_geo`/`geo_to_ssb`; the `orbit`/`sc` arguments must
        reach the LISA-side leg of that composition, not be silently
        dropped in favor of the hard-coded circular orbit.
        """
        lam, beta, pol = 1.0, 0.2, 0.5
        for orbit in (self.taiji_orbit, self.tianqin_orbit):
            for t_ssb in self.times[:5]:
                lisa_params = space.ssb_to_lisa(
                    t_ssb, lam, beta, pol, orbit=orbit)
                t_geo, lam_geo, beta_geo, pol_geo = space.lisa_to_geo(
                    *lisa_params, orbit=orbit, use_astropy=False)
                lisa_rt = space.geo_to_lisa(
                    t_geo, lam_geo, beta_geo, pol_geo, orbit=orbit,
                    use_astropy=False)
                for a, b in zip(lisa_params, lisa_rt):
                    self.assertAlmostEqual(a, b, places=3)


class TestTransformsAcceptOrbitFile(unittest.TestCase):
    """`pycbc.transforms`' SSB<->LISA<->GEO transform plugins
    (`ssb_to_lisa`/`lisa_to_ssb`/`lisa_to_geo`/`geo_to_lisa`) accept an
    optional `orbit-file` config option (an HDF5 file readable by
    `space_orbit.NumericOrbits.from_file`) so a PE config can select any
    constellation orbit (LISA, Taiji, TianQin, numerical, ...), not just the
    default hard-coded circular LISA orbit. This checks that option reaches
    `pycbc.coordinates.space`'s already-generalized functions correctly,
    both when constructing the transform class directly and via
    `from_config`.
    """
    @classmethod
    def setUpClass(cls):
        orbit = AnalyticTianQinOrbit()
        t_grid = numpy.linspace(0.0, 3.15e7, 20000)
        positions = orbit.compute_position(t_grid)
        cls.tmpdir = tempfile.mkdtemp()
        cls.orbit_file = os.path.join(cls.tmpdir, 'tianqin_orbit.hdf5')
        with h5py.File(cls.orbit_file, 'w') as f:
            f.create_dataset('t', data=t_grid)
            f.create_dataset('positions', data=positions)
        cls.numeric_orbit = space_orbit.NumericOrbits.from_file(
            cls.orbit_file)

    def test_ssb_to_lisa_orbit_file_matches_direct_orbit(self):
        t = transforms.SSBToLISA(orbit_file=self.orbit_file)
        out = t.transform({
            'tc': 1.5e7, 'eclipticlongitude': 1.0,
            'eclipticlatitude': 0.2, 'polarization': 0.5})
        expected = space.ssb_to_lisa(
            1.5e7, 1.0, 0.2, 0.5, orbit=self.numeric_orbit)
        self.assertAlmostEqual(out['tc'], expected[0], places=3)
        self.assertAlmostEqual(
            out['eclipticlongitude'], expected[1], places=6)

    def test_lisa_to_geo_orbit_file_matches_direct_orbit(self):
        t = transforms.LISAToGEO(orbit_file=self.orbit_file)
        out = t.transform({
            'tc': 1.5e7, 'eclipticlongitude': 1.0,
            'eclipticlatitude': 0.2, 'polarization': 0.5})
        expected = space.lisa_to_geo(
            1.5e7, 1.0, 0.2, 0.5, orbit=self.numeric_orbit,
            use_astropy=True)
        self.assertAlmostEqual(out['tc'], expected[0], places=3)

    def test_default_orbit_file_none_matches_hardcoded_lisa(self):
        """With no orbit-file, behavior must be unchanged from before this
        option existed.
        """
        t = transforms.SSBToLISA()
        out = t.transform({
            'tc': 1.5e7, 'eclipticlongitude': 1.0,
            'eclipticlatitude': 0.2, 'polarization': 0.5})
        expected = space.ssb_to_lisa(1.5e7, 1.0, 0.2, 0.5)
        self.assertAlmostEqual(out['tc'], expected[0], places=6)

    def test_from_config_parses_orbit_file(self):
        from pycbc.workflow.configuration import WorkflowConfigParser
        ini_text = f"""
[waveform_transforms-tc_lisa+eclipticlongitude_lisa+eclipticlatitude_lisa+polarization_lisa]
name = ssb_to_lisa
tc-ssb = tc
longitude-ssb = eclipticlongitude
latitude-ssb = eclipticlatitude
polarization-ssb = polarization
tc-lisa = tc_lisa
longitude-lisa = eclipticlongitude_lisa
latitude-lisa = eclipticlatitude_lisa
polarization-lisa = polarization_lisa
orbit-file = {self.orbit_file}
"""
        ini_path = os.path.join(self.tmpdir, 'from_config_test.ini')
        with open(ini_path, 'w') as f:
            f.write(ini_text)
        cp = WorkflowConfigParser([ini_path])
        t = transforms.SSBToLISA.from_config(
            cp, 'waveform_transforms',
            'tc_lisa+eclipticlongitude_lisa+eclipticlatitude_lisa+'
            'polarization_lisa')
        self.assertEqual(t.orbit_file, self.orbit_file)
        out = t.transform({
            'tc': 1.5e7, 'eclipticlongitude': 1.0,
            'eclipticlatitude': 0.2, 'polarization': 0.5})
        expected = space.ssb_to_lisa(
            1.5e7, 1.0, 0.2, 0.5, orbit=self.numeric_orbit)
        self.assertAlmostEqual(out['tc_lisa'], expected[0], places=3)


class TestOptionalLisaorbitsDuckTyping(unittest.TestCase):
    """If `lisaorbits` happens to be installed in the test environment,
    verify that a real `lisaorbits.Orbits` instance can be passed directly
    to `constellation_frame` with no adapter code (duck typing), and that
    the result is physically sane (finite, non-degenerate). This test is
    skipped entirely if `lisaorbits` is not installed; `space_orbit` itself
    never imports it.
    """
    def test_lisaorbits_instance_accepted_directly(self):
        try:
            import lisaorbits
        except ImportError:
            self.skipTest('lisaorbits not installed; skipping duck-typing '
                          'cross-check')
        orbit = lisaorbits.EqualArmlengthOrbits()
        times = numpy.array([1e6, 1e7])
        centroid, rotation = space_orbit.constellation_frame(times, orbit)
        self.assertEqual(centroid.shape, (2, 3))
        self.assertEqual(rotation.shape, (2, 3, 3))
        self.assertTrue(numpy.all(numpy.isfinite(centroid)))
        self.assertTrue(numpy.all(numpy.isfinite(rotation)))
        # rotation matrices must be orthogonal
        for r in rotation:
            self.assertLess(
                numpy.max(numpy.abs(r @ r.T - numpy.eye(3))), 1e-8)

    def test_position_velocity_acceleration_match_lisaorbits_exactly(self):
        """`space_orbit`'s analytic position/velocity/acceleration formulas
        (`_equal_arm_orbit_position/_velocity/_acceleration`, underlying
        `LisaAnalyticOrbit`/`TaijiAnalyticOrbit`) are, by construction, the
        exact same closed-form Keplerian expansion (Rubbo, Cornish &
        Poujade 2004) that `lisaorbits.EqualArmlengthOrbits` implements --
        this checks that claim directly, at essentially machine precision,
        rather than merely "close".

        `EqualArmlengthOrbits()`'s defaults give a guiding-center phase
        `mbar(t) = n * t` (`t_init`, `m_init1`, `lambda1` all 0), where `n`
        is its own physically-derived angular frequency
        (`sqrt(GM_SUN/a**3)`) -- not `space_orbit`'s own
        `EARTH_ORBIT_ANGULAR_FREQUENCY` constant (tuned to Earth's actual
        sidereal year, which differs slightly from a bare two-body Kepler
        frequency at 1 AU). Since `_equal_arm_orbit_velocity`/
        `_acceleration` take the angular frequency as an explicit
        parameter (not hardcoded), using `EqualArmlengthOrbits`'s own `n`
        directly here isolates a pure formula-for-formula comparison from
        that unrelated frequency-convention difference.
        """
        try:
            import lisaorbits
        except ImportError:
            self.skipTest('lisaorbits not installed; skipping exact '
                          'position/velocity/acceleration cross-check')
        ref = lisaorbits.EqualArmlengthOrbits(L=ARMLENGTH)
        ref_ecliptic = space_orbit.ICRSOrbitAdapter(ref)

        t = numpy.random.uniform(1e5, 3.1e7, size=20)
        sc = (1, 2, 3)
        alpha = ref.n * t  # matches ref's own mbar(t) for these defaults

        pos_mine = space_orbit._equal_arm_orbit_position(alpha, ARMLENGTH, sc)
        vel_mine = space_orbit._equal_arm_orbit_velocity(
            alpha, ref.n, ARMLENGTH, sc)
        acc_mine = space_orbit._equal_arm_orbit_acceleration(
            alpha, ref.n, ARMLENGTH, sc)

        pos_ref = ref_ecliptic.compute_position(t, sc)
        vel_ref = ref_ecliptic.compute_velocity(t, sc)
        acc_ref = ref_ecliptic.compute_acceleration(t, sc)

        self.assertLess(numpy.max(numpy.abs(pos_mine - pos_ref)), 1e-3)
        self.assertLess(numpy.max(numpy.abs(vel_mine - vel_ref)), 1e-6)
        self.assertLess(numpy.max(numpy.abs(acc_mine - acc_ref)), 1e-9)

    def test_numeric_orbits_matches_interpolated_orbits_on_identical_data(self):
        """`NumericOrbits` and `lisaorbits.InterpolatedOrbits` both build a
        quintic spline (`scipy.interpolate.make_interp_spline`, `k=5` by
        default in both) over the same kind of (times, positions[,
        velocities]) input. This checks that, fed *exactly* the same
        sampled data, the two produce matching position/velocity/
        acceleration at off-grid query times -- the numeric-orbit
        counterpart of `test_position_velocity_acceleration_match_
        lisaorbits_exactly` above (which only covered the analytic-orbit
        path).
        """
        try:
            import lisaorbits
        except ImportError:
            self.skipTest('lisaorbits not installed; skipping numeric '
                          'interpolation cross-check')
        exact_orbit = space_orbit.LisaAnalyticOrbit()
        t_grid = numpy.linspace(0.0, 3.15e7, 400)
        positions = exact_orbit.compute_position(t_grid)

        pycbc_numeric = space_orbit.NumericOrbits(t_grid, positions)
        ref_numeric = lisaorbits.InterpolatedOrbits(
            t_grid, positions, t_init=t_grid[10])

        query_t = numpy.random.uniform(t_grid[20], t_grid[-20], size=30)
        sc = (1, 2, 3)
        for method, tol in [('compute_position', 1e-3),
                             ('compute_velocity', 1e-6),
                             ('compute_acceleration', 1e-6)]:
            mine = getattr(pycbc_numeric, method)(query_t, sc)
            ref = getattr(ref_numeric, method)(query_t, sc)
            self.assertLess(
                numpy.max(numpy.abs(mine - ref)), tol,
                f'{method} disagrees with lisaorbits.InterpolatedOrbits '
                f'by more than {tol}')


class TestNumericOrbitsFileReaders(unittest.TestCase):
    """`NumericOrbits.from_oem_files`/`from_triangle_dat_files` read two
    real-world numeric orbit formats natively (no `lisaorbits`/`oem`
    dependency). These tests write small, synthetic fixture files (a
    known analytic orbit, re-encoded in each format) so they never depend
    on a local, machine-specific data checkout, and check the round trip
    recovers the original orbit.
    """
    def setUp(self):
        self.exact_orbit = space_orbit.LisaAnalyticOrbit()
        self.t_grid = numpy.arange(0.0, 3.15e7, 86400.0)
        self.tmpdir = tempfile.mkdtemp()

    def test_from_triangle_dat_files_round_trip(self):
        positions = self.exact_orbit.compute_position(self.t_grid)
        velocities = self.exact_orbit.compute_velocity(self.t_grid)
        for k, label in enumerate(('1', '2', '3')):
            numpy.savetxt(
                os.path.join(self.tmpdir, f'SCP{label}.dat'),
                positions[:, k, :] / au.value)
            numpy.savetxt(
                os.path.join(self.tmpdir, f'SCV{label}.dat'),
                velocities[:, k, :] / au.value * 86400.0)

        orbit = space_orbit.NumericOrbits.from_triangle_dat_files(
            self.tmpdir, t0=self.t_grid[0], dt=86400.0)
        query_t = numpy.linspace(
            self.t_grid[5], self.t_grid[-5], 30)
        recovered = orbit.compute_position(query_t)
        exact = self.exact_orbit.compute_position(query_t)
        self.assertLess(numpy.max(numpy.abs(recovered - exact)), 1.0)

    def test_from_triangle_dat_files_rejects_mismatched_row_counts(self):
        positions = self.exact_orbit.compute_position(self.t_grid)
        for k, label in enumerate(('1', '2', '3')):
            n_rows = len(self.t_grid) - (1 if label == '3' else 0)
            numpy.savetxt(
                os.path.join(self.tmpdir, f'SCP{label}.dat'),
                positions[:n_rows, k, :] / au.value)
            numpy.savetxt(
                os.path.join(self.tmpdir, f'SCV{label}.dat'),
                positions[:n_rows, k, :] / au.value)
        with self.assertRaises(ValueError):
            space_orbit.NumericOrbits.from_triangle_dat_files(self.tmpdir)

    def _write_oem_file(self, path, positions_icrs_km, velocities_icrs_km_s,
                        epochs_iso, ref_frame='EME2000',
                        time_system='TCB'):
        with open(path, 'w') as f:
            f.write('CCSDS_OEM_VERS = 2.0\n')
            f.write('CREATION_DATE  = 2024-01-01T00:00:00\n')
            f.write('ORIGINATOR     = pycbc-test\n\n')
            f.write('META_START\n')
            f.write('OBJECT_NAME          = TEST\n')
            f.write('OBJECT_ID            = -1\n')
            f.write('CENTER_NAME          = SUN\n')
            f.write(f'REF_FRAME            = {ref_frame}\n')
            f.write(f'TIME_SYSTEM          = {time_system}\n')
            f.write(f'START_TIME           = {epochs_iso[0]}\n')
            f.write(f'STOP_TIME            = {epochs_iso[-1]}\n')
            f.write('INTERPOLATION        = HERMITE\n')
            f.write('INTERPOLATION_DEGREE = 7\n')
            f.write('META_STOP\n\n')
            for epoch, pos, vel in zip(
                    epochs_iso, positions_icrs_km, velocities_icrs_km_s):
                f.write(
                    f'{epoch} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f} '
                    f'{vel[0]:.9f} {vel[1]:.9f} {vel[2]:.9f}\n')

    def test_from_oem_files_round_trip(self):
        from astropy.time import Time

        positions_ecl = self.exact_orbit.compute_position(self.t_grid)
        velocities_ecl = self.exact_orbit.compute_velocity(self.t_grid)
        rotation = space_orbit._icrs_to_ecliptic_rotation_matrix()
        # ecliptic -> ICRS is the inverse (transpose) of the cached
        # ICRS -> ecliptic rotation used throughout this module.
        positions_icrs = (
            positions_ecl.reshape(-1, 3) @ rotation
        ).reshape(positions_ecl.shape) / 1e3  # m -> km
        velocities_icrs = (
            velocities_ecl.reshape(-1, 3) @ rotation
        ).reshape(velocities_ecl.shape) / 1e3  # m/s -> km/s

        epochs_iso = Time(self.t_grid, format='gps').tcb.isot
        oem_paths = []
        for k, label in enumerate(('1', '2', '3')):
            path = os.path.join(self.tmpdir, f'test_sc{label}.oem')
            self._write_oem_file(
                path, positions_icrs[:, k, :], velocities_icrs[:, k, :],
                epochs_iso)
            oem_paths.append(path)

        orbit = space_orbit.NumericOrbits.from_oem_files(*oem_paths)
        query_t = numpy.linspace(self.t_grid[5], self.t_grid[-5], 30)
        recovered = orbit.compute_position(query_t)
        exact = self.exact_orbit.compute_position(query_t)
        # Dominated by ordinary quintic-spline interpolation error at
        # 1-day cadence (consistent with the ~20 m level already measured
        # elsewhere in this test suite for similarly-sampled data), not by
        # the ISO8601/km-precision text round trip, which is far smaller.
        self.assertLess(numpy.max(numpy.abs(recovered - exact)), 20.0)

    def test_from_oem_files_rejects_non_eme2000_frame(self):
        from astropy.time import Time
        positions_ecl = self.exact_orbit.compute_position(self.t_grid)
        epochs_iso = Time(self.t_grid, format='gps').tcb.isot
        oem_paths = []
        for k, label in enumerate(('1', '2', '3')):
            path = os.path.join(self.tmpdir, f'test_sc{label}.oem')
            self._write_oem_file(
                path, positions_ecl[:, k, :] / 1e3,
                positions_ecl[:, k, :] / 1e3, epochs_iso,
                ref_frame='ICRF')
            oem_paths.append(path)
        with self.assertRaises(ValueError):
            space_orbit.NumericOrbits.from_oem_files(*oem_paths)

    def test_from_oem_files_rejects_mismatched_epochs(self):
        from astropy.time import Time
        positions_ecl = self.exact_orbit.compute_position(self.t_grid)
        epochs_iso = Time(self.t_grid, format='gps').tcb.isot
        epochs_iso_shifted = Time(
            self.t_grid + 1.0, format='gps').tcb.isot
        oem_paths = []
        for k, label in enumerate(('1', '2', '3')):
            path = os.path.join(self.tmpdir, f'test_sc{label}.oem')
            these_epochs = epochs_iso_shifted if label == '3' else epochs_iso
            self._write_oem_file(
                path, positions_ecl[:, k, :] / 1e3,
                positions_ecl[:, k, :] / 1e3, these_epochs)
            oem_paths.append(path)
        with self.assertRaises(ValueError):
            space_orbit.NumericOrbits.from_oem_files(*oem_paths)

    def test_from_oem_files_matches_lisaorbits_on_real_esa_data(self):
        """Optional, network-gated real-data cross-check: fetch the
        official ESA lisa-orbit-files trio via lisaorbits' own pooch-based
        download/cache (not a local, machine-specific path), then verify
        the native from_oem_files parser agrees with a real
        lisaorbits.OEMOrbits instance reading the exact same files.
        Skipped entirely if lisaorbits is not installed or the fetch
        fails for any reason.
        """
        try:
            import lisaorbits
            from lisaorbits.oem import OEMOrbits
            import lisaconstants
            from astropy.time import Time
        except ImportError:
            self.skipTest('lisaorbits not installed; skipping real-data '
                          'OEM reader cross-check')
        try:
            paths = OEMOrbits._get_included_paths(
                'esa-trailing', version='2.0.0')
            oem_paths = [
                lisaorbits.oem._ESA_REPOSITORY.fetch(p) for p in paths]
        except Exception as exc:  # pylint: disable=broad-except
            self.skipTest('could not fetch official ESA lisa-orbit-files '
                          f'orbit ({exc!r}); skipping cross-check')

        native = space_orbit.NumericOrbits.from_oem_files(*oem_paths)
        ref = OEMOrbits(*oem_paths)
        ref_ecliptic = space_orbit.ICRSOrbitAdapter(ref)

        query_gps = numpy.linspace(
            native.t_interp[10], native.t_interp[-10], 30)
        query_time = Time(query_gps, format='gps')
        epoch = Time(lisaconstants.LISA_EPOCH_TCB, format='isot',
                    scale='tcb')
        query_tcb_sec = (query_time.tcb - epoch).sec

        pos_native = native.compute_position(query_gps)
        pos_ref = ref_ecliptic.compute_position(query_tcb_sec)
        self.assertLess(numpy.max(numpy.abs(pos_native - pos_ref)), 1.0)


class TestICRSOrbitAdapter(unittest.TestCase):
    """`ICRSOrbitAdapter` wraps an ICRS-frame orbit provider (as produced by
    real spacecraft ephemerides, e.g. CCSDS OEM files) so it can be used as
    a drop-in `OrbitProvider` in pycbc's own SSB-ecliptic convention. These
    tests build a synthetic "ICRS-frame" fixture by rotating one of this
    module's own ecliptic-frame analytic orbits with the *inverse* of the
    adapter's own rotation matrix, then check the adapter recovers the
    original ecliptic-frame values -- a self-contained round trip that
    does not depend on `lisaorbits` or any external orbit file.
    """
    def setUp(self):
        self.ecliptic_orbit = space_orbit.LisaAnalyticOrbit()
        self.rotation = space_orbit._icrs_to_ecliptic_rotation_matrix()
        self.times = numpy.linspace(1e6, 3.15e7, 50)

    def _fake_icrs_orbit(self):
        ecliptic_orbit = self.ecliptic_orbit
        rotation = self.rotation

        class _FakeICRSOrbit:
            def compute_position(self, t, sc=(1, 2, 3)):
                pos_ecliptic = ecliptic_orbit.compute_position(t, sc)
                return pos_ecliptic @ rotation  # ecliptic -> ICRS

        return _FakeICRSOrbit()

    def test_round_trip_recovers_ecliptic_positions(self):
        adapter = space_orbit.ICRSOrbitAdapter(self._fake_icrs_orbit())
        recovered = adapter.compute_position(self.times, (1, 2, 3))
        expected = self.ecliptic_orbit.compute_position(self.times, (1, 2, 3))
        self.assertLess(numpy.max(numpy.abs(recovered - expected)), 1e-3)

    def test_rotation_is_orthogonal(self):
        r = self.rotation
        self.assertLess(numpy.max(numpy.abs(r @ r.T - numpy.eye(3))), 1e-10)

    def test_rotation_matches_independent_astropy_transform(self):
        """A round trip through the adapter alone cannot catch a sign or
        transpose error in the cached rotation matrix (composing a matrix
        with its own transpose is the identity regardless of whether the
        matrix itself is correct). This test instead compares the cached
        matrix's action on an arbitrary, non-basis vector against a fresh,
        independent call to astropy's ICRS -> BarycentricMeanEcliptic
        transform, so a transpose/sign bug in the matrix construction
        cannot hide behind a self-consistent-but-wrong round trip.
        """
        from astropy import units as apy_units
        from astropy.coordinates import ICRS, BarycentricMeanEcliptic

        v_icrs = numpy.array([1.3, -0.7, 2.1])
        icrs = ICRS(x=v_icrs[0] * apy_units.m, y=v_icrs[1] * apy_units.m,
                    z=v_icrs[2] * apy_units.m, representation_type='cartesian')
        ecl = icrs.transform_to(
            BarycentricMeanEcliptic(equinox='J2000')).cartesian
        expected = numpy.array([ecl.x.to(apy_units.m).value,
                                 ecl.y.to(apy_units.m).value,
                                 ecl.z.to(apy_units.m).value])

        result = self.rotation @ v_icrs
        self.assertLess(numpy.max(numpy.abs(result - expected)), 1e-8)

    def test_compute_velocity_rotates_consistently(self):
        t_grid = numpy.linspace(0.0, 3.15e7, 400)
        pos_ecliptic = self.ecliptic_orbit.compute_position(t_grid)
        numeric_ecliptic = space_orbit.NumericOrbits(t_grid, pos_ecliptic)

        rotation = self.rotation

        class _FakeICRSOrbit:
            def compute_position(self, t, sc=(1, 2, 3)):
                return numeric_ecliptic.compute_position(t, sc) @ rotation

            def compute_velocity(self, t, sc=(1, 2, 3)):
                return numeric_ecliptic.compute_velocity(t, sc) @ rotation

        adapter = space_orbit.ICRSOrbitAdapter(_FakeICRSOrbit())
        query_t = numpy.linspace(1e6, 3.1e7, 50)
        recovered_vel = adapter.compute_velocity(query_t, (1, 2, 3))
        expected_vel = numeric_ecliptic.compute_velocity(query_t, (1, 2, 3))
        self.assertLess(
            numpy.max(numpy.abs(recovered_vel - expected_vel)), 1e-3)

    def test_compute_acceleration_rotates_consistently(self):
        t_grid = numpy.linspace(0.0, 3.15e7, 400)
        pos_ecliptic = self.ecliptic_orbit.compute_position(t_grid)
        numeric_ecliptic = space_orbit.NumericOrbits(t_grid, pos_ecliptic)

        rotation = self.rotation

        class _FakeICRSOrbit:
            def compute_position(self, t, sc=(1, 2, 3)):
                return numeric_ecliptic.compute_position(t, sc) @ rotation

            def compute_velocity(self, t, sc=(1, 2, 3)):
                return numeric_ecliptic.compute_velocity(t, sc) @ rotation

            def compute_acceleration(self, t, sc=(1, 2, 3)):
                return numeric_ecliptic.compute_acceleration(t, sc) @ rotation

        adapter = space_orbit.ICRSOrbitAdapter(_FakeICRSOrbit())
        query_t = numpy.linspace(1e6, 3.1e7, 50)
        recovered_acc = adapter.compute_acceleration(query_t, (1, 2, 3))
        expected_acc = numeric_ecliptic.compute_acceleration(query_t, (1, 2, 3))
        self.assertLess(
            numpy.max(numpy.abs(recovered_acc - expected_acc)), 1e-6)

    def test_usable_directly_with_constellation_frame_and_link_vector(self):
        adapter = space_orbit.ICRSOrbitAdapter(self._fake_icrs_orbit())
        centroid, rotation = space_orbit.constellation_frame(
            self.times, adapter)
        expected_centroid, expected_rotation = space_orbit.constellation_frame(
            self.times, self.ecliptic_orbit)
        self.assertLess(
            numpy.max(numpy.abs(centroid - expected_centroid)), 1e-3)
        self.assertLess(
            numpy.max(numpy.abs(rotation - expected_rotation)), 1e-8)

        _, length = space_orbit.link_vector(self.times, adapter, 1, 2)
        _, expected_length = space_orbit.link_vector(
            self.times, self.ecliptic_orbit, 1, 2)
        self.assertLess(numpy.max(numpy.abs(length - expected_length)), 1e-3)


class TestOptionalESAOemOrbitFiles(unittest.TestCase):
    """If `lisaorbits` (and its optional `oem` dependency) are installed and
    network access is available, verify that `space_orbit` correctly
    consumes real ESA LISA orbit files (CCSDS OEM format, published at
    https://github.com/esa/lisa-orbit-files, fetched here via
    `lisaorbits.OEMOrbits.from_included`).

    Unlike `lisaorbits`' own analytic orbits (e.g. `EqualArmlengthOrbits`,
    used in `TestOptionalLisaorbitsDuckTyping` above), OEM files are given in
    the EME2000 (~ICRS) equatorial frame, not the SSB ecliptic frame that
    `space_orbit`/`space` assume. This test therefore also exercises
    `ICRSOrbitAdapter`, which any caller must use (or replicate) before
    treating OEM-derived positions as a `space_orbit` `OrbitProvider`.

    Skipped entirely if `lisaorbits`/`oem` are not installed, or if the
    one-time download of the (~600 kB total) orbit files fails for any
    reason (no network access, hash mismatch, upstream repository changes,
    etc.) -- this test never requires a local, machine-specific data path.
    """
    def test_esa_oem_orbit_matches_design_arm_length(self):
        try:
            import lisaorbits
        except ImportError:
            self.skipTest('lisaorbits not installed; skipping ESA OEM '
                          'orbit-file cross-check')
        try:
            oem_orbit = lisaorbits.OEMOrbits.from_included(
                'esa-trailing', version='2.0.0')
        except Exception as exc:  # pylint: disable=broad-except
            self.skipTest('could not fetch official ESA lisa-orbit-files '
                          f'orbit ({exc!r}); skipping cross-check')

        t_grid = numpy.linspace(
            oem_orbit.t_start + 5 * 86400.0,
            oem_orbit.t_end - 5 * 86400.0, 100)
        esa_orbit = space_orbit.ICRSOrbitAdapter(oem_orbit)
        unit_vec, length = space_orbit.link_vector(t_grid, esa_orbit, 1, 2)
        self.assertTrue(numpy.all(numpy.isfinite(length)))
        self.assertLess(
            numpy.max(numpy.abs(numpy.linalg.norm(unit_vec, axis=-1) - 1.0)),
            1e-8)
        # LISA's design arm length is 2.5 million km; real orbit files
        # deviate from this by up to a few percent due to orbital dynamics.
        self.assertTrue(numpy.all(numpy.abs(length / 1e3 - 2.5e6) < 1e5))

        centroid, rotation = space_orbit.constellation_frame(t_grid, esa_orbit)
        for r in rotation:
            self.assertLess(
                numpy.max(numpy.abs(r @ r.T - numpy.eye(3))), 1e-8)
        centroid_au = numpy.linalg.norm(centroid, axis=-1) / SEMI_MAJOR_AXIS
        self.assertTrue(numpy.all(numpy.abs(centroid_au - 1.0) < 0.05))

    def test_ssb_to_lisa_round_trip_with_real_esa_orbit(self):
        """The lower-level geometric checks above (arm length,
        constellation_frame orthonormality) do not exercise
        `pycbc.coordinates.space`'s actual PE-facing transforms
        (`ssb_to_lisa`/`lisa_to_ssb`, sky-localization + polarization
        angle). This checks those directly against a real numeric orbit,
        not just the analytic toy fixtures used in
        `TestSpaceAcceptsOrbitProvider`.
        """
        try:
            import lisaorbits
        except ImportError:
            self.skipTest('lisaorbits not installed; skipping ESA OEM '
                          'orbit-file cross-check')
        try:
            oem_orbit = lisaorbits.OEMOrbits.from_included(
                'esa-trailing', version='2.0.0')
        except Exception as exc:  # pylint: disable=broad-except
            self.skipTest('could not fetch official ESA lisa-orbit-files '
                          f'orbit ({exc!r}); skipping cross-check')

        esa_orbit = space_orbit.ICRSOrbitAdapter(oem_orbit)
        # Well inside the valid range, away from the interpolation edges.
        t_ssb = (oem_orbit.t_start + oem_orbit.t_end) / 2.0
        lam, beta, pol = 0.9, -0.3, 2.4

        t_det, lam_det, beta_det, pol_det = space.ssb_to_lisa(
            t_ssb, lam, beta, pol, orbit=esa_orbit)
        self.assertTrue(numpy.isfinite(t_det))
        t_rt, lam_rt, beta_rt, pol_rt = space.lisa_to_ssb(
            t_det, lam_det, beta_det, pol_det, orbit=esa_orbit)
        self.assertAlmostEqual(t_rt, t_ssb, places=3)
        self.assertAlmostEqual(lam_rt, lam, places=6)
        self.assertAlmostEqual(beta_rt, beta, places=6)
        self.assertAlmostEqual(pol_rt, pol, places=6)


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestConstellationFrame))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestLinkVector))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestNumericOrbits))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestLisaOrbit))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestTaijiOrbit))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestTianQinOrbit))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestProductionAnalyticOrbits))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestKeplerianOrbits))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestNumericOrbitsGeneralizesBeyondLisa))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestSpaceAcceptsOrbitProvider))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestTransformsAcceptOrbitFile))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestNumericOrbitsFileReaders))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestICRSOrbitAdapter))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestOptionalLisaorbitsDuckTyping))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestOptionalESAOemOrbitFiles))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
