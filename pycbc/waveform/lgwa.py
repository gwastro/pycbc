# Copyright (C) 2026  Shichao Wu, Alex Nitz, Jacopo Tissino, Jan Harms
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
Frequency-domain, stationary-phase-approximation (SPA) response for the
LGWA (Lunar Gravitational Wave Antenna) detector, registered as the
`'LGWA_response'` `pycbc.waveform.fd_det`/`fd_det_sequence` plugin -- this
is what makes `pycbc.inference.models.relbin.Relative` able to use LGWA
data directly (see `relbin.py`'s `still_needs_det_response` dispatch,
driven entirely by `approximant in fd_det_sequence`, with zero LGWA-
specific code needed there).

Why this is a *different* code path from `pycbc.detector.space.
_LGWA_detector.project_wave` (the time-domain response used elsewhere):
relative binning calls this with sparse, non-uniform frequency points
(the bin edges), and realistic LGWA sources (mHz band) can have physical
durations of thousands of days -- both make "generate a time-domain
waveform, call `project_wave`, then FFT" intractable in general. Instead:

1. The intrinsic waveform (`hp`, `hc`) is generated directly at the
   requested frequency points via `get_fd_waveform_sequence`, for
   *any* PyCBC/LAL frequency-domain approximant (`base_approximant`),
   mode-by-mode (`mode_array`) rather than assuming a single dominant
   harmonic -- different (l, m) modes/eccentric harmonics correspond to
   different true emission times at the same frequency, so each mode's
   antenna-pattern factor is evaluated at its own arrival time.
2. Each mode's frequency-to-time relation `t_lm(f)` (stationary-phase
   time-of-arrival, relative to the coalescence reference time) is
   estimated by a genuinely *local* finite difference at each frequency
   point (step size adaptively scaled so the phase change stays well
   under 2*pi), avoiding `numpy.unwrap`, which silently breaks at the low
   frequencies typical of LGWA sources (see the discussion in
   `notebooks/lgwa_response_validation.ipynb`).
3. The detector-orientation antenna pattern is evaluated at each mode's
   own per-frequency arrival time, using the same (cached)
   `_LGWA_detector._detector_frame`/`_antenna_pattern_factors` machinery
   `project_wave` uses -- reusing exactly the same, already bit-exact-
   validated formula (see `test/test_detector_space.py`), just evaluated
   at SPA times instead of uniform time samples.

`ref_position`: the arrival-time reference point. `tc` is interpreted as
the arrival time at this *fixed* (SSB-comoving) point, not necessarily
the SSB origin -- see `pycbc.coordinates.reference_point`'s module
docstring (Tissino et al. 2026 Sec. 3.2) for why a configurable reference
point matters for multiband PE sampling efficiency, and why a fixed point
keeps the conversion to/from SSB a cheap closed-form correction.

Known limitation, not yet resolved (tracked for follow-up rather than
silently ignored): `mode_array` is a spherical-harmonic (l, m) mode
decomposition, the standard way PyCBC/LALSuite expose higher-order-mode
content for quasi-circular approximants (e.g. IMRPhenomXHM). Eccentric
waveform models parameterize their harmonic content differently (and not
always through the same `mode_array` mechanism), so full eccentric-
harmonic support may need model-specific handling investigated separately.
"""

import numpy as np

from pycbc.coordinates import moon as _moon
from pycbc.coordinates import reference_point as _refpt
from pycbc.detector.space import _LGWA_detector
from pycbc.types import Array, FrequencySeries
from pycbc.waveform.waveform import get_fd_waveform_sequence

__all__ = ["lgwa_fd_response"]

# G * M_sun / c**3, in seconds -- used only to pick a safe (adaptive)
# finite-difference step size for the local stationary-phase time
# estimate below, via the leading-order Newtonian chirp-time relation.
# Not assumed to be an accurate estimate of the true (l, m)-mode time of
# arrival itself -- that comes from the real waveform's own phase.
_SUN_MASS_SECONDS = 4.925490947e-6

# Module-level cache of _LGWA_detector instances, keyed by the static
# parameters that fully determine its (expensive-to-build, but here
# reused) detector-orientation grid: (longitude_site, latitude_site,
# cadence). This ensures repeated calls within the same process (e.g.
# once per relbin likelihood evaluation) reuse both the same detector
# object and, in turn, its internal `_frame_cache` -- mirroring bbhx's
# own `cache_generator=True` pattern for the analogous LISA plugin.
_lgwa_detector_cache = {}


def _get_cached_detector(longitude_site, latitude_site, cadence):
    key = (float(longitude_site), float(latitude_site), float(cadence))
    det = _lgwa_detector_cache.get(key)
    if det is None:
        det = _LGWA_detector(
            "LGWA",
            longitude_site=longitude_site,
            latitude_site=latitude_site,
            cadence=cadence,
        )
        _lgwa_detector_cache[key] = det
    return det


def _time_of_f_local(fm, base_approximant, mode, cbc_params):
    """
    Per-frequency-point, unwrap-free stationary-phase time estimate
    `t_lm(f)` for a single (l, m) mode, relative to the model's own
    coalescence-time reference (negative before merger, matching
    `lgwa_response.simple_waveforms.time_to_merger`'s convention) --
    generalized from the technique validated in
    `notebooks/lgwa_response_validation.ipynb`'s full-LGWA-band
    comparison, but batched (one extra `get_fd_waveform_sequence` call
    per +/- offset per retry, not one call per frequency point) and
    generalized to an arbitrary mode rather than assuming the (2,2)
    harmonic.

    The step size is scaled using the leading-order Newtonian chirp-time
    relation (only to pick a step small enough that the phase change
    across it stays well under 2*pi -- not assumed to be an accurate
    estimate of the true, possibly precessing/eccentric, t(f) itself),
    with the mode's own frequency first mapped to an equivalent (2,2)
    frequency (`f * 2/m`, since the (l, m) harmonic's GW frequency is
    approximately m/2 times the orbital-derived (2,2) frequency).
    """
    ell, m = mode
    fm = np.asarray(fm, dtype=float)
    mass1 = float(cbc_params["mass1"])
    mass2 = float(cbc_params["mass2"])
    mchirp_sec = (mass1 * mass2) ** 0.6 / (mass1 + mass2) ** 0.2 * _SUN_MASS_SECONDS

    f_22_equiv = np.abs(fm) * 2.0 / max(abs(m), 1)
    f_22_equiv = np.maximum(f_22_equiv, 1e-8)
    t_est = (
        5.0
        / (256.0 * np.pi ** (8.0 / 3.0))
        / mchirp_sec ** (5.0 / 3.0)
        / f_22_equiv ** (8.0 / 3.0)
    )
    eps = np.minimum(0.01 / (2 * np.pi * t_est), fm * 1e-4)
    eps = np.maximum(eps, fm * 1e-10)

    dphase = None
    for _ in range(15):
        hp_minus, _ = get_fd_waveform_sequence(
            approximant=base_approximant,
            mode_array=[[ell, m]],
            sample_points=Array(fm - eps),
            **cbc_params,
        )
        hp_plus, _ = get_fd_waveform_sequence(
            approximant=base_approximant,
            mode_array=[[ell, m]],
            sample_points=Array(fm + eps),
            **cbc_params,
        )
        dphase = np.angle(hp_plus.numpy() * np.conj(hp_minus.numpy()))
        bad = np.abs(dphase) > np.pi / 2
        if not np.any(bad):
            break
        eps = np.where(bad, eps / 2, eps)

    # LALSimulation/pycbc's h(f) = A(f) * exp(-i*Psi(f)) convention has
    # the opposite phase sign from the "time_to_merger = dphase/df/(2pi)"
    # formula's own derivation (see the sign-convention cross-check in
    # notebooks/lgwa_response_validation.ipynb's IMRPhenomD vs
    # lgwa_response comparison) -- negate so t_rel comes out negative
    # before merger, as documented above.
    return -dphase / (2 * eps) / (2 * np.pi)


def _lgwa_response_core(
    fm,
    base_approximant,
    mode_array,
    ref_position,
    longitude_site,
    latitude_site,
    cadence,
    tc,
    lamb,
    beta,
    polarization_ref,
    cbc_params,
):
    """
    Shared core for both the equal-spaced-grid (`fd_det`) and arbitrary-
    sample-points (`fd_det_sequence`) entry points: returns `(h_X, h_Y)`
    as plain complex numpy arrays evaluated at frequencies `fm`.
    """
    ref_position = np.asarray(ref_position, dtype=float)

    # Ref frame -> SSB: sky localization/polarization pass through
    # unchanged (rotation_matrix_ssb_to_ref is the identity); only the
    # arrival time changes, via a closed-form light-travel correction
    # (see pycbc.coordinates.reference_point's module docstring).
    t_ssb_merger = _refpt.t_ssb_from_t_ref(tc, lamb, beta, ref_position)

    # SSB ecliptic -> ICRS ra/dec/psi for lgwa_response's antenna-pattern
    # formula (no LDC/LAL flip -- see _LGWA_detector.project_wave's own
    # docstring for why). Moon frame == SSB numerically (identity
    # rotation), so lamb/beta/polarization_ref (already SSB-frame values)
    # can be passed directly as the "moon-frame" inputs here.
    _, ra, dec, psi = _moon.moon_to_geo(
        t_moon=0.0,
        longitude_moon=lamb,
        latitude_moon=beta,
        polarization_moon=polarization_ref,
        lal_convention=False,
    )

    # SSB -> LGWA site light-travel correction, evaluated once (not per
    # frequency point) and applied as a near-constant shift -- exactly
    # the same justification already used in
    # _LGWA_detector.project_wave: lunar libration/orbital geometry
    # changes on a timescale of days, far slower than t_lm(f) spreads
    # within a single signal.
    delta_t_site = (
        _moon.t_moon_from_ssb(t_ssb_merger, lamb, beta, longitude_site, latitude_site)
        - t_ssb_merger
    )

    lgwa_det = _get_cached_detector(longitude_site, latitude_site, cadence)
    modes = mode_array if mode_array else [(2, 2)]

    h_x_total = np.zeros(len(fm), dtype=complex)
    h_y_total = np.zeros(len(fm), dtype=complex)

    for ell, m in modes:
        hp_lm, hc_lm = get_fd_waveform_sequence(
            approximant=base_approximant,
            mode_array=[[ell, m]],
            sample_points=Array(fm),
            **cbc_params,
        )

        t_rel = _time_of_f_local(fm, base_approximant, (ell, m), cbc_params)
        t_site = t_ssb_merger + t_rel + delta_t_site

        pad = cadence
        n, x, y = lgwa_det._detector_frame(
            t_site.min() - pad, t_site.max() + pad, t_site
        )
        hpx, hcx, hpy, hcy = lgwa_det._antenna_pattern_factors(n, x, y, ra, dec, psi)

        # get_fd_waveform_sequence's raw output has the coalescence
        # reference at t=0; shift to the actual SSB merger time so the
        # returned waveform is anchored consistently with `tc` (relbin's
        # still_needs_det_response branch applies no further time
        # correction of its own -- see relbin.py's dtc=0./ta[ifo]=0.
        # handling in that branch).
        shift = np.exp(-2j * np.pi * fm * t_ssb_merger)
        hp_lm_arr = hp_lm.numpy() * shift
        hc_lm_arr = hc_lm.numpy() * shift

        h_x_total += hp_lm_arr * hpx + hc_lm_arr * hcx
        h_y_total += hp_lm_arr * hpy + hc_lm_arr * hcy

    return h_x_total, h_y_total


def lgwa_fd_response(
    ifos=None,
    sample_points=None,
    delta_f=None,
    f_final=None,
    f_lower=None,
    approximant=None,
    base_approximant=None,
    mode_array=None,
    ref_position=(0.0, 0.0, 0.0),
    longitude_site=None,
    latitude_site=None,
    cadence=3600.0,
    tc=0.0,
    eclipticlongitude=0.0,
    eclipticlatitude=0.0,
    polarization=0.0,
    **params,
):
    """
    Registered as both `pycbc.waveform.fd_det` (equal-spaced grid, via
    `delta_f`/`f_final`) and `pycbc.waveform.fd_det_sequence` (arbitrary
    `sample_points`) under the approximant name `'LGWA_response'`. See
    the module docstring for the overall approach.

    Parameters
    ----------
    base_approximant : str
        Any PyCBC/LAL frequency-domain approximant name (e.g.
        'IMRPhenomXHM', 'IMRPhenomXPHM', an SEOBNR variant, ...),
        forwarded to `get_fd_waveform_sequence` for the intrinsic
        waveform. Required.
    mode_array : list of (l, m) or None
        Which spherical-harmonic modes to include, each evaluated at its
        own stationary-phase arrival time (see module docstring).
        Default None, meaning [(2, 2)] (dominant mode only).
    ref_position : 3-vector, meters
        Fixed (SSB-comoving) reference point at which `tc` is the
        arrival time -- see `pycbc.coordinates.reference_point`. Default
        (0, 0, 0), the SSB origin.
    longitude_site, latitude_site : float, radians
        Selenodetic coordinates of the LGWA site -- see
        `pycbc.detector.space._LGWA_detector`. Required.
    cadence : float, seconds
        See `_LGWA_detector`. Default 3600.
    tc, eclipticlongitude, eclipticlatitude, polarization : float
        Arrival time (at `ref_position`) and SSB-frame sky
        position/polarization of the source.

    Returns
    -------
    dict
        Keyed by the requested `ifos` ('LGWA_X'/'LGWA_Y'), values are
        `Array` (if `sample_points` was given) or `FrequencySeries` (if
        `delta_f`/`f_final` were given instead).
    """
    # `approximant` (='LGWA_response', the plugin's own registered name)
    # arrives here because callers like relbin forward the full
    # static/fiducial params dict, unchanged, to this function -- it is
    # accepted and discarded rather than left in **params, where it would
    # otherwise collide with the explicit `approximant=base_approximant`
    # passed to get_fd_waveform_sequence below.
    del approximant
    if base_approximant is None:
        raise ValueError(
            "base_approximant is required (any PyCBC/LAL "
            "frequency-domain approximant name)."
        )
    if longitude_site is None or latitude_site is None:
        raise ValueError("longitude_site and latitude_site (radians) are required.")

    is_sequence = sample_points is not None
    if is_sequence:
        fm = (
            sample_points.numpy()
            if isinstance(sample_points, Array)
            else np.asarray(sample_points)
        )
    else:
        if delta_f is None or f_final is None:
            raise ValueError(
                "Either sample_points, or both delta_f and f_final, must be given."
            )
        n_samples = int(f_final / delta_f) + 1
        fm = np.arange(n_samples) * delta_f

    cbc_params = dict(params)
    cbc_params["f_lower"] = f_lower
    cbc_params.setdefault("f_ref", f_lower)

    h_x, h_y = _lgwa_response_core(
        fm,
        base_approximant,
        mode_array,
        ref_position,
        longitude_site,
        latitude_site,
        cadence,
        tc,
        eclipticlongitude,
        eclipticlatitude,
        polarization,
        cbc_params,
    )

    out = {"LGWA_X": h_x, "LGWA_Y": h_y}

    if isinstance(ifos, str):
        ifos = [ifos]
    result = {}
    for ifo in ifos:
        arr = out[ifo]
        if is_sequence:
            result[ifo] = Array(arr)
        else:
            result[ifo] = FrequencySeries(arr, delta_f=delta_f)
    return result


lgwa_fd_response.required = [
    "f_lower",
    "approximant",
    "base_approximant",
    "mass1",
    "mass2",
    "tc",
    "eclipticlongitude",
    "eclipticlatitude",
    "polarization",
    "longitude_site",
    "latitude_site",
]
