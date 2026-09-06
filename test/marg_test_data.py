# Copyright (C) 2026  Alex Nitz
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
"""Shared pieces for the inference validation tests.

These tests check calculations that are approximate, stochastic, or both.
Run against one fixed signal in one fixed noise realization they would
only ever show that that one case works, and the cases that do not are
exactly what they are for: choosing the resolution of the time
marginalization looked settled until a rate that happened to put the peak
between two samples turned up.

So the signal and the noise are drawn afresh on every run, from ranges
wide enough to be worth exploring but narrow enough that the tolerances
still mean something. The seed is reported, and can be set, so anything
that fails can be run again exactly.
"""

import os
import random

import numpy

from pycbc.detector import Detector
from pycbc.filter import sigmasq
from pycbc.noise import noise_from_psd
from pycbc.psd import aLIGOZeroDetHighPower
from pycbc.types import TimeSeries
from pycbc.waveform import get_td_waveform

TC = 1187008882.42840
FLOW, SEGLEN, SRATE = 25., 32, 2048
IFOS = ['H1', 'L1', 'V1']
TARGET_SNR = 25.0


def get_seed(name):
    """The seed for this run, taken from the environment or made up.

    Prints the base seed and returns a per-fixture seed derived from it.
    Set PYCBC_VALIDATION_SEED to the printed base to repeat a run: the
    value printed is exactly the value to feed back, and different
    fixtures still get different draws from the same base. The derived
    seed used to be what was printed while the base was what the
    environment expected, so feeding a printed seed back reproduced a
    different run; that is why past failures looked unreproducible.
    """
    given = os.environ.get('PYCBC_VALIDATION_SEED')
    base = (int(given) if given is not None
            else random.SystemRandom().randrange(2 ** 31))
    print("\n%s: PYCBC_VALIDATION_SEED=%d" % (name, base))
    # keep different fixtures from sharing one realization, reproducibly
    return (base + sum(map(ord, name))) % 2 ** 31


def draw_injection(seed):
    """A binary neutron star signal, drawn from somewhere worth testing.

    The distance is set afterwards to put the network signal to noise
    ratio at a known value, so that how loud the signal is does not vary
    with where it happens to be on the sky.
    """
    rng = numpy.random.RandomState(seed)
    return dict(
        mass1=rng.uniform(1.2, 1.6),
        mass2=rng.uniform(1.2, 1.6),
        inclination=numpy.arccos(rng.uniform(-1, 1)),
        ra=rng.uniform(0, 2 * numpy.pi),
        dec=numpy.arcsin(rng.uniform(-1, 1)),
        polarization=rng.uniform(0, 2 * numpy.pi),
        coa_phase=rng.uniform(0, 2 * numpy.pi),
        distance=100.,
    )


def make_data(seed, ifos=None, target_snr=TARGET_SNR):
    """Simulated noise with a signal in it.

    Returns the data, the power spectra, and the parameters that went in,
    with the distance set to give the network signal to noise ratio asked
    for.
    """
    ifos = list(ifos or IFOS)
    injection = draw_injection(seed)

    flen = int(SRATE * SEGLEN / 2) + 1
    psd = aLIGOZeroDetHighPower(flen, 1. / SEGLEN, FLOW)

    def project(params):
        hp, hc = get_td_waveform(
            approximant='TaylorF2', f_lower=FLOW, delta_t=1. / SRATE,
            **{k: params[k] for k in ('mass1', 'mass2', 'distance',
                                      'inclination', 'coa_phase')})
        # the generator puts coalescence at t=0, so this puts it at tc
        hp.start_time += TC
        hc.start_time += TC
        return {ifo: Detector(ifo).project_wave(
                    hp, hc, params['ra'], params['dec'],
                    params['polarization']) for ifo in ifos}

    def embed(signal):
        """Put a signal into a stretch the length of the analysis."""
        blank = TimeSeries(numpy.zeros(int(SEGLEN * SRATE)),
                           delta_t=1. / SRATE, epoch=TC - SEGLEN / 2)
        return blank.add_into(signal).to_frequencyseries()

    signals = project(injection)
    snr = sum(float(sigmasq(embed(s), psd=psd, low_frequency_cutoff=FLOW))
              for s in signals.values()) ** 0.5
    injection['distance'] *= snr / target_snr
    signals = project(injection)

    data, psds = {}, {}
    for i, ifo in enumerate(ifos):
        noise = noise_from_psd(int(SEGLEN * SRATE), 1. / SRATE, psd,
                               seed=(seed + 977 * (i + 1)) % 2 ** 31)
        noise._epoch = TC - SEGLEN / 2
        data[ifo] = noise.add_into(signals[ifo]).to_frequencyseries()
        psds[ifo] = psd
    return data, psds, injection
