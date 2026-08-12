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
"""How long the likelihoods take, so that changing them can be measured.

Most of the cost of a marginalized likelihood is not the waveform. It is
the bookkeeping around it, which is exactly the part that gets changed
when someone tries to make it faster or more accurate, and exactly the
part nothing was watching. This runs each model over the same signal and
reports the time per call.

Run it to print a table:

    python marg_perf.py

Save what a branch does, then compare another branch against it:

    python marg_perf.py --save before.json
    python marg_perf.py --compare before.json

Configurations whose options a branch does not have are reported as
unsupported rather than failing, so the same benchmark runs either side
of a change that adds one.

The signal is drawn from a fixed seed here, unlike the validation tests:
a benchmark wants the same work every time, not a new sample of it.

Run it on an otherwise idle machine. Each case is the fastest of several
batches, which removes the odd interrupted run, and the comparison
divides out the drift common to the whole board, which removes a machine
that is uniformly busier than it was before. Neither of those saves a
case that spends tens of milliseconds in a Python loop while something
else is contending for the interpreter; those are only meaningful when
the machine is quiet. The relative-binning cases, which are what the
marginalization work is trying to make faster, are short enough to stay
steady through moderate load.
"""

import argparse
import copy
import json
import sys
import time

import numpy

from pycbc.detector import Detector
from pycbc.distributions import (
    CosAngle,
    JointDistribution,
    SinAngle,
    Uniform,
    UniformAngle,
)
from pycbc.filter import sigmasq
from pycbc.inference import models
from pycbc.noise import noise_from_psd
from pycbc.psd import aLIGOZeroDetHighPower
from pycbc.types import TimeSeries
from pycbc.waveform import get_td_waveform

SEED = 20260812
REPEATS = [200]

# The signal every configuration is timed against: a binary neutron star in
# simulated noise, built from a fixed seed so the benchmark does the same
# work each run. Kept self-contained here rather than shared with the
# accuracy tests (which deliberately draw a fresh signal every run); this is
# a copy of just the fixed-signal builder, with the distance set to a known
# network signal to noise ratio.
TC = 1187008882.42840
FLOW, SEGLEN, SRATE = 25., 32, 2048
IFOS = ['H1', 'L1', 'V1']
TARGET_SNR = 25.0


def draw_injection(seed):
    """A binary neutron star signal, at a fixed distance before scaling."""
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
    """Simulated noise with a signal in it, distance set to a known SNR."""
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


def setup():
    """The signal every configuration is timed against."""
    data, psds, inj = make_data(SEED)
    flow = {ifo: FLOW for ifo in data}
    static = dict(mass1=inj['mass1'], mass2=inj['mass2'], f_lower=FLOW,
                  approximant='TaylorF2')
    fid = {'mass1': inj['mass1'], 'tc': TC, 'ra': inj['ra'],
           'dec': inj['dec'], 'polarization': inj['polarization']}
    return data, psds, flow, static, fid, inj


DATA, PSDS, FLOWD, STATIC, FID, INJ = setup()


def distributions(names):
    known = {
        'distance': Uniform(distance=(10, 200)),
        'inclination': SinAngle(inclination=None),
        'tc': Uniform(tc=(TC - 0.1, TC + 0.1)),
        'polarization': UniformAngle(polarization=None),
        'ra': UniformAngle(ra=None),
        'dec': CosAngle(dec=None),
        'coa_phase': UniformAngle(coa_phase=None),
    }
    return [known[n] for n in names]


def build(cls, variable, relbin=True, **kwargs):
    static = {k: v for k, v in STATIC.items() if k not in variable}
    for name in ('ra', 'dec', 'polarization', 'tc'):
        if name not in variable:
            static[name] = FID[name] if name != 'tc' else TC
    if relbin:
        kwargs.setdefault('fiducial_params', FID)
        kwargs.setdefault('epsilon', 0.1)
    return cls(list(variable), copy.deepcopy(DATA),
               low_frequency_cutoff=FLOWD, psds=PSDS, static_params=static,
               prior=JointDistribution(list(variable),
                                       *distributions(variable)), **kwargs)


CASES = [
    # (name, model, variable params, relative binning, options)
    ('gaussian, nothing marginalized',
     models.GaussianNoise, ['distance', 'inclination'], False, {}),
    ('gaussian, phase',
     models.MarginalizedPhaseGaussianNoise,
     ['distance', 'inclination'], False, {}),
    ('gaussian, phase and time',
     models.MarginalizedTime, ['distance', 'inclination', 'tc'], False,
     dict(marginalize_phase=True, marginalize_vector_params='tc',
          sample_rate=4096, marginalize_vector_samples=1000)),
    ('gaussian, phase time and sky',
     models.MarginalizedTime,
     ['distance', 'inclination', 'tc', 'polarization', 'ra', 'dec'], False,
     dict(marginalize_phase=True,
          marginalize_vector_params='tc,ra,dec,polarization',
          sample_rate=4096, marginalize_vector_samples=1000,
          marginalize_sky_initial_samples=1e5)),
    ('relbin, phase',
     models.Relative, ['distance', 'inclination'], True, {}),
    ('relbin, phase and distance',
     models.Relative, ['distance', 'inclination'], True,
     dict(marginalize_distance=True,
          marginalize_distance_param='distance',
          marginalize_distance_interpolator=True,
          marginalize_distance_snr_range=(5, 60),
          marginalize_distance_density=(100, 100),
          marginalize_distance_samples=1000)),
    ('relbin, phase and time',
     models.RelativeTime, ['distance', 'inclination', 'tc'], True,
     dict(marginalize_vector_params='tc', sample_rate=4096,
          marginalize_vector_samples=1000)),
    ('relbin, phase time and sky',
     models.RelativeTime,
     ['distance', 'inclination', 'tc', 'polarization', 'ra', 'dec'], True,
     dict(marginalize_vector_params='tc,ra,dec,polarization',
          sample_rate=4096, marginalize_vector_samples=1000,
          marginalize_sky_initial_samples=1e5)),
    ('relbin, sky with precalculated points',
     models.RelativeTime,
     ['distance', 'inclination', 'tc', 'polarization', 'ra', 'dec'], True,
     dict(marginalize_vector_params='tc,ra,dec,polarization',
          sample_rate=4096, marginalize_vector_samples=1000,
          marginalize_sky_initial_samples=1e5,
          precalculate_marginalization_points=1e4)),
    ('relbin, sky with 4x the points',
     models.RelativeTime,
     ['distance', 'inclination', 'tc', 'polarization', 'ra', 'dec'], True,
     dict(marginalize_vector_params='tc,ra,dec,polarization',
          sample_rate=4096, marginalize_vector_samples=4000,
          marginalize_sky_initial_samples=1e5)),
    ('relbin, finer bins',
     models.Relative, ['distance', 'inclination'], True,
     dict(epsilon=0.01)),
]


def time_case(name, cls, variable, relbin, options, repeats=None):
    repeats = REPEATS[0] if repeats is None else repeats
    """Build one configuration and time its likelihood."""
    numpy.random.seed(SEED)
    try:
        model = build(cls, variable, relbin=relbin, **options)
    except Exception as exc:
        return {'name': name, 'unsupported': str(exc)[:120]}

    point = {'distance': INJ['distance'], 'inclination': INJ['inclination'],
             'tc': TC, 'polarization': INJ['polarization'],
             'ra': INJ['ra'], 'dec': INJ['dec']}
    point = {p: point[p] for p in model.variable_params if p in point}

    model.update(**point)
    assert numpy.isfinite(model.loglr)

    # The minimum over several batches, not the mean over one, because a
    # benchmark wants the cost of the work and not the cost of whatever
    # else the machine was doing. A batch that got interrupted is slow;
    # the fastest batch is the one that was left alone, and that is the
    # number that is stable enough to regress against.
    per_batch = max(1, repeats // 5)
    times = []
    for _ in range(5):
        start = time.perf_counter()
        for i in range(per_batch):
            model.update(**dict(point, inclination=0.5 + 0.001 * (i % 100)))
            assert numpy.isfinite(model.loglr)
        times.append((time.perf_counter() - start) / per_batch * 1e3)
    per_call = min(times)

    nbin = None
    if hasattr(model, 'fedges'):
        nbin = len(model.fedges[list(model.data)[0]])
    return {'name': name, 'ms': per_call, 'nbin': nbin}


def run():
    results = []
    for case in CASES:
        result = time_case(*case)
        results.append(result)
        if 'unsupported' in result:
            print('%-38s %s' % (result['name'], 'unsupported'))
        else:
            print('%-38s %8.3f ms%s'
                  % (result['name'], result['ms'],
                     '' if result['nbin'] is None
                     else '   %d bins' % result['nbin']))
    return results


def compare(results, baseline):
    old = {r['name']: r for r in baseline}
    pairs = [(r, old[r['name']]) for r in results
             if r['name'] in old and 'ms' in r and 'ms' in old[r['name']]]
    ratios = [r['ms'] / was['ms'] for r, was in pairs]
    if not ratios:
        return []

    # everything drifts together when the machine is busy, so the middle
    # of the board is taken as that ambient drift and divided out. What is
    # left is what a change did, not what the load did. This is why the
    # benchmark does not need a perfectly quiet machine, only a change
    # that stands out from the rest of the board.
    ambient = float(numpy.median(ratios))
    print('\nambient drift (median of all cases) %.2fx, divided out below'
          % ambient)
    print('%-38s %10s %10s %8s %8s'
          % ('', 'baseline', 'now', 'raw', 'vs rest'))
    worse = []
    for result, was in pairs:
        ratio = result['ms'] / was['ms']
        relative = ratio / ambient
        flag = ''
        if relative > 1.15:
            flag = '  SLOWER'
            worse.append(result['name'])
        elif relative < 0.87:
            flag = '  faster'
        print('%-38s %8.3f ms %8.3f ms %7.2fx %7.2fx%s'
              % (result['name'], was['ms'], result['ms'], ratio, relative,
                 flag))
    return worse


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--save', metavar='FILE')
    parser.add_argument('--compare', metavar='FILE')
    parser.add_argument('--repeats', type=int, default=REPEATS[0])
    args = parser.parse_args()

    REPEATS[0] = args.repeats

    results = run()

    if args.save:
        with open(args.save, 'w') as handle:
            json.dump(results, handle, indent=1)
        print('\nsaved to %s' % args.save)

    if args.compare:
        with open(args.compare) as handle:
            worse = compare(results, json.load(handle))
        if worse:
            print('\nslower than the baseline: %s' % ', '.join(worse))
            return 1
    return 0


if __name__ == '__main__':
    import logging
    import warnings
    warnings.filterwarnings('ignore')
    logging.disable(logging.CRITICAL)
    from astropy.utils import iers
    iers.conf.auto_download = False
    sys.exit(main())
