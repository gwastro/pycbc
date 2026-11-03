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
"""
A marginalized likelihood is the integral of the likelihood over the prior
of the marginalized parameters. Each test compares one against a grid sum
of the same model, unmarginalized, over the same prior.

Each is checked for every time-marginalizing model, with tc sampled both at
the geocenter and in a detector's frame (tc_ref_frame).
"""

import copy
import unittest

import numpy
from scipy.special import logsumexp
from utils import simple_exit

from pycbc.detector import Detector
from pycbc.distributions import JointDistribution, SinAngle, Uniform
from pycbc.inference import models
from pycbc.noise import noise_from_psd
from pycbc.psd import aLIGOZeroDetHighPower
from pycbc.waveform import get_td_waveform

TC = 1187008882.4284
RA, DEC = 1.7, -0.4
FLOW, SEGLEN, SRATE = 25.0, 32, 2048
POINT = {"distance": 50.0, "inclination": 0.5}
FRAMES = ["geocentric", "H1"]

# The SNR sample rate each model is checked at. RelativeTimeDom and
# SingleTemplate interpolate <s|h> between samples, which at 8192 Hz is off
# by about 0.1 for this signal; at 16384 Hz it is below the test's noise.
MODELS = {
    "RelativeTime": 8192,
    "RelativeTimeDom": 16384,
    "SingleTemplate": 16384,
}


def grid(lo, hi, n):
    """The midpoints of n equal cells between lo and hi."""
    return lo + (numpy.arange(n) + 0.5) * (hi - lo) / n


def frame_delay(frame, ra, dec, t):
    """The time of arrival in frame, less the geocentric time."""
    if frame == "geocentric":
        return 0.0
    return Detector(frame).time_delay_from_earth_center(ra, dec, t)


def frame_tc(frame):
    """The signal's tc in the given frame."""
    return TC + frame_delay(frame, RA, DEC, TC)


class TestMargNormalization(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # a binary neutron star at tc, in Gaussian noise in two detectors
        flen = int(SRATE * SEGLEN / 2) + 1
        psd = aLIGOZeroDetHighPower(flen, 1.0 / SEGLEN, FLOW)
        hp, hc = get_td_waveform(
            approximant="TaylorF2",
            f_lower=FLOW,
            delta_t=1.0 / SRATE,
            mass1=1.42,
            mass2=1.36,
            coa_phase=1.1,
            **POINT,
        )
        hp.start_time += TC
        hc.start_time += TC
        cls.data, cls.psds = {}, {}
        for seed, ifo in enumerate(["H1", "L1"]):
            noise = noise_from_psd(SEGLEN * SRATE, 1.0 / SRATE, psd, seed=seed)
            noise._epoch = TC - SEGLEN / 2
            signal = Detector(ifo).project_wave(hp, hc, RA, DEC, 0.3)
            cls.data[ifo] = noise.add_into(signal).to_frequencyseries()
            cls.psds[ifo] = psd
        cls.static = dict(
            mass1=1.42,
            mass2=1.36,
            f_lower=FLOW,
            approximant="TaylorF2",
            polarization=0.3,
        )

    def model(self, cls, bounds, frame="geocentric", **kwargs):
        """A model with a uniform prior over the given bounds."""
        variable = ["distance", "inclination"] + list(bounds)
        prior = JointDistribution(
            variable,
            SinAngle(inclination=None),
            Uniform(distance=(10, 200)),
            *[Uniform(**{k: v}) for k, v in bounds.items()],
        )
        static = {k: v for k, v in self.static.items()}
        static.update(
            {k: v for k, v in dict(ra=RA, dec=DEC).items() if k not in bounds}
        )
        if frame != "geocentric":
            static["tc_ref_frame"] = frame
        if cls is not models.SingleTemplate:
            # the fiducial tc is in the frame tc is sampled in
            kwargs["fiducial_params"] = dict(
                mass1=1.42, tc=frame_tc(frame), ra=RA, dec=DEC
            )
            kwargs["epsilon"] = 0.1
        return cls(
            list(variable),
            copy.deepcopy(self.data),
            low_frequency_cutoff={ifo: FLOW for ifo in self.data},
            psds=self.psds,
            static_params=static,
            prior=prior,
            **kwargs,
        )

    def integral(self, bounds, cells, frame="geocentric"):
        """
        The log of the mean likelihood ratio over a grid of the prior.

        Each grid point's tc, in the given frame, is moved to the
        geocenter for its sky location, so the sum does not depend on the
        models' own handling of tc_ref_frame.
        """
        model = self.model(models.Relative, bounds)
        axes = [grid(*bounds[k], cells[k]) for k in bounds]
        lr = []
        for values in zip(*[a.ravel() for a in numpy.meshgrid(*axes)]):
            p = dict(zip(bounds, values))
            ra, dec = p.get("ra", RA), p.get("dec", DEC)
            p["tc"] -= frame_delay(frame, ra, dec, p["tc"])
            model.update(**p, **POINT)
            lr.append(model.loglr)
        return logsumexp(lr) - numpy.log(len(lr))

    def marginal(self, name, bounds, frame="geocentric", nseed=6, **kwargs):
        """The marginalized log likelihood ratio and its error."""
        values = []
        for seed in range(nseed):
            numpy.random.seed(seed)
            model = self.model(
                getattr(models, name),
                bounds,
                frame,
                marginalize_vector_params=",".join(bounds),
                marginalize_vector_samples=2000,
                sample_rate=MODELS[name],
                **kwargs,
            )
            model.update(**POINT)
            values.append(model.loglr)
        return numpy.mean(values), numpy.std(values) / nseed**0.5

    def test_time(self):
        for frame in FRAMES:
            # The prior ends at the signal, so about half of the times drawn
            # fall outside it and must carry no weight
            tc = frame_tc(frame)
            bounds = {"tc": (tc - 0.002, tc)}
            integral = self.integral(bounds, {"tc": 4000}, frame)
            for name in MODELS:
                with self.subTest(model=name, frame=frame):
                    mean, error = self.marginal(name, bounds, frame)
                    self.assertAlmostEqual(mean, integral, delta=4 * error)

    def test_sky_and_time(self):
        for frame in FRAMES:
            # A patch of sky small enough that some of the points drawn for
            # it fall outside it, which is when the fraction kept matters
            tc = frame_tc(frame)
            bounds = {
                "tc": (tc - 0.002, tc),
                "ra": (RA - 0.02, RA + 0.02),
                "dec": (DEC - 0.02, DEC + 0.02),
            }
            # a 15 x 15 sky grid sums about 0.07 short of the integral;
            # 31 x 31 is within 0.01 of it
            integral = self.integral(bounds, {"tc": 61, "ra": 31, "dec": 31}, frame)
            # points drawn for each likelihood, and drawn once in advance
            for name in MODELS:
                for kwargs in ({}, {"precalculate_marginalization_points": 4000}):
                    with self.subTest(model=name, frame=frame, **kwargs):
                        mean, error = self.marginal(
                            name,
                            bounds,
                            frame,
                            marginalize_sky_initial_samples=3e4,
                            **kwargs,
                        )
                        self.assertAlmostEqual(
                            mean, integral, delta=max(0.1, 4 * error)
                        )


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMargNormalization))

if __name__ == "__main__":
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
