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
"""Tests the analytic delay-to-sky draw through a model that uses it."""

import copy
import unittest

import numpy
from utils import simple_exit

from pycbc.detector import Detector
from pycbc.distributions import CosAngle, JointDistribution, SinAngle, Uniform
from pycbc.inference import models
from pycbc.noise import noise_from_psd
from pycbc.psd import aLIGOZeroDetHighPower
from pycbc.waveform import get_td_waveform

TC = 1187008882.42
FLOW, SEGLEN, SRATE = 25.0, 32, 2048
# three detectors: with two the second delay is a free azimuth and with
# one there is no delay to invert, so three is where the inversion does
# its full job
IFOS = ["H1", "L1", "V1"]
VARIABLE = ["distance", "inclination", "tc", "ra", "dec"]
POINT = {"distance": 40.0, "inclination": 0.5}
SKY = {"ra": 1.7, "dec": -0.4, "polarization": 0.3}
VSAMPLES = 2000
SRATE_OVERRIDE = 4096


class TestAnalyticSkyDraw(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        hp, hc = get_td_waveform(
            approximant="TaylorF2",
            f_lower=FLOW,
            delta_t=1.0 / SRATE,
            mass1=1.4,
            mass2=1.35,
            distance=40.0,
            inclination=0.5,
            coa_phase=1.1,
        )
        hp.start_time += TC
        hc.start_time += TC
        psd = aLIGOZeroDetHighPower(int(SRATE * SEGLEN / 2) + 1, 1.0 / SEGLEN, FLOW)
        cls.data, cls.psds = {}, {}
        for n, ifo in enumerate(IFOS):
            noise = noise_from_psd(int(SEGLEN * SRATE), 1.0 / SRATE, psd, seed=11 + n)
            noise._epoch = TC - SEGLEN / 2
            sig = Detector(ifo).project_wave(
                hp, hc, SKY["ra"], SKY["dec"], SKY["polarization"]
            )
            cls.data[ifo] = noise.add_into(sig).to_frequencyseries()
            cls.psds[ifo] = psd
        cls.static = {
            "mass1": 1.4,
            "mass2": 1.35,
            "f_lower": FLOW,
            "approximant": "TaylorF2",
            "polarization": SKY["polarization"],
        }
        cls.fiducial = {"mass1": 1.4, "tc": TC, **SKY}

    def model(self, ifos=None, **kwargs):
        ifos = IFOS if ifos is None else ifos
        numpy.random.seed(5)
        prior = JointDistribution(
            list(VARIABLE),
            SinAngle(inclination=None),
            Uniform(distance=(10, 200)),
            Uniform(tc=(TC - 0.1, TC + 0.1)),
            Uniform(ra=(0, 2 * numpy.pi)),
            CosAngle(dec=None),
        )
        return models.RelativeTimeDom(
            list(VARIABLE),
            {i: copy.deepcopy(self.data[i]) for i in ifos},
            low_frequency_cutoff=dict.fromkeys(ifos, FLOW),
            psds={i: self.psds[i] for i in ifos},
            static_params=self.static,
            prior=prior,
            fiducial_params=self.fiducial,
            epsilon=0.1,
            marginalize_phase=True,
            marginalize_vector_params="tc,ra,dec",
            marginalize_vector_samples=VSAMPLES,
            sample_rate=SRATE_OVERRIDE,
            **kwargs,
        )

    def test_it_agrees_with_the_map(self):
        """Both draws estimate the same integral.

        At this signal-to-noise ratio a brute-force average over the sky
        gives 636.6, and the map reads about 0.15 below the analytic draw.
        """
        values = {}
        for analytic in (False, True):
            model = self.model(marginalize_sky_analytic=analytic)
            loglr = []
            for _ in range(5):
                model.update(**POINT)
                loglr.append(model.loglr)
            values[analytic] = numpy.mean(loglr)
            # only the map builds its table of sky points
            self.assertEqual(hasattr(model, "tinfo"), not analytic)
        self.assertAlmostEqual(values[True], values[False], delta=0.6)


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestAnalyticSkyDraw))

if __name__ == "__main__":
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
