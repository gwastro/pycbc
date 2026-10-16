# Copyright (C) 2013 Tito Dal Canton, Josh Willis
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
Unit test for PyCBC's injection module.
"""

import copy
import itertools
import tempfile
import unittest

import lal
from pycbc.types import TimeSeries
from pycbc.detector import Detector, get_available_detectors
from pycbc.inject import InjectionSet
from pycbc.waveform import get_td_waveform
from pycbc.io import FieldArray
import numpy
from igwn_ligolw import ligolw
from igwn_ligolw import lsctables
from igwn_ligolw import utils as ligolw_utils
from utils import parse_args_cpu_only, simple_exit

# Injection tests only need to happen on the CPU
parse_args_cpu_only("Injections")


class MyInjection:
    def fill_sim_inspiral_row(self, row):
        # using dummy values for many fields, should work for our purposes
        row.waveform = "TaylorT4threePointFivePN"
        row.distance = self.distance
        total_mass = self.mass1 + self.mass2
        row.mass1 = self.mass1
        row.mass2 = self.mass2
        row.eta = self.mass1 * self.mass2 / total_mass**2
        row.mchirp = total_mass * row.eta ** (3.0 / 5.0)
        row.latitude = self.latitude
        row.longitude = self.longitude
        row.inclination = self.inclination
        row.polarization = self.polarization
        row.phi0 = 0
        row.f_lower = 20
        row.f_final = lal.C_SI**3 / (
            6.0 ** (3.0 / 2.0) * lal.PI * lal.G_SI * total_mass
        )
        row.spin1x = row.spin1y = row.spin1z = 0
        row.spin2x = row.spin2y = row.spin2z = 0
        row.alpha1 = 0
        row.alpha2 = 0
        row.alpha3 = 0
        row.alpha4 = 0
        row.alpha5 = 0
        row.alpha6 = 0
        row.alpha = 0
        row.beta = 0
        row.theta0 = 0
        row.psi0 = 0
        row.psi3 = 0
        row.geocent_end_time = int(self.end_time)
        row.geocent_end_time_ns = int(1e9 * (self.end_time - row.geocent_end_time))
        row.end_time_gmst = lal.GreenwichMeanSiderealTime(
            lal.LIGOTimeGPS(self.end_time)
        )
        for d in "lhvgt":
            row.__setattr__("eff_dist_" + d, row.distance)
            row.__setattr__(d + "_end_time", row.geocent_end_time)
            row.__setattr__(d + "_end_time_ns", row.geocent_end_time_ns)
        row.amp_order = 0
        row.coa_phase = 0
        row.bandpass = 0
        row.taper = self.taper
        row.numrel_mode_min = 0
        row.numrel_mode_max = 0
        row.numrel_data = None
        row.source = "ANTANI"


class TestInjection(unittest.TestCase):
    def setUp(self):
        available_detectors = get_available_detectors()
        self.assertTrue("H1" in available_detectors)
        self.assertTrue("L1" in available_detectors)
        self.assertTrue("V1" in available_detectors)
        self.detectors = [Detector(d) for d in ["H1", "L1", "V1"]]
        self.sample_rate = 4096.0
        self.earth_time = lal.REARTH_SI / lal.C_SI

        # create a few random injections
        self.injections = []
        start_time = float(lal.GPSTimeNow())
        taper_choices = ("TAPER_NONE", "TAPER_START", "TAPER_END", "TAPER_STARTEND")
        for i, taper in zip(range(20), itertools.cycle(taper_choices)):
            inj = MyInjection()
            inj.end_time = start_time + 40000 * i + numpy.random.normal(scale=3600)
            random = numpy.random.uniform
            inj.mass1 = random(low=1.0, high=20.0)
            inj.mass2 = random(low=1.0, high=20.0)
            inj.distance = random(low=0.9, high=1.1) * 1e6 * lal.PC_SI
            inj.latitude = numpy.arccos(random(low=-1, high=1))
            inj.longitude = random(low=0, high=2 * lal.PI)
            inj.inclination = numpy.arccos(random(low=-1, high=1))
            inj.polarization = random(low=0, high=2 * lal.PI)
            inj.taper = taper
            self.injections.append(inj)

        self.inj_file = self.write_xml(self.injections)

    def write_xml(self, injections):
        # create LIGOLW document
        xmldoc = ligolw.Document()
        xmldoc.appendChild(ligolw.LIGO_LW())

        # create sim inspiral table, link it to document and fill it
        sim_table = lsctables.SimInspiralTable.new()
        xmldoc.childNodes[-1].appendChild(sim_table)
        for i in range(len(injections)):
            row = sim_table.RowType()
            injections[i].fill_sim_inspiral_row(row)
            row.process_id = 0
            row.simulation_id = i
            sim_table.append(row)

        # write document to temp file
        inj_file = tempfile.NamedTemporaryFile(suffix=".xml")
        ligolw_utils.write_fileobj(xmldoc, inj_file)
        return inj_file

    def test_injection_presence(self):
        """Verify presence of signals at expected times"""
        injections = InjectionSet(self.inj_file.name)
        for det in self.detectors:
            for inj in self.injections:
                ts = TimeSeries(
                    numpy.zeros(int(10 * self.sample_rate)),
                    delta_t=1 / self.sample_rate,
                    epoch=lal.LIGOTimeGPS(inj.end_time - 5),
                    dtype=numpy.float64,
                )
                injections.apply(ts, det.name)
                max_amp, max_loc = ts.abs_max_loc()
                # FIXME could test amplitude and time more precisely
                self.assertTrue(max_amp > 0 and max_amp < 1e-10)
                time_error = ts.sample_times.numpy()[max_loc] - inj.end_time
                self.assertTrue(abs(time_error) < 2 * self.earth_time)

    def test_injection_absence(self):
        """Verify absence of signals outside known injection times"""
        clear_times = [
            self.injections[0].end_time - 86400,
            self.injections[-1].end_time + 86400,
        ]
        injections = InjectionSet(self.inj_file.name)
        for det in self.detectors:
            for epoch in clear_times:
                ts = TimeSeries(
                    numpy.zeros(int(10 * self.sample_rate)),
                    delta_t=1 / self.sample_rate,
                    epoch=lal.LIGOTimeGPS(epoch),
                    dtype=numpy.float64,
                )
                injections.apply(ts, det.name)
                max_amp, max_loc = ts.abs_max_loc()
                self.assertEqual(max_amp, 0)

    def test_injection_taper(self):
        """
        Verify that injections are tapered as requested by the injection
        file, for both xml and hdf files
        """
        det = self.detectors[0]
        # use high masses so that the whole signal is in the data
        inj = copy.copy(self.injections[0])
        inj.mass1 = inj.mass2 = 15.0

        def inject(inj_file):
            ts = TimeSeries(
                numpy.zeros(int(16 * self.sample_rate)),
                delta_t=1 / self.sample_rate,
                epoch=lal.LIGOTimeGPS(inj.end_time - 14),
                dtype=numpy.float64,
            )
            InjectionSet(inj_file.name).apply(ts, det.name)
            return ts.numpy()

        def check(tapered, untapered, start, end):
            # only the tapered half (or halves) of the signal should change;
            # adding the signal to the data spreads tiny differences
            # everywhere, so compare to a fraction of the peak
            atol = 1e-3 * abs(untapered).max()
            signal = numpy.nonzero(abs(untapered) > atol)[0]
            mid = (signal[0] + signal[-1]) // 2
            for half, is_tapered in [
                (slice(None, mid), start),
                (slice(mid, None), end),
            ]:
                unchanged = numpy.allclose(
                    tapered[half], untapered[half], rtol=0, atol=atol
                )
                self.assertEqual(unchanged, not is_tapered)

        # xml file
        inj.taper = "TAPER_NONE"
        untapered = inject(self.write_xml([inj]))
        for taper, start, end in [
            ("TAPER_START", True, False),
            ("TAPER_END", False, True),
            ("TAPER_STARTEND", True, True),
        ]:
            inj.taper = taper
            check(inject(self.write_xml([inj])), untapered, start, end)

        # hdf file
        samples = FieldArray.from_kwargs(
            tc=[inj.end_time],
            mass1=[inj.mass1],
            mass2=[inj.mass2],
            distance=[inj.distance / (1e6 * lal.PC_SI)],
            ra=[inj.longitude],
            dec=[inj.latitude],
            inclination=[inj.inclination],
            polarization=[inj.polarization],
        )

        def write_hdf(**taper_args):
            inj_file = tempfile.NamedTemporaryFile(suffix=".hdf")
            InjectionSet.write(
                inj_file.name,
                samples,
                static_args=dict(approximant="TaylorT4", f_lower=20.0, **taper_args),
            )
            return inj_file

        untapered = inject(write_hdf())
        check(inject(write_hdf(taper="start")), untapered, True, False)
        check(
            inject(write_hdf(taper="end", taper_method="constant", taper_window=0.5)),
            untapered,
            False,
            True,
        )


class TestRadiationFrameInjection(unittest.TestCase):
    """Test injecting in the radiation frame (RF) without any sky location or
    polarization parameters.
    """
    def setUp(self):
        self.sample_rate = 2048.
        self.tc = 1126259462.42
        # note: no ra, dec, or polarization
        self.params = {'tc': self.tc, 'mass1': 37., 'mass2': 32.,
                       'inclination': 2.5, 'coa_phase': 1.5,
                       'distance': 100., 'f_ref': 20., 'f_lower': 18.,
                       'approximant': 'IMRPhenomD'}
        samples = FieldArray.from_kwargs(
            **{p: numpy.array([v]) for p, v in self.params.items()})
        self.inj_file = tempfile.NamedTemporaryFile(suffix='.hdf')
        InjectionSet.write(self.inj_file.name, samples)

    def test_rf_injection(self):
        """Verify that an RF injection is the plus polarization at tc"""
        injections = InjectionSet(self.inj_file.name)

        # zero time series over which to inject signal
        ts = TimeSeries(numpy.zeros(int(16 * self.sample_rate)),
                        delta_t=1/self.sample_rate,
                        epoch=lal.LIGOTimeGPS(self.tc - 12),
                        dtype=numpy.float64)
        expected = ts.copy()

        # inject the radiation-frame waveform into zeroes and check it's not
        # silently skipped
        injected = injections.apply(ts, 'RF')
        self.assertEqual(len(injected.table), 1)

        # compare to the plus polarization shifted to tc
        hp, _ = get_td_waveform(delta_t=1/self.sample_rate, **self.params)
        hp.start_time += self.tc

        # inject the plus polarization into zeroes
        expected = expected.inject(hp)

        # check that the injection wasn't silently skipped
        max_amp, max_loc = ts.abs_max_loc()
        self.assertTrue(max_amp > 0)

        # check that the peak of |hp| alone is within a cycle of tc
        time_error = ts.sample_times.numpy()[max_loc] - self.tc
        self.assertTrue(abs(time_error) < 0.01)
        numpy.testing.assert_allclose(ts.numpy(), expected.numpy(),
                                      rtol=0, atol=1e-3 * max_amp)

    def test_detector_injection_needs_sky_location(self):
        """Verify that injecting in a real detector still requires the sky
        location"""
        injections = InjectionSet(self.inj_file.name)
        zero_ts = TimeSeries(numpy.zeros(int(16 * self.sample_rate)),
                        delta_t=1/self.sample_rate,
                        epoch=lal.LIGOTimeGPS(self.tc - 12),
                        dtype=numpy.float64)
        with self.assertRaises(Exception):
            injections.apply(zero_ts, 'H1')


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestInjection))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestRadiationFrameInjection))

if __name__ == "__main__":
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
