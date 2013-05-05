# Copyright (C) 2013 Tito Dal Canton
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

import sys
import tempfile
import lal
import lalsimulation
import pycbc
from pycbc.types import TimeSeries
from pycbc.detector import Detector
from pycbc.inject import InjectionSet
import unittest
import numpy
import optparse
from glue.ligolw import ligolw
from glue.ligolw import lsctables
from glue.ligolw import utils


class MyInjection(object):
    def fill_sim_inspiral_row(self, row):
        # using dummy values for many fields, should work for our purposes
        row.waveform = 'TaylorT4'
        row.distance = self.distance
        total_mass = self.mass1 + self.mass2
        row.mass1 = self.mass1
        row.mass2 = self.mass2
        row.eta = self.mass1 * self.mass2 / total_mass ** 2
        row.mchirp = total_mass * row.eta ** (3. / 5.)
        row.latitude = self.latitude
        row.longitude = self.longitude
        row.inclination = self.inclination
        row.polarization = self.polarization
        row.phi0 = 0
        row.f_lower = 20
        row.f_final = lal.LAL_C_SI ** 3 / \
                (6. ** (3. / 2.) * lal.LAL_PI * lal.LAL_G_SI * total_mass)
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
                lal.LIGOTimeGPS(self.end_time))
        for d in 'lhvgt':
            row.__setattr__('eff_dist_' + d, row.distance)
            row.__setattr__(d + '_end_time', row.geocent_end_time)
            row.__setattr__(d + '_end_time_ns', row.geocent_end_time_ns)
        row.amp_order = 0
        row.coa_phase = 0
        row.bandpass = 0
        row.taper = lalsimulation.LAL_SIM_INSPIRAL_TAPER_NONE
        row.numrel_mode_min = 0
        row.numrel_mode_max = 0
        row.numrel_data = 0
        row.source = 'ANTANI'

_parser = optparse.OptionParser()

def _check_scheme(option, opt_str, scheme, parser):
    if scheme == 'cuda' and not pycbc.HAVE_CUDA:
        raise optparse.OptionValueError("CUDA not found")
    if scheme == 'opencl' and not pycbc.HAVE_OPENCL:
        raise optparse.OptionValueError("OpenCL not found")
    setattr(parser.values, option.dest, scheme)

_parser.add_option('--scheme', '-s', action='callback', type='choice',
    choices=('cpu', 'cuda', 'opencl'), default='cpu', dest='scheme',
    callback=_check_scheme,
    help='specifies processing scheme, can be cpu [default], cuda, or opencl')

_parser.add_option('--device-num', '-d', action='store', type='int',
    dest='devicenum', default=0,
    help='specifies a GPU device to use for CUDA or OpenCL, 0 by default')

(_options, _args) = _parser.parse_args()

if _options.scheme == 'cuda':
    _context = pycbc.scheme.CUDAScheme(device_num=_options.devicenum)
elif _options.scheme == 'opencl':
    _context = pycbc.scheme.OpenCLScheme(device_num=_options.devicenum)
elif _options.scheme == 'cpu':
    _context = pycbc.scheme.CPUScheme()

class TestInjection(unittest.TestCase):
    def setUp(self):
        self.detectors = [Detector(d) for d in ['H1', 'L1', 'V1']]
        self.sample_rate = 4096.
        self.earth_time = lal.LAL_REARTH_SI / lal.LAL_C_SI

        # create a few random injections
        self.injections = []
        start_time = float(lal.GPSTimeNow())
        for i in xrange(10):
            inj = MyInjection()
            inj.end_time = start_time + 40000 * i + \
                    numpy.random.normal(scale=3600)
            random = numpy.random.uniform
            inj.mass1 = random(low=1., high=20.)
            inj.mass2 = random(low=1., high=20.)
            inj.distance = random(low=0.9, high=1.1) * 1e6 * lal.LAL_PC_SI
            inj.latitude = numpy.arccos(random(low=-1, high=1))
            inj.longitude = random(low=0, high=2 * lal.LAL_PI)
            inj.inclination = numpy.arccos(random(low=-1, high=1))
            inj.polarization = random(low=0, high=2 * lal.LAL_PI)
            self.injections.append(inj)

        # create LIGOLW document
        xmldoc = ligolw.Document()
        xmldoc.appendChild(ligolw.LIGO_LW())

        # create sim inspiral table, link it to document and fill it
        sim_table = lsctables.New(lsctables.SimInspiralTable)
        xmldoc.childNodes[-1].appendChild(sim_table)
        for i in xrange(len(self.injections)):
            row = sim_table.RowType()
            self.injections[i].fill_sim_inspiral_row(row)
            row.process_id = 'process:process_id:0'
            row.simulation_id = 'sim_inspiral:simulation_id:%d' % i
            sim_table.append(row)

        # write document to temp file
        self.inj_file = tempfile.NamedTemporaryFile(suffix='.xml')
        utils.write_fileobj(xmldoc, self.inj_file)
    
    def test_injection_presence(self):
        """Verify presence of signals at expected times"""
        injections = InjectionSet(self.inj_file.name)
        for det in self.detectors:
            for inj in self.injections:
                ts = TimeSeries(numpy.zeros(10 * self.sample_rate),
                                delta_t=1/self.sample_rate,
                                epoch=lal.LIGOTimeGPS(inj.end_time - 5),
                                dtype=numpy.float64)
                injections.apply(ts, det.name)
                max_amp, max_loc = ts.abs_max_loc()
                # FIXME could test amplitude and time more precisely
                self.assertTrue(max_amp > 0 and max_amp < 1e-10)
                time_error = ts.sample_times.numpy()[max_loc] - inj.end_time
                self.assertTrue(abs(time_error) < 1.5 * self.earth_time)

    def test_injection_absence(self):
        """Verify absence of signals outside known injection times"""
        clear_times = [
            self.injections[0].end_time - 86400,
            self.injections[-1].end_time + 86400
        ]
        injections = InjectionSet(self.inj_file.name)
        for det in self.detectors:
            for epoch in clear_times:
                ts = TimeSeries(numpy.zeros(10 * self.sample_rate),
                                delta_t=1/self.sample_rate,
                                epoch=lal.LIGOTimeGPS(epoch),
                                dtype=numpy.float64)
                injections.apply(ts, det.name)
                max_amp, max_loc = ts.abs_max_loc()
                self.assertEqual(max_amp, 0)

suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestInjection))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    
    NotImpErrors = 0
    for error in results.errors:
        for errormsg in error:
            if type(errormsg) is str:
                if 'NotImplemented' in errormsg:
                    NotImpErrors += 1
                    break
    if results.wasSuccessful():
        sys.exit(0)
    elif len(results.failures) == 0 and len(results.errors) == NotImpErrors:
        sys.exit(1)
    else:
        sys.exit(2)
