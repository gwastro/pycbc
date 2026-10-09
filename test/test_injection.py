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

import tempfile
import lal
from pycbc.types import TimeSeries
from pycbc.detector import Detector, get_available_detectors
from pycbc.inject import (InjectionSet, read_injection_table,
                          read_hdf_param_table, get_table_column)
from pycbc.io import HFile, WaveformArray
import unittest
import numpy
import itertools
from igwn_ligolw import ligolw
from igwn_ligolw import lsctables
from igwn_ligolw import utils as ligolw_utils
from utils import parse_args_cpu_only, simple_exit

# Injection tests only need to happen on the CPU
parse_args_cpu_only("Injections")

class MyInjection(object):
    def fill_sim_inspiral_row(self, row):
        # using dummy values for many fields, should work for our purposes
        row.waveform = 'TaylorT4threePointFivePN'
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
        row.f_final = lal.C_SI ** 3 / \
                (6. ** (3. / 2.) * lal.PI * lal.G_SI * total_mass)
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
        row.taper = self.taper
        row.numrel_mode_min = 0
        row.numrel_mode_max = 0
        row.numrel_data = None
        row.source = 'ANTANI'

class TestInjection(unittest.TestCase):
    def setUp(self):
        available_detectors = get_available_detectors()
        self.assertTrue('H1' in available_detectors)
        self.assertTrue('L1' in available_detectors)
        self.assertTrue('V1' in available_detectors)
        self.detectors = [Detector(d) for d in ['H1', 'L1', 'V1']]
        self.sample_rate = 4096.
        self.earth_time = lal.REARTH_SI / lal.C_SI

        # create a few random injections
        self.injections = []
        start_time = float(lal.GPSTimeNow())
        taper_choices = ('TAPER_NONE', 'TAPER_START', 'TAPER_END', 'TAPER_STARTEND')
        for i, taper in zip(range(20), itertools.cycle(taper_choices)):
            inj = MyInjection()
            inj.end_time = start_time + 40000 * i + \
                    numpy.random.normal(scale=3600)
            random = numpy.random.uniform
            inj.mass1 = random(low=1., high=20.)
            inj.mass2 = random(low=1., high=20.)
            inj.distance = random(low=0.9, high=1.1) * 1e6 * lal.PC_SI
            inj.latitude = numpy.arccos(random(low=-1, high=1))
            inj.longitude = random(low=0, high=2 * lal.PI)
            inj.inclination = numpy.arccos(random(low=-1, high=1))
            inj.polarization = random(low=0, high=2 * lal.PI)
            inj.taper = taper
            self.injections.append(inj)

        # create LIGOLW document
        xmldoc = ligolw.Document()
        xmldoc.appendChild(ligolw.LIGO_LW())

        # create sim inspiral table, link it to document and fill it
        sim_table = lsctables.SimInspiralTable.new()
        xmldoc.childNodes[-1].appendChild(sim_table)
        for i in range(len(self.injections)):
            row = sim_table.RowType()
            self.injections[i].fill_sim_inspiral_row(row)
            row.process_id = 0
            row.simulation_id = i
            sim_table.append(row)

        # write document to temp file
        self.inj_file = tempfile.NamedTemporaryFile(suffix='.xml')
        ligolw_utils.write_fileobj(xmldoc, self.inj_file)

    def test_injection_presence(self):
        """Verify presence of signals at expected times"""
        injections = InjectionSet(self.inj_file.name)
        for det in self.detectors:
            for inj in self.injections:
                ts = TimeSeries(numpy.zeros(int(10 * self.sample_rate)),
                                delta_t=1/self.sample_rate,
                                epoch=lal.LIGOTimeGPS(inj.end_time - 5),
                                dtype=numpy.float64)
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
            self.injections[-1].end_time + 86400
        ]
        injections = InjectionSet(self.inj_file.name)
        for det in self.detectors:
            for epoch in clear_times:
                ts = TimeSeries(numpy.zeros(int(10 * self.sample_rate)),
                                delta_t=1/self.sample_rate,
                                epoch=lal.LIGOTimeGPS(epoch),
                                dtype=numpy.float64)
                injections.apply(ts, det.name)
                max_amp, max_loc = ts.abs_max_loc()
                self.assertEqual(max_amp, 0)

def _fill_sngl_inspiral_row(row, mass1, mass2, spin1z, spin2z):
    # Every column must be set or writing the document fails.
    for col_key, col_type in lsctables.SnglInspiralTable.validcolumns.items():
        col = col_key.split(':')[-1]
        setattr(row, col, 0 if 'int' in col_type else
               0.0 if 'real' in col_type else '')
    row.mass1 = mass1
    row.mass2 = mass2
    row.spin1z = spin1z
    row.spin2z = spin2z
    row.ifo = 'H1'


class TestReadInjectionTable(unittest.TestCase):
    """Tests for pycbc.inject's shared HDF/XML table-reading helpers:
    read_injection_table, read_hdf_param_table, get_table_column.
    """

    def setUp(self):
        self.mass1 = numpy.array([1.4, 30.0])
        self.mass2 = numpy.array([1.3, 25.0])
        self.spin1z = numpy.array([0.1, -0.2])
        self.spin2z = numpy.array([0.0, 0.3])

        # sim_inspiral XML (injection set)
        xmldoc = ligolw.Document()
        xmldoc.appendChild(ligolw.LIGO_LW())
        sim_table = lsctables.SimInspiralTable.new()
        xmldoc.childNodes[-1].appendChild(sim_table)
        for i in range(len(self.mass1)):
            inj = MyInjection()
            inj.end_time = 1000000000 + i
            inj.mass1 = self.mass1[i]
            inj.mass2 = self.mass2[i]
            inj.distance = 1e6 * lal.PC_SI
            inj.latitude = 0.1
            inj.longitude = 0.2
            inj.inclination = 0.3
            inj.polarization = 0.4
            inj.taper = 'TAPER_NONE'
            row = sim_table.RowType()
            inj.fill_sim_inspiral_row(row)
            row.spin1z = self.spin1z[i]
            row.spin2z = self.spin2z[i]
            row.process_id = 0
            row.simulation_id = i
            sim_table.append(row)
        self.sim_xml_file = tempfile.NamedTemporaryFile(suffix='.xml')
        ligolw_utils.write_fileobj(xmldoc, self.sim_xml_file)

        # sngl_inspiral-only XML (template bank)
        xmldoc = ligolw.Document()
        xmldoc.appendChild(ligolw.LIGO_LW())
        sngl_table = lsctables.SnglInspiralTable.new()
        xmldoc.childNodes[-1].appendChild(sngl_table)
        for i in range(len(self.mass1)):
            row = sngl_table.RowType()
            _fill_sngl_inspiral_row(row, self.mass1[i], self.mass2[i],
                                    self.spin1z[i], self.spin2z[i])
            row.event_id = i
            sngl_table.append(row)
        self.sngl_xml_file = tempfile.NamedTemporaryFile(suffix='.xml')
        ligolw_utils.write_fileobj(xmldoc, self.sngl_xml_file)

        # HDF injection set (pycbc_create_injections style: varying
        # params as datasets, constants via 'static_args', requires 'tc')
        self.hdf_injset_file = tempfile.NamedTemporaryFile(suffix='.hdf')
        with HFile(self.hdf_injset_file.name, 'w') as f:
            f['mass1'] = self.mass1
            f['mass2'] = self.mass2
            f['spin1z'] = self.spin1z
            f['spin2z'] = self.spin2z
            f['tc'] = numpy.array([1000000000.0, 1000000001.0])
            f.attrs['static_args'] = ['approximant', 'f_lower']
            f.attrs['approximant'] = 'IMRPhenomD'
            f.attrs['f_lower'] = 20.0

        # Plain HDF parameter table (e.g. an HDF template bank): no 'tc',
        # so cannot be read as an injection set.
        self.hdf_table_file = tempfile.NamedTemporaryFile(suffix='.hdf')
        with HFile(self.hdf_table_file.name, 'w') as f:
            f.attrs['parameters'] = ['mass1', 'mass2', 'spin1z', 'spin2z']
            f['mass1'] = self.mass1
            f['mass2'] = self.mass2
            f['spin1z'] = self.spin1z
            f['spin2z'] = self.spin2z

    def test_read_injection_table_xml(self):
        """XML files are read as sim_inspiral by default, falling back to
        sngl_inspiral (e.g. a template bank) if that's what's present.
        """
        for xml_file in (self.sim_xml_file, self.sngl_xml_file):
            table = read_injection_table(xml_file.name)
            self.assertEqual(len(table), len(self.mass1))
            numpy.testing.assert_array_almost_equal(
                [row.mass1 for row in table], self.mass1)
            numpy.testing.assert_array_almost_equal(
                [row.spin1z for row in table], self.spin1z)

    def test_read_injection_table_xml_tables_order(self):
        """xml_tables controls which table is preferred; a file missing
        every listed table raises.
        """
        table = read_injection_table(
            self.sngl_xml_file.name, xml_tables=('sngl_inspiral',)
        )
        self.assertEqual(len(table), len(self.mass1))
        with self.assertRaises(ValueError):
            read_injection_table(
                self.sngl_xml_file.name, xml_tables=('sim_inspiral',)
            )

    def test_read_injection_table_hdf(self):
        """HDF files with 'tc' are read via InjectionSet (including
        static parameters); without 'tc', falls back to a plain table.
        """
        table = read_injection_table(self.hdf_injset_file.name)
        self.assertEqual(len(table), len(self.mass1))
        numpy.testing.assert_array_almost_equal(table['mass1'], self.mass1)
        numpy.testing.assert_array_almost_equal(
            table['f_lower'], numpy.full(len(self.mass1), 20.0))

        table = read_injection_table(self.hdf_table_file.name)
        self.assertEqual(len(table), len(self.mass1))
        numpy.testing.assert_array_almost_equal(table['mass1'], self.mass1)

        # read_hdf_param_table backs the fallback above; check directly too.
        table = read_hdf_param_table(self.hdf_table_file.name)
        self.assertIsInstance(table, WaveformArray)
        numpy.testing.assert_array_almost_equal(table['spin2z'],
                                                self.spin2z)

    def test_get_table_column(self):
        """get_table_column reads a column the same way from HDF and XML
        tables, defaulting a field missing from an HDF table (e.g.
        spin1x for aligned-spin-only) to zero like a ligolw row would.
        """
        hdf_table = read_hdf_param_table(self.hdf_table_file.name)
        numpy.testing.assert_array_almost_equal(
            get_table_column(hdf_table, 'mass1'), self.mass1)
        self.assertNotIn('spin1x', hdf_table.fieldnames)
        numpy.testing.assert_array_equal(
            get_table_column(hdf_table, 'spin1x'),
            numpy.zeros(len(self.mass1)))

        xml_table = read_injection_table(self.sim_xml_file.name)
        numpy.testing.assert_array_almost_equal(
            get_table_column(xml_table, 'mass1'), self.mass1)


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestInjection))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestReadInjectionTable))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
