# Copyright (C) 2018 Tito Dal Canton
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

import unittest
import os
import shutil
import random
import tempfile
import itertools
import numpy as np
from utils import parse_args_cpu_only, simple_exit
from pycbc.types import TimeSeries, FrequencySeries
from pycbc.io.live import SingleCoincForGraceDB
from glue.ligolw import ligolw
from glue.ligolw import lsctables
from glue.ligolw import table
from glue.ligolw import utils as ligolw_utils
from lal import series as lalseries

# if we have the GraceDb module then we can do deeper tests,
# otherwise just fall back to quicker ones
try:
    from ligo.gracedb.rest import GraceDb
except ImportError:
    GraceDb = None


parse_args_cpu_only("io.live")

class ContentHandler(ligolw.LIGOLWContentHandler):
    pass
lsctables.use_in(ContentHandler)

class TestIOLive(unittest.TestCase):
    def setUp(self):
        self.template = {'template_id': 0,
                         'mass1': 10,
                         'mass2': 11,
                         'spin1x': 0,
                         'spin1y': 0,
                         'spin1z': 0,
                         'spin2x': 0,
                         'spin2y': 0,
                         'spin2z': 0}

        self.possible_ifos = 'H1 L1 V1 K1 I1'.split()

    def do_test(self, n_ifos, n_ifos_followup):
        # choose a random selection of interferometers
        # n_ifos will be used to generate the simulated trigger
        # n_ifos_followup will be used as followup-only
        all_ifos = random.sample(self.possible_ifos, n_ifos + n_ifos_followup)
        trig_ifos = all_ifos[0:n_ifos]

        results = {'foreground/stat': np.random.uniform(4, 20),
                   'foreground/ifar': np.random.uniform(0.01, 1000)}
        followup_data = {}
        for ifo in all_ifos:
            offset = 10000 + np.random.uniform(-0.02, 0.02)
            amplitude = np.random.uniform(4, 20)

            # generate a mock SNR time series with a peak
            n = 201
            dt = 1. / 2048.
            t = np.arange(n) * dt
            t_peak = dt * n / 2
            snr = np.exp(-(t - t_peak) ** 2 * 3e-3 ** -2) * amplitude
            snr_series = TimeSeries((snr + 1j * 0).astype(np.complex64),
                                    delta_t=dt, epoch=offset)

            # generate a mock PSD
            psd_samples = np.random.exponential(size=1024)
            psd = FrequencySeries(psd_samples, delta_f=1.)

            # fill in the various fields
            if ifo in trig_ifos:
                base = 'foreground/' + ifo + '/'
                results[base + 'end_time'] = t_peak + offset
                results[base + 'snr'] = amplitude
                results[base + 'sigmasq'] = np.random.uniform(1e6, 2e6)
            followup_data[ifo] = {'snr_series': snr_series,
                                  'psd': psd}

        for ifo, k in itertools.product(trig_ifos, self.template):
            results['foreground/' + ifo + '/' + k] = self.template[k]

        channel_names = {ifo: 'TEST' for ifo in all_ifos}
        kwargs = {'psds': {ifo: followup_data[ifo]['psd'] for ifo in all_ifos},
                  'low_frequency_cutoff': 20.,
                  'followup_data': followup_data,
                  'channel_names': channel_names}
        coinc = SingleCoincForGraceDB(trig_ifos, results, **kwargs)

        tempdir = tempfile.mkdtemp()

        coinc_file_name = os.path.join(tempdir, 'coinc.xml.gz')

        if GraceDb is not None:
            # pretend to upload the event to GraceDB.
            # The upload will fail, but it should not raise an exception
            # and it should still leave the event file around
            coinc.upload(coinc_file_name, gracedb_server='localhost',
                         testing=True)
        else:
            # no GraceDb module, so just save the coinc file
            coinc.save(coinc_file_name)

        # read back and check the coinc document
        read_coinc = ligolw_utils.load_filename(
                coinc_file_name, verbose=False, contenthandler=ContentHandler)
        single_table = table.get_table(
                read_coinc, lsctables.SnglInspiralTable.tableName)
        self.assertEqual(len(single_table), len(all_ifos))
        coinc_table = table.get_table(
                read_coinc, lsctables.CoincInspiralTable.tableName)
        self.assertEqual(len(coinc_table), 1)

        # make sure lalseries can read the PSDs
        psd_doc = ligolw_utils.load_filename(
                coinc_file_name, verbose=False,
                contenthandler=lalseries.PSDContentHandler)
        psd_dict = lalseries.read_psd_xmldoc(psd_doc)
        self.assertEqual(set(psd_dict.keys()), set(all_ifos))

        shutil.rmtree(tempdir)

    def test_2_ifos_no_followup(self):
        self.do_test(2, 0)

    def test_3_ifos_no_followup(self):
        self.do_test(3, 0)

    def test_4_ifos_no_followup(self):
        self.do_test(4, 0)

    def test_5_ifos_no_followup(self):
        self.do_test(5, 0)

    def test_2_ifos_1_followup(self):
        self.do_test(2, 1)

    def test_2_ifos_2_followup(self):
        self.do_test(2, 2)

    def test_2_ifos_3_followup(self):
        self.do_test(2, 3)

    def test_3_ifos_1_followup(self):
        self.do_test(3, 1)

    def test_3_ifos_2_followup(self):
        self.do_test(3, 2)

    def test_4_ifos_1_followup(self):
        self.do_test(4, 1)


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestIOLive))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
