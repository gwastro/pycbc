"""Exercise offline localization XML with and without follow-up-only detectors."""
import pathlib
import subprocess
import sys
import tempfile
import unittest

import h5py
import lal
import numpy as np
from igwn_ligolw import ligolw, lsctables, utils
from ligo.skymap.io.events.ligolw import LigoLWEventSource
from ligo.skymap.bayestar import condition

import pycbc
from pycbc.io.ligolw import (
    LIGOLWContentHandler, create_process_table, return_empty_sngl,
)
from pycbc.types import TimeSeries


SCRIPT = pathlib.Path(__file__).resolve().parents[1] / 'bin/all_sky_search/pycbc_prepare_xml_for_gracedb'


class TestOfflineLocalizationXML(unittest.TestCase):
    def run_case(self, ifos):
        with tempfile.TemporaryDirectory() as dirname:
            root = pathlib.Path(dirname)
            doc = ligolw.Document()
            doc.appendChild(ligolw.LIGO_LW())
            proc = create_process_table(doc, program_name='pycbc', detectors=['H1', 'L1'])
            tables = {}
            for cls in (lsctables.SnglInspiralTable, lsctables.CoincTable,
                        lsctables.CoincInspiralTable, lsctables.CoincMapTable,
                        lsctables.TimeSlideTable):
                tables[cls] = table = cls.new()
                doc.childNodes[0].appendChild(table)
            for index, ifo in enumerate(['H1', 'L1']):
                single = return_empty_sngl()
                single.event_id = index
                single.process_id = proc.process_id
                single.ifo = ifo
                single.mass1 = single.mass2 = 1.4
                single.f_final = 512
                single.end = lal.LIGOTimeGPS(100)
                single.snr = 8 + index
                tables[lsctables.SnglInspiralTable].append(single)
                mapping = lsctables.CoincMap()
                mapping.coinc_event_id = 0
                mapping.table_name = 'sngl_inspiral'
                mapping.event_id = index
                tables[lsctables.CoincMapTable].append(mapping)
                slide = lsctables.TimeSlide()
                slide.process_id = proc.process_id
                slide.time_slide_id = 0
                slide.instrument = ifo
                slide.offset = 0
                tables[lsctables.TimeSlideTable].append(slide)
            event = lsctables.Coinc()
            event.coinc_event_id = event.coinc_def_id = event.time_slide_id = 0
            event.process_id = proc.process_id
            event.nevents = 2
            event.instruments = 'H1,L1'
            event.likelihood = 12
            tables[lsctables.CoincTable].append(event)
            inspiral = lsctables.CoincInspiral()
            for column in lsctables.CoincInspiralTable.validcolumns:
                setattr(inspiral, ligolw.Column.ColumnName(column), None)
            inspiral.coinc_event_id = 0
            inspiral.end = lal.LIGOTimeGPS(100)
            inspiral.instruments = ['H1', 'L1']
            inspiral.snr = 12
            inspiral.combined_far = 1e-8
            tables[lsctables.CoincInspiralTable].append(inspiral)
            source = root / 'input.xml'
            output = root / 'output.xml'
            utils.write_filename(doc, str(source))
            snr_files, psd_files = [], []
            for ifo in ifos:
                snr_file = root / (ifo + '-snr.hdf')
                series = TimeSeries(np.ones(2049, dtype=np.complex64),
                                    delta_t=1/4096, epoch=99.75)
                series.save(str(snr_file), group='snr')
                with h5py.File(snr_file, 'a') as f:
                    f.attrs['ifo'] = ifo
                    f.attrs['event_time'] = 100.
                psd_file = root / (ifo + '-psd.hdf')
                with h5py.File(psd_file, 'w') as f:
                    f.attrs['low_frequency_cutoff'] = 20.
                    group = f.create_group(ifo)
                    group['start_time'] = [99.]
                    group['end_time'] = [101.]
                    psd = group.create_dataset('psds/0', data=np.full(4097, pycbc.DYN_RANGE_FAC**2 * 1e-46))
                    psd.attrs['delta_f'] = .25
                    psd.attrs['epoch'] = 99
                snr_files.append(str(snr_file))
                psd_files.append(str(psd_file))
            subprocess.run([sys.executable, str(SCRIPT), '--input-file', str(source),
                            '--output-file', str(output), '--event-id', '0',
                            '--snr-timeseries', *snr_files, '--psd-files', *psd_files,
                            '--channel-name', *[ifo + ':STRAIN' for ifo in ifos],
                            '--snr-timeseries-plot', str(root / 'snr.png')],
                           check=True, capture_output=True, text=True)
            result = utils.load_filename(str(output), contenthandler=LIGOLWContentHandler)
            singles = lsctables.SnglInspiralTable.get_table(result)
            self.assertEqual({row.ifo for row in singles}, set(ifos))
            self.assertEqual(len({row.event_id for row in singles}), len(ifos))
            self.assertEqual(len(lsctables.CoincMapTable.get_table(result)), len(ifos))
            self.assertEqual(lsctables.CoincTable.get_table(result)[0].nevents, len(ifos))
            ci = lsctables.CoincInspiralTable.get_table(result)[0]
            self.assertEqual(ci.snr, 12)
            self.assertEqual(ci.combined_far, 1e-8)
            parsed = LigoLWEventSource(str(output), psd_file=str(output), coinc_def=None)[0]
            self.assertEqual({single.detector for single in parsed.singles}, set(ifos))
            for single in parsed.singles:
                self.assertIsNotNone(single.snr_series)
                self.assertGreater(len(single.psd.data.data), 0)
                self.assertEqual(single.time, 100.)
                if single.detector == 'V1':
                    self.assertIsNone(single.snr)
                    self.assertIsNone(single.phase)
            self.assertEqual([row.snr for row in singles if row.ifo in ('H1', 'L1')], [8, 9])
            # Exercise the same detector and time-series preparation as localization.
            conditioned = condition(parsed, waveform='TaylorF2', f_low=20)
            self.assertEqual(conditioned[3].shape[0], len(ifos))

    def test_two_triggered_detectors(self):
        self.run_case(['H1', 'L1'])

    def test_subthreshold_detector(self):
        self.run_case(['H1', 'L1', 'V1'])


if __name__ == '__main__':
    unittest.main()
