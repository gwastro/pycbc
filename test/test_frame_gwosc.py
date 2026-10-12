# Copyright (C) 2026 Luca Cirfeta
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

"""Network-free regression tests for GWOSC frame rate selection."""

import json
import unittest
from unittest import mock

from pycbc.frame import gwosc, query_and_read_frame


class GWOSCFrameTest(unittest.TestCase):
    def test_run_name_preserves_defaults_and_selects_4khz(self):
        cases = [
            (815726592, 'H1', None, 'S5'),
            (930960000, 'H1', None, 'S6'),
            (1126259462, 'H1', None, 'O1'),
            (1126259462, 'H1', 16384, 'O1_16KHZ'),
            (1170000000, 'H1', None, 'O2_16KHZ_R1'),
            (1170000000, 'H1', 4096, 'O2_4KHZ_R1'),
            (1238166018, 'L1', 4096, 'O3a_4KHZ_R1'),
            (1368195220, 'H1', 4096, 'O4a_4KHZ_R1'),
            (1180922494, 'H1', None, 'BKGW170608_16KHZ_R1'),
            (1180922494, 'L1', None, 'O2_16KHZ_R1'),
        ]
        for time, ifo, rate, expected in cases:
            with self.subTest(time=time, ifo=ifo, rate=rate):
                self.assertEqual(gwosc.get_run(time, ifo, rate), expected)

    def test_unpublished_and_invalid_rates(self):
        for time, ifo, rate in [
            (815726592, 'H1', 16384),
            (930960000, 'H1', 16384),
            (1180922494, 'H1', 4096),
        ]:
            with self.subTest(time=time, ifo=ifo, rate=rate):
                with self.assertRaisesRegex(ValueError, 'not published'):
                    gwosc.get_run(time, ifo, rate)
        with self.assertRaisesRegex(ValueError, 'Unsupported GWOSC sample'):
            gwosc.get_run(1170000000, 'H1', 8192)
        with self.assertRaisesRegex(ValueError, 'not available'):
            gwosc.get_run(1000000000, 'H1')

    @mock.patch('pycbc.frame.gwosc.get_file', return_value='metadata.json')
    def test_json_uses_pycbc_mirrored_downloader(self, get_file):
        expected = {'strain': [{'format': 'gwf', 'sampling_rate': 4096,
                                'url': 'https://gwosc.org/strain.gwf'}]}
        with mock.patch('builtins.open', mock.mock_open(
            read_data=json.dumps(expected)
        )) as open_file:
            actual = gwosc.gwosc_frame_json(
                'H1', 1170000000, 1170000004, sample_rate=4096
            )

        self.assertEqual(actual, expected)
        get_file.assert_called_once_with(
            'https://www.gwosc.org/archive/links/'
            'O2_4KHZ_R1/H1/1170000000/1170000004/json/', cache=False
        )
        open_file.assert_called_once_with('metadata.json', 'r')

    def test_json_rejects_multiple_runs_before_network(self):
        with mock.patch('pycbc.frame.gwosc.get_file') as get_file:
            with self.assertRaisesRegex(ValueError, 'Spanning multiple runs'):
                gwosc.gwosc_frame_json('L1', 1187733618, 1238166018)
        get_file.assert_not_called()

    @mock.patch('pycbc.frame.gwosc.get_file', side_effect=OSError('offline'))
    def test_json_preserves_download_error(self, get_file):
        with self.assertRaisesRegex(ValueError, 'Failed to find gwf files') as ctx:
            gwosc.gwosc_frame_json('L1', 1238166018, 1238166020)
        self.assertIsInstance(ctx.exception.__cause__, OSError)
        get_file.assert_called_once()

    @mock.patch('pycbc.frame.gwosc.gwosc_frame_json')
    def test_urls_filter_exact_rate_and_format(self, frame_json):
        frame_json.return_value = {'strain': [
            {'format': 'hdf5', 'sampling_rate': 4096, 'url': 'data.hdf5'},
            {'format': 'gwf', 'sampling_rate': 16384, 'url': '16k.gwf'},
            {'format': 'gwf', 'sampling_rate': 4096, 'url': '4k.gwf'},
        ]}
        self.assertEqual(
            gwosc.gwosc_frame_urls('H1', 1170000000, 1170000004, 4096),
            ['4k.gwf'],
        )
        frame_json.assert_called_once_with('H1', 1170000000, 1170000004,
                                           4096)

    @mock.patch('pycbc.frame.gwosc.read_frame')
    @mock.patch('pycbc.frame.gwosc.get_file', return_value='data.gwf')
    @mock.patch('pycbc.frame.gwosc.gwosc_frame_urls',
                return_value=['https://gwosc.org/data.gwf'])
    def test_read_frame_uses_pycbc_downloader(
        self, frame_urls, get_file, read_frame
    ):
        result = gwosc.read_frame_gwosc(
            'H1:GWOSC-4KHZ_R1_STRAIN', 1170000000, 1170000004, 4096
        )
        self.assertIs(result, read_frame.return_value)
        frame_urls.assert_called_once_with('H1', 1170000000, 1170000004,
                                           4096)
        get_file.assert_called_once_with('https://gwosc.org/data.gwf',
                                         cache=True)
        read_frame.assert_called_once_with(
            ['data.gwf'], 'H1:GWOSC-4KHZ_R1_STRAIN',
            start_time=1170000000, end_time=1170000004
        )

    @mock.patch('pycbc.frame.gwosc.gwosc_frame_urls', return_value=[])
    def test_read_frame_reports_missing_data(self, frame_urls):
        with self.assertRaisesRegex(ValueError, 'No data found for H1'):
            gwosc.read_frame_gwosc('H1:GWOSC-4KHZ_R1_STRAIN',
                                   1170000000, 1170000004, 4096)

    @mock.patch('pycbc.frame.gwosc.read_frame')
    @mock.patch('pycbc.frame.gwosc.get_file', side_effect=lambda url, **_: url)
    @mock.patch('pycbc.frame.gwosc.gwosc_frame_urls',
                side_effect=lambda ifo, *_: [f'{ifo}.gwf'])
    def test_read_frame_multiple_channels(self, frame_urls, get_file,
                                          read_frame):
        channels = ['H1:GWOSC-4KHZ_R1_STRAIN',
                    'L1:GWOSC-4KHZ_R1_STRAIN']
        read_frame.side_effect = ['H1 strain', 'L1 strain']
        self.assertEqual(
            gwosc.read_frame_gwosc(channels, 1170000000, 1170000004, 4096),
            ['H1 strain', 'L1 strain'],
        )
        self.assertEqual(frame_urls.call_count, 2)

    @mock.patch('pycbc.frame.gwosc.read_frame_gwosc')
    def test_read_strain_channel_and_default(self, read_frame_gwosc):
        cases = [
            (1170000000, None, 'GWOSC-16KHZ_R1_STRAIN'),
            (1170000000, 4096, 'GWOSC-4KHZ_R1_STRAIN'),
            (1126259462, None, 'LOSC-STRAIN'),
            (1126259462, 16384, 'GWOSC-16KHZ_R1_STRAIN'),
        ]
        for time, rate, channel in cases:
            with self.subTest(time=time, rate=rate):
                gwosc.read_strain_gwosc('H1', time, time + 4, rate)
                read_frame_gwosc.assert_called_with(
                    f'H1:{channel}', time, time + 4, rate
                )
        self.assertEqual(read_frame_gwosc.call_count, len(cases))

    def test_read_strain_rejects_invalid_rate(self):
        with self.assertRaisesRegex(ValueError, 'Unsupported GWOSC sample'):
            gwosc.read_strain_gwosc('H1', 1170000000, 1170000004, 8192)

    @mock.patch('pycbc.frame.gwosc.read_strain_gwosc')
    def test_query_forwards_rate_to_strain_reader(self, read_strain_gwosc):
        query_and_read_frame('GWOSC_STRAIN', 'H1:GWOSC-4KHZ_R1_STRAIN',
                             1170000000, 1170000004, sample_rate=4096)
        read_strain_gwosc.assert_called_once_with(
            'H1', 1170000000, 1170000004, 4096
        )

    @mock.patch('pycbc.frame.gwosc.read_frame_gwosc')
    def test_query_forwards_rate_to_frame_reader(self, read_frame_gwosc):
        query_and_read_frame('GWOSC', 'H1:GWOSC-4KHZ_R1_STRAIN',
                             1170000000, 1170000004, sample_rate=4096)
        read_frame_gwosc.assert_called_once_with(
            'H1:GWOSC-4KHZ_R1_STRAIN', 1170000000, 1170000004, 4096
        )


if __name__ == '__main__':
    unittest.main()
