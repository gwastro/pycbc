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

import unittest
from unittest import mock

from pycbc.frame import gwosc


class GWOSCFrameTest(unittest.TestCase):
    @mock.patch('pycbc.frame.gwosc.fetch_run_json')
    def test_frame_json_uses_gwosc_client(self, fetch_run_json):
        expected = {'strain': []}
        fetch_run_json.return_value = expected

        self.assertIs(
            gwosc.gwosc_frame_json('H1', 1180922494, 1180922495),
            expected,
        )
        fetch_run_json.assert_called_once_with(
            'BKGW170608_16KHZ_R1',
            'H1',
            gpsstart=1180922494,
            gpsend=1180922495,
        )

    def test_frame_json_rejects_multiple_runs(self):
        with self.assertRaisesRegex(ValueError, 'Spanning multiple runs'):
            gwosc.gwosc_frame_json('L1', 1187733618, 1238166018)

    @mock.patch('pycbc.frame.gwosc.fetch_run_json')
    def test_frame_json_wraps_client_error(self, fetch_run_json):
        fetch_run_json.side_effect = RuntimeError('request failed')

        with self.assertRaisesRegex(ValueError, 'Failed to find gwf files'):
            gwosc.gwosc_frame_json('L1', 1238166018, 1238166020)

    @mock.patch('pycbc.frame.gwosc.get_urls')
    @mock.patch('pycbc.frame.gwosc.find_datasets')
    def test_frame_urls_sample_rate(self, find_datasets, get_urls):
        expected = ['https://gwosc.org/data.gwf']
        find_datasets.return_value = ['O3a']
        get_urls.return_value = expected

        self.assertEqual(
            gwosc.gwosc_frame_urls('H1', 1234.5, 1240.5, sample_rate=4096),
            expected,
        )
        find_datasets.assert_called_once_with(
            detector='H1', type='run', segment=(1234, 1240)
        )
        get_urls.assert_called_once_with(
            'H1',
            1234,
            1240,
            dataset='O3a',
            sample_rate=4096,
            format='gwf',
        )

    @mock.patch('pycbc.frame.gwosc.get_urls')
    @mock.patch('pycbc.frame.gwosc.find_datasets')
    def test_frame_urls_default_sample_rate(self, find_datasets, get_urls):
        find_datasets.return_value = []
        get_urls.return_value = []

        gwosc.gwosc_frame_urls('L1', 1238166018, 1238166020)

        get_urls.assert_called_once_with(
            'L1',
            1238166018,
            1238166020,
            sample_rate=16384,
            format='gwf',
        )

    @mock.patch('pycbc.frame.gwosc.get_urls')
    @mock.patch('pycbc.frame.gwosc.find_datasets')
    def test_frame_urls_legacy_default_is_4khz(
        self, find_datasets, get_urls
    ):
        find_datasets.return_value = ['O1']
        get_urls.return_value = []

        gwosc.gwosc_frame_urls('H1', 1126259462, 1126259466)

        get_urls.assert_called_once_with(
            'H1',
            1126259462,
            1126259466,
            dataset='O1',
            sample_rate=4096,
            format='gwf',
        )

    @mock.patch('pycbc.frame.gwosc.get_urls')
    @mock.patch('pycbc.frame.gwosc.find_datasets')
    def test_frame_urls_require_full_run_coverage(
        self, find_datasets, get_urls
    ):
        find_datasets.return_value = ['O2', 'O3a']
        first_error = ValueError('O2 does not cover the requested interval')
        expected = ['https://gwosc.org/data.gwf']
        get_urls.side_effect = [first_error, expected]

        self.assertEqual(
            gwosc.gwosc_frame_urls('H1', 1234, 1240),
            expected,
        )
        self.assertEqual(get_urls.call_count, 2)
        self.assertEqual(get_urls.call_args_list[0].kwargs['dataset'], 'O2')
        self.assertEqual(get_urls.call_args_list[1].kwargs['dataset'], 'O3a')

    @mock.patch('pycbc.frame.gwosc.get_urls')
    @mock.patch('pycbc.frame.gwosc.find_datasets')
    def test_frame_urls_report_missing_run_coverage(
        self, find_datasets, get_urls
    ):
        find_datasets.return_value = ['O2']
        expected = ValueError('run does not cover the requested interval')
        fallback_error = ValueError('no event data cover the interval')
        get_urls.side_effect = [expected, fallback_error]

        with self.assertRaises(ValueError) as context:
            gwosc.gwosc_frame_urls('H1', 1234, 1240)

        self.assertIs(context.exception, expected)
        self.assertIs(context.exception.__cause__, fallback_error)

    @mock.patch('pycbc.frame.gwosc.gwosc_frame_json')
    @mock.patch('pycbc.frame.gwosc.get_urls')
    @mock.patch('pycbc.frame.gwosc.find_datasets')
    def test_gw170608_keeps_background_release(
        self, find_datasets, get_urls, frame_json
    ):
        find_datasets.return_value = ['O2']
        get_urls.side_effect = ValueError('O2 has no H1 data here')
        url = 'https://gwosc.org/archive/data/BKGW170608.gwf'
        frame_json.return_value = {
            'strain': [{'format': 'gwf', 'url': url}],
        }

        self.assertEqual(
            gwosc.gwosc_frame_urls('H1', 1180922494, 1180922498),
            [url],
        )
        frame_json.assert_called_once_with('H1', 1180922494, 1180922498)

    @mock.patch('pycbc.frame.gwosc.get_urls')
    @mock.patch('pycbc.frame.gwosc.find_datasets')
    def test_gw170608_4khz_selects_short_event_file(
        self, find_datasets, get_urls
    ):
        find_datasets.return_value = ['O2']
        base = 'https://gwosc.org/eventapi/json/GW170608/'
        long_url = base + 'H-H1_GWOSC_4KHZ_R1-1180920447-4096.gwf'
        short_url = base + 'H-H1_GWOSC_4KHZ_R1-1180922479-32.gwf'
        get_urls.side_effect = [
            ValueError('O2 has no H1 data here'),
            [long_url, short_url],
        ]

        self.assertEqual(
            gwosc.gwosc_frame_urls(
                'H1', 1180922494, 1180922498, sample_rate=4096
            ),
            [short_url],
        )

    @mock.patch('pycbc.frame.gwosc.read_frame')
    @mock.patch('pycbc.frame.gwosc.get_file')
    @mock.patch('pycbc.frame.gwosc.gwosc_frame_urls')
    def test_read_frame_uses_pycbc_downloader(
        self, frame_urls, get_file, read_frame
    ):
        frame_urls.return_value = [
            'https://gwosc.org/first.gwf',
            'https://gwosc.org/second.gwf',
        ]
        get_file.side_effect = ['first.gwf', 'second.gwf']
        read_frame.return_value = object()

        result = gwosc.read_frame_gwosc(
            'H1:GWOSC-4KHZ_R1_STRAIN', 1234, 1240, sample_rate=4096
        )

        self.assertIs(result, read_frame.return_value)
        frame_urls.assert_called_once_with(
            'H1', 1234, 1240, sample_rate=4096
        )
        get_file.assert_has_calls([
            mock.call('https://gwosc.org/first.gwf', cache=True),
            mock.call('https://gwosc.org/second.gwf', cache=True),
        ])
        read_frame.assert_called_once_with(
            ['first.gwf', 'second.gwf'],
            'H1:GWOSC-4KHZ_R1_STRAIN',
            start_time=1234,
            end_time=1240,
        )

    @mock.patch('pycbc.frame.gwosc.gwosc_frame_urls', return_value=[])
    def test_read_frame_reports_missing_data(self, frame_urls):
        with self.assertRaisesRegex(ValueError, 'No data found for H1'):
            gwosc.read_frame_gwosc(
                'H1:GWOSC-4KHZ_R1_STRAIN',
                1238166018,
                1238166020,
                sample_rate=4096,
            )

    @mock.patch('pycbc.frame.gwosc.read_frame')
    @mock.patch('pycbc.frame.gwosc.get_file')
    @mock.patch('pycbc.frame.gwosc.gwosc_frame_urls')
    def test_read_frame_multiple_channels(
        self, frame_urls, get_file, read_frame
    ):
        frame_urls.side_effect = lambda ifo, *args, **kwargs: [f'{ifo}.gwf']
        get_file.side_effect = lambda url, **kwargs: url
        read_frame.side_effect = ['H1 strain', 'L1 strain']
        channels = [
            'H1:GWOSC-4KHZ_R1_STRAIN',
            'L1:GWOSC-4KHZ_R1_STRAIN',
        ]

        result = gwosc.read_frame_gwosc(
            channels, 1238166018, 1238166020, sample_rate=4096
        )

        self.assertEqual(result, ['H1 strain', 'L1 strain'])
        self.assertEqual(frame_urls.call_count, 2)

    @mock.patch('pycbc.frame.gwosc.read_frame_gwosc')
    def test_read_strain_selects_4khz_channel(self, read_frame_gwosc):
        gwosc.read_strain_gwosc(
            'H1', 1238166018, 1238166020, sample_rate=4096
        )

        read_frame_gwosc.assert_called_once_with(
            'H1:GWOSC-4KHZ_R1_STRAIN',
            1238166018,
            1238166020,
            sample_rate=4096,
        )

    @mock.patch('pycbc.frame.gwosc.read_frame_gwosc')
    def test_read_strain_default_is_16khz(self, read_frame_gwosc):
        gwosc.read_strain_gwosc('L1', 1238166018, 1238166020)

        read_frame_gwosc.assert_called_once_with(
            'L1:GWOSC-16KHZ_R1_STRAIN',
            1238166018,
            1238166020,
            sample_rate=16384,
        )

    @mock.patch('pycbc.frame.gwosc.read_frame_gwosc')
    def test_read_strain_legacy_default_is_4khz(self, read_frame_gwosc):
        gwosc.read_strain_gwosc('H1', 1126259462, 1126259466)

        read_frame_gwosc.assert_called_once_with(
            'H1:LOSC-STRAIN',
            1126259462,
            1126259466,
            sample_rate=4096,
        )

    def test_invalid_sample_rate(self):
        message = 'Unsupported GWOSC sample rate'
        with self.assertRaisesRegex(ValueError, message):
            gwosc.read_strain_gwosc(
                'H1', 1238166018, 1238166020, sample_rate=8192
            )


if __name__ == '__main__':
    unittest.main()
