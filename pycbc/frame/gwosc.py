# Copyright (C) 2017 Alex Nitz
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Generals
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""
This modules contains functions for getting data from the Gravitational Wave
Open Science Center (GWOSC).
"""
import logging

from gwosc.api import fetch_run_json
from gwosc.datasets import find_datasets
from gwosc.locate import get_urls
from gwosc.utils import url_segment

from pycbc.io import get_file
from pycbc.frame import read_frame

logger = logging.getLogger('pycbc.frame.gwosc')

_GWOSC_SAMPLE_RATES = {
    4096: 'GWOSC-4KHZ_R1_STRAIN',
    16384: 'GWOSC-16KHZ_R1_STRAIN',
}


def get_run(time, ifo=None):
    """Return the run name for a given time.

    Parameters
    ----------
    time: int
        The GPS time.
    ifo: str
        The interferometer prefix string. Optional and normally unused,
        except for some special times where data releases were made for a
        single detector under unusual circumstances. For example, to get
        the data around GW170608 in the Hanford detector.
    """
    cases = [
        (
            # ifo is only needed in this special case, otherwise,
            # the run name is the same for all ifos
            1180911618 <= time <= 1180982427 and ifo == 'H1',
            'BKGW170608_16KHZ_R1'
        ),
        (1396417050 <= time <= 1422118818, 'O4b_16KHZ_R1'),
        (1368195220 <= time <= 1389456018, 'O4a_16KHZ_R1'),
        (1253977219 <= time <= 1320363336, 'O3b_16KHZ_R1'),
        (1238166018 <= time <= 1253977218, 'O3a_16KHZ_R1'),
        (1164556817 <= time <= 1187733618, 'O2_16KHZ_R1'),
        (1126051217 <= time <= 1137254417, 'O1'),
        (815011213 <= time <= 875318414, 'S5'),
        (930787215 <= time <= 971568015, 'S6')
    ]
    for condition, name in cases:
        if condition:
            return name
    raise ValueError(f'Time {time} not available in a public dataset')


def _get_sample_rate(time, sample_rate):
    if sample_rate is None:
        sample_rate = 4096 if time < 1164556817 else 16384

    if sample_rate not in _GWOSC_SAMPLE_RATES:
        rates = ', '.join(str(rate) for rate in _GWOSC_SAMPLE_RATES)
        raise ValueError(
            f'Unsupported GWOSC sample rate {sample_rate}; choose {rates}'
        )
    return sample_rate


def _get_channel(time, sample_rate=None):
    sample_rate = _get_sample_rate(time, sample_rate)
    if time < 1164556817:
        return 'LOSC-STRAIN'
    return _GWOSC_SAMPLE_RATES[sample_rate]


def gwosc_frame_json(ifo, start_time, end_time):
    """Get the information about the public data files in a duration of time.

    Parameters
    ----------
    ifo: str
        The name of the interferometer to find the information about.
    start_time: int
        The start time in GPS seconds.
    end_time: int
        The end time in GPS seconds.

    Returns
    -------
    info: dict
        A dictionary containing information about the files that span the
        requested times.
    """
    run = get_run(start_time, ifo)
    run2 = get_run(end_time, ifo)
    if run != run2:
        raise ValueError(
            'Spanning multiple runs is not currently supported. '
            f'You have requested data that uses both {run} and {run2}'
        )

    try:
        return fetch_run_json(
            run,
            ifo,
            gpsstart=int(start_time),
            gpsend=int(end_time),
        )
    except Exception as exc:
        msg = ('Failed to find gwf files for '
               f'ifo={ifo}, run={run}, between {start_time}-{end_time}')
        raise ValueError(msg) from exc


def gwosc_frame_urls(ifo, start_time, end_time, sample_rate=None):
    """Get a list of URLs to GWOSC frame files.

    Parameters
    ----------
    ifo: str
        The name of the interferometer to find the information about.
    start_time: int
        The start time in GPS seconds.
    end_time: int
        The end time in GPS seconds.
    sample_rate: int, optional
        Sample rate of the requested frame files in Hz. GWOSC strain data are
        available at 4096 Hz and 16384 Hz. By default this preserves the
        previous behavior: 4096 Hz before O2 and 16384 Hz from O2 onward.

    Returns
    -------
    frame_files: list
        URLs of frame files that span the requested times.
    """
    start_time = int(start_time)
    end_time = int(end_time)
    sample_rate = _get_sample_rate(start_time, sample_rate)
    datasets = find_datasets(
        detector=ifo,
        type='run',
        segment=(start_time, end_time),
    )

    # Restrict discovery to observing-run datasets where possible. Without
    # this, an interval around an event can return overlapping event files of
    # several durations instead of the contiguous run frames used here.
    error = None
    for dataset in datasets:
        try:
            return get_urls(
                ifo,
                start_time,
                end_time,
                dataset=dataset,
                sample_rate=sample_rate,
                format='gwf',
            )
        except ValueError as exc:
            error = exc

    # H1 data around GW170608 were published outside the normal O2 release.
    # Keep using that release for the historical 16 kHz default.
    if error is not None and sample_rate == 16384:
        try:
            if get_run(start_time, ifo) == 'BKGW170608_16KHZ_R1':
                data = gwosc_frame_json(ifo, start_time, end_time)['strain']
                return [item['url'] for item in data
                        if item['format'] == 'gwf']
        except ValueError:
            pass

    try:
        urls = get_urls(
            ifo,
            start_time,
            end_time,
            sample_rate=sample_rate,
            format='gwf',
        )
    except ValueError as exc:
        if error is not None:
            raise error from exc
        raise

    # Event releases may include both short and long files covering the same
    # interval. Use the shortest single file that spans the request.
    covering = [url for url in urls
                if url_segment(url)[0] <= start_time
                and url_segment(url)[1] >= end_time]
    if covering:
        return [min(covering, key=lambda url: url_segment(url)[1]
                    - url_segment(url)[0])]
    return urls


def read_frame_gwosc(channels, start_time, end_time, sample_rate=None):
    """Read channels from GWOSC data.

    Parameters
    ----------
    channels: str or list
        The channel name to read or list of channel names.
    start_time: int
        The start time in GPS seconds.
    end_time: int
        The end time in GPS seconds.
    sample_rate: int, optional
        Sample rate of the requested frame files in Hz. By default this is
        4096 Hz before O2 and 16384 Hz from O2 onward, matching the previous
        behavior.

    Returns
    -------
    ts: TimeSeries
        Returns a timeseries or list of timeseries with the requested data.
    """
    if not isinstance(channels, list):
        channels = [channels]
    sample_rate = _get_sample_rate(start_time, sample_rate)
    ifos = [c[0:2] for c in channels]
    urls = {}
    for ifo in ifos:
        urls[ifo] = gwosc_frame_urls(
            ifo,
            start_time,
            end_time,
            sample_rate=sample_rate,
        )
        if len(urls[ifo]) == 0:
            raise ValueError("No data found for %s so we "
                             "can't produce a time series" % ifo)

    fnames = {ifo: [] for ifo in ifos}
    for ifo in ifos:
        for url in urls[ifo]:
            fname = get_file(url, cache=True)
            fnames[ifo].append(fname)

    ts_list = [read_frame(fnames[channel[0:2]], channel,
                          start_time=start_time, end_time=end_time)
               for channel in channels]
    if len(ts_list) == 1:
        return ts_list[0]
    return ts_list


def read_strain_gwosc(ifo, start_time, end_time, sample_rate=None):
    """Get the strain data from the GWOSC data.

    Parameters
    ----------
    ifo: str
        The name of the interferometer to read data for. Ex. 'H1', 'L1', 'V1'.
    start_time: int
        The start time in GPS seconds.
    end_time: int
        The end time in GPS seconds.
    sample_rate: int, optional
        Sample rate of the requested strain in Hz. Supported values are 4096
        and 16384. By default this is 4096 Hz before O2 and 16384 Hz from O2
        onward, matching the previous behavior.

    Returns
    -------
    ts: TimeSeries
        Returns a timeseries with the strain data.
    """
    sample_rate = _get_sample_rate(start_time, sample_rate)
    channel = _get_channel(start_time, sample_rate)
    return read_frame_gwosc(
        f'{ifo}:{channel}',
        start_time,
        end_time,
        sample_rate=sample_rate,
    )
