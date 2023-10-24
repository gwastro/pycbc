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
import json
from urllib.request import urlopen
from pycbc.io import get_file
from pycbc.frame import read_frame


_GWOSC_URL = "https://www.gwosc.org/archive/links/%s/%s/%s/%s/json/"


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


def _get_channel(time):
    if time < 1164556817:
        return 'LOSC-STRAIN'
    return 'GWOSC-16KHZ_R1_STRAIN'


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
    run = get_run(start_time)
    run2 = get_run(end_time)
    if run != run2:
        raise ValueError(
            'Spanning multiple runs is not currently supported. '
            f'You have requested data that uses both {run} and {run2}'
        )

    url = _GWOSC_URL % (run, ifo, int(start_time), int(end_time))

    try:
        return json.loads(urlopen(url).read().decode())
    except Exception as exc:
        msg = ('Failed to find gwf files for '
               f'ifo={ifo}, run={run}, between {start_time}-{end_time}')
        raise ValueError(msg) from exc


def gwosc_frame_urls(ifo, start_time, end_time):
    """Get a list of URLs to GWOSC frame files.

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
    frame_files: list
        A dictionary containing information about the files that span the
        requested times.
    """
    data = gwosc_frame_json(ifo, start_time, end_time)['strain']
    return [d['url'] for d in data if d['format'] == 'gwf']


def read_frame_gwosc(channels, start_time, end_time):
    """Read channels from GWOSC data.

    Parameters
    ----------
    channels: str or list
        The channel name to read or list of channel names.
    start_time: int
        The start time in GPS seconds.
    end_time: int
        The end time in GPS seconds.

    Returns
    -------
    ts: TimeSeries
        Returns a timeseries or list of timeseries with the requested data.
    """
    if not isinstance(channels, list):
        channels = [channels]
    ifos = [c[0:2] for c in channels]
    urls = {}
    for ifo in ifos:
        urls[ifo] = gwosc_frame_urls(ifo, start_time, end_time)
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


def read_strain_gwosc(ifo, start_time, end_time):
    """Get the strain data from the GWOSC data.

    Parameters
    ----------
    ifo: str
        The name of the interferometer to read data for. Ex. 'H1', 'L1', 'V1'.
    start_time: int
        The start time in GPS seconds.
    end_time: int
        The end time in GPS seconds.

    Returns
    -------
    ts: TimeSeries
        Returns a timeseries with the strain data.
    """
    channel = _get_channel(start_time)
    return read_frame_gwosc(f'{ifo}:{channel}', start_time, end_time)
