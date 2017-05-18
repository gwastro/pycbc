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
This modules contains functions for getting data from the LOSC
"""

_losc_url = "https://losc.ligo.org/archive/links/%s/%s/%s/%s/json/"

def _get_run(time):
    if 815011213 < time < 875318414:
        return 'S5'
    elif 930787215 < time < 971568015:
        return 'S6'
    else:
        raise ValueError('Time %s not available in a public dataset' % time)
    

def losc_frame_json(ifo, start_time, end_time):
    """ Get the information about the public data files in a duration of time
    
    Parameters
    ----------
    ifo: str
        The name of the IFO to find the information about.
    start_time: int
        The gps time in GPS seconds
    end_time: int
        The end time in GPS seconds
    
    Returns
    -------
    info: dict
        A dictionary containing information about the files that span the
        requested times.
    """
    import urllib, json
    run = _get_run(start_time)
    run2 = _get_run(end_time)
    if run != run2:
        raise ValueError('Spanning multiple runs is not currently supported.'
                         'You have requested data that uses '
                         'both %s and %s' % (run, run2))
    
    url = _losc_url % (run, ifo, int(start_time), int(end_time))
    
    try:
        return json.loads(urllib.urlopen(url).read())
    except Exception:
        raise ValueError('Failed to find gwf files for '
            'ifo=%s, run=%s, between %s-%s' % (ifo, run, start_time, end_time))
           
def losc_frame_urls(ifo, start_time, end_time):
    """ Get a list of urls to losc frame files
    
    Parameters
    ----------
    ifo: str
        The name of the IFO to find the information about.
    start_time: int
        The gps time in GPS seconds
    end_time: int
        The end time in GPS seconds
    
    Returns
    -------
    frame_files: list
        A dictionary containing information about the files that span the
        requested times.
    """
    data = losc_frame_json(ifo, start_time, end_time)['strain']
    return [d['url'] for d in data if d['format'] == 'gwf']
    
def read_frame_losc(channels, start_time, end_time):
    """ Read channels from losc data

    Parameters
    ----------
    channels: str or list
        The channel name to read or list of channel names.
    start_time: int
        The gps time in GPS seconds
    end_time: int
        The end time in GPS seconds
      
    Returns
    -------
    ts: TimeSeries
        Returns a timeseries or list of timeseries with the requested data.
    """    
    import urllib
    from pycbc.frame import read_frame
    if not isinstance(channels, list):
        channels = [channels]
    ifos = [c[0:2] for c in channels]
    urls = {}
    for ifo in ifos:
        urls[ifo] = losc_frame_urls(ifo, start_time, end_time)
        if len(urls[ifo]) == 0:
            raise ValueError("No data found for %s so we "
                             "can't produce a time series" % ifo)
                             
    fnames = {ifo:[] for ifo in ifos}
    for ifo in ifos:
        for url in urls[ifo]:
            fname, _ = urllib.urlretrieve(url)
            fnames[ifo].append(fname)
    
    ts = [read_frame(fnames[channel[0:2]], channel,
           start_time=start_time, end_time=end_time) for channel in channels]
    if len(ts) == 1:
        return ts[0]
    else:
        return ts
           
def read_strain_losc(ifo, start_time, end_time):
    """ Get the strain data from the LOSC data

    Parameters
    ----------
    ifo: str
        The name of the IFO to read data for. Ex. 'H1', 'L1', 'V1'
    start_time: int
        The gps time in GPS seconds
    end_time: int
        The end time in GPS seconds
      
    Returns
    -------
    ts: TimeSeries
        Returns a timeseries with the strain data.
    """   
    return read_frame_losc('%s:LOSC-STRAIN' % ifo, start_time, end_time)
    
