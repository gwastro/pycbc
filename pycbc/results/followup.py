# Copyright (C) 2014 Alex Nitz
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

#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
""" This module provides functions to generate followup plots and trigger
time series.
"""
import h5py, numpy, matplotlib
# Only if a backend is not already set ... This should really *not* be done
# here, but in the executables you should set matplotlib.use()
# This matches the check that matplotlib does internally, but this *may* be
# version dependenant. If this is a problem then remove this and control from
# the executables directly.
import sys
if 'matplotlib.backends' not in sys.modules:
    matplotlib.use('agg')
import pylab, mpld3, mpld3.plugins
from ligo.segments import segment

def columns_from_file_list(file_list, columns, ifo, start, end):
    """ Return columns of information stored in single detector trigger
    files.

    Parameters
    ----------
    file_list_file : string
        pickle file containing the list of single detector
    triggers.
    ifo : string
        The ifo to return triggers for.
    columns : list of strings
        The list of columns to read from the trigger files.
    start : int
        The start time to get triggers from
    end : int
        The end time to get triggers from

    Returns
    -------
    trigger_dict : dict
        A dictionary of column vectors with column names as keys.
    """
    file_list = file_list.find_output_with_ifo(ifo)
    file_list = file_list.find_all_output_in_range(ifo, segment(start, end))

    trig_dict = {}
    for trig_file in file_list:
        f = h5py.File(trig_file.storage_path, 'r')

        time = f['end_time'][:]
        pick = numpy.logical_and(time < end, time > start)
        pick_loc = numpy.where(pick)[0]

        for col in columns:
            if col not in trig_dict:
                trig_dict[col] = []
            trig_dict[col] = numpy.concatenate([trig_dict[col], f[col][:][pick_loc]])

    return trig_dict

ifo_color = {'H1': 'blue', 'L1':'red', 'V1':'green'}

def coinc_timeseries_plot(coinc_file, start, end):
    fig = pylab.figure()
    f = h5py.File(coinc_file, 'r')

    stat1 = f['foreground/stat1']
    stat2 = f['foreground/stat2']
    time1 = f['foreground/time1']
    time2 = f['foreground/time2']
    ifo1 = f.attrs['detector_1']
    ifo2 = f.attrs['detector_2']

    pylab.scatter(time1, stat1, label=ifo1, color=ifo_color[ifo1])
    pylab.scatter(time2, stat2, label=ifo2, color=ifo_color[ifo2])

    fmt = '.12g'
    mpld3.plugins.connect(fig, mpld3.plugins.MousePosition(fmt=fmt))
    pylab.legend()
    pylab.xlabel('Time (s)')
    pylab.ylabel('NewSNR')
    pylab.grid()
    return mpld3.fig_to_html(fig)

def trigger_timeseries_plot(file_list, ifos, start, end):

    fig = pylab.figure()
    for ifo in ifos:
        trigs = columns_from_file_list(file_list,
                                       ['snr', 'end_time'],
                                       ifo, start, end)
        print(trigs)
        pylab.scatter(trigs['end_time'], trigs['snr'], label=ifo,
                      color=ifo_color[ifo])

        fmt = '.12g'
        mpld3.plugins.connect(fig, mpld3.plugins.MousePosition(fmt=fmt))
    pylab.legend()
    pylab.xlabel('Time (s)')
    pylab.ylabel('SNR')
    pylab.grid()
    return mpld3.fig_to_html(fig)

def times_to_urls(times, window, tag):
    base = '/../followup/%s/%s/%s'
    return times_to_links(times, window, tag, base=base)

def times_to_links(times, window, tag, base=None):
    if base is None:
        base = "<a href='/../followup/%s/%s/%s' target='_blank'>followup</a>"

    urls = []
    for time in times:
        start = time - window
        end = time + window
        urls.append(base % (tag, start, end))
    return urls

