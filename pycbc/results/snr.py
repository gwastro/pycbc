# Copyright (C) 2022
# Tito Dal Canton, Gareth Cabourn Davies, Ian Harry
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
"""
Module to generate PSD figures
"""
from pycbc.results import ifo_color
from pycbc import DYN_RANGE_FAC
import matplotlib
matplotlib.use('Agg')
import pylab as pl
import numpy as np

def generate_snr_plot(snrdict, output_filename, sngl_table):
    """
    Generate an SNR timeseries plot as used for upload to GraceDB.

    Parameters
    ----------

    snrdict: dictionary
        A dictionary keyed on ifo containing the SNRs as
        TimeSeries objects

    output_filename: string
        The filename for the plot to be saved to

    Returns
    -------
        None
    """
    events = {sngl.ifo: sngl for sngl in sngl_table}
    pl.figure()
    ref_time = np.floor(np.mean([evnt.end_time + 1e-9 * evnt.end_time_ns
                                 for evnt in events.values()]))
    for ifo in sorted(snrdict.keys()):
        curr_snrs = snrdict[ifo]

        pl.plot(curr_snrs.sample_times - ref_time, abs(curr_snrs),
                c=ifo_color(ifo), label=ifo)
        pl.plot(events[ifo].end_time + 1e-9 * events[ifo].end_time_ns - ref_time,
                events[ifo].snr, marker='x', c=ifo_color(ifo))


    pl.legend()
    pl.xlabel(f'GPS time from {ref_time} (s)')
    pl.ylabel('SNR')
    pl.savefig(output_filename)
    pl.close()

__all__ = ["generate_snr_plot"]
