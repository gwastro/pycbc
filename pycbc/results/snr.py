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
Module to generate SNR figures
"""
import pylab as pl
from pycbc.results import ifo_color


def generate_snr_plot(snrdict, output_filename, triggers, ref_time):
    """
    Generate an SNR timeseries plot as used for upload to GraceDB.

    Parameters
    ----------

    snrdict: dictionary
        A dictionary keyed on ifo containing the SNR
        TimeSeries objects
    output_filename: string
        The filename for the plot to be saved to
    triggers : dictionary of tuples
        A dictionary keyed on IFO, containing (trigger time, trigger snr)
    ref_time : number, GPS seconds
        Reference time which will be used as the zero point of the plot
        This should be an integer value, but doesn't need to be an integer

    Returns
    -------
        None
    """
    pl.figure()
    ref_time = int(ref_time)
    for ifo in sorted(snrdict):
        curr_snrs = snrdict[ifo]

        pl.plot(curr_snrs.sample_times - ref_time, abs(curr_snrs),
                c=ifo_color(ifo), label=ifo)
        if ifo in triggers:
            pl.plot(triggers[ifo][0] - ref_time,
                    triggers[ifo][1], marker='x', c=ifo_color(ifo))

    pl.legend()
    pl.xlabel(f'GPS time from {ref_time:d} (s)')
    pl.ylabel('SNR')
    pl.savefig(output_filename)
    pl.close()


__all__ = ["generate_snr_plot"]
