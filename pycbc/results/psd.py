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


def generate_asd_plot(psddict, output_filename, f_min=10.):
    """
    Generate an ASD plot as used for upload to GraceDB.

    Parameters
    ----------
    psddict: dictionary
        A dictionary keyed on ifo containing the PSDs as
        FrequencySeries objects

    output_filename: string
        The filename for the plot to be saved to

    f_min: float
        Minimum frequency at which anything should be plotted

    Returns
    -------
        None
    """
    from matplotlib import pyplot as plt
    asd_fig, asd_ax = plt.subplots(1)
    asd_min = [1E-24]  # Default minimum to plot

    for ifo in sorted(psddict.keys()):
        curr_psd = psddict[ifo]
        freqs = curr_psd.sample_frequencies
        physical = (freqs >= f_min)  # Ignore lower frequencies
        asd_to_plot = curr_psd[physical] ** 0.5 / DYN_RANGE_FAC
        asd_min.append(min(asd_to_plot))
        asd_ax.loglog(freqs[physical],
                      asd_to_plot,
                      c=ifo_color(ifo),
                      label=ifo)

    asd_ax.grid(True)
    asd_ax.legend()
    asd_ax.set_xlim([f_min, 1300])
    asd_ax.set_ylim([min(asd_min), 1E-20])
    asd_ax.set_xlabel('Frequency (Hz)')
    asd_ax.set_ylabel('ASD')
    asd_fig.savefig(output_filename)


__all__ = ["generate_asd_plot"]
