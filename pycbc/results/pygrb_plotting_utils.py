# Copyright (C) 2019 Francesco Pannarale, Gino Contestabile, Cameron Mills
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


# =============================================================================
# Preamble
# =============================================================================

"""
Module to generate PyGRB figures: scatter plots and timeseries.
"""

import copy
import numpy
from pycbc.results import save_fig_with_metadata


#
# Used locally
#

# =============================================================================
# Function to calculate chi-square weight for the reweighted SNR
# =============================================================================
def new_snr_chisq(snr, new_snr, chisq_dof, chisq_index=4.0, chisq_nhigh=3.0):
    """Returns the chi-square value needed to weight SNR into new SNR"""

    chisqnorm = (snr/new_snr)**chisq_index
    if chisqnorm <= 1:
        return 1E-20

    return chisq_dof * (2*chisqnorm - 1)**(chisq_nhigh/chisq_index)


# =============================================================================
# Plot contours in a scatter plot where SNR is on the horizontal axis
# =============================================================================
def contour_plotter(axis, snr_vals, contours, colors, vert_spike=False):
    """Plot contours in a scatter plot where SNR is on the horizontal axis"""

    for i, _ in enumerate(contours):
        plot_vals_x = []
        plot_vals_y = []
        if vert_spike:
            for j, _ in enumerate(snr_vals):
                # Workaround to ensure vertical spike is shown on veto plots
                if contours[i][j] > 1E-15 and not plot_vals_x:
                    plot_vals_x.append(snr_vals[j])
                    plot_vals_y.append(0.1)
                if contours[i][j] > 1E-15 and plot_vals_x:
                    plot_vals_x.append(snr_vals[j])
                    plot_vals_y.append(contours[i][j])
        else:
            plot_vals_x = snr_vals
            plot_vals_y = contours[i]
        axis.plot(plot_vals_x, plot_vals_y, colors[i])


#
# Used (also) in executables
#

# =============================================================================
# Given the trigger and injection values of a quantity, determine the maximum
# =============================================================================
def axis_max_value(trig_values, inj_values, inj_file):
    """Deterime the maximum of a quantity in the trigger and injection data"""

    axis_max = trig_values.max()
    if inj_file and inj_values.size and inj_values.max() > axis_max:
        axis_max = inj_values.max()

    return axis_max


# =============================================================================
# Master plotting function: fits all plotting needs in for PyGRB results
# =============================================================================
def pygrb_plotter(trigs, injs, xlabel, ylabel, opts,
                  snr_vals=None, conts=None, shade_cont_value=None,
                  colors=None, vert_spike=False, cmd=None):
    """Master function to plot PyGRB results"""
    from matplotlib import pyplot as plt

    # Set up plot
    fig = plt.figure()
    cax = fig.gca()
    # Plot trigger-related and (if present) injection-related quantities
    cax_plotter = cax.loglog if opts.use_logs else cax.plot
    cax_plotter(trigs[0], trigs[1], 'bx')
    if not (injs[0] is None and injs[1] is None):
        cax_plotter(injs[0], injs[1], 'r+')
    cax.grid()
    # Plot contours
    if conts is not None:
        contour_plotter(cax, snr_vals, conts, colors, vert_spike=vert_spike)
    # Add shading above a specific contour (typically used for vetoed area)
    if shade_cont_value is not None:
        limy = cax.get_ylim()[1]
        polyx = copy.deepcopy(snr_vals)
        polyy = copy.deepcopy(conts[shade_cont_value])
        polyx = numpy.append(polyx, [max(snr_vals), min(snr_vals)])
        polyy = numpy.append(polyy, [limy, limy])
        cax.fill(polyx, polyy, color='#dddddd')
    # Axes: labels and limits
    cax.set_xlabel(xlabel)
    cax.set_ylabel(ylabel)
    if opts.x_lims:
        x_lims = map(float, opts.x_lims.split(','))
        cax.set_xlim(x_lims)
    if opts.y_lims:
        y_lims = map(float, opts.y_lims.split(','))
        cax.set_ylim(y_lims)
    # Wrap up
    plt.tight_layout()
    save_fig_with_metadata(fig, opts.output_file, cmd=cmd,
                           title=opts.plot_title,
                           caption=opts.plot_caption)
    plt.close()
