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

import sys
import os
import logging
import copy
import numpy
from pycbc.results import save_fig_with_metadata
# Only if a backend is not already set ... This should really *not* be done
# here, but in the executables you should set matplotlib.use()
# This matches the check that matplotlib does internally, but this *may* be
# version dependent. If this is a problem then remove this and control from
# the executables directly.
import matplotlib
if 'matplotlib.backends' not in sys.modules:  # nopep8
    matplotlib.use('agg')
from matplotlib import rc
from matplotlib import pyplot as plt

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
# Calculate all chi-square contours for diagnostic plots
# =============================================================================
# TODO: tailor to plot_null_stat and plot_chisq_veto
def calculate_contours(trig_data, opts, new_snrs=None):
    """Generate the plot contours for veto plots"""

    if new_snrs is None:
        new_snrs = [5.5, 6, 6.5, 7, 8, 9, 10, 11]
    chisq_index = opts.chisq_index
    chisq_nhigh = opts.chisq_nhigh
    new_snr_thresh = opts.newsnr_threshold
    null_thresh = []
    for val in map(float, opts.null_snr_threshold.split(',')):
        null_thresh.append(val)
    null_thresh = null_thresh[-1]
    null_grad_snr = opts.null_grad_thresh
    null_grad_val = opts.null_grad_val
    chisq_dof = trig_data['chisq_dof'][0]
    bank_chisq_dof = trig_data['bank_chisq_dof'][0]
    cont_chisq_dof = trig_data['cont_chisq_dof'][0]

    # Add the new SNR threshold contour to the list if necessary
    # and keep track of where it is
    cont_value = None
    try:
        cont_value = new_snrs.index(new_snr_thresh)
    except ValueError:
        new_snrs.append(new_snr_thresh)
        cont_value = -1

    # Initialise chisq contour values and colours
    colors = ["k-" if snr == new_snr_thresh else
              "y-" if snr == int(snr) else
              "y--" for snr in new_snrs]

    # Get SNR values for contours
    snr_low_vals = numpy.arange(4, 30, 0.1)
    snr_high_vals = numpy.arange(30, 500, 1)
    snr_vals = numpy.asarray(list(snr_low_vals) + list(snr_high_vals))

    # Initialise contours
    bank_conts = numpy.zeros([len(new_snrs), len(snr_vals)],
                             dtype=numpy.float64)
    auto_conts = numpy.zeros([len(new_snrs), len(snr_vals)],
                             dtype=numpy.float64)
    chi_conts = numpy.zeros([len(new_snrs), len(snr_vals)],
                            dtype=numpy.float64)
    null_cont = []

    # Loop over each and calculate chisq variable needed for SNR contour
    for j, snr in enumerate(snr_vals):
        for i, new_snr in enumerate(new_snrs):
            bank_conts[i][j] = new_snr_chisq(snr, new_snr, bank_chisq_dof,
                                             chisq_index, chisq_nhigh)
            auto_conts[i][j] = new_snr_chisq(snr, new_snr, cont_chisq_dof,
                                             chisq_index, chisq_nhigh)
            chi_conts[i][j] = new_snr_chisq(snr, new_snr, chisq_dof,
                                            chisq_index, chisq_nhigh)

        if snr > null_grad_snr:
            null_cont.append(null_thresh + (snr-null_grad_snr)*null_grad_val)
        else:
            null_cont.append(null_thresh)
    null_cont = numpy.asarray(null_cont)

    return bank_conts, auto_conts, chi_conts, null_cont, snr_vals, \
        cont_value, colors


# =============================================================================
# Contains plotting setups shared by PyGRB plots
# =============================================================================
def pygrb_shared_plot_setups():
    """Master function to plot PyGRB results"""

    # Get rcParams
    rc('font', size=14)
    # Set color for out-of-range values
    #plt.cm.spring.set_over('g')


# =============================================================================
# Master plotting function: fits all plotting needs in for PyGRB results
# =============================================================================
def pygrb_plotter(trig_x, trig_y, inj_x, inj_y, inj_file, xlabel, ylabel,
                  fig_path, snr_vals=None, conts=None,
                  shade_cont_value=None, colors=None, vert_spike=False,
                  xlims=None, ylims=None, use_logs=True,
                  cmd=None, plot_title=None, plot_caption=None):
    """Master function to plot PyGRB results"""

    fig_name = os.path.split(os.path.abspath(fig_path))[1]
    logging.info(" * %s (%s vs %s)...", fig_name, xlabel, ylabel)
    # Set up plot
    fig = plt.figure()
    cax = fig.gca()
    # Plot trigger-related quantities
    if use_logs:
        cax.loglog(trig_x, trig_y, 'bx')
    else:
        cax.plot(trig_x, trig_y, 'bx')
    cax.grid()
    # Plot injection-related quantities
    if inj_file:
        if use_logs:
            cax.loglog(inj_x, inj_y, 'r+')
        else:
            cax.plot(inj_x, inj_y, 'r+')
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
    if xlims:
        cax.set_xlim(xlims)
    if ylims:
        cax.set_ylim(ylims)
    # Wrap up
    plt.tight_layout()
    save_fig_with_metadata(fig, fig_path, cmd=cmd, title=plot_title,
                           caption=plot_caption)
    # fig_kwds=fig_kwds,
    plt.close()
