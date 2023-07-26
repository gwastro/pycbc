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
from ligo import segments
from pycbc.results import save_fig_with_metadata


# =============================================================================
# Used locally: plot contours in a scatter plot with SNR as horizontal axis
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
# Functions used in executables
#

# =============================================================================
# Plot trigger time and offsource extent over segments
# Courtesy of Alex Dietz
# =============================================================================
def make_grb_segments_plot(wkflow, science_segs, trigger_time, trigger_name,
                           out_dir, coherent_seg=None, fail_criterion=None):
    """Plot trigger time and offsource extent over segments"""

    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from matplotlib.lines import Line2D
    from pycbc.results.color import ifo_color

    ifos = wkflow.ifos
    if len(science_segs.keys()) == 0:
        extent = segments.segment(int(wkflow.cp.get("workflow", "start-time")),
                                  int(wkflow.cp.get("workflow", "end-time")))
    else:
        pltpad = [science_segs.extent_all()[1] - trigger_time,
                  trigger_time - science_segs.extent_all()[0]]
        extent = segments.segmentlist([science_segs.extent_all(),
                                       segments.segment(trigger_time
                                                        - pltpad[0],
                                                        trigger_time
                                                        + pltpad[1])]).extent()

    ifo_colors = {}
    for ifo in ifos:
        ifo_colors[ifo] = ifo_color(ifo)
        if ifo not in science_segs.keys():
            science_segs[ifo] = segments.segmentlist([])

    # Make plot
    fig, subs = plt.subplots(len(ifos), sharey=True)
    if len(ifos) == 1:
        subs = [subs]
    plt.xticks(rotation=20, ha='right')
    for sub, ifo in zip(subs, ifos):
        for seg in science_segs[ifo]:
            sub.add_patch(Rectangle((seg[0], 0.1), abs(seg), 0.8,
                                    facecolor=ifo_colors[ifo],
                                    edgecolor='none'))
        if coherent_seg:
            if len(science_segs[ifo]) > 0 and \
                    coherent_seg in science_segs[ifo]:
                sub.plot([trigger_time, trigger_time], [0, 1], '-',
                         c='orange')
                sub.add_patch(Rectangle((coherent_seg[0], 0),
                                        abs(coherent_seg), 1, alpha=0.5,
                                        facecolor='orange', edgecolor='none'))
            else:
                sub.plot([trigger_time, trigger_time], [0, 1], ':',
                         c='orange')
                sub.plot([coherent_seg[0], coherent_seg[0]], [0, 1], '--',
                         c='orange', alpha=0.5)
                sub.plot([coherent_seg[1], coherent_seg[1]], [0, 1], '--',
                         c='orange', alpha=0.5)
        else:
            sub.plot([trigger_time, trigger_time], [0, 1], ':k')
        if fail_criterion:
            if len(science_segs[ifo]) > 0:
                style_str = '--'
            else:
                style_str = '-'
            sub.plot([fail_criterion[0], fail_criterion[0]], [0, 1], style_str,
                     c='black', alpha=0.5)
            sub.plot([fail_criterion[1], fail_criterion[1]], [0, 1], style_str,
                     c='black', alpha=0.5)

        sub.set_frame_on(False)
        sub.set_yticks([])
        sub.set_ylabel(ifo, rotation=45)
        sub.set_ylim([0, 1])
        sub.set_xlim([float(extent[0]), float(extent[1])])
        sub.get_xaxis().get_major_formatter().set_useOffset(False)
        sub.get_xaxis().get_major_formatter().set_scientific(False)
        sub.get_xaxis().tick_bottom()
        if sub is subs[-1]:
            sub.tick_params(labelsize=10, pad=1)
        else:
            sub.get_xaxis().set_ticks([])
            sub.get_xaxis().set_ticklabels([])

    xmin, xmax = fig.axes[-1].get_xaxis().get_view_interval()
    ymin, _ = fig.axes[-1].get_yaxis().get_view_interval()
    fig.axes[-1].add_artist(Line2D((xmin, xmax), (ymin, ymin), color='black',
                                   linewidth=2))
    fig.axes[-1].set_xlabel('GPS Time')

    fig.axes[0].set_title('Science Segments for GRB%s' % trigger_name)
    plt.tight_layout()
    fig.subplots_adjust(hspace=0)

    plot_name = 'GRB%s_segments.png' % trigger_name
    plot_url = 'file://localhost%s/%s' % (out_dir, plot_name)
    fig.savefig('%s/%s' % (out_dir, plot_name))

    return [ifos, plot_name, extent, plot_url]


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
