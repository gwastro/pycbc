#!/usr/bin/env python

# Copyright (C) 2015 Andrew R. Williamson
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

from ligo import segments

def make_grb_segments_plot(wkflow, science_segs, trigger_time, trigger_name,
                           out_dir, coherent_seg=None, fail_criterion=None):

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
            segments.segment(trigger_time - pltpad[0],
                             trigger_time + pltpad[1])]).extent()

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
                                    facecolor=ifo_colors[ifo], edgecolor='none'))
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
