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

# =============================================================================
# Preamble
# =============================================================================

from __future__ import division

import re
import os
from argparse import ArgumentParser
from ligo import segments
from glue import markup

def initialize_page(title, style, script, header=None):
    """
    A function that returns a markup.py page object with the required html
    header.
    """

    page = markup.page(mode="strict_html")
    page._escape = False
    page.init(title=title, css=style, script=script, header=header)

    return page


def write_banner(title, text='&nbsp;'):
    """
    Write html <title> tag into markup.page object
    """

    page = markup.page(mode="strict_html")
    page._escape = False

    page.div(id="header")
    page.h1()
    page.add(title)
    page.h1.close()
    page.h3()
    page.add(text)
    page.h3.close()

    page.hr(class_="short")
    page.hr(class_="long")

    page.div.close()

    page.div(id="container")

    return page


def write_table(page, headers, data, cl=''):

    """
    Write table in html
    """

    page.table(class_=cl)

    # list
    if cl=='list':
        for i in range(len(headers)):

            page.tr()
            page.th()
            page.add('%s' % headers[i])
            page.th.close()
            page.td()
            page.add('%s' % data[i])
            page.td.close()
            page.tr.close()

    else:
        page.tr()
        for n in headers:
            page.th()
            page.add('%s' % n)
            page.th.close()
        page.tr.close()

        if data and not re.search('list',str(type(data[0]))):
            data = [data]

        for row in data:
            page.tr()
            for item in row:
                page.td()
                page.add('%s' % item)
                page.td.close()
            page.tr.close()

    page.table.close()

    return page


def write_summary(page, args, ifos, skyError=None, ipn=False, ipnError=False):

    """
        Write summary of information to markup.page object page
    """
    from pylal import antenna
    from lal.gpstime import gps_to_utc, LIGOTimeGPS

    gps = args.start_time
    grbdate = gps_to_utc(LIGOTimeGPS(gps))\
                                .strftime("%B %d %Y, %H:%M:%S %ZUTC")
    page.h3()
    page.add('Basic information')
    page.h3.close()

    if ipn:
        ra = []
        dec = []
        td1 = []
        td2 = []
        td3 = []
        timedelay = {}
        search_file = '../../../S5IPN_GRB%s_search_180deg.txt' % args.grb_name
        for line in open(search_file):
            ra.append(line.split()[0])
            dec.append(line.split()[1])
        th1 = [ 'GPS', 'Date', 'Error Box (sq.deg.)', 'IFOs' ]
        td1 = [ gps, grbdate, ipnError, ifos ]
        th2 = [ 'RA', 'DEC' ]
        th3 = ['Timedelays (ms)', '', '' ]
        for ra_i,dec_i in zip(ra,dec):
            td_i = [ ra_i, dec_i ]
            td2.append(td_i)
        ifo_list = [ ifos[i*2:(i*2)+2] for i in range(int(len(ifos)/2)) ]
        for j in td2:
            for p in range(0, len(ifo_list)):
                for q in range(0, len(ifo_list)):
                    pairs = [ifo_list[p], ifo_list[q]]
                    ifo_pairs = "".join(pairs)
                    timedelay[ifo_pairs] = antenna.timeDelay(int(gps),
                            float(j[0]), float(j[1]), 'degree', ifo_list[p],
                            ifo_list[q])
                    timedelay[ifo_pairs]="%.4f" % timedelay[ifo_pairs]
            if ifos == 'H1H2L1':
                td3.append(['H1L1: %f' % float(timedelay['H1L1'])])
            if ifos == 'H1H2L1V1':
                td3.append(['H1L1: %f' % float(timedelay['H1L1']),
                            'H1V1: %f' % float(timedelay['H1V1']),
                            'L1V1: %f' % float(timedelay['L1V1'])])
            if ifos == 'L1V1':
                td3.append(['L1V1: %f' % float(timedelay['L1V1'])])
        page = write_table(page, th1, td1)
        page = write_table(page, th2, td2)
        page = write_table(page, th3, td3)

    else:
        ra = args.ra
        dec = args.dec
        if skyError:
            th = [ 'GPS', 'Date', 'RA', 'DEC', 'Sky Error', 'IFOs' ]
            td = [ gps, grbdate, ra, dec, skyError, ifos ]
        else:
            th = [ 'GPS', 'Date', 'RA', 'DEC', 'IFOs' ]
            td = [ gps, grbdate, ra, dec, ifos ]

        page = write_table(page, th, td)

    return page


def write_antenna(page, args, seg_plot=None, grid=False, ipn=False):

    """
    Write antenna factors to merkup.page object page and generate John's
    detector response plot.
    """
    from pylal import antenna

    page.h3()
    page.add('Antenna factors and sky locations')
    page.h3.close()

    th = []
    td = []
    th2 = []
    td2 = []

    ifos = [args.ifo_tag[i:i+2] for i in range(0, len(args.ifo_tag), 2)]

    if ipn:
        antenna_ifo = {}
        ra = []
        dec = []
        # FIXME: Remove hardcoding here and show this in all cases
        search_file = open('../../../S5IPN_GRB%s_search_180deg.txt'
                           % args.grb_name)
        for line in search_file:
            ra.append(line.split()[0])
            dec.append(line.split()[1])
        for ifo in ifos:
            antenna_ifo[ifo] = []
            for k, l in zip(ra, dec):
                _, _, _, f_q = antenna.response(args.start_time, float(k),
                                                float(l), 0.0, 0.0, 'degree',
                                                ifo)
                antenna_ifo[ifo].append(round(f_q,3))
        dectKeys = antenna_ifo.keys()

        for elements in range(len(antenna_ifo.values()[0])):
            newDict={}
            for detectors in range(len(antenna_ifo.keys())):
                newDict[dectKeys[detectors]] = antenna_ifo[\
                                               dectKeys[detectors]][elements]
            for key in newDict.keys():
                th.append(key)
            td.append(newDict.values())
        page = write_table(page, list(set(th)), td)
    for ifo in ifos:
        _, _, _, f_q = antenna.response(args.start_time, args.ra, args.dec,
                                        0.0, 0.0, 'degree',ifo)
        th.append(ifo)
        td.append(round(f_q, 3))

    #FIXME: Work out a way to make these external calls safely
    #cmmnd = 'projectedDetectorTensor --gps-sec %d --ra-deg %f --dec-deg %f' \
    #         % (args.start_time,args.ra, args.dec)
    #for ifo in ifos:
    #    if ifo == 'H1':
    #        cmmnd += ' --display-lho'
    #    elif ifo == 'L1':
    #        cmmnd += ' --display-llo'
    #    elif ifo == 'V1':
    #        cmmnd += ' --display-virgo'
    #status = make_external_call(cmmnd)

    page = write_table(page, th, td)

#    plot = markup.page()
#    p = "projtens.png"
#    plot.a(href=p, title="Detector response and polarization")
#    plot.img(src=p)
#    plot.a.close()
#    th2 = ['Response Diagram']
#    td2 = [plot() ]

# FIXME: Add these in!!
#    plot = markup.page()
#    p = "ALL_TIMES/plots_clustered/GRB%s_search.png"\
#        % args.grb_name
#    plot.a(href=p, title="Error Box Search")
#    plot.img(src=p)
#    plot.a.close()
#    th2.append('Error Box Search')
#    td2.append(plot())

#    plot = markup.page()
#    p = "ALL_TIMES/plots_clustered/GRB%s_simulations.png"\
#        % args.grb_name
#    plot.a(href=p, title="Error Box Simulations")
#    plot.img(src=p)
#    plot.a.close()
#    th2.append('Error Box Simulations')
#    td2.append(plot())

    if seg_plot is not None:
        plot = markup.page()
        p = os.path.basename(seg_plot)
        plot.a(href=p, title="Science Segments")
        plot.img(src=p)
        plot.a.close()
        th2.append('Science Segments')
        td2.append(plot())

    plot = markup.page()
    p = "ALL_TIMES/plots_clustered/GRB%s_sky_grid.png"\
            % args.grb_name
    plot.a(href=p, title="Sky Grid")
    plot.img(src=p)
    plot.a.close()
    th2.append('Sky Grid')
    td2.append(plot())

#    plot = markup.page()
#    p = "GRB%s_inspiral_horizon_distance.png"\
#            % args.grb_name
#    plot.a(href=p, title="Inspiral Horizon Distance")
#    plot.img(src=p)
#    plot.a.close()
#    th2.append('Inspiral Horizon Distance')
#    td2.append(plot())

    page = write_table(page, th2, td2)

    return page


def write_offsource(page, args, grbtag, onsource=False):

    """
        Write offsource SNR versus time plots to markup.page object page
    """

    th = ['Re-weighted SNR', 'Coherent SNR']

    if args.time_slides:
        if onsource:
            out_dir = 'ZEROLAG_ALL'
        else:
            out_dir = 'ZEROLAG_OFF'
    else:
        if onsource:
            out_dir = 'ALL_TIMES'
        else:
            out_dir = 'OFFSOURCE'

    plot = markup.page()
    p = "%s/plots_clustered/GRB%s_bestnr_vs_time_noinj.png" % (out_dir, grbtag)
    plot.a(href=p, title="Detection statistic versus time")
    plot.img(src=p)
    plot.a.close()
    td = [ plot() ]

    plot = markup.page()
    p = "%s/plots_clustered/GRB%s_triggers_vs_time_noinj.png" % (out_dir, grbtag)
    plot.a(href=p, title="Coherent SNR versus time")
    plot.img(src=p)
    plot.a.close()
    td.append(plot())

    ifos = [args.ifo_tag[i:i+2] for i in range(0, len(args.ifo_tag), 2)]
    for ifo in ifos:
        th.append('%s SNR' % ifo)
        plot = markup.page()
        p = "%s/plots_clustered/GRB%s_%s_triggers_vs_time_noinj.png"\
            % (out_dir, grbtag, ifo)
        plot.a(href=p, title="%s SNR versus time" % ifo)
        plot.img(src=p)
        plot.a.close()
        td.append(plot())

    page = write_table(page, th, td)

    return page


def write_chisq(page, injList, grbtag):

    """
        Write injection chisq plots to markup.page object page
    """

    if injList:
        th = ['']+injList + ['OFFSOURCE']
    else:
        th= ['','OFFSOURCE']
        injList = ['OFFSOURCE']
    td = []

    plots = ['bank_veto','auto_veto','chi_square', 'mchirp']

    for test in plots:
        pTag = test.replace('_',' ').title()
        d = [pTag]
        for inj in injList + ['OFFSOURCE']:
            plot = markup.page()
            p = "%s/plots_clustered/GRB%s_%s_vs_snr_zoom.png" % (inj, grbtag,
                                                                 test)
            plot.a(href=p, title="%s %s versus SNR" % (inj, pTag))
            plot.img(src=p)
            plot.a.close()
            d.append(plot())

        td.append(d)

    page = write_table(page, th, td)

    return page


def write_inj_snrs(page, ifos, injList, grbtag):

    """
        Write injection chisq plots to markup.page object page
    """

    if injList:
        th = ['']+injList + ['OFFSOURCE']
    else:
        th= ['','OFFSOURCE']
        injList = ['OFFSOURCE']
    td = []

    ifos = [ifos[i:i+2] for i in range(0, len(ifos), 2)]
    plots = ['null_stat2']+['%s_snr' % ifo for ifo in ifos]

    for row in plots:
        pTag = row.replace('_',' ').title()
        d = [pTag]
        for inj in injList + ['OFFSOURCE']:
            plot = markup.page()
            p = "%s/plots_clustered/GRB%s_%s_vs_snr_zoom.png" % (inj, grbtag,
                                                                 row)
            plot.a(href=p, title="%s %s versus SNR" % (inj, pTag))
            plot.img(src=p)
            plot.a.close()
            d.append(plot())
        td.append(d)

    page = write_table(page, th, td)

    return page


def write_found_missed(page, args, injList):

    """
        Write injection found/missed plots to markup.page object page
    """

    th = ['']+injList
    td = []

    #FIXME: Work out a way to make externals calls safely
    #d = ['Number of injections']
    #for inj in injList:
    #    cmmnd = 'lwtprint ../*' + inj + '*MISSED*xml -t sim_inspiral | wc -l'
    #    output,status = make_external_call(cmmnd, shell=True)
    #    numInjs = int(output)
    #    cmmnd = 'lwtprint ../*' + inj + '*FOUND*xml -t sim_inspiral | wc -l'
    #    output,status = make_external_call(cmmnd, shell=True)
    #    numInjs += int(output)
    #    d.append(str(numInjs))
    #td.append(d)

    plots = []
    text  = {}
    ifos = [args.ifo_tag[i:i+2] for i in range(0, len(args.ifo_tag), 2)]
    plots.extend(['dist', 'dist_time'])
    text['dist'] = 'Dist vs Mchirp'
    text['dist_time'] = 'Dist vs Time'
    for ifo in ifos:
        plots.extend(['effdist_%s' % ifo[0].lower(),\
                      'effdist_time_%s' % ifo[0].lower()])
        text['effdist_%s' % ifo[0].lower()] = 'Eff. dist. %s vs Mchirp' % ifo
        text['effdist_time_%s' % ifo[0].lower()] = 'Eff. dist %s vs Time' % ifo

    for row in plots:
        pTag = text[row]
        d = [pTag]
        for inj in injList:
            plot = markup.page()
            p = "%s/efficiency_OFFTRIAL_1/found_missed_injections_%s.png"\
                % (inj, row)
            plot.a(href=p, title=pTag)
            plot.img(src=p)
            plot.a.close()
            d.append(plot())
        td.append(d)

    td.append(['Close injections without FAP = 0']+\
              ['<a href="%s/efficiency_OFFTRIAL_1/quiet_found_triggers.html"> '
               'here</a>' % inj for inj in injList])

    page = write_table(page, th, td)

    return page


def write_recovery(page, injList):

    """
        Write injection recovery plots to markup.page object page
    """

    th = ['']+injList
    td = []

    plots = ['sky_error_time','sky_error_mchirp','sky_error_distance']
    text = { 'sky_error_time':'Sky error vs time',\
                      'sky_error_mchirp':'Sky error vs mchirp',\
                      'sky_error_distance':'Sky error vs distance' }

    for row in plots:
        pTag = text[row]
        d = [pTag]
        for inj in injList:
            plot = markup.page()
            plot = markup.page()
            p = "%s/efficiency_OFFTRIAL_1/found_%s.png" % (inj, row)
            plot.a(href=p, title=pTag)
            plot.img(src=p)
            plot.a.close()
            d.append(plot())
        td.append(d)

    page = write_table(page, th, td)

    return page


def write_loudest_events(page, bins, onsource=False):

    """
        Write injection chisq plots to markup.page object page
    """

    th = ['']+['Mchirp %s - %s' % tuple(bin) for bin in bins]
    td = []

    plots = ['BestNR','SNR']

    if onsource:
        trial = 'ONSOURCE'
    else:
        trial = 'OFFTRIAL_1'

    for pTag in plots:
        row = pTag.lower()
        d = [pTag]
        for bin in bins:
            b = '%s_%s' % tuple(bin)
            plot = markup.page()
            p = "%s/efficiency/%s_vs_fap_%s.png" % (trial, row, b)
            plot.a(href=p, title="FAP versus %s" % pTag)
            plot.img(src=p)
            plot.a.close()
            d.append(plot())
        td.append(d)

    row = 'snruncut'
    d = ['SNR after cuts <br> have been applied']
    for bin in bins:
        b = '%s_%s' % tuple(bin)
        plot = markup.page()
        p = "%s/efficiency/%s_vs_fap_%s.png" % (trial, row, b)
        plot.a(href=p, title="FAP versus %s" % pTag)
        plot.img(src=p)
        plot.a.close()
        d.append(plot())
    td.append(d)

    page = write_table(page, th, td)

    page.add('For more details on the loudest offsource events see')
    page.a(href='%s/efficiency/loudest_offsource_trigs.html' % (trial))
    page.add('here.')
    page.a.close()


    return page


def write_exclusion_distances(page , trial, injList, massbins, reduced=False,
                              onsource=False):
    file = open('%s/efficiency/loud_numbers.txt' % (trial), 'r')
    FAPS = []
    for line in file:
        line = line.replace('\n','')
        if line == "-2":
            FAPS.append('No event')
        else:
            FAPS.append(line)

    file.close()

    th = ['']+['Mchirp %s - %s' % tuple(bin) for bin in massbins]
    td = ['FAP']+FAPS
    page = write_table(page, th, td)
    page.add('For more details on the loudest onsource events see')
    page.a(href='%s/efficiency/loudest_events.html' % (trial))
    page.add('here.')
    page.a.close()

    if reduced or not injList:
        return page

    page.h3()
    page.add('Detection efficiency plots - injections louder than loudest '
             'background trigger')
    page.h3.close()

    th = injList
    td = []
    d = []
    for inj in injList:
        plot = markup.page()
        p = "%s/efficiency_%s/BestNR_max_efficiency.png" % (inj, trial)
        plot.a(href=p, title="Detection efficiency")
        plot.img(src=p)
        plot.a.close()
        d.append(plot())
    td.append(d)

    page = write_table(page, th, td)

    page.h3()
    page.add('Exclusion distance plots - injections louder than loudest '
             'foreground trigger')
    page.h3.close()

    th = injList
    td = []
    d = []
    for inj in injList:
        plot = markup.page()
        p = "%s/efficiency_%s/BestNR_on_efficiency.png" % (inj, trial)
        plot.a(href=p, title="Exclusion efficiency")
        plot.img(src=p)
        plot.a.close()
        d.append(plot())
    td.append(d)

    page = write_table(page, th, td)

    for percentile in [90, 50]:
        page.h3()
        page.add('%d%% confidence exclusion distances (Mpc)' % percentile)
        th = injList
        td = []
        d = []
        for inj in injList:
            file = open('%s/efficiency_%s/exclusion_distance_%d.txt'
                        % (inj, trial, percentile), 'r')
            for line in file:
                line = line.replace('\n','')
                excl_dist = float(line)
            d.append(excl_dist)
            file.close()
        td.append(d)
        page = write_table(page, th, td)

    page.h3.close()

    return page


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
