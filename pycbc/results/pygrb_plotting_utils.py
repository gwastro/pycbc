# Copyright (C) 2019 Francesco Pannarale, Gino Contestabile
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
import argparse
import copy
import numpy
from pycbc.results import save_fig_with_metadata
# TODO: get rid of this (or move get_bestnr from pylal)
try:
    from pylal.coh_PTF_pyutils import get_bestnr, get_det_response
except ImportError:
    import logging
    MSG = "Running the PyGRB post-processing requires a pylal installation"
    logging.error(MSG)

import matplotlib
from matplotlib import rc
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Get rcParams
rc('font', size=14)

plt.cm.spring.set_over('g')


# =============================================================================
# Parse command line
# =============================================================================

def pygrb_plot_opts_parser(usage='', description=None, version=None):
    """Parses options for PyGRB plotting scripts"""
    parser = argparse.ArgumentParser(usage=usage, description=description)

    parser.add_argument("--version", action="version", version=version)

    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="verbose output")

    parser.add_argument("-t", "--trig-file", action="store",
                        default=None, required=True,
                        help="The location of the trigger file")

    parser.add_argument("-I", "--inj-file", action="store", default=None,
                        help="The location of the injection file")

    parser.add_argument("-a", "--segment-dir", action="store",
                        required=True, help="directory holding buffer, on " +
                        "and off source segment files.")

    parser.add_argument("-o", "--output-file", default=None, required=True,
                        help="Output file.")

    parser.add_argument("-O", "--zoomed-output-file", default=None,
                        required=False, help="Output file for a zoomed in " +
                        "version of the plot.")

    parser.add_argument("-Q", "--chisq-index", action="store", type=float,
                        default=4.0, help="chisq_index for newSNR calculation")

    parser.add_argument("-N", "--chisq-nhigh", action="store", type=float,
                        default=3.0, help="nhigh for newSNR calculation")

    parser.add_argument("-B", "--sngl-snr-threshold", action="store",
                        type=float, default=4.0, help="Single detector SNR " +
                        "threshold, the two most sensitive detectors " +
                        "should have SNR above this")

    parser.add_argument("-d", "--snr-threshold", action="store", type=float,
                        default=6.0, help="SNR threshold for recording " +
                        "triggers")

    parser.add_argument("-c", "--newsnr-threshold", action="store", type=float,
                        default=None, help="NewSNR threshold for " +
                        "calculating the chisq of triggers (based on value " +
                        "of auto and bank chisq  values. By default will " +
                        "take the same value as snr-threshold")

    parser.add_argument("-A", "--null-snr-threshold", action="store",
                        default="4.25,6",
                        help="comma separated lower,higher null SNR " +
                        "threshold for null SNR cut")

    parser.add_argument("-C", "--null-grad-thresh", action="store", type=float,
                        default=20., help="Threshold above which to " +
                        "increase the values of the null SNR cut")

    parser.add_argument("-D", "--null-grad-val", action="store", type=float,
                        default=0.2, help="Rate the null SNR cut will " +
                        "increase above the threshold")

    parser.add_argument("-l", "--veto-directory", action="store", default=None,
                        help="The location of the CATX veto files")

    parser.add_argument("-b", "--veto-category", action="store", type=int,
                        default=None, help="Apply vetoes up to this level " +
                        "inclusive")

    parser.add_argument("--ifo", default=None, help="IFO to be used")

    return parser.parse_args()


# =============================================================================
# Format single detector chi-square data as numpy array and floor at 0.005
# =============================================================================

def format_single_chisqs(trig_ifo_cs, ifos):
    """Format single IFO chi-square data as numpy array and floor at 0.005"""
    for ifo in ifos:
        trig_ifo_cs[ifo] = numpy.asarray(trig_ifo_cs[ifo])
        numpy.putmask(trig_ifo_cs[ifo], trig_ifo_cs[ifo] == 0, 0.005)

    return trig_ifo_cs


# =============================================================================
# Extract trigger/injection data produced by PyGRB
# =============================================================================
class PygrbFilterOutput(object):
    """Extract trigger/injection data produced by PyGRB search"""
    def __init__(self, trigs_or_injs, ifos, columns, output_type, chisq_index,
                 chisq_nhigh, null_thresh, snrThresh, snglSnrThresh,
                 newSnrThresh, nullGradThresh, nullGradVal, verbose=False):
        if verbose:
            sys.stdout.write("Extracting data related to the %s just " +
                             "loaded...\n" % output_type)
        # Initialize all content of self
        self.time = None
        self.snr = numpy.array(None)
        self.reweighted_snr = None
        self.null_snr = None
        self.null_stat = None
        self.trace_snr = None
        self.chi_square = numpy.array(None)
        self.bank_veto = None
        self.auto_veto = None
        self.coinc_snr = None
        self.ifo_snr = dict((ifo, None) for ifo in ifos)
        self.first_snr = None
        self.second_snr = None
        self.third_snr = None
        self.ifo_bank_cs = dict((ifo, None) for ifo in ifos)
        self.ifo_auto_cs = dict((ifo, None) for ifo in ifos)
        self.ifo_stan_cs = dict((ifo, None) for ifo in ifos)
        self.rel_amp_1 = None
        self.norm_3 = None
        self.rel_amp_2 = None
        self.inclination = None
        # Exctract data and fill in content of self
        if trigs_or_injs is not None:
            # Work out if using sngl chisqs
            ifo_att = {'G1': 'g', 'H1': 'h1', 'H2': 'h2', 'L1': 'l', 'V1': 'v',
                       'T1': 't'}
            i = ifo_att[ifos[0]]

            self.sngl_chisq = 'chisq_%s' % i in columns
            self.sngl_bank_chisq = 'bank_chisq_%s' % i in columns
            self.sngl_cont_chisq = 'cont_chisq_%s' % i in columns

            # Set basic data
            self.time = numpy.asarray(trigs_or_injs.get_end())
            self.snr = numpy.asarray(trigs_or_injs.get_column('snr'))
            self.reweighted_snr = [get_bestnr(t, q=chisq_index, n=chisq_nhigh,
                                              null_thresh=null_thresh,
                                              snr_threshold=snrThresh,
                                              sngl_snr_threshold=snglSnrThresh,
                                              chisq_threshold=newSnrThresh,
                                              null_grad_thresh=nullGradThresh,
                                              null_grad_val=nullGradVal)
                                   for t in trigs_or_injs]
            self.reweighted_snr = numpy.array(self.reweighted_snr)
            self.null_snr = numpy.asarray(trigs_or_injs.get_null_snr())
            self.null_stat = numpy.asarray(trigs_or_injs.get_column(
                'null_statistic'))
            self.trace_snr = numpy.asarray(trigs_or_injs.get_column(
                'null_stat_degen'))

            # Get chisq data
            self.chi_square = numpy.asarray(trigs_or_injs.get_column('chisq'))
            self.bank_veto = numpy.asarray(trigs_or_injs.get_column(
                'bank_chisq'))
            self.auto_veto = numpy.asarray(trigs_or_injs.get_column(
                'cont_chisq'))
            numpy.putmask(self.chi_square, self.chi_square == 0, 0.005)
            numpy.putmask(self.bank_veto, self.bank_veto == 0, 0.005)
            numpy.putmask(self.auto_veto, self.auto_veto == 0, 0.005)

            # Get single detector data
            self.coinc_snr = (trigs_or_injs.get_column('coinc_snr'))
            self.ifo_snr = dict((ifo, trigs_or_injs.get_sngl_snr(ifo))
                                for ifo in ifos)
            tmp = numpy.sort(self.ifo_snr.values(), 0)
            if ifos:
                self.first_snr = tmp[-1, :]
            if len(ifos) > 1:
                self.second_snr = tmp[-2, :]
            if len(ifos) > 2:
                self.third_snr = tmp[-3, :]
            if self.sngl_bank_chisq:
                self.ifo_bank_cs = trigs_or_injs.get_sngl_bank_chisqs(ifos)
                self.ifo_bank_cs = format_single_chisqs(self.ifo_bank_cs, ifos)
            if self.sngl_cont_chisq:
                self.ifo_auto_cs = trigs_or_injs.get_sngl_cont_chisqs(ifos)
                self.ifo_auto_cs = format_single_chisqs(self.ifo_auto_cs, ifos)
            if self.sngl_chisq:
                self.ifo_stan_cs = trigs_or_injs.get_sngl_chisqs(ifos)
                self.ifo_stan_cs = format_single_chisqs(self.ifo_stan_cs, ifos)

            # Initiate amplitude generator
            num_amp = 4
            amplitudes = xrange(1, num_amp+1)

            # Get amplitude terms
            amp = dict((amplitude,
                        numpy.asarray(trigs_or_injs.get_column(
                            'amp_term_%d' % amplitude)))
                       for amplitude in amplitudes)
            # print numpy.count_nonzero(amp[1]) # All trig_amp's are 0
            # print numpy.count_nonzero(amp[2]) # Hence the 3 warnings
            # print numpy.count_nonzero(amp[3])
            # print numpy.count_nonzero(amp[4])
            self.rel_amp_1 = numpy.sqrt((amp[1]**2 + amp[2]**2) /
                                        (amp[3]**2 + amp[4]**2))
            gamma_r = amp[1] - amp[4]
            gamma_i = amp[2] + amp[3]
            delta_r = amp[1] + amp[4]
            delta_i = amp[3] - amp[2]
            norm_1 = delta_r*delta_r + delta_i*delta_i
            norm_2 = gamma_r*gamma_r + gamma_i*gamma_i
            self.norm_3 = ((norm_1**0.25) + (norm_2**0.25))**2
            amp_plus = (norm_1)**0.5 + (norm_2)**0.5
            amp_cross = abs((norm_1)**0.5 - (norm_2)**0.5)
            self.rel_amp_2 = amp_plus/amp_cross
            self.inclination = amp_cross/self.norm_3

            num_trigs_or_injs = len(trigs_or_injs)
            if num_trigs_or_injs < 1:
                sys.stderr.write("WARNING: No %s found." % output_type)
            elif num_trigs_or_injs >= 1 and verbose:
                sys.stdout.write("%d %s foundd.\n"
                                 % (num_trigs_or_injs, output_type))
            if output_type == "triggers":
                sigma = trigs_or_injs.get_sigmasqs()
                self.sigma_tot = numpy.zeros(num_trigs_or_injs)
                # Get antenna response based parameters
                self.longitude = numpy.degrees(trigs_or_injs.get_column('ra'))
                self.latitude = numpy.degrees(trigs_or_injs.get_column('dec'))
                self.f_resp = dict((ifo, numpy.empty(num_trigs_or_injs))
                                   for ifo in ifos)
                for i in xrange(num_trigs_or_injs):
                    # Calculate f_resp for each IFO if we haven't done so yet
                    f_plus, f_cross = get_det_response(self.longitude[i],
                                                       self.latitude[i],
                                                       self.time[i])
                    for ifo in ifos:
                        self.f_resp[ifo][i] = sum(numpy.array([f_plus[ifo],
                                                               f_cross[ifo]]
                                                             )**2)
                        self.sigma_tot[i] += sigma[ifo][i] * \
                                             self.f_resp[ifo][i]

                for ifo in ifos:
                    self.f_resp[ifo] = self.f_resp[ifo].mean()

                # Normalise trig_sigma
                self.sigma_tot = numpy.array(self.sigma_tot)
                for ifo in ifos:
                    sigma[ifo] = numpy.asarray(sigma[ifo]) / self.sigma_tot

                self.sigma_mean = {}
                self.sigma_max = {}
                self.sigma_min = {}
                for ifo in ifos:
                    try:
                        self.sigma_mean[ifo] = sigma[ifo].mean()
                        self.sigma_max[ifo] = sigma[ifo].max()
                        self.sigma_min[ifo] = sigma[ifo].min()
                    except ValueError:
                        self.sigma_mean[ifo] = 0
                        self.sigma_max[ifo] = 0
                        self.sigma_min[ifo] = 0
        if verbose:
            sys.stdout.write("%s parameters extractedd\n"
                             % (output_type))


# =============================================================================
# Plot contours in a scatter plot where SNR is on the horizontal axis
# =============================================================================

def contour_plotter(axis, snr_vals, contours, colors, vert_spike=False):
    """Plot contours in a scatter plot where SNR is on the horizontal axis"""
    for i in contours:
        plot_vals_x = []
        plot_vals_y = []
        if vert_spike:
            for j in snr_vals:
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


# =============================================================================
# Master plotting fucntion: fits all plotting needs in for PyGRB results
# =============================================================================

def pygrb_plotter(trig_x, trig_y, inj_x, inj_y, inj_file, xlabel, ylabel,
                  fig_path, snr_vals=None, conts=None,
                  shade_cont_value=None, colors=None, vert_spike=False,
                  xlims=None, ylims=None, use_logs=True, verbose=False):
    """Master function to plot PyGRB results"""
    if verbose:
        fig_name = os.path.split(os.path.abspath(fig_path))[1]
        sys.stdout.write(" * %s (%s vs %s)...\n" % (fig_name, xlabel, ylabel))
    # Plot trigger-related quantities
    fig = plt.figure()
    ax = fig.gca()
    if use_logs:
        ax.loglog(trig_x, trig_y, 'bx')
    else:
        ax.plot(trig_x, trig_y, 'bx')
    ax.grid()
    # Plot injection-related quantities
    if inj_file:
        if use_logs:
            ax.loglog(inj_x, inj_y, 'r+')
        else:
            ax.plot(inj_x, inj_y, 'r+')
    # Plot contours
    if conts is not None:
        contour_plotter(ax, snr_vals, conts, colors, vert_spike=vert_spike)
    # Add shading above a specific contour (typically used for vetoed area)
    if shade_cont_value is not None:
        limy = ax.get_ylim()[1]
        polyx = copy.deepcopy(snr_vals)
        polyy = copy.deepcopy(conts[shade_cont_value])
        polyx = numpy.append(polyx, [max(snr_vals), min(snr_vals)])
        polyy = numpy.append(polyy, [limy, limy])
        ax.fill(polyx, polyy, color='#dddddd')
    # Axes: labels and limits
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if xlims:
        ax.set_xlim(xlims)
    if ylims:
        ax.set_ylim(ylims)
    # Wrap up
    save_fig_with_metadata(fig, fig_path)
    # fig_kwds=fig_kwds, title=plot_title,
    # cmd=' '.join(sys.argv),
    # caption=plot_caption)
    plt.close()
