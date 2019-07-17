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
import logging
import argparse
import copy
import numpy
from pycbc.results import save_fig_with_metadata
# TODO: imports to fix/remove
try:
    from glue import segments
    from glue.ligolw import utils, lsctables, ligolw, table
except ImportError:
    pass
try:
    from pylal import MultiInspiralUtils
    from pylal.coh_PTF_pyutils import get_bestnr, get_det_response
    from pylal.coh_PTF_pyutils import readSegFiles
    from pylal.dq import dqSegmentUtils
except ImportError:
    pass
# Only if a backend is not already set ... This should really *not* be done
# here, but in the executables you should set matplotlib.use()
# This matches the check that matplotlib does internally, but this *may* be
# version dependenant. If this is a problem then remove this and control from
# the executables directly.
import matplotlib
if 'matplotlib.backends' not in sys.modules:  # nopep8
    matplotlib.use('agg')
from matplotlib import rc
from matplotlib import pyplot as plt


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

    parser.add_argument("-i", "--ifo", default=None, help="IFO used for IFO " +
                        "specific plots")

    parser.add_argument("--use-sngl-ifo-snr", default=False,
                        action="store_true", help="Plots are vs single IFO " +
                        "SNR, rather than coherent SNR")

    parser.add_argument("--variable", default=None, help="Quantity to plot " +
                        "the vertical axis. Supported choices are: " +
                        "coherent, single, reweighted, or null (for " +
                        "timeeries plots), standard, bank, or auto (for " +
                        "chi-square veto plots), coincident, nullstat, " +
                        "or overwhitened (for null statistics plots)")

    parser.add_argument('--plot-title',
                        help="If given, use this as the plot caption")

    parser.add_argument('--plot-caption',
                        help="If given, use this as the plot caption")

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
# Reset times so that t=0 is corresponds to the GRB trigger time
# =============================================================================

def reset_times(seg_dir, trig_data, inj_data, inj_file):
    """Reset times so that t=0 is corresponds to the GRB trigger time"""
    segs = readSegFiles(seg_dir)
    grb_time = segs['on'][1] - 1
    start = int(min(trig_data.time)) - grb_time
    end = int(max(trig_data.time)) - grb_time
    duration = end-start
    start -= duration*0.05
    end += duration*0.05
    trig_data.time = [t-grb_time for t in trig_data.time]

    if inj_file:
        inj_data.time = [t-grb_time for t in inj_data.time]

    return grb_time, start, end, trig_data, inj_data


# =============================================================================
# Extract trigger/injection data produced by PyGRB
# =============================================================================

class PygrbFilterOutput(object):
    """Extract trigger/injection data produced by PyGRB search"""
    def __init__(self, trigs_or_injs, ifos, columns, output_type, opts):
        logging.info("Extracting data from the %s just loaded...", output_type)
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
        self.ifo_bank_cs = dict((ifo, None) for ifo in ifos)
        self.ifo_auto_cs = dict((ifo, None) for ifo in ifos)
        self.ifo_stan_cs = dict((ifo, None) for ifo in ifos)
        self.rel_amp_1 = None
        self.norm_3 = None
        self.rel_amp_2 = None
        self.inclination = None
        # Exctract data and fill in content of self
        null_thresh = map(float, opts.null_snr_threshold.split(','))
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
            self.reweighted_snr = [get_bestnr(t, q=opts.chisq_index,
                                              n=opts.chisq_nhigh,
                                              null_thresh=null_thresh,
                                              snr_threshold=opts.snr_threshold,
                                              sngl_snr_threshold=opts.sngl_snr_threshold,
                                              chisq_threshold=opts.newsnr_threshold,
                                              null_grad_thresh=opts.null_grad_thresh,
                                              null_grad_val=opts.null_grad_val)
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
            amplitudes = range(1, num_amp+1)

            # Get amplitude terms
            amp = dict((amplitude,
                        numpy.asarray(trigs_or_injs.get_column(
                            'amp_term_%d' % amplitude)))
                       for amplitude in amplitudes)
            #
            # All 0, hence the 3 warnings
            # for i in amplitudes:
            #     print numpy.count_nonzero(amp[amplitudes])
            #
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
                logging.warning("No %s found.", output_type)
            elif num_trigs_or_injs >= 1:
                logging.info("%d %s found.", num_trigs_or_injs, output_type)
            # Deal with the sigma-squares (historically called sigmas here)
            if output_type == "triggers":
                sigma = trigs_or_injs.get_sigmasqs()
                self.sigma_tot = numpy.zeros(num_trigs_or_injs)
                # Get antenna response based parameters
                self.longitude = numpy.degrees(trigs_or_injs.get_column('ra'))
                self.latitude = numpy.degrees(trigs_or_injs.get_column('dec'))
                self.f_resp = dict((ifo, numpy.empty(num_trigs_or_injs))
                                   for ifo in ifos)
                for i in range(num_trigs_or_injs):
                    # Calculate f_resp for each IFO if we haven't done so yet
                    f_plus, f_cross = get_det_response(self.longitude[i],
                                                       self.latitude[i],
                                                       self.time[i])
                    for ifo in ifos:
                        self.f_resp[ifo][i] = sum(numpy.array([f_plus[ifo],
                                                               f_cross[ifo]]
                                                              )**2)
                        self.sigma_tot[i] += (sigma[ifo][i] *
                                              self.f_resp[ifo][i])

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
        logging.info("%s parameters extracted", output_type)


# =============================================================================
# Function to open trigger and injection xml files
# =============================================================================

def load_xml_file(filename):
    """Wrapper to ligolw's utils.load_filename"""

    xml_doc = utils.load_filename(filename, gz=filename.endswith("gz"),
                                  contenthandler=lsctables.use_in(
                                      ligolw.LIGOLWContentHandler))

    return xml_doc


# =============================================================================
# Function to extract ifos
# =============================================================================

def extract_ifos(trig_file):
    """Extracts IFOs from search summary table"""

    # Load search summary
    xml_doc = load_xml_file(trig_file)
    search_summ = table.get_table(xml_doc,
                                  lsctables.SearchSummaryTable.tableName)

    # Extract IFOs
    ifos = sorted(map(str, search_summ[0].get_ifos()))

    return ifos


# =============================================================================
# Function to extract vetoes
# =============================================================================

def extract_vetoes(veto_files, ifos):
    """Extracts vetoes from veto filelist"""

    # Initialize vetoe containers
    vetoes = segments.segmentlistdict()
    for ifo in ifos:
        vetoes[ifo] = segments.segmentlist()

    # Construct veto list from veto filelist
    if veto_files:
        for file in veto_files:
            ifo = os.path.basename(file)[:2]
            if ifo in ifos:
                # This returns a coalesced list of the vetoes
                tmp_veto_segs = dqSegmentUtils.fromsegmentxml(open(file, 'r'))
                for entry in tmp_veto_segs:
                    vetoes[ifo].append(entry)
    for ifo in ifos:
        vetoes[ifo].coalesce()

    return vetoes


# =============================================================================
# Function to load triggers
# =============================================================================

def load_triggers(trig_file, vetoes, ifos):
    """"Loads triggers from PyGRB output file"""
    logging.info("Loading triggers...")

    # Extract time-slides
    multis, slide_dict, _ = \
        MultiInspiralUtils.ReadMultiInspiralTimeSlidesFromFiles([trig_file])
    num_slides = len(slide_dict)
    lsctables.MultiInspiralTable.loadcolumns =\
        [slot for slot in multis[0].__slots__ if hasattr(multis[0], slot)]

    # Extract triggers
    trigs = lsctables.New(lsctables.MultiInspiralTable,
                          columns=lsctables.MultiInspiralTable.loadcolumns)
    logging.info("%d triggers found.", len(trigs))

    # Time-slid vetoes
    for slide_id in range(num_slides):
        slid_vetoes = copy.deepcopy(vetoes)
        for ifo in ifos:
            slid_vetoes[ifo].shift(-slide_dict[slide_id][ifo])

        # Add time-slid triggers
        vets = slid_vetoes.union(slid_vetoes.keys())
        trigs.extend(t for t in multis.veto(vets)
                     if int(t.time_slide_id) == slide_id)

    logging.info("%d triggers found when including timeslides.", len(trigs))

    return trigs


# =============================================================================
# Function to load injections
# =============================================================================

def load_injections(inj_file, vetoes):
    """"Loads injections from PyGRB output file"""
    logging.info("Loading injections...")

    # Load injection file
    xml_doc = load_xml_file(inj_file)
    multis = table.get_table(xml_doc, lsctables.MultiInspiralTable.tableName)

    # Extract injections
    injs = lsctables.New(lsctables.MultiInspiralTable,
                         columns=lsctables.MultiInspiralTable.loadcolumns)

    # Injections in time-slid non-vetoed data
    injs.extend(t for t in multis if t.get_end() not in vetoes)

    logging.info("%d injections found.", len(injs))

    return injs


# =============================================================================
# Function to load injections
# =============================================================================

def new_snr_chisq(snr, new_snr, chisq_dof, chisq_index=4.0, chisq_nhigh=3.0):
    """Returns the chi-square value needed to weight SNR into new SNR"""

    chisqnorm = (snr/new_snr)**chisq_index
    if chisqnorm <= 1:
        return 1E-20

    return chisq_dof * (2*chisqnorm - 1)**(chisq_nhigh/chisq_index)


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

def calculate_contours(trigs, opts, new_snrs=None):
    """Generate the plot contours for chisq variable plots"""

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
    chisq_dof = trigs[0].chisq_dof
    bank_chisq_dof = trigs[0].bank_chisq_dof
    cont_chisq_dof = trigs[0].cont_chisq_dof

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


# =============================================================================
# Contains plotting setups shared by PyGRB plots
# =============================================================================

def pygrb_shared_plot_setups():
    """Master function to plot PyGRB results"""

    # Get rcParams
    rc('font', size=14)
    # Set color for out-of-range values
    plt.cm.spring.set_over('g')


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
