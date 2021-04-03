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
import glob
import os
import logging
import argparse
import copy
import numpy
from scipy import stats
import logging
from pycbc.detector import Detector
from pycbc.results import save_fig_with_metadata
# TODO: imports to fix/remove
try:
    from glue import segments
    from glue.ligolw import utils, lsctables, ligolw, table
    from glue.segmentsUtils import fromsegwizard
except ImportError:
    pass
try:
    #TODO: Look at pycbc/io/hdf.py
    from pylal import MultiInspiralUtils
    from pylal.coh_PTF_pyutils import get_bestnr as pylal_get_bestnr
except ImportError:
    pass
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


# =============================================================================
# Function to calculate the antenna factors F+ and Fx
# =============================================================================
# TODO: the call, e.g., Detector("H1", reference_time=None) will not always
# work because we are on Python 2.7 and therefore on an old version of astropy
# which cannot download recent enough IERS tables. TEMPORARILY use the
# default time (GW150914) as reference, thus approximating the sidereal time.
def get_antenna_factors(antenna, ra, dec, geocent_time):
    """Returns the antenna responses F+ and Fx of an IFO (passed as pycbc """
    """Detector type) at a given sky location and time."""

    f_plus, f_cross = antenna.antenna_pattern(ra, dec, 0, geocent_time)

    return f_plus, f_cross


# =============================================================================
# Function to calculate the antenna response F+^2 + Fx^2
# =============================================================================
def get_antenna_single_response(antenna, ra, dec, geocent_time):
    """Returns the antenna response F+^2 + Fx^2 of an IFO (passed as pycbc """
    """Detector type) at a given sky location and time."""

    fp, fc = get_antenna_factors(antenna, ra, dec, geocent_time)

    return fp**2 + fc**2

# Vectorize the function above on all but the first argument
get_antenna_responses = numpy.vectorize(get_antenna_single_response,\
                                    otypes=[float])
get_antenna_responses.excluded.add(0)


# =============================================================================
# Wrapper to avoid import pylal in executables: switch here from this pylal
# function to PyCBC functions when ready.
# 
# Function to calculate the detection statistic
# =============================================================================
def get_bestnr(trig, q=4.0, n=3.0, null_thresh=(4.25,6), snr_threshold=6.,\
               sngl_snr_threshold=4., chisq_threshold = None,\
               null_grad_thresh=20., null_grad_val = 0.2):

    bestNR = pylal_get_bestnr(trig, q=q, n=n, null_thresh=null_thresh,\
                              snr_threshold=snr_threshold,\
                              sngl_snr_threshold=sngl_snr_threshold,\
                              chisq_threshold=chisq_threshold,\
                              null_grad_thresh=null_grad_thresh,\
                              null_grad_val=null_grad_val)

    return bestNR


#def get_bestnr( trig, q=4.0, n=3.0, null_thresh=(4.25,6), snr_threshold=6.,\
#                sngl_snr_threshold=4., chisq_threshold = None,\
#                null_grad_thresh=20., null_grad_val = 1./5.):
#    """
#    Calculate BestNR (coh_PTF detection statistic) through signal based vetoes:
#    The signal based vetoes are as follows:
#      * Coherent SNR < 6
#      * Bank chi-squared reduced (new) SNR < 6
#      * Auto veto reduced (new) SNR < 6
#      * Single-detector SNR (from two most sensitive IFOs) < 4
#      * Null SNR (CoincSNR^2 - CohSNR^2)^(1/2) < nullthresh
#    Returns BestNR as float
#    """
#
#    snr = trig.snr
#    if not chisq_threshold:
#        chisq_threshold = snr_threshold
#
#    # coherent SNR and null SNR cut
#    if (snr < snr_threshold) \
#         or (trig.get_new_snr(index=q, nhigh=n, column='bank_chisq')\
#                          < chisq_threshold) \
#         or (trig.get_new_snr(index=q, nhigh=n, column='cont_chisq')\
#                          < chisq_threshold):
#        return 0
#
#    # define IFOs for sngl cut
#    ifos = map(str,trig.get_ifos())
#
#    # single detector SNR cut
#    sens = {}
#    fPlus, fCross = get_det_response(numpy.degrees(trig.ra),\
#                                     numpy.degrees(trig.dec),\
#                                     trig.get_end())
#    for ifo in ifos:
#        if ifo.lower()[0] == 'h':
#            i = ifo.lower()
#        else:
#            i = ifo[0].lower()
##        sens[ifo] = getattr(trig, 'sigmasq_%s' % i.lower()) * \
##                        sum(numpy.array([fPlus[ifo], fCross[ifo]])**2)
#        sens[ifo] = getattr(trig, 'sigmasq_%s' % i.lower()) * \
#                        get_antenna_single_response(ifo, self.ra[i],
#                                                    self.dec[i], self.time[i])
#    ifos.sort(key=lambda ifo: sens[ifo], reverse=True)
#    if len(ifos) > 1:
#        for i in xrange(0, 2):
#            if ifos[i].lower()[0] == 'h':
#                i = ifos[i].lower()
#            else:
#                i = ifos[i][0].lower()
#            if getattr(trig, 'snr_%s' % i) < sngl_snr_threshold:
#                return 0
#
#    # get chisq reduced (new) SNR
#    bestNR = trig.get_bestnr(index=q, nhigh=n, \
#             null_snr_threshold=null_thresh[0], \
#             null_grad_thresh=null_grad_thresh, null_grad_val=null_grad_val)
#
#    # If we got this far, the bestNR is non-zero. Verify that chisq actually
#    # was calculated for the trigger
#    if trig.chisq == 0:
#        # Some stuff for debugging
#        print >> sys.stderr,\
#            "Chisq not calculated for trigger with end time and snr:"
#        print >> sys.stderr,  trig.get_end(),trig.snr
#        raise ValueError("Chisq has not been calculated for trigger.")
#
#    return bestNR


# =============================================================================
# Parse command line
# =============================================================================

# TODO: regroup options that are now all in this unique parser
def pygrb_plot_opts_parser(usage='', description=None, version=None):
    """Parses options for PyGRB post-processing scripts"""
    parser = argparse.ArgumentParser(usage=usage, description=description,
                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--version", action="version", version=version)

    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="verbose output")

    parser.add_argument("-t", "--trig-file", action="store",
                        default=None, #required=True,
                        help="The location of the trigger file")

    parser.add_argument("-I", "--inj-file", action="store", default=None,
                        help="The location of the injection file")

    parser.add_argument("-a", "--segment-dir", action="store",
                        required=True, help="directory holding buffer, on " +
                        "and off source segment files.")

    parser.add_argument("-o", "--output-file", default=None, #required=True,
                        help="Output file.")

    parser.add_argument("-O", "--zoomed-output-file", default=None,
                        required=False, help="Output file for a zoomed in " +
                        "version of the plot.")

    parser.add_argument("--x-lims", action="store", default=None,
                        help="Comma separated minimum and maximum values "+
                        "for the horizontal axis. When using negative values "+
                        "an equal sign after --x-lims is necessary.")

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

    parser.add_argument("-T", "--null-grad-thresh", action="store", type=float,
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

    # These 3 are for found/missed injections plots
    parser.add_argument("--found-missed-output-file", default=None, #required=True,
                        help="Found-missed injections output file.")

    parser.add_argument("--missed-found-output-file", default=None, #required=True,
                        help="Missed-found injections output file.")

    parser.add_argument("-x", "--x-variable", default=None, help="Quantity " +
                        "to plot on the horizontal axis. Supported choices " +
                        "are: distance, mchirp, time (for sky error plots).")

    # Html and h5 output files of pycbc_pygrb_page_tables
    parser.add_argument("--quiet-found-injs-output-file", default=None, #required=True,
                        help="Quiet-found injections html output file.")

    parser.add_argument("--missed-found-injs-output-file", default=None, #required=True,
                        help="Missed-found injections html output file.")

    parser.add_argument("--quiet-found-injs-h5-output-file", default=None, #required=True,
                        help="Quiet-found injections h5 output file.")

    parser.add_argument("--loudest-offsource-trigs-output-file", default=None, #required=True,
                        help="Loudest offsource triggers html output file.")

    parser.add_argument("--loudest-offsource-trigs-h5-output-file", default=None, #required=True,
                        help="Loudest offsource triggers h5 output file.")
 
    parser.add_argument("--loudest-onsource-trig-output-file", default=None, #required=True,
                        help="Loudest onsource trigger html output file.")

    parser.add_argument("--loudest-onsource-trig-h5-output-file", default=None, #required=True,
                        help="Loudest onsource trigger h5 output file.")

    parser.add_argument("--num-loudest-off-trigs", action="store",
                        type=int, default=30, help="Number of loudest " +
                        "offsouce triggers to output details about.")


    # This is originally for SNR and chi-square veto plots 
    parser.add_argument("-y", "--y-variable", default=None, help="Quantity " +
                        "to plot on the vertical axis. Supported choices " +
                        "are: coherent, single, reweighted, or null (for " +
                        "timeseries plots), standard, bank, or auto (for " +
                        "chi-square veto plots), coincident, nullstat, " +
                        "or overwhitened (for null statistics plots)")

    parser.add_argument('--plot-title', default=None,
                        help="If given, use this as the plot caption")

    parser.add_argument('--plot-caption', default=None,
                        help="If given, use this as the plot caption")
    
    # This is originally for SNR timeseries plots 
    parser.add_argument('--central-time', type=float, default=None,
                        help="Center plot on the given GPS time. If omitted, "+
                        "use the GRB trigger time")

    # pygrb_efficiency only options: start here
    # Does this differ from trig-file? --> it doesn't contain the onsource.
    # NB: for now removed the requirement to specify trig_file.
    # instead added a requirement at the end of the parser function
    # to specify either offsource-file or trig-file.
    parser.add_argument("-F", "--offsource-file", action="store",
                        default=None, help="Location of off-source trigger file")

    # As opposed to offsource-file and trig-file, this only contains onsource
    parser.add_argument("--onsource-file", action="store",
                        default=None, help="Location of on-source trigger file")

    parser.add_argument("-f", "--found-file", action="store",
                        default=None,
                        help="Location of the found injections file")

    parser.add_argument("-m", "--missed-file",action="store",
                        default=None,
                        help="Location of the missed injections file")

    parser.add_argument("--background-output-file", default=None, #required=True,
                        help="Detection efficiency output file.")

    parser.add_argument("--onsource-output-file", default=None, #required=True,
                        help="Exclusion distance output file.")

    # FIXME: eventually remove below argument and require output-file
    # be specified. 
    parser.add_argument("--output-path", default=os.getcwd(), 
                        help="Output path for plots")

    parser.add_argument("-s", "--segment-length", action="store", type=float,
                        default=None, help="The length of analysis segments.")

    parser.add_argument("-e", "--pad-data", action="store", type=float,
                        default=None,
                        help="The length of padding around analysis chunk.")

    parser.add_argument("-g", "--glitch-check-factor", action="store",
                        type=float,default=1.0, help="When deciding " +
                        "exclusion efficiencies this value is multiplied " +
                        "to the offsource around the injection trigger to " +
                        "determine if it is just a loud glitch.")

    parser.add_argument("-C", "--cluster-window", action="store", type=float,
                        default=0.1,help="The cluster window used " +
                        "to cluster triggers in time.")

    parser.add_argument("-U", "--upper-inj-dist", action="store",
                        type=float,default=1000,help="The upper distance " +
                        "of the injections in Mpc, if used.")

    parser.add_argument("-L", "--lower-inj-dist", action="store",
                        type=float,default=0,help="The lower distance of " +
                        "the injections in Mpc, if used.")

    parser.add_argument("-n", "--num-bins", action="store", type=int,
                        default=0,help="The number of bins used to " +
                        "calculate injection efficiency.")

    parser.add_argument("-M", "--num-mc-injections", action="store",
                        type=int, default=100, help="Number of Monte " +
                        "Carlo injection simulations to perform.")

    parser.add_argument("-S", "--seed", action="store", type=int,
                        default=1234, help="Seed to initialize Monte Carlo.")

    parser.add_argument("-w", "--waveform-error", action="store",
                        type=float, default=0, help="The standard " +
                        "deviation to use when calculating the waveform error.")

    parser.add_argument("--h1-cal-error", action="store", type=float,
                        default=0, help="The standard deviation to use when " +
                        "calculating the H1 calibration amplitude error.")

    parser.add_argument("--k1-cal-error", action="store", type=float,
                        default=0, help="The standard deviation to use when " +
                        "calculating the K1 calibration amplitude error.")

    parser.add_argument("--l1-cal-error", action="store", type=float,
                        default=0, help="The standard deviation to use when " +
                        "calculating the L1 calibration amplitude error.")

    parser.add_argument("--v1-cal-error", action="store", type=float,
                        default=0, help="The standard deviation to use when " +
                        "calculating the V1 calibration amplitude error.")

    parser.add_argument("--h1-dc-cal-error", action="store",
                        type=float, default=1.0, help="The scaling factor " +
                        "to use when calculating the H1 calibration " +
                        "amplitude error.")

    parser.add_argument("--k1-dc-cal-error", action="store",
                        type=float, default=1.0, help="The scaling factor " +
                        "to use when calculating the K1 calibration " +
                        "amplitude error.")

    parser.add_argument("--l1-dc-cal-error", action="store",
                        type=float, default=1.0, help="The scaling factor " +
                        "to use when calculating the L1 calibration " +
                        "amplitude error.")

    parser.add_argument("--v1-dc-cal-error", action="store",
                        type=float, default=1.0, help="The scaling factor " +
                        "to use when calculating the V1 calibration " +
                        "amplitude error.")
    # pygrb_efficiency only options end here

    # pycbc_pygrb_slice_inj_followup uses this
    parser.add_argument("--followup-inj-idx", action="store", type=int, 
                        default=None, help="Index of missed/quiet injeciton " + 
                        "to followup.")

    args = parser.parse_args()
    if not (args.trig_file or args.offsource_file):
        parser.error('Must specify either trig-file or offsource-file.')

    return args


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
# Wrapper to read segments files
# =============================================================================

def read_seg_files(seg_dir):
    """Given the segments directory, read segments files"""

    times = {}
    keys = ["buffer", "off", "on"]
    file_names = ["bufferSeg.txt","offSourceSeg.txt","onSourceSeg.txt"]

    for key,file_name in zip(keys, file_names):
        segs = fromsegwizard(open(os.path.join(seg_dir,file_name), 'r'))
        if len(segs)>1:
            logging.error('More than one segment, an error has occured.')
            sys.exit()
        times[key] = segs[0]

    return times


# =============================================================================
# Find GRB trigger time
# =============================================================================

def get_grb_time(seg_dir):
    """Determine GRB trigger time"""
    segs = read_seg_files(seg_dir)
    grb_time = segs['on'][1] - 1

    return grb_time


# =============================================================================
# Find start and end times of trigger/injecton data relative to a given time
# =============================================================================
def get_start_end_times(data, central_time):
    """Determine start and end times of data relative to central_time"""
    start = int(min(data.time)) - central_time
    end = int(max(data.time)) - central_time
    duration = end-start
    start -= duration*0.05
    end += duration*0.05

    return start, end


# =============================================================================
# Reset times so that t=0 is corresponds to the given trigger time
# =============================================================================

def reset_times(data, trig_time):
    """Reset times in data so that t=0 corresponds to the trigger time provided"""
    data.time = [t-trig_time for t in data.time]

    return data

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
                # Get antenna response based parameters
                self.ra = trigs_or_injs.get_column('ra')
                self.longitude = numpy.degrees(self.ra)
                self.dec = trigs_or_injs.get_column('dec')
                self.latitude = numpy.degrees(self.dec)
                self.f_resp = {}
                sigma = trigs_or_injs.get_sigmasqs()
                self.sigma_tot = numpy.zeros(num_trigs_or_injs)
                for ifo in ifos:
                    antenna = Detector(ifo)
                    self.f_resp[ifo] = get_antenna_responses(antenna,
                                                             self.ra,
                                                             self.dec,
                                                             self.time)
                    self.sigma_tot += (sigma[ifo] * self.f_resp[ifo])
                    # After the detailed calculations, this only stores the mean responses
                    self.f_resp[ifo] = self.f_resp[ifo].mean()

                # Normalise trig_sigma
                for ifo in ifos:
                    sigma[ifo] /= self.sigma_tot

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

def load_xml_file(file_name):
    """Wrapper to ligolw's utils.load_filename"""

    xml_doc = utils.load_filename(file_name, gz=file_name.endswith("gz"),
                                  contenthandler=lsctables.use_in(
                                      ligolw.LIGOLWContentHandler))

    return xml_doc


# =============================================================================
# Function to load a table from an xml file
# =============================================================================

def load_xml_table(file_name, table_name):
    """Load xml table from file."""

    xml_doc = load_xml_file(file_name)

    return table.get_table(xml_doc, table_name)


# =============================================================================
# Function to extract ifos
# =============================================================================

def extract_ifos(trig_file):
    """Extracts IFOs from search summary table"""

    # Load search summary
    search_summ = load_xml_table(trig_file,
                                 lsctables.SearchSummaryTable.tableName)

    # Extract IFOs
    ifos = sorted(map(str, search_summ[0].get_ifos()))

    return ifos


# ==============================================================================
# Function to load segments from an xml file
# ==============================================================================
def load_segments_from_xml(xml_doc, return_dict=False, select_id=None):
    """Read a glue.segments.segmentlist from the file object file containing an
    xml segment table.

    Parameters
    ----------
        xml_doc: name of segment xml file

        Keyword Arguments:
            return_dict : [ True | False ]
                return a glue.segments.segmentlistdict containing coalesced
                glue.segments.segmentlists keyed by seg_def.name for each entry
                in the contained segment_def_table. Default False
            select_id : int
                return a glue.segments.segmentlist object containing only
                those segments matching the given segment_def_id integer

    """

    # Load SegmentDefTable and SegmentTable
    seg_def_table = load_xml_table(xml_doc, lsctables.SegmentDefTable.tableName)
    seg_table = load_xml_table(xml_doc, lsctables.SegmentTable.tableName)

    if return_dict:
        segs = segments.segmentlistdict()
    else:
        segs = segments.segmentlist()

    seg_id = {}
    for seg_def in seg_def_table:
        seg_id[int(seg_def.segment_def_id)] = str(seg_def.name)
        if return_dict:
            segs[str(seg_def.name)] = segments.segmentlist()

    for seg in seg_table:
        if return_dict:
            segs[seg_id[int(seg.segment_def_id)]]\
                .append(segments.segment(seg.start_time, seg.end_time))
            continue
        if select_id and int(seg.segment_def_id)==select_id:
            segs.append(segments.segment(seg.start_time, seg.end_time))
            continue
        segs.append(segments.segment(seg.start_time, seg.end_time))

    if return_dict:
       for seg_name in seg_id.values():
           segs[seg_name] = segs[seg_name].coalesce()
    else:
        segs = segs.coalesce()

    return segs


# =============================================================================
# Function to extract vetoes
# =============================================================================

def extract_vetoes(veto_files, ifos):
    """Extracts vetoes from veto filelist"""

    # Initialize veto containers
    vetoes = segments.segmentlistdict()
    for ifo in ifos:
        vetoes[ifo] = segments.segmentlist()

    # Construct veto list from veto filelist
    if veto_files:
        for veto_file in veto_files:
            ifo = os.path.basename(veto_file)[:2]
            if ifo in ifos:
                # This returns a coalesced list of the vetoes
                tmp_veto_segs = load_segments_from_xml(veto_file)
                for entry in tmp_veto_segs:
                    vetoes[ifo].append(entry)
    for ifo in ifos:
        vetoes[ifo].coalesce()

    return vetoes


# =============================================================================
# Function to extract IFOs and vetoes
# =============================================================================

def extract_ifos_and_vetoes(trig_file, veto_dir, veto_cat):
    """Extracts IFOs from search summary table and vetoes from a directory"""

    # Extract IFOs 
    ifos = extract_ifos(trig_file)

    # Extract vetoes
    veto_files = []
    if veto_dir:
        veto_string = ','.join([str(i) for i in range(2,veto_cat+1)])
        veto_files = glob.glob(veto_dir +'/*CAT[%s]*.xml' %(veto_string))
    vetoes = extract_vetoes(veto_files, ifos)

    return ifos, vetoes


# =============================================================================
# Function to load triggers
# =============================================================================

def load_triggers(trig_file, vetoes, ifos):
    """Loads triggers from PyGRB output file"""
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
    """Loads injections from PyGRB output file"""
    logging.info("Loading injections...")

    # Load injections in injection file
    multis = load_xml_table(inj_file, lsctables.MultiInspiralTable.tableName)

    # Extract injections in time-slid non-vetoed data
    injs = lsctables.New(lsctables.MultiInspiralTable,
                         columns=lsctables.MultiInspiralTable.loadcolumns)
    injs.extend(t for t in multis if t.get_end() not in vetoes)

    logging.info("%d injections found.", len(injs))

    return injs


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
# Function to get the ID numbers from a LIGO-LW table
# =============================================================================

def get_id_numbers(ligolw_table, column):
    """Grab the IDs of a LIGO-LW table"""

    ids = [int(getattr(row, column)) for row in ligolw_table]

    return ids


# =============================================================================
# Function to load timeslides
# =============================================================================

def load_time_slides(xml_file):
    """Loads timeslides from PyGRB output file"""

    time_slide = load_xml_table(xml_file, lsctables.TimeSlideTable.tableName)
    time_slide_unsorted = [dict(i) for i in time_slide.as_dict().values()]
    # NB: sorting not necessary if using python 3
    sort_idx = numpy.argsort(numpy.array([
            int(time_slide.get_time_slide_id(ov)) for ov in time_slide_unsorted
    ]))
    time_slide_list = numpy.array(time_slide_unsorted)[sort_idx]
    # Check time_slide_ids are ordered correctly. 
    ids = get_id_numbers(time_slide, "time_slide_id")[::len(time_slide_list[0].keys())]
    if not (numpy.all(ids[1:] == numpy.array(ids[:-1])+1) and ids[0]==0):
        err_msg = "time_slide_ids list should start at zero and increase by "
        err_msg += "one for every element"
        logging.err(err_msg)
        sys.exit()
    # Check that the zero-lag slide has time_slide_id == 0.
    if not numpy.all(numpy.array(list(time_slide_list[0].values())) == 0):
        err_msg = "The zero-lag slide should have time_slide_id == 0 "
        err_msg += "but the first element of time_slide_list is "
        err_msg += "%s \n" % time_slide_list[0]
        logging.err(err_msg)
        sys.exit()

    return time_slide_list


# =============================================================================
# Function to determine the id of the zero-lag timeslide
# =============================================================================

def find_zero_lag_slide_id(slide_dict):
    """Loads timeslides from PyGRB output file"""

    zero_lag_slide_id = None
    num_slides = len(slide_dict)
    for i in range(num_slides):
        # Is this slide the zero-lag?
        if max(slide_dict[i].values()) == 0:
            if min(slide_dict[i].values()) == 0:
                if zero_lag_slide_id is None:
                    zero_lag_slide_id = i
                else:
                    err_msg = 'zero_lag_slide_id was already assigned: there'
                    err_msg += 'seems to be more than one zero-lag slide!'
                    logging.error(err_msg)
                    sys.exit()
    
    if zero_lag_slide_id is None:
        err_msg = 'Unable to assign zero_lag_slide_id: '
        err_msg += 'there seems to be no zero-lag slide!'
        logging.error(err_msg)
        sys.exit()

    return zero_lag_slide_id


# =============================================================================
# Function to calculate the error bars and fraction of recovered injections
# (used for efficiency/distance plots)
# =============================================================================

def efficiency_with_errs(found_bestnr, num_injections, num_mc_injs=0):
    """Calculates the fraction of recovered injections and its error bars"""

    if not isinstance(num_mc_injs, int):
        err_msg = "The parameter num_mc_injs is the number of Monte-Carlo "
        err_msg += "injections.  It must be an integer."
        logging.error(err_msg)
        sys.exit()

    only_found_injs = found_bestnr[:-1]
    all_injs = num_injections[:-1]
    fraction = only_found_injs / all_injs

    # Divide by Monte-Carlo iterations
    if num_mc_injs:
        only_found_injs = only_found_injs / num_mc_injs
        all_injs = all_injs / num_mc_injs

    # TODO: optimize
    err_common = all_injs * (2 * only_found_injs + 1)
    err_denom = 2 * all_injs * (all_injs + 1)
    err_vary = 4 * all_injs * only_found_injs * (all_injs - only_found_injs) \
                + all_injs**2
    err_vary = err_vary**0.5
    err_low = (err_common - err_vary)/err_denom
    err_low_mc = fraction - err_low
    err_high = (err_common + err_vary)/err_denom
    err_high_mc = err_high - fraction

    return err_low_mc, err_high_mc, fraction


# =============================================================================
# Function to load the segment dicitonary
# =============================================================================

def load_segment_dict(xml_file):
    """Loads the segment dictionary """

    # Get the mapping table
    # TODO: unclear whether this step is necessary (seems the 
    # segment_def_id and time_slide_id are always identical)
    time_slide_map_table = load_xml_table(xml_file, lsctables.TimeSlideSegmentMapTable.tableName)
    segment_map = {
            int(entry.segment_def_id): int(entry.time_slide_id) 
            for entry in time_slide_map_table
    }
    # Extract the segment table
    segment_table = load_xml_table(
        xml_file,lsctables.SegmentTable.tableName
        )
    segmentDict = {}
    for entry in segment_table:
        currSlidId = segment_map[int(entry.segment_def_id)]
        currSeg = entry.get()
        if not currSlidId in segmentDict.keys():
            segmentDict[currSlidId] = segments.segmentlist()
        segmentDict[currSlidId].append(currSeg)
        segmentDict[currSlidId].coalesce()

    return segmentDict


# =============================================================================
# Construct the trials from the timeslides, segments, and vetoes 
# =============================================================================

def construct_trials(num_slides, segs, segment_dict, ifos, slide_dict, vetoes):
    """Constructs trials from triggers, timeslides, segments and vetoes"""

    trial_dict = {}

    # Separate segments
    trial_time = abs(segs['on'])

    for slide_id in range(num_slides):
        # These can only *reduce* the analysis time
        curr_seg_list = segment_dict[slide_id]
    
        # Construct the buffer segment list
        seg_buffer = segments.segmentlist()
        for ifo in ifos:
            slide_offset = slide_dict[slide_id][ifo]
            seg_buffer.append(segments.segment(segs['buffer'][0] - slide_offset,\
                                               segs['buffer'][1] - slide_offset))
        seg_buffer.coalesce()
    
        # Construct the ifo list
        slid_vetoes = copy.deepcopy(vetoes)
        for ifo in ifos:
            slid_vetoes[ifo].shift(-slide_dict[slide_id][ifo])
    
        # Construct trial list and check against buffer
        trial_dict[slide_id] = segments.segmentlist()
        for curr_seg in curr_seg_list:
            iter_int = 0
            while 1:
                if (curr_seg[0] + trial_time*(iter_int+1)) > curr_seg[1]:
                    break
                curr_trial = segments.segment(curr_seg[0] + trial_time*iter_int,\
                                              curr_seg[0] + trial_time*(iter_int+1))
                if not seg_buffer.intersects_segment(curr_trial):
                    for ifo in ifos:
                        if slid_vetoes[ifo].intersects_segment(curr_trial):
                            break
                    else:
                        trial_dict[slide_id].append(curr_trial)
                iter_int += 1
    
    return trial_dict


# =============================================================================
# Construct the sorted triggers from the trials
# =============================================================================

def sort_trigs(trial_dict, trigs, num_slides, segment_dict):
    """Constructs sorted triggers"""

    sorted_trigs = {}

    # Begin by sorting the triggers into each slide
    # New seems pretty slow, so run it once and then use deepcopy
    tmp_table = lsctables.New(lsctables.MultiInspiralTable)
    for slide_id in range(num_slides):
        sorted_trigs[slide_id] = copy.deepcopy(tmp_table)
    for trig in trigs:
        sorted_trigs[int(trig.time_slide_id)].append(trig)

    for slide_id in range(num_slides):
        # These can only *reduce* the analysis time
        curr_seg_list = segment_dict[slide_id]
    
        ###### TODO: below is a check we can possibly remove #####
        # Check the triggers are all in the analysed segment lists
        for trig in sorted_trigs[slide_id]:
            if trig.end_time not in curr_seg_list:
                # This can be raised if the trigger is on the segment boundary, so
                # check if the trigger is within 1/100 of a second within the list
                if trig.get_end() + 0.01 in curr_seg_list:
                    continue
                if trig.get_end() - 0.01 in curr_seg_list:
                    continue
                err_msg = "Triggers found in input files not in the list of "
                err_msg += "analysed segments. This should not happen."
                logging.error(err_msg)
                sys.exit()
        ###### end of check #####
    
        # The below line works like the inverse of .veto and only returns trigs
        # that are within the segment specified by trial_dict[slide_id]
        sorted_trigs[slide_id] = sorted_trigs[slide_id].vetoed(trial_dict[slide_id])
    
    return sorted_trigs


# =============================================================================
# Find max and median of loudest SNRs or BestNRs
# =============================================================================

def sort_stat(time_veto_max_stat):
    """Sort a dictionary of loudest SNRs/BestNRs"""

    full_time_veto_max_stat = numpy.concatenate(time_veto_max_stat.values())
    full_time_veto_max_stat.sort()

    return full_time_veto_max_stat


# =============================================================================
# Find max and median of loudest SNRs or BestNRs
# =============================================================================

def max_median_stat(num_slides, time_veto_max_stat, trig_stat, total_trials):
    """Deterime the maximum and median of the loudest SNRs/BestNRs"""

    max_stat = max([trig_stat[slide_id].max() if trig_stat[slide_id].size \
                   else 0 for slide_id in range(num_slides)])
    
    full_time_veto_max_stat = sort_stat(time_veto_max_stat)

    if total_trials % 2:
        median_stat = full_time_veto_max_stat[(total_trials - 1) // 2]
    else:
        median_stat = numpy.mean((full_time_veto_max_stat)\
                              [total_trials//2 - 1 : total_trials//2 + 1])

    return max_stat, median_stat, full_time_veto_max_stat 


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



# =============================================================================
# Incorporate calibration and waveform errors for efficiency plots
# =============================================================================
# Used by pycbc_pygrb_page_tables

def mc_cal_wf_errs(num_mc_injs, num_injs, inj_dists, cal_err, wf_err, max_dc_cal_err):
    """Includes calibration and waveform errors by running an MC"""

    # The numbers for the efficiency plots include calibration and waveform
    # errors incorporated by running over each injection num_mc_injs times,
    # where each time we draw a random value of distance.

    # num_injs = len(num_mc_injs) or len(inj_dists)

    inj_dist_mc = numpy.ndarray((num_mc_injs+1, num_injs))
    inj_dist_mc[0, :] = inj_dists
    for i in range(num_mc_injs):
        cal_dist_red = stats.norm.rvs(size=num_injs) * cal_err
        wf_dist_red = numpy.abs(stats.norm.rvs(size=num_injs) * wf_err)
        inj_dist_mc[i+1, :] = inj_dists / (max_dc_cal_err * \
                                     (1 + cal_dist_red) * (1 + wf_dist_red))

    return inj_dist_mc


# =============================================================================
# Dumping material from pycbc_pygrb_efficiency and pycbc_pygrb_page_tables
# that is not strictly needed but that we may want to turn into functions,
# additional tools, unit tests
# =============================================================================
# Calculate and print fraction of trials with an event (highest FAP value)
#num_events = 0
#for slide_id in range(num_slides):
#    for trial in range(len(trial_dict[slide_id])):
#        if time_veto_max_bestnr[slide_id][trial] > 0:
#            num_events += 1
#quietest_file = open('%s/quiet_fap_val.txt' % outdir, 'w')
#quietest_file.write('%s' % (num_events/total_trials))
#quietest_file.close()


# Write inclination recovery to file
# GRB start time
#grb_time = segs['on'][1] - 1
#f_incl_txt = open('%s/found_inclinations.txt' % outdir, 'w')
#f_incl_txt.write('GPS time\tTime since %d\tInclination\n\n' % grb_time)
#stacked = np.column_stack([found_injs["time"][zero_fap],
#                           found_injs["time"][zero_fap] - grb_time,
#                           found_injs["inclination"][zero_fap]])
#np.savetxt(f_incl_txt, stacked, delimiter='\t')
#f_incl_txt.close()
#t_incl_txt = open('%s/total_inclinations.txt' % outdir, 'w')
#t_incl_txt.write('GPS time\tTime since %d\tInclination\n\n' % grb_time)
#stacked = np.column_stack([np.concatenate((found_injs['time'], missed_injs['time'])),
#                           np.concatenate((found_injs['time'] - grb_time,
#                                           missed_injs['time'] - grb_time)),
#                           np.concatenate((found_injs['inclination'], missed_injs['inclination']))])
#np.savetxt(t_incl_txt, stacked, delimiter='\t')
#t_incl_txt.close()


# Save the efficiency values to disk
#efficiency_txt = open('%s/efficiency_numbers.txt' % outdir, 'w')
#efficiency_txt.write('distance (Mpc) \tfraction\tyerr_low\tyerr_high\n\n')
#stacked = np.column_stack([dist_plot_vals, fraction_mc,
#                           yerr_low_mc, yerr_high_mc])
#np.savetxt(efficiency_txt, stacked, fmt='%.8e', delimiter='\t')
#efficiency_txt.close()

# Print efficiency curve to file
#efficiency_curve_txt = open("%s/efficiency_curve.txt" % outdir, "w")
#efficiency_curve_txt.write('Distance (Mpc)\tEfficiency including counting errors\n\n')
#stacked = np.column_stack([dist_plot_vals, red_efficiency])
#np.savetxt(efficiency_curve_txt, stacked, delimiter='\t')
#efficiency_curve_txt.close()

# Useful while developping: exploit as unit test?
#np.savetxt('%s/found_maxbestnr.txt' % outdir, found_max_bestnr.T)
#np.savetxt('%s/found_maxbestnrnomc.txt' % outdir, found_max_bestnr_no_mc.T)
#np.savetxt('%s/foundonbestnr.txt' % outdir, found_on_bestnr.T)
#np.savetxt('%s/foundonbestnrnomc.txt' % outdir, found_on_bestnr_no_mc.T)
#np.savetxt('%s/numinjections.txt' % outdir, num_injections.T)
#np.savetxt('%s/numinjectionsnomc.txt' % outdir, num_injections_no_mc.T)
#inj_rec_file = open("%s/injection_recovery.html" % outdir, "w")
#inj_rec_file.write("Total injections found louder than all background using %s "\
#                   "is: %s<br>\n" % ('BestNR', found_max_bestnr[-1]))
#inj_rec_file.write("Total injections found louder than all background (and "\
#                   "nearby triggers in the offsource) using %s is: %s<br>"\
#                   % ('BestNR', found_on_bestnr[-1]))
#inj_rec_file.close()
