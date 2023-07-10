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

import os
import logging
import argparse
import copy
import numpy
import h5py
from scipy import stats
from pycbc.detector import Detector
# All/most of these final imports will become obsolete with hdf5 switch
try:
    from ligo import segments
    from ligo.lw import utils, lsctables
    from ligo.lw.table import Table
    from ligo.segments.utils import fromsegwizard
    # Handle MultiInspiral xml-tables with glue,
    # as ligo.lw no longer supports them
    from glue.ligolw import lsctables as glsctables
    # from glue.ligolw.ilwd import ilwdchar as gilwdchar
    from glue.ligolw.ligolw import LIGOLWContentHandler
except ImportError:
    pass


# =============================================================================
# Arguments functions:
# * Initialize a parser object with arguments shared by all plotting scripts
# * Add to the parser object the arguments used for Monte-Carlo on distance
# * Add to the parser object the arguments used for BestNR calculation
# * Add to the parser object the arguments for found/missed injection files
# =============================================================================
def pygrb_initialize_plot_parser(description=None, version=None):
    """Sets up a basic argument parser object for PyGRB plotting scripts"""

    formatter_class = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=formatter_class)
    parser.add_argument("--version", action="version", version=version)
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="Verbose output")
    parser.add_argument("-o", "--output-file", default=None,
                        help="Output file.")
    parser.add_argument("--x-lims", action="store", default=None,
                        help="Comma separated minimum and maximum values " +
                        "for the horizontal axis. When using negative " +
                        "values an equal sign after --x-lims is necessary.")
    parser.add_argument("--y-lims", action="store", default=None,
                        help="Comma separated minimum and maximum values " +
                        "for the vertical axis. When using negative values " +
                        "an equal sign after --y-lims is necessary.")
    parser.add_argument("--use-logs", default=False, action="store_true",
                        help="Produce a log-log plot")
    parser.add_argument("-i", "--ifo", default=None, help="IFO used for IFO " +
                        "specific plots")
    parser.add_argument("-a", "--seg-files", nargs="+", action="store",
                        default=None, help="The location of the buffer, " +
                        "onsource and offsource segment files.")
    parser.add_argument("-V", "--veto-files", nargs="+", action="store",
                        default=None, help="The location of the CATX veto " +
                        "files provided as a list of space-separated values.")
    parser.add_argument("-b", "--veto-category", action="store", type=int,
                        default=None, help="Apply vetoes up to this level " +
                        "inclusive.")
    parser.add_argument('--plot-title', default=None,
                        help="If provided, use the given string as the plot " +
                        "title.")
    parser.add_argument('--plot-caption', default=None,
                        help="If provided, use the given string as the plot " +
                        "caption")

    return parser


def pygrb_add_injmc_opts(parser):
    """Add to parser object the arguments used for Monte-Carlo on distance."""
    if parser is None:
        parser = argparse.ArgumentParser()
    parser.add_argument("-M", "--num-mc-injections", action="store",
                        type=int, default=100, help="Number of Monte " +
                        "Carlo injection simulations to perform.")
    parser.add_argument("-S", "--seed", action="store", type=int,
                        default=1234, help="Seed to initialize Monte Carlo.")
    parser.add_argument("-U", "--upper-inj-dist", action="store",
                        type=float, default=1000, help="The upper distance " +
                        "of the injections in Mpc, if used.")
    parser.add_argument("-L", "--lower-inj-dist", action="store",
                        type=float, default=0, help="The lower distance of " +
                        "the injections in Mpc, if used.")
    parser.add_argument("-n", "--num-bins", action="store", type=int,
                        default=0, help="The number of bins used to " +
                        "calculate injection efficiency.")
    parser.add_argument("-w", "--waveform-error", action="store",
                        type=float, default=0, help="The standard deviation " +
                        "to use when calculating the waveform error.")
    for ifo in ["g1", "h1", "k1", "l1", "v1"]:
        parser.add_argument(f"--{ifo}-cal-error", action="store", type=float,
                            default=0, help="The standard deviation to use " +
                            f"when calculating the {ifo.upper()} " +
                            "calibration amplitude error.")
        parser.add_argument(f"--{ifo}-dc-cal-error", action="store",
                            type=float, default=1.0, help="The scaling " +
                            "factor to use when calculating the " +
                            f"{ifo.upper()} calibration amplitude error.")


def pygrb_add_bestnr_opts(parser):
    """Add to the parser object the arguments used for BestNR calculation."""
    if parser is None:
        parser = argparse.ArgumentParser()
    parser.add_argument("-Q", "--chisq-index", action="store", type=float,
                        default=6.0, help="chisq_index for newSNR calculation")
    parser.add_argument("-N", "--chisq-nhigh", action="store", type=float,
                        default=2.0, help="nhigh for newSNR calculation")
    parser.add_argument("-B", "--sngl-snr-threshold", action="store",
                        type=float, default=4.0, help="Single detector SNR " +
                        "threshold, the two most sensitive detectors " +
                        "should have SNR above this.")
    parser.add_argument("-d", "--snr-threshold", action="store", type=float,
                        default=6.0, help="SNR threshold for recording " +
                        "triggers.")
    parser.add_argument("-c", "--newsnr-threshold", action="store", type=float,
                        default=6.0, help="NewSNR threshold for " +
                        "calculating the chisq of triggers (based on value " +
                        "of auto and bank chisq  values. By default will " +
                        "take the same value as snr-threshold.")
    parser.add_argument("-A", "--null-snr-threshold", action="store",
                        default="3.5,5.25",
                        help="Comma separated lower,higher null SNR " +
                        "threshold for null SNR cut")
    parser.add_argument("-T", "--null-grad-thresh", action="store", type=float,
                        default=20., help="Threshold above which to " +
                        "increase the values of the null SNR cut")
    parser.add_argument("-D", "--null-grad-val", action="store", type=float,
                        default=0.2, help="Rate the null SNR cut will " +
                        "increase above the threshold")


# =============================================================================
# Wrapper to read segments files
# =============================================================================
def _read_seg_files(seg_files):
    """Read segments txt files"""

    if len(seg_files) != 3:
        err_msg = "The location of three segment files is necessary."
        err_msg += "[bufferSeg.txt, offSourceSeg.txt, onSourceSeg.txt]"
        raise RuntimeError(err_msg)

    times = {}
    keys = ["buffer", "off", "on"]

    for key, seg_file in zip(keys, seg_files):
        segs = fromsegwizard(open(seg_file, 'r'))
        if len(segs) > 1:
            err_msg = 'More than one segment, an error has occured.'
            raise RuntimeError(err_msg)
        times[key] = segs[0]

    return times


# =============================================================================
# Function to load a table from an xml file
# =============================================================================
def load_xml_table(file_name, table_name):
    """Load xml table from file."""

    xml_doc = utils.load_filename(
        file_name,
        compress='auto',
        contenthandler=glsctables.use_in(LIGOLWContentHandler)
    )
    return Table.get_table(xml_doc, table_name)


# ==============================================================================
# Function to load segments from an xml file
# ==============================================================================
def _load_segments_from_xml(xml_doc, return_dict=False, select_id=None):
    """Read a ligo.segments.segmentlist from the file object file containing an
    xml segment table.

    Parameters
    ----------
        xml_doc: name of segment xml file

        Keyword Arguments:
            return_dict : [ True | False ]
                return a ligo.segments.segmentlistdict containing coalesced
                ligo.segments.segmentlists keyed by seg_def.name for each entry
                in the contained segment_def_table. Default False
            select_id : int
                return a ligo.segments.segmentlist object containing only
                those segments matching the given segment_def_id integer

    """

    # Load SegmentDefTable and SegmentTable
    seg_def_table = load_xml_table(xml_doc,
                                   glsctables.SegmentDefTable.tableName)
    seg_table = load_xml_table(xml_doc, glsctables.SegmentTable.tableName)

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
        if select_id and int(seg.segment_def_id) == select_id:
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
def _extract_vetoes(all_veto_files, ifos, veto_cat):
    """Extracts vetoes from veto filelist"""

    if all_veto_files and (veto_cat is None):
        err_msg = "Must supply veto category to apply vetoes."
        raise RuntimeError(err_msg)

    # Initialize veto containers
    vetoes = segments.segmentlistdict()
    for ifo in ifos:
        vetoes[ifo] = segments.segmentlist()

    veto_files = []
    veto_cats = range(2, veto_cat+1)
    for cat in veto_cats:
        veto_files += [vf for vf in all_veto_files if "CAT"+str(cat) in vf]
    n_found = len(veto_files)
    n_expected = len(ifos)*len(veto_cats)
    if n_found != n_expected:
        err_msg = f"Found {n_found} veto files instead of the expected "
        err_msg += f"{n_expected}; check the options."
        raise RuntimeError(err_msg)

    # Construct veto list from veto filelist
    if veto_files:
        for veto_file in veto_files:
            ifo = os.path.basename(veto_file)[:2]
            if ifo in ifos:
                # This returns a coalesced list of the vetoes
                tmp_veto_segs = _load_segments_from_xml(veto_file)
                for entry in tmp_veto_segs:
                    vetoes[ifo].append(entry)
    for ifo in ifos:
        vetoes[ifo].coalesce()

    return vetoes


# =============================================================================
# Function to get the ID numbers from a LIGO-LW table
# =============================================================================
def _get_id_numbers(ligolw_table, column):
    """Grab the IDs of a LIGO-LW table"""

    ids = [int(getattr(row, column)) for row in ligolw_table]

    return ids


# =============================================================================
# Function to build a dictionary (indexed by ifo) of time-slid vetoes
# =============================================================================
def _slide_vetoes(vetoes, slide_dict_or_list, slide_id):
    """Build a dictionary (indexed by ifo) of time-slid vetoes"""

    # Copy vetoes
    slid_vetoes = copy.deepcopy(vetoes)

    # Slide them
    ifos = vetoes.keys()
    for ifo in ifos:
        slid_vetoes[ifo].shift(-slide_dict_or_list[slide_id][ifo])

    return slid_vetoes


#
# Used (also) in executables
#

# =============================================================================
# Function to load triggers
# =============================================================================
def load_triggers(input_file, vetoes):
    """Loads triggers from PyGRB output file"""

    trigs = h5py.File(input_file, 'r')

    if vetoes is not None:
        # Developers: see PR 3972 for previous implementation
        raise NotImplementedError

    return trigs


# =============================================================================
# Detector utils:
# * Function to calculate the antenna response F+^2 + Fx^2
# * Function to calculate the antenna distance factor
# =============================================================================
def _get_antenna_single_response(antenna, ra, dec, geocent_time):
    """Returns the antenna response F+^2 + Fx^2 of an IFO (passed as pycbc
    Detector type) at a given sky location and time."""

    fp, fc = antenna.antenna_pattern(ra, dec, 0, geocent_time)

    return fp**2 + fc**2


# Vectorize the function above on all but the first argument
get_antenna_responses = numpy.vectorize(_get_antenna_single_response,
                                        otypes=[float])
get_antenna_responses.excluded.add(0)


def get_antenna_dist_factor(antenna, ra, dec, geocent_time, inc=0.0):
    """Returns the antenna factors (defined as eq. 4.3 on page 57 of
    Duncan Brown's Ph.D.) for an IFO (passed as pycbc Detector type) at
    a given sky location and time."""

    fp, fc = antenna.antenna_pattern(ra, dec, 0, geocent_time)

    return numpy.sqrt(fp ** 2 * (1 + numpy.cos(inc)) ** 2 / 4 + fc ** 2)


# =============================================================================
# Function to calculate the detection statistic of a list of triggers
# =============================================================================
def get_bestnrs(trigs, q=4.0, n=3.0, null_thresh=(4.25, 6), snr_threshold=6.,
                sngl_snr_threshold=4., chisq_threshold=None,
                null_grad_thresh=20., null_grad_val=0.2):
    """Calculate BestNR (coh_PTF detection statistic) of triggers through
    signal based vetoes.  The (default) signal based vetoes are:
    * Coherent SNR < 6
    * Bank chi-squared reduced (new) SNR < 6
    * Auto veto reduced (new) SNR < 6
    * Single-detector SNR (from two most sensitive IFOs) < 4
    * Null SNR (CoincSNR^2 - CohSNR^2)^(1/2) < null_thresh

    Returns Numpy array of BestNR values.
    """

    if not trigs:
        return numpy.array([])

    # Grab sky position and timing
    ra = trigs.get_column('ra')
    dec = trigs.get_column('dec')
    time = trigs.get_end()

    # Initialize BestNRs
    snr = trigs.get_column('snr')
    bestnr = numpy.ones(len(snr))

    # Coherent SNR cut
    bestnr[numpy.asarray(snr) < snr_threshold] = 0

    # Bank and auto chi-squared cuts
    if not chisq_threshold:
        chisq_threshold = snr_threshold
    for chisq in ['bank_chisq', 'cont_chisq']:
        bestnr[numpy.asarray(trigs.get_new_snr(index=q, nhigh=n,
                                               column=chisq))
               < chisq_threshold] = 0

    # Define IFOs for sngl cut
    ifos = list(map(str, trigs[0].get_ifos()))

    # Single detector SNR cut
    sens = {}
    sigmasqs = trigs.get_sigmasqs()
    ifo_snr = dict((ifo, trigs.get_sngl_snr(ifo)) for ifo in ifos)
    for ifo in ifos:
        antenna = Detector(ifo)
        sens[ifo] = sigmasqs[ifo] * get_antenna_responses(antenna, ra,
                                                          dec, time)
    # Apply this cut only if there is more than 1 IFO
    if len(ifos) > 1:
        for i_trig, _ in enumerate(trigs):
            # Apply only to triggers that were not already cut previously
            if bestnr[i_trig] != 0:
                ifos.sort(key=lambda ifo, j=i_trig: sens[ifo][j], reverse=True)
                if (ifo_snr[ifos[0]][i_trig] < sngl_snr_threshold or
                        ifo_snr[ifos[1]][i_trig] < sngl_snr_threshold):
                    bestnr[i_trig] = 0
    for i_trig, trig in enumerate(trigs):
        # Get chisq reduced (new) SNR for triggers that were not cut so far
        # NOTE: .get_bestnr is in glue.ligolw.lsctables.MultiInspiralTable
        if bestnr[i_trig] != 0:
            bestnr[i_trig] = trig.get_bestnr(index=q, nhigh=n,
                                             null_snr_threshold=null_thresh[0],
                                             null_grad_thresh=null_grad_thresh,
                                             null_grad_val=null_grad_val)
            # If we got this far and the bestNR is non-zero, verify that chisq
            # was actually calculated for the trigger, otherwise raise an
            # error with info useful to figure out why this happened.
            if bestnr[i_trig] != 0 and trig.chisq == 0:
                err_msg = "Chisq not calculated for trigger with end time "
                err_msg += f"{trig.get_end()} and SNR {trig.snr}."
                raise RuntimeError(err_msg)

    return bestnr


# =============================================================================
# Construct sorted triggers from trials
# =============================================================================
def sort_trigs(trial_dict, trigs, slide_dict, seg_dict):
    """Constructs sorted triggers from a trials dictionary"""

    sorted_trigs = {}

    # Begin by sorting the triggers into each slide
    # New seems pretty slow, so run it once and then use deepcopy
    tmp_table = glsctables.New(glsctables.MultiInspiralTable)
    for slide_id in slide_dict:
        sorted_trigs[slide_id] = copy.deepcopy(tmp_table)
    for trig in trigs:
        sorted_trigs[int(trig.time_slide_id)].append(trig)

    for slide_id in slide_dict:
        # These can only *reduce* the analysis time
        curr_seg_list = seg_dict[slide_id]

        # Check the triggers are all in the analysed segment lists
        for trig in sorted_trigs[slide_id]:
            if trig.end_time not in curr_seg_list:
                # This can be raised if the trigger is on the segment boundary,
                # so check if the trigger is within 1/100 of a second within
                # the list
                if trig.get_end() + 0.01 in curr_seg_list:
                    continue
                if trig.get_end() - 0.01 in curr_seg_list:
                    continue
                err_msg = "Triggers found in input files not in the list of "
                err_msg += "analysed segments. This should not happen."
                raise RuntimeError(err_msg)
        # END OF CHECK #

        # The below line works like the inverse of .veto and only returns trigs
        # that are within the segment specified by trial_dict[slide_id]
        sorted_trigs[slide_id] = \
            sorted_trigs[slide_id].vetoed(trial_dict[slide_id])

    return sorted_trigs


# =============================================================================
# Extract basic trigger properties and store them as dictionaries
# =============================================================================
def extract_basic_trig_properties(trial_dict, trigs, slide_dict, seg_dict,
                                  opts):
    """Extract and store as dictionaries time, SNR, and BestNR of
    time-slid triggers"""

    # Sort the triggers into each slide
    sorted_trigs = sort_trigs(trial_dict, trigs, slide_dict, seg_dict)
    logging.info("Triggers sorted.")

    # Local copies of variables entering the BestNR definition
    chisq_index = opts.chisq_index
    chisq_nhigh = opts.chisq_nhigh
    null_thresh = list(map(float, opts.null_snr_threshold.split(',')))
    snr_thresh = opts.snr_threshold
    sngl_snr_thresh = opts.sngl_snr_threshold
    new_snr_thresh = opts.newsnr_threshold
    null_grad_thresh = opts.null_grad_thresh
    null_grad_val = opts.null_grad_val

    # Build the 3 dictionaries
    trig_time = {}
    trig_snr = {}
    trig_bestnr = {}
    for slide_id in slide_dict:
        slide_trigs = sorted_trigs[slide_id]
        if slide_trigs:
            trig_time[slide_id] = numpy.asarray(slide_trigs.get_end()).\
                                  astype(float)
            trig_snr[slide_id] = numpy.asarray(slide_trigs.get_column('snr'))
        else:
            trig_time[slide_id] = numpy.asarray([])
            trig_snr[slide_id] = numpy.asarray([])
        trig_bestnr[slide_id] = get_bestnrs(slide_trigs,
                                            q=chisq_index,
                                            n=chisq_nhigh,
                                            null_thresh=null_thresh,
                                            snr_threshold=snr_thresh,
                                            sngl_snr_threshold=sngl_snr_thresh,
                                            chisq_threshold=new_snr_thresh,
                                            null_grad_thresh=null_grad_thresh,
                                            null_grad_val=null_grad_val)
    logging.info("Time, SNR, and BestNR of triggers extracted.")

    return trig_time, trig_snr, trig_bestnr


# =============================================================================
# Function to extract ifos from hdfs
# =============================================================================
def extract_ifos(trig_file):
    """Extracts IFOs from hdf file"""

    # Load hdf file
    hdf_file = h5py.File(trig_file, 'r')

    # Extract IFOs
    ifos = sorted(list(hdf_file.keys()))

    # Remove 'network' key from list of ifos
    if 'network' in ifos:
        ifos.remove('network')

    return ifos


# =============================================================================
# Function to extract IFOs and vetoes
# =============================================================================
def extract_ifos_and_vetoes(trig_file, veto_files, veto_cat):
    """Extracts IFOs from HDF files and vetoes from a directory"""

    logging.info("Extracting IFOs and vetoes.")

    # Extract IFOs
    ifos = extract_ifos(trig_file)

    # Extract vetoes
    if veto_files is not None:
        vetoes = _extract_vetoes(veto_files, ifos, veto_cat)
    else:
        vetoes = None

    return ifos, vetoes


# =============================================================================
# Function to load injections
# =============================================================================
def load_injections(inj_file, vetoes, sim_table=False, label=None):
    """Loads injections from PyGRB output file"""

    if label is None:
        logging.info("Loading injections...")
    else:
        logging.info("Loading %s...", label)

    insp_table = glsctables.MultiInspiralTable
    if sim_table:
        insp_table = glsctables.SimInspiralTable

    # Load injections in injection file
    inj_table = load_xml_table(inj_file, insp_table.tableName)

    # Extract injections in time-slid non-vetoed data
    injs = lsctables.New(insp_table, columns=insp_table.loadcolumns)
    injs.extend(inj for inj in inj_table if inj.get_end() not in vetoes)

    if label is None:
        logging.info("%d injections found.", len(injs))
    else:
        logging.info("%d %s found.", len(injs), label)

    return injs


# =============================================================================
# Function to load timeslides
# =============================================================================
def load_time_slides(xml_file):
    """Loads timeslides from PyGRB output file as a dictionary"""

    # Get all timeslides: these are number_of_ifos * number_of_timeslides
    time_slide = load_xml_table(xml_file, glsctables.TimeSlideTable.tableName)
    # Get a list of unique timeslide dictionaries
    time_slide_list = [dict(i) for i in time_slide.as_dict().values()]
    # Turn it into a dictionary indexed by the timeslide ID
    time_slide_dict = {int(time_slide.get_time_slide_id(ov)): ov
                       for ov in time_slide_list}
    # Check time_slide_ids are ordered correctly.
    ids = _get_id_numbers(time_slide,
                          "time_slide_id")[::len(time_slide_dict[0].keys())]
    if not (numpy.all(ids[1:] == numpy.array(ids[:-1])+1) and ids[0] == 0):
        err_msg = "time_slide_ids list should start at zero and increase by "
        err_msg += "one for every element"
        raise RuntimeError(err_msg)
    # Check that the zero-lag slide has time_slide_id == 0.
    if not numpy.all(numpy.array(list(time_slide_dict[0].values())) == 0):
        err_msg = "The slide with time_slide_id == 0 should be the "
        err_msg += "zero-lag-slide but it has non-zero slide values: "
        err_msg += f"{time_slide_dict[0]}."
        raise RuntimeError(err_msg)

    return time_slide_dict


# =============================================================================
# Function to load the segment dicitonary
# =============================================================================
def load_segment_dict(xml_file):
    """Loads the segment dictionary """

    # Get the mapping table
    time_slide_map_table = \
        load_xml_table(xml_file, glsctables.TimeSlideSegmentMapTable.tableName)
    # Perhaps unnecessary as segment_def_id and time_slide_id seem to always
    # be identical identical
    segment_map = {
        int(entry.segment_def_id): int(entry.time_slide_id)
        for entry in time_slide_map_table
    }
    # Extract the segment table
    segment_table = load_xml_table(
        xml_file, glsctables.SegmentTable.tableName)
    segment_dict = {}
    for entry in segment_table:
        curr_slid_id = segment_map[int(entry.segment_def_id)]
        curr_seg = entry.get()
        if curr_slid_id not in segment_dict:
            segment_dict[curr_slid_id] = segments.segmentlist()
        segment_dict[curr_slid_id].append(curr_seg)
        segment_dict[curr_slid_id].coalesce()

    return segment_dict


# =============================================================================
# Construct the trials from the timeslides, segments, and vetoes
# =============================================================================
def construct_trials(seg_files, seg_dict, ifos, slide_dict, vetoes):
    """Constructs trials from triggers, timeslides, segments and vetoes"""

    trial_dict = {}

    # Get segments
    segs = _read_seg_files(seg_files)

    # Separate segments
    trial_time = abs(segs['on'])

    for slide_id in slide_dict:
        # These can only *reduce* the analysis time
        curr_seg_list = seg_dict[slide_id]

        # Construct the buffer segment list
        seg_buffer = segments.segmentlist()
        for ifo in ifos:
            slide_offset = slide_dict[slide_id][ifo]
            seg_buffer.append(segments.segment(segs['buffer'][0] -
                                               slide_offset,
                                               segs['buffer'][1] -
                                               slide_offset))
        seg_buffer.coalesce()

        # Construct the ifo-indexed dictionary of slid veteoes
        slid_vetoes = _slide_vetoes(vetoes, slide_dict, slide_id)

        # Construct trial list and check against buffer
        trial_dict[slide_id] = segments.segmentlist()
        for curr_seg in curr_seg_list:
            iter_int = 1
            while 1:
                trial_end = curr_seg[0] + trial_time*iter_int
                if trial_end > curr_seg[1]:
                    break
                curr_trial = segments.segment(trial_end - trial_time,
                                              trial_end)
                if not seg_buffer.intersects_segment(curr_trial):
                    intersect = numpy.any([slid_vetoes[ifo].
                                           intersects_segment(curr_trial)
                                           for ifo in ifos])
                    if not intersect:
                        trial_dict[slide_id].append(curr_trial)

                iter_int += 1

    return trial_dict


# =============================================================================
# Find max and median of loudest SNRs or BestNRs
# =============================================================================
def sort_stat(time_veto_max_stat):
    """Sort a dictionary of loudest SNRs/BestNRs"""

    full_time_veto_max_stat = list(time_veto_max_stat.values())
    full_time_veto_max_stat = numpy.concatenate(full_time_veto_max_stat)
    full_time_veto_max_stat.sort()

    return full_time_veto_max_stat


# =============================================================================
# Find max and median of loudest SNRs or BestNRs
# =============================================================================
def max_median_stat(slide_dict, time_veto_max_stat, trig_stat, total_trials):
    """Deterine the maximum and median of the loudest SNRs/BestNRs"""

    max_stat = max([trig_stat[slide_id].max() if trig_stat[slide_id].size
                   else 0 for slide_id in slide_dict])

    full_time_veto_max_stat = sort_stat(time_veto_max_stat)

    if total_trials % 2:
        median_stat = full_time_veto_max_stat[(total_trials - 1) // 2]
    else:
        median_stat = numpy.mean((full_time_veto_max_stat)
                                 [total_trials//2 - 1: total_trials//2 + 1])

    return max_stat, median_stat, full_time_veto_max_stat


# =============================================================================
# Function to determine calibration and waveform errors for injection sets
# =============================================================================
def mc_cal_wf_errs(num_mc_injs, inj_dists, cal_err, wf_err, max_dc_cal_err):
    """Includes calibration and waveform errors by running an MC"""

    # The efficiency calculations include calibration and waveform
    # errors incorporated by running over each injection num_mc_injs times,
    # where each time we draw a random value of distance.

    num_injs = len(inj_dists)

    inj_dist_mc = numpy.ndarray((num_mc_injs+1, num_injs))
    inj_dist_mc[0, :] = inj_dists
    for i in range(num_mc_injs):
        cal_dist_red = stats.norm.rvs(size=num_injs) * cal_err
        wf_dist_red = numpy.abs(stats.norm.rvs(size=num_injs) * wf_err)
        inj_dist_mc[i+1, :] = inj_dists / (max_dc_cal_err *
                                           (1 + cal_dist_red) *
                                           (1 + wf_dist_red))

    return inj_dist_mc


# =============================================================================
# Function to calculate the coincident SNR
# =============================================================================
def get_coinc_snr(trigs_or_injs, ifos):
    """ Calculate coincident SNR using single IFO SNRs"""

    num_trigs_or_injs = len(trigs_or_injs['network/end_time_gc'][:])

    # Calculate coincident SNR
    single_snr_sq = dict((ifo, None) for ifo in ifos)
    snr_sum_square = numpy.zeros(num_trigs_or_injs)
    for ifo in ifos:
        key = ifo + '/snr_' + ifo.lower()
        if ifo.lower() != 'h1':
            key = key[:-1]
        # Square the individual SNRs
        single_snr_sq[ifo] = numpy.square(
            trigs_or_injs[key][:])
        # Add them
        snr_sum_square = numpy.add(snr_sum_square, single_snr_sq[ifo])
    # Obtain the square root
    coinc_snr = numpy.sqrt(snr_sum_square)

    return coinc_snr
