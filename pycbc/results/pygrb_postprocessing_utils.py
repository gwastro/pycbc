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
import ligo.segments as segments
from pycbc.events.coherent import reweightedsnr_cut
# All/most of these final imports will become obsolete with hdf5 switch
try:
    from ligo.lw import utils
    from ligo.lw.table import Table
    from ligo.segments.utils import fromsegwizard
    from glue.ligolw import lsctables as glsctables
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
    """Add to the parser object the arguments used for BestNR calculation"""
    if parser is None:
        parser = argparse.ArgumentParser()
    parser.add_argument("-Q", "--chisq-index", action="store", type=float,
                        default=6.0, help="chisq_index for newSNR " +
                        "calculation (default: 6)")
    parser.add_argument("-N", "--chisq-nhigh", action="store", type=float,
                        default=2.0, help="chisq_nhigh for newSNR " +
                        "calculation (default: 2")


def pygrb_add_null_snr_opts(parser):
    """Add to the parser object the arguments used for null SNR calculation
    and null SNR cut."""
    parser.add_argument("-A", "--null-snr-threshold", action="store",
                        default=5.25,
                        type=float,
                        help="Null SNR threshold for null SNR cut "
                        "(default: 5.25)")
    parser.add_argument("-T", "--null-grad-thresh", action="store", type=float,
                        default=20., help="Threshold above which to " +
                        "increase the values of the null SNR cut")
    parser.add_argument("-D", "--null-grad-val", action="store", type=float,
                        default=0.2, help="Rate the null SNR cut will " +
                        "increase above the threshold")


def pygrb_add_single_snr_cut_opt(parser):
    """Add to the parser object an argument to place a threshold on single
    detector SNR."""
    parser.add_argument("-B", "--sngl-snr-threshold", action="store",
                        type=float, default=4.0, help="Single detector SNR " +
                        "threshold, the two most sensitive detectors " +
                        "should have SNR above this.")


def pygrb_add_bestnr_cut_opt(parser):
    """Add to the parser object an argument to place a threshold on BestNR."""
    if parser is None:
        parser = argparse.ArgumentParser()
    parser.add_argument("--newsnr-threshold", type=float, metavar='THRESHOLD',
                        default=0.,
                        help="Cut triggers with NewSNR less than THRESHOLD" +
                        "Default 0: all events are considered.")


# =============================================================================
# Wrapper to read segments files
# =============================================================================
def _read_seg_files(seg_files):
    """Read segments txt files"""

    if len(seg_files) != 3 or seg_files is None:
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


# =============================================================================
# Function to load segments from an xml file
# =============================================================================
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
def _slide_vetoes(vetoes, slide_dict_or_list, slide_id, ifos):
    """Build a dictionary (indexed by ifo) of time-slid vetoes"""

    # Copy vetoes
    if vetoes is not None:
        slid_vetoes = copy.deepcopy(vetoes)
        # Slide them
        for ifo in ifos:
            slid_vetoes[ifo].shift(-slide_dict_or_list[slide_id][ifo])
    else:
        slid_vetoes = {ifo: segments.segmentlist() for ifo in ifos}

    return slid_vetoes


#
# Used (also) in executables
#

# =============================================================================
# Functions to load triggers
# =============================================================================
def dataset_iterator(g, prefix=''):
    """Reach all datasets in and HDF file"""

    for key, item in g.items():
        # Avoid slash as first character
        path = prefix[1:] + '/' + key
        if isinstance(item, h5py.Dataset):
            yield (path, item)
        elif isinstance(item, h5py.Group):
            yield from dataset_iterator(item, path)


def load_triggers(input_file, ifos, vetoes, rw_snr_threshold=None):
    """Loads triggers from PyGRB output file, returning a dictionary"""

    trigs = h5py.File(input_file, 'r')
    rw_snr = trigs['network/reweighted_snr'][:]
    net_ids = trigs['network/event_id'][:]
    ifo_ids = {}
    for ifo in ifos:
        ifo_ids[ifo] = trigs[ifo+'/event_id'][:]
    trigs.close()

    if vetoes is not None:
        # Developers: see PR 3972 for previous implementation
        raise NotImplementedError

    # Apply the reweighted SNR cut on the reweighted SNR
    if rw_snr_threshold is not None:
        rw_snr = reweightedsnr_cut(rw_snr, rw_snr_threshold)

    # Establish the indices of data not surviving the cut
    above_thresh = rw_snr > 0
    num_orig_pts = len(above_thresh)

    # Do not assume that IFO and network datasets are sorted the same way:
    # find where each surviving network/event_id is placed in the IFO/event_id
    ifo_ids_above_thresh_locations = {}
    for ifo in ifos:
        ifo_ids_above_thresh_locations[ifo] = \
            numpy.array([numpy.where(ifo_ids[ifo] == net_id)[0][0]
                         for net_id in net_ids[above_thresh]])

    # Apply the cut on all the data by remove points with reweighted SNR = 0
    trigs_dict = {}
    with h5py.File(input_file, "r") as trigs:
        for (path, dset) in dataset_iterator(trigs):
            # The dataset contains information other than trig/inj properties:
            # just copy it
            if len(dset) != num_orig_pts:
                trigs_dict[path] = dset[:]
            # The dataset is relative to an IFO: cut with the correct index
            elif path[:2] in ifos:
                ifo = path[:2]
                if ifo_ids_above_thresh_locations[ifo].size != 0:
                    trigs_dict[path] = \
                        dset[:][ifo_ids_above_thresh_locations[ifo]]
                else:
                    trigs_dict[path] = numpy.array([])
            # The dataset is relative to the network: cut it before copying
            else:
                trigs_dict[path] = dset[above_thresh]

    return trigs_dict


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
# Construct sorted triggers from trials
# =============================================================================
def sort_trigs(trial_dict, trigs, slide_dict, seg_dict):
    """Constructs sorted triggers from a trials dictionary"""

    sorted_trigs = {}

    # Begin by sorting the triggers into each slide
    for slide_id in slide_dict:
        sorted_trigs[slide_id] = []
    for slide_id, event_id in zip(trigs['network/slide_id'],
                                  trigs['network/event_id']):
        sorted_trigs[slide_id].append(event_id)

    for slide_id in slide_dict:
        # These can only *reduce* the analysis time
        curr_seg_list = seg_dict[slide_id]

        # Check the triggers are all in the analysed segment lists
        for event_id in sorted_trigs[slide_id]:
            index = numpy.flatnonzero(trigs['network/event_id'] == event_id)[0]
            end_time = trigs['network/end_time_gc'][index]
            if end_time not in curr_seg_list:
                # This can be raised if the trigger is on the segment boundary,
                # so check if the trigger is within 1/100 of a second within
                # the list
                if end_time + 0.01 in curr_seg_list:
                    continue
                if end_time - 0.01 in curr_seg_list:
                    continue
                err_msg = "Triggers found in input files not in the list of "
                err_msg += "analysed segments. This should not happen."
                raise RuntimeError(err_msg)
        # END OF CHECK #

        # Keep triggers that are in trial_dict
        sorted_trigs[slide_id] = [event_id for event_id in
                                  sorted_trigs[slide_id]
                                  if trigs['network/end_time_gc'][
                                      trigs['network/event_id'] == event_id][0]
                                  in trial_dict[slide_id]]

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

    # Build the 3 dictionaries
    trig_time = {}
    trig_snr = {}
    trig_bestnr = {}
    for slide_id in slide_dict:
        slide_trigs = sorted_trigs[slide_id]
        indices = numpy.nonzero(
            numpy.isin(trigs['network/event_id'], slide_trigs))[0]
        if slide_trigs:
            trig_time[slide_id] = trigs['network/end_time_gc'][
                indices]
            trig_snr[slide_id] = trigs['network/coherent_snr'][
                indices]
        else:
            trig_time[slide_id] = numpy.asarray([])
            trig_snr[slide_id] = numpy.asarray([])
        trig_bestnr[slide_id] = reweightedsnr_cut(
            trigs['network/reweighted_snr'][indices],
            opts.newsnr_threshold)

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
# Function to load timeslides
# =============================================================================
def load_time_slides(hdf_file_path):
    """Loads timeslides from PyGRB output file as a dictionary"""
    hdf_file = h5py.File(hdf_file_path, 'r')
    ifos = extract_ifos(hdf_file_path)
    ids = numpy.arange(len(hdf_file[f'{ifos[0]}/search/time_slides']))
    time_slide_dict = {
        slide_id: {
            ifo: hdf_file[f'{ifo}/search/time_slides'][slide_id]
            for ifo in ifos}
        for slide_id in ids}

    # Check time_slide_ids are ordered correctly.
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
def load_segment_dict(hdf_file_path):
    """
    Loads the segment dictionary with the format
    {slide_id: segmentlist(segments analyzed)}
    """

    # Long time slides will require mapping between slides and segments
    hdf_file = h5py.File(hdf_file_path, 'r')
    ifos = extract_ifos(hdf_file_path)
    # Get slide IDs
    slide_ids = numpy.arange(len(hdf_file[f'{ifos[0]}/search/time_slides']))
    # Get segment start/end times
    seg_starts = hdf_file['network/search/segments/start_times'][:]
    seg_ends = hdf_file['network/search/segments/end_times'][:]
    # Write list of segments
    seg_list = segments.segmentlist([segments.segment(seg_start, seg_ends[i])
                                    for i, seg_start in enumerate(seg_starts)])

    # Write segment_dict in proper format
    # At the moment of this comment, there is only one segment
    segment_dict = {slide: seg_list.coalesce() for slide in slide_ids}

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
        slid_vetoes = _slide_vetoes(vetoes, slide_dict, slide_id, ifos)

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
def get_coinc_snr(trigs_or_injs):
    """ Calculate coincident SNR using coherent and null SNRs"""

    coh_snr_sq = numpy.square(trigs_or_injs['network/coherent_snr'][:])
    null_snr_sq = numpy.square(trigs_or_injs['network/null_snr'][:])
    coinc_snr = numpy.sqrt(coh_snr_sq + null_snr_sq)

    return coinc_snr


def template_hash_to_id(trigger_file, bank_path):
    """
    This function converts the template hashes from a trigger file
    into 'template_id's that represent indices of the
    templates within the bank.
    Parameters
    ----------
    trigger_file: h5py File object for trigger file
    bank_file: filepath for template bank
    """
    with h5py.File(bank_path, "r") as bank:
        hashes = bank['template_hash'][:]
    ifos = [k for k in trigger_file.keys() if k != 'network']
    trig_hashes = trigger_file[f'{ifos[0]}/template_hash'][:]
    trig_ids = numpy.zeros(trig_hashes.shape[0], dtype=int)
    for idx, t_hash in enumerate(hashes):
        matches = numpy.where(trig_hashes == t_hash)
        trig_ids[matches] = idx
    return trig_ids
