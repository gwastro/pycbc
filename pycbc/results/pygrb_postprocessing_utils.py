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

import logging
import argparse
import copy
import numpy
import h5py
import igwn_segments as segments

from scipy import stats
from igwn_segments.utils import fromsegwizard
from pycbc.events.coherent import reweightedsnr_cut
from pycbc.events import veto
from pycbc import add_common_pycbc_options
from pycbc.io.hdf import HFile

logger = logging.getLogger('pycbc.results.pygrb_postprocessing_utils')


# =============================================================================
# Arguments functions:
# * Initialize a parser object with arguments shared by all plotting scripts
# * Add to the parser object the arguments used for Monte-Carlo on distance
# * Add to the parser object the arguments used for BestNR calculation
# * Add to the parser object the arguments for found/missed injection files
# =============================================================================
def pygrb_initialize_plot_parser(description=None):
    """Sets up a basic argument parser object for PyGRB plotting scripts"""

    formatter_class = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=formatter_class)
    add_common_pycbc_options(parser)
    parser.add_argument("-o", "--output-file", default=None,
                        help="Output file.")
    parser.add_argument("--x-lims", default=None,
                        help="Comma separated minimum and maximum values "
                        "for the horizontal axis. When using negative "
                        "values an equal sign after --x-lims is necessary.")
    parser.add_argument("--y-lims", default=None,
                        help="Comma separated minimum and maximum values "
                        "for the vertical axis. When using negative values "
                        "an equal sign after --y-lims is necessary.")
    parser.add_argument("--use-logs", default=False, action="store_true",
                        help="Produce a log-log plot")
    parser.add_argument("-i", "--ifo", default=None, help="IFO used for IFO "
                        "specific plots")
    parser.add_argument("-a", "--seg-files", nargs="+",
                        default=[], help="The location of the buffer, "
                        "onsource and offsource txt segment files.")
    parser.add_argument("-V", "--veto-file",
                        help="The location of the xml veto file.")
    parser.add_argument('--plot-title', default=None,
                        help="If provided, use the given string as the plot "
                        "title.")
    parser.add_argument('--plot-caption', default=None,
                        help="If provided, use the given string as the plot "
                        "caption")
    return parser


def pygrb_add_slide_opts(parser):
    """Add to parser object arguments related to short timeslides"""
    parser.add_argument("--slide-id", type=str, default='0',
                        help="Select a specific slide or set to all to plot "
                        "results from all short slides.")


def slide_opts_helper(args):
    """
       This function overwrites the types of input slide_id information
       when loading data in postprocessing scripts.
    """
    if args.slide_id.isdigit():
        args.slide_id = int(args.slide_id)
    elif args.slide_id.lower() == "all":
        args.slide_id = None
    else:
        raise ValueError("--slide-id must be the string all or an int")


def pygrb_add_injmc_opts(parser):
    """Add to parser object the arguments used for Monte-Carlo on distance."""
    if parser is None:
        parser = argparse.ArgumentParser()
    parser.add_argument("-M", "--num-mc-injections",
                        type=int, default=100, help="Number of Monte "
                        "Carlo injection simulations to perform.")
    parser.add_argument("-S", "--seed", type=int,
                        default=1234, help="Seed to initialize Monte Carlo.")
    parser.add_argument("-U", "--upper-inj-dist",
                        type=float, default=1000, help="The upper distance "
                        "of the injections in Mpc, if used.")
    parser.add_argument("-L", "--lower-inj-dist",
                        type=float, default=0, help="The lower distance of "
                        "the injections in Mpc, if used.")
    parser.add_argument("-n", "--num-bins", type=int,
                        default=0, help="The number of bins used to "
                        "calculate injection efficiency.")
    parser.add_argument("-w", "--waveform-error",
                        type=float, default=0, help="The standard deviation "
                        "to use when calculating the waveform error.")
    for ifo in ["g1", "h1", "k1", "l1", "v1"]:
        parser.add_argument(f"--{ifo}-cal-error", type=float,
                            default=0, help="The standard deviation to use "
                            f"when calculating the {ifo.upper()} "
                            "calibration amplitude error.")
        parser.add_argument(f"--{ifo}-dc-cal-error",
                            type=float, default=1.0, help="The scaling "
                            "factor to use when calculating the "
                            f"{ifo.upper()} calibration amplitude error.")


def pygrb_add_bestnr_opts(parser):
    """Add to the parser object the arguments used for BestNR calculation"""
    if parser is None:
        parser = argparse.ArgumentParser()
    parser.add_argument("-Q", "--chisq-index", type=float,
                        default=6.0, help="chisq_index for newSNR "
                        "calculation (default: 6)")
    parser.add_argument("-N", "--chisq-nhigh", type=float,
                        default=2.0, help="chisq_nhigh for newSNR "
                        "calculation (default: 2")


def pygrb_add_null_snr_opts(parser):
    """Add to the parser object the arguments used for null SNR calculation
    and null SNR cut."""
    parser.add_argument("-A", "--null-snr-threshold",
                        default=5.25,
                        type=float,
                        help="Null SNR threshold for null SNR cut "
                        "(default: 5.25)")
    parser.add_argument("-T", "--null-grad-thresh", type=float,
                        default=20., help="Threshold above which to "
                        "increase the values of the null SNR cut")
    parser.add_argument("-D", "--null-grad-val", type=float,
                        default=0.2, help="Rate the null SNR cut will "
                        "increase above the threshold")


def pygrb_add_bestnr_cut_opt(parser):
    """Add to the parser object an argument to place a threshold on BestNR."""
    if parser is None:
        parser = argparse.ArgumentParser()
    parser.add_argument("--newsnr-threshold", type=float, metavar='THRESHOLD',
                        default=0.,
                        help="Cut triggers with NewSNR less than THRESHOLD. "
                        "Default 0: all events are considered.")


# =============================================================================
# Wrapper to pick triggers with a given slide_id
# =============================================================================
# Underscore starts name of functions not called outside this file
def _slide_filter(trig_file, data, slide_id=None):
    """
    This function adds the capability to select triggers with specific
    slide_ids during the postprocessing stage of PyGRB.
    """
    if slide_id is None:
        return data
    mask = numpy.where(trig_file['network/slide_id'][:] == slide_id)[0]
    return data[mask]


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
    # Needs to be in this order for consistency with build_segment_filelist
    keys = ["buffer", "off", "on"]

    for key, seg_file in zip(keys, seg_files):
        segs = fromsegwizard(open(seg_file, 'r'))
        if len(segs) > 1:
            err_msg = 'More than one segment, an error has occured.'
            raise RuntimeError(err_msg)
        times[key] = segs[0]

    return times


# =============================================================================
# Function to extract vetoes
# =============================================================================
def _extract_vetoes(veto_file, ifos, offsource):
    """Extracts the veto segments from the veto File"""

    clean_segs = {}
    vetoes = segments.segmentlistdict()

    if veto_file:
        for ifo in ifos:
            segs = veto.select_segments_by_definer(veto_file, ifo=ifo)
            segs.coalesce()
            if len(segs) > 0:
                clean_segs[ifo] = segs
            else:
                clean_segs[ifo] = segments.segmentlist([offsource])

    if clean_segs:
        for ifo in ifos:
            vetoes[ifo] = segments.segmentlist([offsource]) - clean_segs[ifo]
        vetoes.coalesce()
        for ifo in ifos:
            for v in vetoes[ifo]:
                v_span = v[1] - v[0]
                logging.info("%ds of data vetoed at GPS time %d",
                             v_span, v[0])

    return vetoes


# =============================================================================
# Function to build a dictionary (indexed by ifo) of time-slid vetoes
# =============================================================================
def _slide_vetoes(vetoes, slide_dict_or_list, slide_id, ifos):
    """Build a dictionary (indexed by ifo) of time-slid vetoes"""

    # Copy vetoes
    if vetoes:
        slid_vetoes = copy.deepcopy(vetoes)
        # Slide them
        for ifo in ifos:
            slid_vetoes[ifo].shift(-slide_dict_or_list[slide_id][ifo])
    else:
        slid_vetoes = {ifo: segments.segmentlist() for ifo in ifos}

    return slid_vetoes


# =============================================================================
# Recursive function to reach all datasets in an HDF file handle
# =============================================================================
def _dataset_iterator(g, prefix=''):
    """Reach all datasets in an HDF file handle"""

    for key, item in g.items():
        # Avoid slash as first character
        pref = prefix[1:] if prefix.startswith('/') else prefix
        path = pref + '/' + key
        if isinstance(item, h5py.Dataset):
            yield (path, item)
        elif isinstance(item, h5py.Group):
            yield from _dataset_iterator(item, path)


# =============================================================================
# Function to load trigger/injection data
# =============================================================================
def load_data(input_file, ifos, rw_snr_threshold=None, data_tag=None,
              slide_id=None):
    """Load data from a trigger/injection PyGRB output file, returning a
    dictionary. If the input_file is None, None is returned. data_tag enables
    logging information about the number of triggers/injections found, so the
    user should not set it to 'trigs'/'injs' when processing the onsource."""

    if not input_file:
        return None

    trigs = HFile(input_file, 'r')
    rw_snr = trigs['network/reweighted_snr'][:] \
        if 'network/reweighted_snr' in trigs.keys() else numpy.array([])
    net_ids = trigs['network/event_id'][:] \
        if 'network/event_id' in trigs.keys() \
        else numpy.array([], dtype=numpy.int64)

    # Output the number of items loaded only upon a request by the user who is
    # expected not to set data_tag to 'trigs'or 'injs' when processing the
    # onsource
    if data_tag == 'trigs':
        logging.info("%d triggers loaded.", len(rw_snr))
    elif data_tag == 'injs':
        logging.info("%d injections loaded.", len(rw_snr))
    else:
        logging.info("Loading triggers.")
    ifo_ids = {}
    for ifo in ifos:
        ifo_ids[ifo] = trigs[ifo+'/event_id'][:] \
            if ifo+'/event_id' in trigs.keys() \
            else numpy.array([], dtype=numpy.int64)
    trigs.close()

    # Apply the reweighted SNR cut on the reweighted SNR
    if rw_snr_threshold is not None:
        rw_snr = reweightedsnr_cut(rw_snr, rw_snr_threshold)

    # Establish the indices of data not surviving the cut
    above_thresh = rw_snr > 0

    # Output the number of items surviging vetoes with the same logic as above
    msg = ""
    if data_tag == 'trigs':
        msg += f"{sum(above_thresh)} triggers "
    elif data_tag == 'injs':
        msg = f"{sum(above_thresh)} injections "
    if msg:
        msg += f"surviving reweighted SNR cut at {rw_snr_threshold}."
        logging.info(msg)

    # Do not assume that IFO and network datasets are sorted the same way:
    # find where each surviving network/event_id is placed in the IFO/event_id
    ifo_ids_above_thresh_locations = {}
    for ifo in ifos:
        ifo_ids_above_thresh_locations[ifo] = \
            numpy.array([numpy.where(ifo_ids[ifo] == net_id)[0][0]
                         for net_id in net_ids[above_thresh]])

    # Apply the cut on all the data by removing points with reweighted SNR = 0
    trigs_dict = {}
    with HFile(input_file, "r") as trigs:
        for (path, dset) in _dataset_iterator(trigs):
            # The dataset contains search information or missed injections
            # information, not properties of triggers or found injections:
            # just copy it
            if 'search' in path or 'missed' in path or 'gating' in path:
                trigs_dict[path] = dset[:]
            # The dataset is trig/inj info at an IFO:
            # cut with the correct index
            elif path[:2] in ifos:
                ifo = path[:2]
                if ifo_ids_above_thresh_locations[ifo].size != 0:
                    trigs_dict[path] = \
                        dset[:][ifo_ids_above_thresh_locations[ifo]]
                else:
                    trigs_dict[path] = numpy.array([])
            # The dataset is trig/inj network info: cut it before copying
            else:
                trigs_dict[path] = dset[above_thresh]

            if 'network/slide_id' in trigs.keys():
                if trigs_dict[path].size == trigs['network/slide_id'][:].size:
                    trigs_dict[path] = _slide_filter(trigs, trigs_dict[path],
                                                     slide_id=slide_id)

    return trigs_dict


# =============================================================================
# Function to apply vetoes to found injections
# =============================================================================
def apply_vetoes_to_found_injs(found_missed_file, found_injs, ifos,
                               veto_file=None, keys=None):
    """Separate injections surviving vetoes from vetoed injections.

    Parameters
    ----------
        found_missed_file: injections results File

        found_injs: dictionary of found injections

        ifos: list of interferometers to use in vetoing

        veto_file: vetoed segments File (optional)

        keys: list of desired dataset names (optional)

    Return
    ------
        found_after_vetoes: dictionary of injections surviving vetoes

        missed_after_vetoes: dictionary of vetoed injections

        found_idx: numpy.array of indices of surviving injections

        veto_idx: numpy.array of indices of vetoed injections
    """

    keep_keys = keys if keys else found_injs.keys()

    if not found_missed_file or ifos[0]+'/end_time' not in found_injs.keys():
        t_id_key = 'network/template_id'
        if t_id_key not in keep_keys:
            keep_keys = list(keep_keys+[t_id_key])
        empty_dict = dict.fromkeys(keep_keys, numpy.array([]))
        empty_dict[t_id_key] = \
            empty_dict[t_id_key].astype(dtype=numpy.int64)
        return (empty_dict, empty_dict, None, None)

    found_idx = numpy.arange(len(found_injs[ifos[0]+'/end_time'][:]))
    veto_idx = numpy.array([], dtype=numpy.int64)

    if veto_file:
        logging.info("Applying data vetoes to found injections...")
        for ifo in ifos:
            inj_time = found_injs[ifo+'/end_time'][:]
            segs = veto.select_segments_by_definer(veto_file,
                                                   segment_name=None,
                                                   ifo=ifo)
            if len(segs) > 0:
                idx, _ = veto.indices_outside_segments(inj_time,
                                                       [veto_file],
                                                       ifo,
                                                       None)
                veto_idx = numpy.append(veto_idx, idx)
                logging.info("%d injections vetoed due to %s.", len(idx), ifo)
                idx, _ = veto.indices_within_segments(inj_time,
                                                      [veto_file],
                                                      ifo,
                                                      None)
                found_idx = numpy.intersect1d(found_idx, idx)

        veto_idx = numpy.unique(veto_idx)
        logging.info("%d injections vetoed.", len(veto_idx))
        logging.info("%d injections surviving vetoes.", len(found_idx))

    found_after_vetoes = {}
    missed_after_vetoes = {}
    for key in keep_keys:
        if key == 'network/coincident_snr':
            found_injs[key] = get_coinc_snr(found_injs)
        if isinstance(found_injs[key], numpy.ndarray):
            found_after_vetoes[key] = found_injs[key][found_idx]
            missed_after_vetoes[key] = found_injs[key][veto_idx]

    return found_after_vetoes, missed_after_vetoes, found_idx, veto_idx


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
    """Constructs sorted triggers from a trials dictionary for the slides
    requested via slide_dict."""

    sorted_trigs = {}

    # Begin by sorting the triggers into each slide
    for slide_id in slide_dict:
        sorted_trigs[slide_id] = []
    for slide_id, event_id in zip(trigs['network/slide_id'],
                                  trigs['network/event_id']):
        if slide_id in slide_dict:
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
        num_trigs_before = len(sorted_trigs[slide_id])
        sorted_trigs[slide_id] = [event_id for event_id in
                                  sorted_trigs[slide_id]
                                  if trigs['network/end_time_gc'][
                                      trigs['network/event_id'] == event_id][0]
                                  in trial_dict[slide_id]]

        # Check that the number of triggers has not increased after vetoes
        assert len(sorted_trigs[slide_id]) <= num_trigs_before, \
            f"Slide {slide_id} has {num_trigs_before} triggers before the "\
            f"trials dictionary was used and {len(sorted_trigs[slide_id])} "\
            "after. This should not happen."
        # END OF CHECK #

    return sorted_trigs


# =============================================================================
# Extract trigger properties and store them as dictionaries
# =============================================================================
def extract_trig_properties(trial_dict, trigs, slide_dict, seg_dict, keys):
    """Extract and store as dictionaries specific keys of time-slid
    triggers (trigs) compatibly with the trials dictionary (trial_dict)"""

    # Sort the triggers into each slide
    sorted_trigs = sort_trigs(trial_dict, trigs, slide_dict, seg_dict)
    n_surviving_trigs = sum(len(i) for i in sorted_trigs.values())
    msg = f"{n_surviving_trigs} triggers found within the trials dictionary "
    msg += "and sorted."
    logging.info(msg)

    found_trigs = {}
    for key in keys:
        found_trigs[key] = {}

    for slide_id in slide_dict:
        slide_trigs = sorted_trigs[slide_id]
        indices = numpy.nonzero(
            numpy.isin(trigs['network/event_id'], slide_trigs))[0]
        for key in keys:
            if slide_trigs:
                found_trigs[key][slide_id] = get_coinc_snr(trigs)[indices] \
                    if key == 'network/coincident_snr' else trigs[key][indices]
            else:
                found_trigs[key][slide_id] = numpy.asarray([])

    logging.info("Triggers information extracted.")

    return found_trigs


# =============================================================================
# Function to extract ifos from hdfs
# =============================================================================
def extract_ifos(trig_file, ifo=None):
    """Extracts IFOs from hdf file and checks for presence of a specific IFO"""

    # Load hdf file
    hdf_file = HFile(trig_file, 'r')

    # Extract IFOs
    ifos = sorted(list(hdf_file.keys()))

    # Remove unwanted keys from key list to reduce it to the ifos
    for key in ['network', 'found', 'missed']:
        if key in ifos:
            ifos.remove(key)

    # Exit gracefully if the requested IFO is not available
    if ifo and ifo not in ifos:
        err_msg = "The IFO selected with --ifo is unavailable in the data."
        raise RuntimeError(err_msg)

    return ifos


# =============================================================================
# Function to load timeslides
# =============================================================================
def load_time_slides(hdf_file_path):
    """Loads timeslides from PyGRB output file as a dictionary"""
    logging.info("Loading timeslides.")
    hdf_file = HFile(hdf_file_path, 'r')
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

    logging.info("Loading segments.")

    # Long time slides will require mapping between slides and segments
    hdf_file = HFile(hdf_file_path, 'r')
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
def construct_trials(seg_files, seg_dict, ifos, slide_dict, veto_file,
                     hide_onsource=True):
    """Constructs trials from segments, timeslides, and vetoes"""

    logging.info("Constructing trials.")
    trial_dict = {}

    # Get segments
    segs = _read_seg_files(seg_files)

    # Separate segments
    trial_time = abs(segs['on'])

    # Determine the veto segments
    vetoes = _extract_vetoes(veto_file, ifos, segs['off'])

    # Slide vetoes over trials: this can only *reduce* the analysis time
    for slide_id in slide_dict:
        curr_seg_list = seg_dict[slide_id]

        # Fill in the buffer segment list if the onsource data must be hidden
        seg_buffer = segments.segmentlist()
        if hide_onsource:
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

    total_trials = sum(len(trial_dict[slide_id]) for slide_id in slide_dict)
    logging.info("%d trials generated.", total_trials)

    return trial_dict, total_trials


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
    """Return maximum and median of trig_stat and sorted time_veto_max_stat"""

    max_stat = max(trig_stat[slide_id].max() if trig_stat[slide_id].size
                   else 0 for slide_id in slide_dict)

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

    coinc_snr = numpy.array([])
    if 'network/coherent_snr' in trigs_or_injs.keys() and \
            'network/null_snr' in trigs_or_injs.keys():
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
    trigger_file: HFile object for trigger file
    bank_file: filepath for template bank
    """
    ifos = [k for k in trigger_file.keys() if k != 'network']
    if ifos[0]+'/template_hash' not in trigger_file.keys():
        return numpy.array([], dtype=int)
    with HFile(bank_path, "r") as bank:
        hashes = bank['template_hash'][:]
    trig_hashes = trigger_file[f'{ifos[0]}/template_hash'][:]
    trig_ids = numpy.zeros(trig_hashes.shape[0], dtype=int)
    for idx, t_hash in enumerate(hashes):
        matches = numpy.where(trig_hashes == t_hash)
        trig_ids[matches] = idx
    return trig_ids
