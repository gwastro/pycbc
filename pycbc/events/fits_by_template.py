"""
Functions and utils to fit single-detector trigger density in each template
according to a density model
"""

import numpy as np
import logging

import pycbc
from pycbc.events import veto, trigger_fits as trstats
from pycbc import bin_utils

logger = logging.getLogger('pycbc.events.fits_by_template')


def get_stat(trigs, rank_method, threshold, additional_datasets=None):
    """
    Get statisic values from the trigger file, do so efficiently
    """
    logger.info('Calculating stat values')

    # For now this is using the single detector ranking. If we want, this
    # could use the Stat classes in stat.py using similar code as in hdf/io.py
    # This requires additional options, so only change this if it's useful!

    chunk_size = 2 ** 23
    stat = []
    select = []
    size = trigs['end_time'].size
    chunk_start = 0
    additional_datasets = additional_datasets \
        if additional_datasets is not None else []
    ad_dict = {adk: [] for adk in additional_datasets}
    while chunk_start < size:
        chunk_end = chunk_start + chunk_size \
            if (chunk_start + chunk_size) <= size else size

        # read and format chunk of data so it can be read by key
        # as the stat classes expect.
        chunk = {k: trigs[k][chunk_start:chunk_end]
                 for k in trigs if len(trigs[k]) == size}

        chunk_stat = rank_method.get_sngl_ranking(chunk)

        above = chunk_stat >= threshold
        stat.append(chunk_stat[above])
        select.append(above)
        for adk in ad_dict.keys():
            ad_dict[adk].append(chunk[adk][above])
        chunk_start += chunk_size

    ad_dict = {adk: np.concatenate(adv) for adk, adv in ad_dict.items()}
    # Return boolean area that selects the triggers above threshold
    # along with the stat values above threshold
    return np.concatenate(select), np.concatenate(stat), ad_dict


def count_in_template(trigger_set, template_f):
    """
    Function to return the number of triggers in each template.

    Lots of this is index juggling, as the triggers are stored in
    template hash order in the HDF_TRIGGER_MERGE file.

    Parameters
    ----------
    trigger_set: h5py Group
        The group in the trigger file containing triggers to be considered.
        Needed fop
    template_f: h5py File
        The file containing template information - this is needed for the
        template hash information
    """
    logger.info('Counting number of triggers in each template')
    # template boundaries dataset is in order of template_id
    template_boundaries = trigger_set['template_boundaries'][:]
    template_id = np.arange(len(template_boundaries))
    # template boundary values ascend in the same order as template hash
    # hence sort by hash
    hash_sort = np.argsort(template_f['template_hash'][:])
    template_boundaries_hashorder = template_boundaries[hash_sort]
    # reorder template IDs in parallel to the boundary values
    template_id_hashorder = template_id[hash_sort]

    # Calculate the differences between the boundary indices to get the
    # number in each template
    # adding on total number at the end to get number in the last template
    total_number = trigger_set['template_id'].size
    count_in_template_hashorder = \
        np.diff(np.append(template_boundaries_hashorder, total_number))

    # re-reorder values from hash order to template_id order
    template_id_sort = np.argsort(template_id_hashorder)
    return count_in_template_hashorder[template_id_sort]


def retained_segments(trigger_file, ifo, veto_file=None,
                      veto_segment_name=None, gating_veto_windows=None):
    logger.info("Getting segment information")
    # Calculate total time being analysed from segments
    # use set() to eliminate duplicates
    trigger_group = trigger_file[ifo]
    segment_starts = sorted(set(trigger_group['search/start_time'][:]))
    segment_ends = sorted(set(trigger_group['search/end_time'][:]))
    all_segments = veto.start_end_to_segments(
        segment_starts,
        segment_ends
    )

    logger.debug("%.2fs total time", abs(all_segments))
    veto_file = veto_file if veto_file is not None else []
    veto_segment_name = veto_segment_name \
        if veto_segment_name is not None else []

    # now do vetoing
    for vfile, vsegment_name in zip(veto_file, veto_segment_name):
        vetoed_segments = veto.select_segments_by_definer(
            vfile,
            ifo=ifo,
            segment_name=vsegment_name
        )
        all_segments -= vetoed_segments
        logger.info(
           '%2.fs removed when vetoing with %s in %s',
           abs(vetoed_segments),
           vsegment_name,
           vfile,
        )

    # Include gating vetoes
    if gating_veto_windows:
        logger.info("Removing triggers near gates")
        gating_veto = gating_veto_windows[ifo].split(',')
        gveto_before = float(gating_veto[0])
        gveto_after = float(gating_veto[1])
        if gveto_before > 0 or gveto_after < 0:
            raise ValueError("Gating veto window values must be negative "
                             "before gates and positive after gates.")
        if not (gveto_before == 0 and gveto_after == 0):
            autogate_times = np.unique(trigger_group['gating/auto/time'][:])
            if 'gating/file' in trigger_group:
                detgate_times = trigger_group['gating/file/time'][:]
            else:
                detgate_times = []
            gate_times = np.concatenate((autogate_times, detgate_times))
            gveto_segs = veto.start_end_to_segments(
                gate_times + gveto_before,
                gate_times + gveto_after).coalesce()
            all_segments -= gveto_segs
            logger.info('%.2f removed near gates', abs(gveto_segs))

    return all_segments


def segment_mask(time, keep_segments):
    """
    Given a set of segments, which triggers should be kept
    (note that this ensures we retain the correct order)
    """
    vmask = np.zeros(time.size, dtype=bool)
    start, end = veto.segments_to_start_end(keep_segments)
    idx_keep = veto.indices_within_times(
        time,
        start,
        end,
    )
    vmask[idx_keep] = True
    return vmask


def prune_mask(stat, time, prune_param, prune_number=2, prune_bins=2, prune_window=0.1, log=False):
    # retain_trigs are the triggers to be retained after pruning
    retain_trigs = np.ones(stat.size, dtype=bool)
    # test_trigs are the triggers to be tested in the pruning loop
    test_trigs = np.ones(stat.size, dtype=bool)

    # initialize bin storage
    prunedtimes = {}
    for i in range(prune_bins):
        prunedtimes[i] = []

    # Set up the bins for pruning
    binclass = bin_utils.LogarithmicBins if log else bin_utils.LinearBins
    bins = binclass(prune_param.min(), prune_param.max(), prune_bins)

    # many trials may be required to prune in 'quieter' bins
    pruning_done = False
    skipped_prunes = 0
    for j in range(1000):
        # are all the bins full already?
        numpruned = sum([len(prunedtimes[i]) for i in range(prune_bins)])
        if numpruned == prune_bins * prune_number:
            logging.info('Finished pruning!')
            pruning_done = True
            break
        if numpruned > prune_bins * prune_number:
            logging.error(
                'Uh-oh, we pruned too many things .. %i, to be precise',
                numpruned
            )
            raise RuntimeError
        if not any(test_trigs):
            logging.error(
                "Too many things have been removed in the pruning loop, "
                "and there is nothing left"
            )
            raise RuntimeError
        loudest = np.argmax(stat[test_trigs])
        lstat = stat[test_trigs][loudest]
        ltime = time[test_trigs][loudest]
        lbin = bins[prune_param[test_trigs][loudest]]

        # is the bin where the loudest trigger lives full already?
        if len(prunedtimes[lbin]) == prune_number:
            logging.debug(
                '%i - Bin %i full, not pruning event with stat %.3f at '
                'time %.3f',
                j,
                lbin,
                lstat,
                ltime
            )
            # prune the reference trigger array
            remove = abs(time - ltime) <= prune_window
            test_trigs[remove] = False
            skipped_prunes += 1
        else:
            if skipped_prunes > 0:
                logging.info(
                    '%d events skipped due to full bins',
                    skipped_prunes
                )
            logging.info(
                'Pruning event with stat %.3f at %.3f in bin %i',
                lstat,
                ltime,
                lbin
            )
            # now do the pruning
            remove = abs(time - ltime) <= prune_window
            retain_trigs[remove] = False
            # also for the reference trig arrays
            test_trigs[remove] = False
            prunedtimes[lbin].append(ltime)
        del remove

    if not pruning_done:
        logging.warning(
            "Not enough triggers have been pruned after the loop over "
            "1000 triggers - something could be wrong!"
        )
    logging.info(
        '%i trigs remain after pruning loop',
        np.count_nonzero(retain_trigs)
    )
    return retain_trigs


def fit_triggers(stat, template_id, fit_function, stat_threshold, template_ids):
    logging.info("Fitting in each template")
    # sorting into template_id order
    tsort = template_id.argsort()
    template_id = template_id[tsort]
    stat = stat[tsort]

    logging.info("Getting trigger ranges for each template id")
    left = np.searchsorted(template_id, template_ids, side='left')
    right = np.searchsorted(template_id, template_ids, side='right')

    counts_above = right - left
    fits = np.ones(len(template_ids), dtype=np.float64) * -100

    logging.info("Calculating fit coefficients")
    for j in irange(len(template_ids)):
        if not j % 10000:
            logging.info(
                "Fitting template %d out of %d",
                j, len(template_ids),
            )
        elif not j % 1000:
            logging.debug(
                "Fitting template %d out of %d",
                j, len(template_ids),
            )

        if counts_above[j] == 0:
            continue
        stat_in_template = stat[left[j]:right[j]]
        alpha, _ = trstats.fit_above_thresh(
            fit_function,
            stat_in_template,
            stat_threshold
        )

        fits[j] = alpha

    return counts_above, fits
