"""
Functions for combining mutiple fit coefficient files
"""

import numpy as np
import logging
import copy

import pycbc
from pycbc.events import veto, trigger_fits as trstats
from pycbc import bin_utils

logger = logging.getLogger('pycbc.events.fits_combination')

def dist(i1, i2, parvals, smoothing_width, sort=None):
    """
    Computes the vector of parameter values at index/indices i1 and
    index/indices i2, and gives the Euclidean distance between
    the two with a metric of 1/(smoothing width^2)
    """
    dsq = 0
    for v, s in zip(parvals, smoothing_width):
        if sort is not None:
            v = v[sort]
        dsq += (v[i2] - v[i1]) ** 2.0 / s ** 2.0
    return dsq ** 0.5


def smooth_templates(nabove, invalphan, ntotal, template_idx,
                     weights=None):
    """
    Find the smoothed values according to the specified templates,
    weighted appropriately.
    The max likelihood fit for 1/alpha is linear in the trigger
    statistic values, so we perform a possibly-weighted average
    of (n_above / alpha) over templates and then invert this
    and multiply by (smoothed) nabove to obtain smoothed alpha.

    Parameters
    ----------
    nabove: ndarray
        The array of counts of triggers above threshold for all templates
    invalphan: ndarray
        The array of n_above / alpha values for all templates
    ntotal: ndarray
        The array of count of triggers in the template, both above and
        below threshold
    template_idx: ndarray of ints
        The indices of the templates to be used for the smoothing

    Optional Parameters
    -------------------
    weights: ndarray
        Weighting factor to apply to the templates specified by template_idx

    Returns
    -------
    tuple: 3 floats
        First float: the smoothed count above threshold value
        Second float: the smoothed fit coefficient (alpha) value
        Third float: the smoothed total count in template value

    """
    if weights is None: weights = np.ones_like(template_idx)
    nabove_t_smoothed = np.average(nabove[template_idx], weights=weights)
    ntotal_t_smoothed = np.average(ntotal[template_idx], weights=weights)
    invalphan_mean = np.average(invalphan[template_idx], weights=weights)

    return_tuple = (nabove_t_smoothed,
                    nabove_t_smoothed / invalphan_mean,
                    ntotal_t_smoothed)
    return return_tuple


def smooth_tophat(nabove, invalphan, ntotal, dists, **kwargs): # pylint:disable=unused-argument
    """
    Smooth templates using a tophat function with templates within unit
    dists
    """
    idx_within_area = np.flatnonzero(dists < 1.)
    return smooth_templates(
        nabove,
        invalphan,
        ntotal,
        idx_within_area
    )


def smooth_n_closest(nabove, invalphan, ntotal, dists, total_trigs=500, **kwargs): # pylint:disable=unused-argument
    """
    Smooth templates according to the closest N templates
    No weighting is applied
    """
    dist_sort = np.argsort(dists)
    templates_required = 0
    n_triggers = 0
    # Count number of templates required to gather total_trigs templates,
    # start at closest
    # total_trigs, if supplied on command line will be a str so convert to int
    while n_triggers < int(total_trigs):
        n_triggers += nabove[dist_sort[n_triggers]]
        templates_required += 1
    logger.debug(
        "%d templates required to obtain %d triggers",
        templates_required,
        n_triggers
    )
    idx_to_smooth = dist_sort[:templates_required]
    return smooth_templates(nabove, invalphan, ntotal, idx_to_smooth)


def smooth_distance_weighted(nabove, invalphan, ntotal, dists, **kwargs): # pylint:disable=unused-argument
    """
    Smooth templates weighted according to dists in a unit-width normal
    distribution, truncated at three sigma
    """
    idx_within_area = np.flatnonzero(dists < 3.)
    weights = norm.pdf(dists[idx_within_area])
    return smooth_templates(nabove, invalphan, ntotal,
                            idx_within_area, weights=weights)

_smooth_dist_func = {
    'smooth_tophat': smooth_tophat,
    'n_closest': smooth_n_closest,
    'distance_weighted': smooth_distance_weighted
}


def smooth(nabove, invalphan, ntotal, dists, **kwargs):
    """
    Wrapper for smoothing according to a function defined by smoothing_method

    nabove, invalphan, ntotal are as defined in the above smooth_templates
    function docstring

    dists is an array of the distances of the templates from the
    template of interest
    """
    return _smooth_dist_func[kwargs['method']](nabove, invalphan,
                                               ntotal, dists, **kwargs)


# Number of smoothing lengths around the current template where
# distances will be calculated
# n_closest has no limit as it needs to contain enough
# templates to contain n triggers, which we cannot know beforehand

_smooth_cut = {
    'smooth_tophat': np.inf, #1,
    'n_closest': np.inf,
    'distance_weighted': 3,
}


def smoothing_kwargs(args):
    kwarg_dict = {}

    # Input sanitisation
    assert len(args.log_param) == len(args.fit_param) == len(args.smoothing_width)

    kwarg_dict['fit_param'] = args.fit_param
    kwarg_dict['method'] = args.smoothing_method
    kwarg_dict['width'] = args.smoothing_width
    kwarg_dict['parnames'] = []
    # work out the function needed to 
    for param, slog in zip(args.fit_param, args.log_param):
        log_param = eval(slog[:5].title())
        par_func = np.log if log_param else lambda x: x
        kwarg_dict[f'{param}_func'] = par_func
        parname = param if not log_param else f'log({param})'
        kwarg_dict['parnames'].append(parname)
        if not log_param:
            logger.info('Using param: %s', param)
        else:
            logger.info('Using log param: %s', param)

    if args.smoothing_keywords:
        smooth_kwargs = args.smoothing_keywords
    else:
        smooth_kwargs = []

    for inputstr in smooth_kwargs:
        try:
            key, value = inputstr.split(':')
            kwarg_dict[key] = value
        except ValueError:
                err_txt = "--smoothing-keywords must take input in the " \
                          "form KWARG1:VALUE1 KWARG2:VALUE2 KWARG3:VALUE3 ... " \
                          "Received {}".format(' '.join(args.smoothing_keywords))
                raise ValueError(err_txt)
    return kwarg_dict

def report_percentage(i, length):
    """
    Convenience function - report how long through the loop we are.
    Every ten percent
    Parameters
    ----------
    i: integer
        index being looped through
    length : integer
        number of loops we will go through in total
    """
    pc = int(np.floor(i / length * 100))
    pc_last = int(np.floor((i - 1) / length * 100))
    if not pc % 10 and pc_last % 10:
        logger.info(f"Template {i} out of {length} ({pc:.0f}%)")

def oned_tophat(parvals, smooth_kwargs, nabove, invalphan, ntotal):
    # Handle the one-dimensional case of tophat smoothing separately
    # as it is easier to optimize computational performance.
    logger.info("Using efficient 1D tophat smoothing")
    sort = parvals[0].argsort()
    parvals_0 = parvals[0][sort]
    # For each template, find the range of nearby templates which fall within
    # the chosen window.
    left = np.searchsorted(
        parvals_0,
        parvals[0] - smooth_kwargs['width'][0]
    )
    right = np.searchsorted(
        parvals_0,
        parvals[0] + smooth_kwargs['width'][0]
    ) - 1

    del parvals_0
    # Precompute the sums so we can quickly look up differences between
    # templates
    nasum = nabove.cumsum()
    invsum = invalphan.cumsum()
    ntsum = ntotal.cumsum()
    num = right - left

    logger.info("Smoothing ...")
    nabove_smoothed = (nasum[right] - nasum[left]) / num
    invmean = (invsum[right] - invsum[left]) / num
    alpha_smoothed = nabove_smoothed / invmean
    ntotal_smoothed = (ntsum[right] - ntsum[left]) / num

    return nabove_smoothed, alpha_smoothed, ntotal_smoothed

def cut_templates(parvals, smoothing_method, smoothing_width, val_names):
    c = _smooth_cut[smoothing_method]
    if not np.isfinite(c):
        # This is a non-cutting smoothing method, return values so that
        # nothing happens
        slices = [slice(0, parvals[0].size) for _ in parvals[0]]
        return None, slices

    cut_lengths = [s * c for s in smoothing_width]
    # Find the "longest" dimension in cut lengths
    sort_dim = np.argmax([(v.max() - v.min()) / c
                          for v, c in zip(parvals, cut_lengths)])
    logger.info("Sorting / Cutting on %s", val_names[sort_dim])

    # Sort parvals by the sort dimension
    par_sort = np.argsort(parvals[sort_dim])
    pvals = copy.deepcopy(parvals[sort_dim])[par_sort]

    # For each template, find the range of nearby templates which fall within
    # the chosen window.
    lefts = np.searchsorted(
        pvals,
        pvals - cut_lengths[sort_dim]
    )
    rights = np.searchsorted(
        pvals,
        pvals + cut_lengths[sort_dim]
    )

    slices = [slice(l,r) for l, r in zip(lefts, rights)]

    n_removed = len(pvals) - rights + lefts
    logger.info(
        "Cutting between %d and %d templates for each smoothing",
        n_removed.min(),
        n_removed.max(),
    )

    return par_sort, slices

def smooth_samples(parvals, smooth_kwargs, nabove, invalphan, ntotal):
    par_sort, slices = cut_templates(
        parvals,
        smooth_kwargs['method'],
        smooth_kwargs['width'],
        smooth_kwargs['parnames']
    )

    if par_sort is not None:
        # Sort parvals by the sort dimension
        parvals = [p[par_sort] for p in parvals]

        # Sort the values to be smoothed by parameter value
        nabove = nabove[par_sort]
        invalphan = invalphan[par_sort]
        ntotal = ntotal[par_sort]

    logging.info("Smoothing ...")

    nabove_smoothed = []
    alpha_smoothed = []
    ntotal_smoothed = []
    for i in np.arange(0, len(nabove)):
        report_percentage(i, len(nabove))
        slc = slices[i]
        d = dist(i, slc, parvals, smooth_kwargs['width'])

        smoothed_tuple = smooth(
            nabove[slc],
            invalphan[slc],
            ntotal[slc],
            d,
            **smooth_kwargs
        )
        nabove_smoothed.append(smoothed_tuple[0])
        alpha_smoothed.append(smoothed_tuple[1])
        ntotal_smoothed.append(smoothed_tuple[2])

    if par_sort is not None:
        # Undo the sorts
        unsort = np.argsort(par_sort)
        parvals = [p[unsort] for p in parvals]
        nabove_smoothed = np.array(nabove_smoothed)[unsort]
        alpha_smoothed = np.array(alpha_smoothed)[unsort]
        ntotal_smoothed = np.array(ntotal_smoothed)[unsort]

    return nabove_smoothed, alpha_smoothed, ntotal_smoothed


#def smooth_samples(parvals, smoothing_method, smoothing_width, nabove, invalphan):
#
#    if len(parvals) == 1 and smoothing_method == 'smooth_tophat':
#        return oned_tophat(parvals, smoothing_width, nabove, invalphan)
#    c = _smooth_cut[smoothing_method]
#    cut_lengths = [s * c for s in smoothing_width]
#    # Find the "longest" dimension in cut lengths
#    sort_dim = np.argmax([(v.max() - v.min()) / c
#                              for v, c in zip(parvals, cut_lengths)])
#    logger.info("Sorting / Cutting on dimension %s", parnames[sort_dim])
#
#    # Sort parvals by the sort dimension
#    par_sort = np.argsort(parvals[sort_dim])
#    parvals = [p[par_sort] for p in parvals]
#
#    # For each template, find the range of nearby templates which fall within
#    # the chosen window.
#    lefts = np.searchsorted(parvals[sort_dim],
#            parvals[sort_dim] - cut_lengths[sort_dim])
#    rights = np.searchsorted(parvals[sort_dim],
#            parvals[sort_dim] + cut_lengths[sort_dim])
#    n_removed = len(parvals[0]) - rights + lefts
#    logger.info(
#        "Cutting between %d and %d templates for each smoothing",
#        n_removed.min(),
#       n_removed.max()
#    )
#    # Sort the values to be smoothed by parameter value
#    nabove = nabove[par_sort]
#    invalphan = invalphan[par_sort]
#    ntotal = ntotal[par_sort]
#    logger.info("Smoothing ...")
#    slices = [slice(l,r) for l, r in zip(lefts, rights)]
#    for i in rang:
#        report_percentage(i, rang.max())
#        slc = slices[i]
#        d = dist(i, slc, parvals, smoothing_width)
#
#        smoothed_tuple = smooth(
#            nabove[slc],
#            invalphan[slc],
#            ntotal[slc],
#            d,
#            smoothing_method,
#            **kwarg_dict
#        )
#        nabove_smoothed.append(smoothed_tuple[0])
#        alpha_smoothed.append(smoothed_tuple[1])
#        ntotal_smoothed.append(smoothed_tuple[2])
#
#    # Undo the sorts
#    unsort = np.argsort(par_sort)
#    parvals = [p[unsort] for p in parvals]
#    nabove_smoothed = np.array(nabove_smoothed)[unsort]
#    alpha_smoothed = np.array(alpha_smoothed)[unsort]
#    ntotal_smoothed = np.array(ntotal_smoothed)[unsort]
#else:
#    logger.info("Smoothing ...")
#    for i in rang:
#        report_percentage(i, rang.max())
#        d = dist(i, rang, parvals, args.smoothing_width)
#        smoothed_tuple = smooth(
#            nabove,
#            invalphan,
#            ntotal,
#            d,
#            args.smoothing_method,
#            **kwarg_dict
#        )
#        nabove_smoothed.append(smoothed_tuple[0])
#        alpha_smoothed.append(smoothed_tuple[1])
#        ntotal_smoothed.append(smoothed_tuple[2])

