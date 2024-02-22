"""
Functions for combining mutiple fit coefficient files
"""

import logging
import copy
import numpy as np
from scipy.stats import norm

logger = logging.getLogger('pycbc.events.fits_combination')


def dist(idx_1, idx_2, parvals, smoothing_width, sort=None):
    """
    Computes the vector of parameter values at index/indices idx_1 and
    index/indices idx_2, and gives the Euclidean distance between
    the two with a metric of 1/(smoothing width^2)

    If idx_1 and idx_2 are both more than one index, then they must be
    the same size.
    """
    dsq = 0
    for pval, swidth in zip(parvals, smoothing_width):
        if sort is not None:
            pval = pval[sort]
        dsq += (pval[idx_2] - pval[idx_1]) ** 2.0 / swidth ** 2.0
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
    if weights is None:
        # Equal weighting to each template
        weights = np.ones_like(template_idx)

    nabove_t_smoothed = np.average(nabove[template_idx], weights=weights)
    ntotal_t_smoothed = np.average(ntotal[template_idx], weights=weights)
    invalphan_mean = np.average(invalphan[template_idx], weights=weights)

    return_tuple = (nabove_t_smoothed,
                    nabove_t_smoothed / invalphan_mean,
                    ntotal_t_smoothed)
    return return_tuple


def smooth_tophat(nabove, invalphan, ntotal, dists, **kwargs):  # pylint:disable=unused-argument
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


def smooth_n_closest(nabove, invalphan, ntotal, dists, total_trigs=500, **kwargs):  # pylint:disable=unused-argument
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


def smooth_distance_weighted(nabove, invalphan, ntotal, dists, **kwargs):  # pylint:disable=unused-argument
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


def digest_smoothing_kwargs(args, parser):
    """
    Check inputs and pack the arguments into a standard dictionary.

    Parameters
    ----------
    args : argparse.ArgumentParser.parse_args() output
       parser command line arguments, must contain certain values
    parser : argparse.ArgumentParser
       parser, so that we can raise errors appropriately

    Returns
    -------
    kwarg_dict : dictionary
        Dictionary containing the required information to control the
        template fits smoothing
    """
    kwarg_dict = {}

    # Input sanitisation
    if not len(args.log_param) == len(args.fit_param) or \
            not len(args.log_param) == len(args.smoothing_width):
        parser.error("--log-param, --fit-param and --smoothing-width "
                     "must all have the same number of inputs")

    kwarg_dict['fit_param'] = args.fit_param
    kwarg_dict['method'] = args.smoothing_method
    kwarg_dict['width'] = args.smoothing_width
    kwarg_dict['parnames'] = []

    # work out whether we need to apply a log function to the parameter
    for param, slog in zip(args.fit_param, args.log_param):
        log_param = slog[:5].title()
        if not log_param in ["True", "False"]:
            logger.error(
                "--log-param value cannot b eavluated to True or False, "
                "provide this correctly"
            )
            raise ValueError
        # Provide a function which can calculate the log if needed,
        # or is transparent if not
        par_func = np.log if log_param == 'True' else lambda x: x
        kwarg_dict[f'{param}_func'] = par_func

        # For use when reporting the parameters
        parname = param if not log_param == 'True' else f'log({param})'
        kwarg_dict['parnames'].append(parname)

        logger.info('Using param: %s', parname)

    # Don't fail if not given
    if args.smoothing_keywords is not None:
        smooth_kwargs = args.smoothing_keywords
    else:
        smooth_kwargs = []

    for inputstr in smooth_kwargs:
        # Pack these into a dictionary
        try:
            key, value = inputstr.split(':')
            kwarg_dict[key] = value
        except ValueError as ve:
            logger.error(
                "--smoothing-keywords must take input in the "
                "form KWARG1:VALUE1 KWARG2:VALUE2 KWARG3:VALUE3 ... "
                "Received %s", ' '.join(args.smoothing_keywords)
            )
            raise ve
    return kwarg_dict


def report_percentage(i, length, pc_report=10, log_func=logger.info):
    """
    Convenience function - report how long through the loop we are.

    Parameters
    ----------
    i: integer
        index being looped through
    length : integer
        number of loops we will go through in total
    pc_report : integer
        When to report, default every 10 percent
    log_func : function
        The logging function used to report the value, this can be
        changed in order to use more or less verbosity if required
    """
    pc_now = int(np.floor(i / length * 100))
    # This bit stops getting loads of logging each time
    pc_last = int(np.floor((i - 1) / length * 100))
    if not pc_now % pc_report and not pc_last == pc_now:
        log_func("Template %d out of %d (%.0f%%)", i, length, pc_now)


def oned_tophat(parvals, smooth_kwargs, nabove, invalphan, ntotal):
    """
    Handle the one-dimensional case of tophat smoothing separately
    as it is easier to optimize computational performance.

    Parameters
    ----------
    parvals: tuple of numpy.ndarray (should be only one)
        Parameter values being used to smooth over
    smooth_kwargs: dict
        Dictionary defining how the smoothing is to be completed
        This function only uses the 'width' value
        Made e.g. using smoothing_kwargs above
    nabove: numpy.ndarray
        The count of triggers above the fit threshold in each template
    invalphan: numpy.ndarray
        The array of n_above / alpha values for all templates
    ntotal: numpy.ndarray
        The total number of triggers in each template

    Returns
    -------
    nabove_smoothed: numpy.ndarray
        The smoothed count above threshold value
    alpha_smoothed: numpy.ndarray
        The smoothed fit coefficient (alpha) value
    ntotal_smoothed: numpy.ndarray
        The smoothed total count in template value
    """
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

    # Number of templates in the region considered
    num = right - left

    logger.info("Smoothing ...")
    nabove_smoothed = (nasum[right] - nasum[left]) / num
    invmean = (invsum[right] - invsum[left]) / num
    alpha_smoothed = nabove_smoothed / invmean
    ntotal_smoothed = (ntsum[right] - ntsum[left]) / num

    return nabove_smoothed, alpha_smoothed, ntotal_smoothed


# Number of smoothing lengths around the current template where
# distances will be calculated
# n_closest has no limit as it needs to contain enough
# templates to contain n triggers, which we cannot know beforehand

_smooth_cut = {
    'smooth_tophat': 1,
    'n_closest': np.inf,
    'distance_weighted': 3,
}


def cut_templates(parvals, smoothing_method, smoothing_width, val_names):
    """
    Work out which templates are not going to be used in the smoothing
    and remove them from the list we are checking

    Parameters
    ----------
    parvals: tuple of ndarrays
       Values of parameters of the different dimensions being considered
       for the distance
    smoothing_method: string
       String defining which method is being used for the smoothing, this
       is used to determine how many 'widths' should be cut.
    smoothing_width: list of floats
       Same size as the tuple parvals
       The re-scaling value in each dimension 
    val_names: list of strings
       (Just for logging) List of names corresponding to the parvals
       Same size as the tuple parvals

    Returns
    -------
    par_sort: numpy ndarray of int
       The ordering of the templates required to convert into the same
       as used by the slice objects

    slices: list of slice objects
       Once the templates have been ordered by par_sort, these slices
       define the region which will be used in the smoothing for that
       template
    """

    cnum = _smooth_cut[smoothing_method]
    if not np.isfinite(cnum):
        # This is a non-cutting smoothing method, return values so that
        # nothing happens
        slices = [slice(0, parvals[0].size) for _ in parvals[0]]
        return None, slices

    # We only cut along the dimension which will work best - maybe this can
    # be done in multiple dimensions, but it may be a case of diminishing
    # returns
    cut_lengths = [s * cnum for s in smoothing_width]
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

    slices = [slice(left, right) for left, right in zip(lefts, rights)]

    n_removed = len(pvals) - rights + lefts
    logger.info(
        "Cutting between %d and %d templates for each smoothing",
        n_removed.min(),
        n_removed.max(),
    )

    return par_sort, slices


def smooth_samples(parvals, smooth_kwargs, nabove, invalphan, ntotal):
    """
    Parameters
    ----------
    parvals: tuple of numpy.ndarray
        Parameter values being used to smooth over
    smooth_kwargs: dict
        Dictionary defining how the smoothing is to be completed
        Made e.g. using smoothing_kwargs above
    nabove: numpy.ndarray
        The count of triggers above the fit threshold in each template
    invalphan: numpy.ndarray
        The array of n_above / alpha values for all templates
    ntotal: numpy.ndarray
        The total number of triggers in each template

    Returns
    -------
    nabove_smoothed: numpy.ndarray
        The smoothed count above threshold value
    alpha_smoothed: numpy.ndarray
        The smoothed fit coefficient (alpha) value
    ntotal_smoothed: numpy.ndarray
        The smoothed total count in template value

    """
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

    logger.info("Smoothing ...")

    nabove_smoothed = []
    alpha_smoothed = []
    ntotal_smoothed = []
    for i in np.arange(0, len(nabove)):
        # Info logging every 10 percent (defaults)
        report_percentage(i, len(nabove))
        # Debug logging every percent
        report_percentage(i, len(nabove), pc_report=1, log_func=logger.debug)
        slc = slices[i]
        temp_dists = dist(i, slc, parvals, smooth_kwargs['width'])

        smoothed_tuple = smooth(
            nabove[slc],
            invalphan[slc],
            ntotal[slc],
            temp_dists,
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
