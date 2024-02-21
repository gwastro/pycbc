# Copyright (C) 2022 Gareth Cabourn Davies
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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#

"""
This module contains functions to calculate the significance
through different estimation methods of the background, and functions that
read in the associated options to do so.
"""
import logging
import copy
import lal
import numpy as np
from pycbc.events import trigger_fits as trstats

logger = logging.getLogger('pycbc.events.significance')


def count_n_louder(bstat, fstat, dec,
                   **kwargs):  # pylint:disable=unused-argument
    """ Calculate for each foreground event the number of background events
    that are louder than it.

    Parameters
    ----------
    bstat: numpy.ndarray
        Array of the background statistic values
    fstat: numpy.ndarray or scalar
        Array of the foreground statistic values or single value
    dec: numpy.ndarray
        Array of the decimation factors for the background statistics

    Returns
    -------
    cum_back_num: numpy.ndarray
        The cumulative array of background triggers.
    fore_n_louder: numpy.ndarray
        The number of background triggers above each foreground trigger
    """
    sort = bstat.argsort()
    bstat = bstat[sort]
    dec = dec[sort]

    # calculate cumulative number of triggers louder than the trigger in
    # a given index. We need to subtract the decimation factor, as the cumsum
    # includes itself in the first sum (it is inclusive of the first value)
    n_louder = dec[::-1].cumsum()[::-1] - dec

    # Determine how many values are louder than the foreground ones
    # We need to subtract one from the index, to be consistent with definition
    # of n_louder, as here we do want to include the background value at the
    # found index
    idx = np.searchsorted(bstat, fstat, side='left') - 1

    # If the foreground are *quieter* than the background or at the same value
    # then the search sorted algorithm will choose position -1, which does not
    # exist. We force it back to zero.
    if isinstance(idx, np.ndarray):  # Case where our input is an array
        idx[idx < 0] = 0
    else:  # Case where our input is just a scalar value
        if idx < 0:
            idx = 0

    fore_n_louder = n_louder[idx]

    unsort = sort.argsort()
    back_cum_num = n_louder[unsort]
    return back_cum_num, fore_n_louder


def n_louder_from_fit(back_stat, fore_stat, dec_facs,
                      fit_function='exponential', fit_threshold=0,
                      **kwargs):  # pylint:disable=unused-argument
    """
    Use a fit to events in back_stat in order to estimate the
    distribution for use in recovering the estimate count of louder
    background events. Below the fit threshold, use the n_louder
    method for these triggers

    back_stat: numpy.ndarray
        Array of the background statistic values
    fore_stat: numpy.ndarray or scalar
        Array of the foreground statistic values or single value
    dec_facs: numpy.ndarray
        Array of the decimation factors for the background statistics
    fit_function: str
        Name of the function to be used for the fit to background
        statistic values
    fit_threshold: float
        Threshold above which triggers use the fitted value, below this
        the counted number of louder events will be used

    Returns
    -------
    bg_n_louder: numpy.ndarray
        The estimated number of background events louder than each
        background event
    fg_n_louder: numpy.ndarray
        The estimated number of background events louder than each
        foreground event
    """

    # Calculate the fitting factor of the ranking statistic distribution
    alpha, _ = trstats.fit_above_thresh(fit_function, back_stat,
                                        thresh=fit_threshold,
                                        weights=dec_facs)

    # Count background events above threshold as the cum_fit is
    # normalised to 1
    bg_above = back_stat > fit_threshold
    n_above = np.sum(dec_facs[bg_above])
    fg_above = fore_stat > fit_threshold

    # These will be overwritten, but just to silence a warning
    # in the case where trstats.cum_fit returns zero
    bg_n_louder = np.zeros_like(back_stat)
    fg_n_louder = np.zeros_like(fore_stat)

    # Ue the fit above the threshold
    bg_n_louder[bg_above] = n_above * trstats.cum_fit(fit_function,
                                                      back_stat[bg_above],
                                                      alpha,
                                                      fit_threshold)
    fg_n_louder[fg_above] = n_above * trstats.cum_fit(fit_function,
                                                      fore_stat[fg_above],
                                                      alpha,
                                                      fit_threshold)

    # Below the fit threshold, we expect there to be sufficient events
    # to use the count_n_louder method, and the distribution may deviate
    # from the fit function
    fg_below = np.logical_not(fg_above)
    bg_below = np.logical_not(bg_above)

    # Count the number of below-threshold background events louder than the
    # bg and foreground
    bg_n_louder[bg_below], fg_n_louder[fg_below] = \
        count_n_louder(back_stat[bg_below], fore_stat[fg_below], dec_facs)

    # As we have only counted the louder below-threshold events, need to
    # add the above threshold events, which by definition are louder than
    # all the below-threshold events
    bg_n_louder[bg_below] += n_above
    fg_n_louder[fg_below] += n_above

    return bg_n_louder, fg_n_louder


_significance_meth_dict = {
    'trigger_fit': n_louder_from_fit,
    'n_louder': count_n_louder
}

_default_opt_dict = {
    'method': 'n_louder',
    'fit_threshold': None,
    'fit_function': None,
    'far_limit': 0.}


def get_n_louder(back_stat, fore_stat, dec_facs,
                 method=_default_opt_dict['method'],
                 **kwargs):  # pylint:disable=unused-argument
    """
    Wrapper to find the correct n_louder calculation method using standard
    inputs
    """
    return _significance_meth_dict[method](
        back_stat,
        fore_stat,
        dec_facs,
        **kwargs)


def get_far(back_stat, fore_stat, dec_facs,
            background_time,
            method=_default_opt_dict['method'],
            **kwargs):  # pylint:disable=unused-argument
    """
    Return the appropriate FAR given the significance calculation method

    If the n_louder method is used, find the IFAR according to Eq.17-18
    of Usman et al., arXiv:1508.02357. The p-value of a candidate in a
    search of duration T, with n_bg louder time shifted events over a
    total background time T_bg is
    `p = 1 - exp(-T * (n_bg + 1) / T_bg)`
    corresponding to an effective false alarm rate of (n_bg + 1) / T_bg.

    If the trigger_fit method is used, we are extrapolating the background
    for the specific aim of FAR not being limited to 1 / T_bg, and so we
    do not add 1 to n_bg

    Parameters
    ----------
    See description in get_n_louder for most parameters

    background_time: float
       The amount of time to convert the number of louder events into
       a FAR

    """
    bg_n_louder, fg_n_louder = get_n_louder(
        back_stat,
        fore_stat,
        dec_facs,
        method=method,
        **kwargs
    )

    # If we are counting the number of louder events in the background,
    # we add one. This is part of the p-value calculation in Usman 2015.
    # If we are doing trigger fit extrapolation, this is not needed
    if method == 'n_louder':
        bg_n_louder += 1
        fg_n_louder += 1

    bg_far = bg_n_louder / background_time
    fg_far = fg_n_louder / background_time

    return bg_far, fg_far


def insert_significance_option_group(parser):
    """
    Add some options for use when a significance is being estimated from
    events or event distributions.
    """
    parser.add_argument('--far-calculation-method', nargs='+',
                        default=[],
                        help="Method used for FAR calculation in each "
                             "detector combination, given as "
                             "combination:method pairs, i.e. "
                             "H1:trigger_fit H1L1:n_louder H1L1V1:n_louder "
                             "etc. Method options are ["
                             + ",".join(_significance_meth_dict.keys()) +
                             "]. Default = n_louder for all not given")
    parser.add_argument('--fit-threshold', nargs='+', default=[],
                        help="Trigger statistic fit thresholds for FAN "
                             "estimation, given as combination-value pairs "
                             "ex. H1:0 L1:0 V1:-4 for all combinations with "
                             "--far-calculation-method = trigger_fit")
    parser.add_argument("--fit-function", nargs='+', default=[],
                        help="Functional form for the statistic slope fit if "
                             "--far-calculation-method is 'trigger_fit'. "
                             "Given as combination:function pairs, i.e. "
                             "H1:exponential H1L1:n_louder H1L1V1:n_louder. "
                             "Options: ["
                             + ",".join(trstats.fitalpha_dict.keys()) + "]. "
                             "Default = exponential for all")
    parser.add_argument('--limit-ifar', nargs='+', default=[],
                        help="Impose upper limits on IFAR values (years)"
                             ". Given as combination:value pairs, eg "
                             "H1L1:10000 L1:1000. Used to avoid under/"
                             "overflows for loud signals and injections "
                             "using the fit extrapolation method. A value"
                             " 0 or no value means unlimited IFAR")


def positive_float(inp):
    """
    Wrapper around float conversion which ensures that the float must be
    positive or zero
    """
    fl_in = float(inp)
    if fl_in < 0:
        logger.warning("Value provided to positive_float is less than zero, "
                       "this is not allowed")
        raise ValueError
    return fl_in


def check_significance_options(args, parser):
    """
    Check the significance group options
    """
    # Check that the combo:method/function/threshold are in the
    # right format, and are in allowed combinations
    lists_to_check = [(args.far_calculation_method, str,
                       _significance_meth_dict.keys()),
                      (args.fit_function, str,
                       trstats.fitalpha_dict.keys()),
                      (args.fit_threshold, float, None),
                      (args.limit_ifar, positive_float, None)]

    for list_to_check, type_to_convert, allowed_values in lists_to_check:
        combo_list = []
        for combo_value in list_to_check:
            try:
                combo, value = tuple(combo_value.split(':'))
            except ValueError:
                parser.error("Need combo:value format, got %s" % combo_value)

            if combo in combo_list:
                parser.error("Duplicate combo %s in a significance "
                             "option" % combo)
            combo_list.append(combo)

            try:
                type_to_convert(value)
            except ValueError:
                err_fmat = "Value {} of combo {} can't be converted"
                parser.error(err_fmat.format(value, combo))

            if allowed_values is not None and \
                    type_to_convert(value) not in allowed_values:
                err_fmat = "Value {} of combo {} is not in allowed values: {}"
                parser.error(err_fmat.format(value, combo, allowed_values))

    # Are the functions/thresholds appropriate for the methods given?
    methods = {}
    # A method has been specified
    for combo_value in args.far_calculation_method:
        combo, value = tuple(combo_value.split(':'))
        methods[combo] = value

    # A function or threshold has been specified
    function_or_thresh_given = []
    for combo_value in args.fit_function + args.fit_threshold:
        combo, _ = tuple(combo_value.split(':'))
        if combo not in methods:
            # Assign the default method for use in further tests
            methods[combo] = _default_opt_dict['method']
        function_or_thresh_given.append(combo)

    for combo, value in methods.items():
        if value != 'trigger_fit' and combo in function_or_thresh_given:
            # Function/Threshold given for combo not using trigger_fit method
            parser.error("--fit-function and/or --fit-threshold given for "
                         + combo + " which has method " + value)
        elif value == 'trigger_fit' and combo not in function_or_thresh_given:
            # Threshold not given for trigger_fit combo
            parser.error("Threshold required for combo " + combo)


def ifar_opt_to_far_limit(ifar_str):
    """
    Convert the string of an IFAR limit in years into a
    float FAR limit in Hz.

    Parameters
    ----------
    ifar_str: string
        Upper limit on IFAR in years. Zero indicates no upper limit

    """
    ifar_float = positive_float(ifar_str)

    far_hz = 0. if (ifar_float == 0.) else 1. / (lal.YRJUL_SI * ifar_float)

    return far_hz


def digest_significance_options(combo_keys, args):
    """
    Read in information from the significance option group and ensure
    it makes sense before putting into a dictionary

    Parameters
    ----------

    combo_keys: list of strings
        list of detector combinations for which options are needed

    args: parsed arguments
        from argparse ArgumentParser parse_args()

    Returns
    -------
    significance_dict: dictionary
        Dictionary containing method, threshold and function for trigger fits
        as appropriate, and any limit on FAR (Hz)
    """

    lists_to_unpack = [('method', args.far_calculation_method, str),
                       ('fit_function', args.fit_function, str),
                       ('fit_threshold', args.fit_threshold, float),
                       ('far_limit', args.limit_ifar, ifar_opt_to_far_limit)]

    significance_dict = {}
    # Set everything as a default to start with:
    for combo in combo_keys:
        significance_dict[combo] = copy.deepcopy(_default_opt_dict)

    # Unpack everything from the arguments into the dictionary
    for argument_key, arg_to_unpack, conv_func in lists_to_unpack:
        for combo_value in arg_to_unpack:
            combo, value = tuple(combo_value.split(':'))
            if combo not in significance_dict:
                # Allow options for detector combos that are not actually
                # used/required for a given job. Such options have
                # no effect, but emit a warning for (e.g.) diagnostic checks
                logger.warning("Key %s not used by this code, uses %s",
                               combo, combo_keys)
                significance_dict[combo] = copy.deepcopy(_default_opt_dict)
            significance_dict[combo][argument_key] = conv_func(value)

    return significance_dict


def apply_far_limit(far, significance_dict, combo=None):
    """
    Apply a FAR limit to events according to command line options.

    If far_limit in significance_dict is zero, no limit is applied.

    Parameters
    ----------
    far: numpy array
    significance_dict: dictionary
        Dictionary containing any limit on FAR (Hz), made by
        digest_significance_options
    active_combination: numpy.array or string
        Array of IFO combinations given as utf-8 encoded strings, or a string
        which defines the IFO combination for all events

    Returns
    -------
    far_out: numpy array
        FARs with far limit applied as appropriate

    """
    far_out = copy.deepcopy(far)
    if isinstance(combo, str):
        # Single IFO combo used
        if significance_dict[combo]['far_limit'] == 0:
            return far_out
        far_limit_str = f"{significance_dict[combo]['far_limit']:.3e}"
        logger.info("Applying FAR limit of %s to %s events",
                    far_limit_str, combo)
        far_out = np.maximum(far, significance_dict[combo]['far_limit'])
    else:
        # IFO combo supplied as an array, by e.g. pycbc_add_statmap
        # Need to check which events are in which IFO combo in order to
        # apply the right limit to each
        for ifo_combo in significance_dict:
            if significance_dict[ifo_combo]['far_limit'] == 0:
                continue
            far_limit_str = f"{significance_dict[ifo_combo]['far_limit']:.3e}"
            logger.info("Applying FAR limit of %s to %s events",
                        far_limit_str, ifo_combo)
            this_combo_idx = combo == ifo_combo.encode('utf-8')
            far_out[this_combo_idx] = np.maximum(
                far[this_combo_idx],
                significance_dict[ifo_combo]['far_limit']
            )
    return far_out
