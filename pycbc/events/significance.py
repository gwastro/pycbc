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
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

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
import sys
import numpy as np
from pycbc.events import coinc
from pycbc.events import trigger_fits as trstats


def n_louder(bstat, fstat, dec, skip_background=False, **kwargs):  # pylint:disable=unused-argument
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
    skip_background: optional, {boolean, False}
        Skip calculating cumulative numbers for background triggers

    Returns
    -------
    cum_back_num: numpy.ndarray
        The cumulative array of background triggers. Does not return this
        argument if skip_background == True
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
    # We need to subtract one from the index, to be consistent with the definition
    # of n_louder, as here we do want to include the background value at the
    # found index
    idx = np.searchsorted(bstat, fstat, side='left') - 1

    # If the foreground are *quieter* than the background or at the same value
    # then the search sorted algorithm will choose position -1, which does not exist
    # We force it back to zero.
    if isinstance(idx, np.ndarray):  # Case where our input is an array
        idx[idx < 0] = 0
    else:  # Case where our input is just a scalar value
        if idx < 0:
            idx = 0

    fore_n_louder = n_louder[idx]

    if not skip_background:
        unsort = sort.argsort()
        back_cum_num = n_louder[unsort]
        return back_cum_num, fore_n_louder
    else:
        return fore_n_louder


def trig_fit(back_stat, fore_stat, dec_facs, fit_func='exponential',
             fit_thresh=0):
    """
    Use a fit to events in back_stat in order to estimate the
    distribution for use in recovering the ifars. Below the
    fit threshold, use the n_louder method for these triggers
    """

    # Calculate the fitting factor of the ranking statistic distribution
    alpha, _ = trstats.fit_above_thresh(fit_func, back_stat,
                                        thresh=fit_thresh,
                                        weights=dec_facs)

    # Count background events above threshold as the cum_fit is
    # normalised to 1
    bg_above = back_stat > fit_thresh
    bg_above_thresh = np.sum(dec_facs[bg_above])
    fg_above = fore_stat > fit_thresh

    # These will be overwritten, but just to silence a warning
    # in the case where trstats.cum_fit returns zero
    back_cnum = np.zeros_like(back_stat)
    fnlouder = np.zeros_like(fore_stat)

    # Ue the fit above the threshold
    back_cnum[bg_above] = trstats.cum_fit(fit_func, back_stat[bg_above],
                                          alpha, fit_thresh) * bg_above_thresh
    fnlouder[fg_above] = trstats.cum_fit(fit_func, fore_stat[fg_above],
                                         alpha, fit_thresh) * bg_above_thresh

    # below the fit threshold, we count the number of louder events,
    # as things get complicated by clustering below this point
    fg_below = np.logical_not(fg_above)
    bg_below = np.logical_not(bg_above)

    back_cnum[bg_below], fnlouder[fg_below] = \
        n_louder(back_stat, fore_stat, dec_facs)

    return back_cnum, fnlouder


_significance_meth_dict = {
    'trigger_fit': trig_fit,
    'n_louder': n_louder
}


def get_n_louder(method_dict, back_stat, fore_stat, dec_facs):
    """
    Wrapper to find the correct n_louder calculation method using standard
    inputs
    """
    calculation_method = method_dict['method']
    fit_func = method_dict['function']
    fit_thresh = method_dict['threshold']
    return _significance_meth_dict[calculation_method](back_stat, fore_stat,
                                                       dec_facs,
                                                       fit_func=fit_func,
                                                       fit_thresh=fit_thresh)


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
                      (args.fit_threshold, float,
                       None)]

    for list_to_check, type_to_convert, allowed_values in lists_to_check:
        combo_list = []
        for combo_value in list_to_check:
            try:
                combo, value = tuple(combo_value.split(':'))
            except ValueError:
                parser.error("Need combo:value format, got %s" % combo_value)

            if combo in combo_list:
                parser.error("Duplicate combo %s in a significance option" % combo)
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
            methods[combo] = 'n_louder'
        function_or_thresh_given.append(combo)

    for combo, value in methods.items():
        if value != 'trigger_fit' and combo in function_or_thresh_given:
            # Function/Threshold given for combo not using trigger_fit method
            parser.error("--fit-function and/or --fit-threshold given for "
                         + combo + " which has method " + value)
        elif value == 'trigger_fit' and combo not in function_or_thresh_given:
            # Threshold not given for trigger_fit combo
            parser.error("Threshold required for combo " + combo)


def digest_significance_options(combo_keys, args):
    """
    Read in information from the significance option group and ensure
    it makes sense before putting into a dictionary

    Parameters
    ----------

    combo_keys: list of strings
        list of combinations of detectors which could be used as keys
        for the options

    args: parsed arguments
        from argparse ArgumentParser parse_args()

    Returns
    -------
    significance_dict: dictionary
        Dictionary containing method, threshold and function for trigger fits
        as appropriate
    """

    lists_to_unpack = [('method', args.far_calculation_method, str),
                       ('function', args.fit_function, str),
                       ('threshold', args.fit_threshold, float)]

    # Set up the defaults
    default_dict = {'method': 'n_louder', 'threshold': None, 'function': None}

    significance_dict = {}
    # Set everything as a default to start with:
    for combo in combo_keys:
        significance_dict[combo] = copy.deepcopy(default_dict)

    # Unpack everything from the arguments into the dictionary
    for argument_key, arg_to_unpack, conv_func in lists_to_unpack:
        for combo_value in arg_to_unpack:
            combo, value = tuple(combo_value.split(':'))
            if combo not in significance_dict:
                # This is a newly added combo, not actually used by the code
                # This is a warning not an exit, so we can reuse the same
                # settings for multiple jobs in a workflow. However we don't
                # just want to accept this silently
                logging.warning("Key %s not used by this code, uses %s",
                                combo, combo_keys)
                significance_dict[combo] = copy.deepcopy(default_dict)
            significance_dict[combo][argument_key] = conv_func(value)

    return significance_dict

