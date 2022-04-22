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
import numpy as np
from pycbc.events import coinc
from pycbc.events import trigger_fits as trstats
from pycbc.types.optparse import MultiDetOptionAction

def n_louder(bstat, fstat, dec, **kwargs):  # pylint:disable=unused-argument
    """
    Count the number of louder background events, taking decimation
    into account
    """
    return coinc.calculate_n_louder(bstat, fstat, dec)


def trig_fit(back_stat, fore_stat, dec_facs, fit_func=None,
             fit_thresh=None, **kwargs):  # pylint:disable=unused-argument
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

    # estimate number of louder events according to the fit
    back_cnum[bg_above] = trstats.cum_fit(fit_func, back_stat[bg_above],
                                        alpha, fit_thresh) * bg_above_thresh
    fnlouder[fg_above] = trstats.cum_fit(fit_func, fore_stat[fg_above],
                                        alpha, fit_thresh) * bg_above_thresh

    # below the fit threshold, we count the number of louder events,
    # as things get complicated by clustering below this point
    below_back_cnum, below_fnlouder = n_louder(back_stat, fore_stat, dec_facs)

    back_cnum[np.logical_not(bg_above)] = below_back_cnum[np.logical_not(bg_above)]
    fnlouder[np.logical_not(fg_above)] = below_fnlouder[np.logical_not(fg_above)]

    return back_cnum, fnlouder


_significance_meth_dict = {
    'trigger_fit': trig_fit,
    'n_louder': n_louder
}


def calculate_n_louder(method_dict, back_stat, fore_stat, dec_facs):
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
                        help="Method used for FAR calculation in each "
                             "detector combination, given as "
                             "combination:method pairs, i.e. "
                             "H1:trigger_fit H1L1:n_louder H1L1V1:n_louder "
                             "etc. n_louder counts the rate of louder "
                             "events, trigger_fit fits the triggers to a "
                             "distribution and extrapolates. "
                             "Default = n_louder for all not given")
    parser.add_argument('--fit-threshold', nargs='+',
                        help="Thresholds for the fits to statistic "
                             "values for FAN approximation if "
                             "--far-calculation-method is 'trigger_fit'. "
                             "Given as combination-values pairs, e.g. "
                             "H1:0 L1:0 V1:-4, for all combinations which "
                             "have --far-calculation-method "
                             "of trigger_fit")
    parser.add_argument("--fit-function", nargs='+',
                        help="Functional form for the statistic slope fit if "
                             "--far-calculation-method is 'trigger_fit'. "
                             "Given as combination:function pairs, i.e. "
                             "H1:exponential H1L1:n_louder H1L1V1:n_louder "
                             "for all combinations with "
                             "--far-calculation-method of 'trigger_fit'"
                             "Default = exponential for all")


def digest_significance_options(combo_keys, args, parser):
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

    parser: argparse ArgumentParser instance
        just for the parser.error method

    Returns
    -------
    significance_dict: dictionary
        Dictionary containing method, threshold and function for trigger fits
        as appropriate
    """
    # First: check that the arguments can be read as lists:
    if args.far_calculation_method:
        calc_methods = args.far_calculation_method
    else:
        calc_methods = []
    fit_threshes = args.fit_threshold if args.fit_threshold else []
    fit_functions = args.fit_function if args.fit_function else []

    # Second: Check that the key:method/function/threshold are in the
    # right format, and are in allowed combinations
    for list_to_check in [calc_methods, fit_threshes, fit_functions]:
        for key_value in list_to_check:
            try:
                key, value = tuple(key_value.split(':'))
            except ValueError:
                parser.error("Need exactly one colon in argument %s" % key_value)
            if key not in combo_keys:
                # This is a warning not an exit, so we can reuse the same
                # settings for multiple jobs in workflow
                logging.warn("Key %s not in allowed list: %s", key, combo_keys)


    # Third: Unpack the arguments into a standard-format dictionary
    significance_dict = {}
    # Go through --far-calculation-method list
    for key_method in calc_methods:
        key, method = tuple(key_method.split(':'))
        if key in significance_dict:
            parser.error("--far-calculation-method key %s given already" % key)
        if method not in _significance_meth_dict:
            parser.error("Method %s is not possible" % method)
        significance_dict[key] = {}
        significance_dict[key]['method'] = method

    # Apply the default values to combinations not already filled:
    for key in combo_keys:
        if key in significance_dict: continue
        significance_dict[key] = {}
        significance_dict[key]['method'] = 'n_louder'
        significance_dict[key]['threshold'] = None
        significance_dict[key]['function'] = None

    # Grab the fit threshold for each key:
    for key_thresh in fit_threshes:
        key, thresh = tuple(key_thresh.split(':'))
        if significance_dict[key]['method'] == 'n_louder':
            parser.error("Fit threshold given for detector "
                         "combination %s which has method "
                         "other than trigger_fit" % key)
        significance_dict[key]['threshold'] = float(thresh)

    # Grab the fit function for each key:
    for key_function in fit_functions:
        key, function = tuple(key_function.split(':'))
        if key not in combo_keys:
            parser.error("key %s not in allowed list" % key)
        if not significance_dict[key]['method'] == 'trigger_fit':
            parser.error("Fit function given for detector "
                         "combination %s which has method "
                         "other than trigger_fit. IGNORING" % key)
        if 'threshold' not in significance_dict[key]:
            parser.error("Function given without threshold for %s" % key)
        significance_dict[key]['function'] = function

    # Use the default exponential trigger fit where function not given
    for key in significance_dict.keys():
        if not significance_dict[key]['method'] == 'trigger_fit': continue
        if 'function' not in significance_dict[key]:
            significance_dict[key]['function'] = 'exponential'

    # Check that the threshold and function were given for all keys where the
    # method is trigger_fit
    for key in significance_dict.keys():
        if not significance_dict[key]['method'] == 'trigger_fit': continue
        if 'threshold' not in significance_dict[key]:
            parser.error( "No threshold given for %s" % key)

    # If n_louder is specified, then set the default theshold/function values
    for key in significance_dict.keys():
        if 'threshold' not in significance_dict[key]:
            significance_dict[key]['threshold'] = None
        if 'function' not in significance_dict[key]:
            significance_dict[key]['function'] = None

    return significance_dict

