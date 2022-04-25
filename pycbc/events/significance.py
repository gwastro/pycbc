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


def n_louder(bstat, fstat, dec, **kwargs):  # pylint:disable=unused-argument
    """
    Count the number of louder background events, taking decimation
    into account
    """
    return coinc.calculate_n_louder(bstat, fstat, dec)


def trig_fit(back_stat, fore_stat, dec_facs, fit_func=None,
             fit_thresh=None):
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
                        default=[],
                        help="Method used for FAR calculation in each "
                             "detector combination, given as "
                             "combination:method pairs, i.e. "
                             "H1:trigger_fit H1L1:n_louder H1L1V1:n_louder "
                             "etc. Method options are ["
                             + ",".join(_significance_meth_dict.keys()) + 
                             "]. Default = n_louder for all not given")
    parser.add_argument('--fit-threshold', nargs='+', default=[],
                        help="Thresholds for the fits to statistic "
                             "values for FAN approximation if "
                             "--far-calculation-method is 'trigger_fit'. "
                             "Given as combination-values pairs, e.g. "
                             "H1:0 L1:0 V1:-4, for all combinations which "
                             "have --far-calculation-method "
                             "of trigger_fit")
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
    Check that the arguments given for the insert_significance_option_group
    options make sense
    """
    # Check that the key:method/function/threshold are in the
    # right format, and are in allowed combinations
    lists_to_check = [args.far_calculation_method,
                      args.fit_threshold,
                      args.fit_function]

    for list_to_check in lists_to_check:
        key_list = []
        for key_value in list_to_check:
            try:
                key, value = tuple(key_value.split(':'))
            except ValueError:
                parser.error("Need one colon in argument, got %s" % key_value)
            if key in key_list:
                parser.error("Key %s duplicated in significance option" % key)
            key_list.append(key)

    # Keep track of the methods used for later checks
    methods = {}
    for key_value in args.far_calculation_method:
        key, value = tuple(key_value.split(':'))
        if value not in _significance_meth_dict.keys():
            parser.error(("--far-calculation-method value %s for key %s "
                          "is not valid, choose from [" + 
                          ','.join(_significance_meth_dict.keys()) +
                          "].") % (value, key))
        methods[key] = value

    function_given = []
    for key_value in args.fit_function:
        key, value = tuple(key_value.split(':'))
        if value not in trstats.fitalpha_dict.keys():
            parser.error(("--fit-function value %s for key %s "
                          "is not valid, choose from [" +
                          ','.join(trstats.fitalpha_dict.keys()) +
                          "].") % (value, key))
        function_given.append(key)

    thresh_given = []
    for key_value in args.fit_threshold:
        # must be able to convert to a float
        key, value = tuple(key_value.split(':'))
        try:
            float(value)
        except ValueError:
            parser.error("--fit-threshold value %s for key %s "
                         "cannot be converted to a float" % (value, key))
        thresh_given.append(key)

    # Get places where the threshold or function are given
    function_or_thresh_given = set(function_given + thresh_given)

    for key in function_or_thresh_given:
        if key not in methods:
            methods[key] = 'n_louder'

    for key, value in methods.items():
        if value != 'trigger_fit' and key in function_or_thresh_given:
            # Function/Threshold given for key not using trigger_fit method
            parser.error("--fit-function or --fit-threshold given for key "
                         + key + " which has method " + value)
        elif value == 'trigger_fit' and key not in thresh_given:
            parser.error("Threshold required for key " + key)

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

    # Unpack the string arguments into a standard-format dictionary
    significance_dict = {}
    for unpack_key, arg_to_unpack, conv_func in lists_to_unpack:
        for key_value in arg_to_unpack:
            key, value = tuple(key_value.split(':'))
            if key not in combo_keys:
                # This is a warning not an exit, so we can reuse the same
                # settings for multiple jobs in a workflow. However we don't
                # just want to accept this silently
                logging.warning("Key %s not used by this code, uses %s",
                                key, combo_keys)
            if key not in significance_dict:
                significance_dict[key] = {}
            significance_dict[key][unpack_key] = conv_func(value)

    # Apply the default values to combinations not already filled:
    for key in list(significance_dict.keys()) + combo_keys:
        if key not in significance_dict:
            significance_dict[key] = {}
        # If method not given, then default is n_louder
        if 'method' not in significance_dict[key]:
            significance_dict[key]['method'] = 'n_louder'

        if significance_dict[key]['method'] == 'n_louder':
            # If method is n_louder, dont need threshold or function,
            # but they need to be given
            significance_dict[key]['threshold'] = None
            significance_dict[key]['function'] = None
        elif significance_dict[key]['method'] == 'trigger_fit':
            # Default function with trigger_fit method is exponential
            if 'function' not in significance_dict[key]:
                significance_dict[key]['function'] = 'exponential'

    return significance_dict
