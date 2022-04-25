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
import logging, copy
import numpy as np
from pycbc import bin_utils
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
    # Manually set some defaults:
    fit_func = fit_func if fit_func else 'exponential'
    fit_thresh = fit_thresh if fit_thresh else 0

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
    lists_to_check = [(args.far_calculation_method, str,
                       _significance_meth_dict.keys()),
                      (args.fit_function, str,
                       trstats.fitalpha_dict.keys()),
                      (args.fit_threshold, float,
                       range(-10000, 10000))]
    # (10,000 is suitable for current use but can be changed if wanted)

    for list_to_check, type_to_convert, allowed_values in lists_to_check:
        key_list = []
        for key_value in list_to_check:
            try:
                key, value = tuple(key_value.split(':'))
            except ValueError:
                parser.error("Need key:value format, got %s" % key_value)

            if key in key_list:
                parser.error("Key %s duplicated in a significance option" % key)
            key_list.append(key)

            try:
                type_to_convert(value)
            except ValueError:
                parser.error("Value %s of key %s "
                             "cannot be converted as appropriate" % (value, key))

            if type_to_convert(value) not in allowed_values:
                parser.error("Value %s of key %s is not in allowed values: %s" %
                             (value, key, allowed_values))

    # Are the functions/thresholds appropriate for the methods given?
    methods = {}
    # A method has been specified
    for key_value in args.far_calculation_method:
        key, value = tuple(key_value.split(':'))
        methods[key] = value

    # A function or threshold has been specified
    function_or_thresh_given = []
    for key_value in args.fit_function + args.fit_threshold:
        key, _ = tuple(key_value.split(':'))
        if key not in methods:
            methods[key] = 'n_louder'
        function_or_thresh_given.append(key)

    for key, value in methods.items():
        if value != 'trigger_fit' and key in function_or_thresh_given:
            # Function/Threshold given for key not using trigger_fit method
            parser.error("--fit-function and/or --fit-threshold given for "
                         + key + " which has method " + value)
        elif value == 'trigger_fit' and key not in function_or_thresh_given:
            # Threshold not given for trigger_fit key
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

    # Set up the defaults
    default_dict = {'method': 'n_louder', 'threshold': None, 'function': None}

    # Unpack the string arguments into a standard-format dictionary
    significance_dict = {}

    # Set everything as a default to start with:
    for key in combo_keys:
        significance_dict[key] = copy.copy(default_dict)

    # Unpack everything from the arguments into the dictionary
    for unpack_key, arg_to_unpack, conv_func in lists_to_unpack:
        for key_value in arg_to_unpack:
            key, value = tuple(key_value.split(':'))
            if key not in significance_dict:
                # This is a newly added key, not actually used by the code
                # This is a warning not an exit, so we can reuse the same
                # settings for multiple jobs in a workflow. However we don't
                # just want to accept this silently
                logging.warning("Key %s not used by this code, uses %s",
                                key, combo_keys)
                significance_dict[key] = copy.copy(default_dict)
            significance_dict[key][unpack_key] = conv_func(value)

    return significance_dict
