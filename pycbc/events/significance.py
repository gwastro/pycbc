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

import numpy as np
from pycbc.events import trigger_fits as trstats
from pycbc.types.optparse import MultiDetOptionAction


def n_louder(back_stat, fore_stat, dec_facs, **kwargs):  # pylint:disable=unused-argument
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
        The cumulative array of background triggers
    fore_n_louder: numpy.ndarray
        The number of background triggers above each foreground trigger
    """
    sort = bstat.argsort()
    bstat = bstat[sort]
    dec = dec[sort]

    # calculate cumulative number of triggers louder than the trigger in
    # a given index. We need to subtract the decimation factor, as the cumsum
    # includes itself in the first sum (it is inclusive of the first value)
    n_bg_louder = dec[::-1].cumsum()[::-1] - dec

    # Determine how many values are louder than the foreground ones
    # We need to subtract one from the index, to be consistent with the definition
    # of n_bg_louder, as here we do want to include the background value at the
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

    fore_n_louder = n_bg_louder[idx]

    # For background, need to get the number of louder events back into
    # the original order
    unsort = sort.argsort()
    back_cum_num = n_bg_louder[unsort]

    return back_cum_num, fore_n_louder


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

    back_cnum[np.logical_not(fg_above)] = below_back_cnum[np.logical_not(fg_above)]
    fnlouder[np.logical_not(bg_above)] = below_fnlouder[np.logical_not(bg_above)]

    return back_cnum, fnlouder


_significance_meth_dict = {
    'trigger_fit': trig_fit,
    'n_louder': n_louder
}


def calculate_n_louder(method_dict, back_stat, fore_stat, dec_facs,
                       fit_func=None, fit_thresh=None):
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
                        default={},
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
    parser.add_argument("--fit-function", default='exponential',
                        choices=["exponential", "rayleigh", "power"],
                        help="Functional form for the statistic slope fit if "
                             "--far-calculation-method is 'trigger_fit'. "
                             "Given as combination:function pairs, i.e. "
                             "H1:exponential H1L1:n_louder H1L1V1:n_louder "
                             "for all combinations with "
                             "--far-calculation-method of 'trigger_fit'"
                             "Default = exponential for all")


def digest_significance_option_group(combo_keys, args, parser):
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
        This is just for being able to use the parser error

    Returns
    -------
    significance_dict: dictionary
        Dictionary containing method, threshold and function for trigger fits
        ass appropriate
    """

    significance_dict = {}
    # There should be an entry into the --far-calculation-method list for
    # each detector combination:
    for key_method in args.far_calculation_method:
        key, method = key_method.split(':')
        if key in significance_dict:
            raise parser.error("key %s is duplicated", key)
        significance_dict[key] = {}
        significance_dict[key]['method'] = method

    # Apply the default n_louder to combinations not already filled:
    for combo_key in combo_keys:
        if combo_key in significance_dict: continue
        significance_dict[key] = {}
        significance_dict[key]['method'] = 'n_louder'
        significance_dict[key]['threshold'] = None
        significance_dict[key]['function'] = None

    # Grab the fit threshold for each key:
    for key_thresh in args.fit_threshold:
        key, thresh = key_thresh.split(':')
        if significance_dict[key]['method'] == 'n_louder':
            logging.warn("Fit threshold given for detector "
                         "combination %s which has method "
                         "other than trigger_fit. IGNORING",
                         key)
            continue
        significance_dict[key]['threshold'] = thresh

    # Grab the fit function for each key:
    for key_function in args.fit_function:
        key, function = key_function.split(':')
        if not significance_dict[key]['method'] == 'trigger_fit':
            logging.warn("Fit function given for detector "
                         "combination %s which has method "
                         "other than trigger_fit. IGNORING",
                          key)
            continue
        if 'threshold' not in significance_dict[key]:
            parser.error("Function given without threshold for %s", key)
        significance_dict[key]['function'] = function

    # Use the default exponential trigger fit where function not given
    for combo_key in combo_keys:
        if not significance_dict[key]['method'] == 'trigger_fit': continue
        if 'function' not in significance_dict[key]:
            significance_dict[key]['function'] = 'exponential'

    # Check that the threshold and function were given for all keys where the
    # method is trigger_fit
    for combo_key in combo_keys:
        if not significance_dict[key]['method'] == 'trigger_fit': continue
        if 'threshold' not in significance_dict[key]:
            parser.error("No threshold given for %s", key)

    return significance_dict
