# Copyright (C) 2017 Christopher M. Biwer
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
""" This modules contains functions for reading single and coincident triggers
from the command line.
"""

import h5py
import numpy
from pycbc import types
from pycbc.events import coinc
from pycbc.io import hdf

def insert_bank_bins_option_group(parser):
    """ Add options to the optparser object for selecting templates in bins.

    Parameters
    -----------
    parser : object
        OptionParser instance.
    """
    bins_group = parser.add_argument_group(
                                 "Options for selecting templates in bins.")
    bins_group.add_argument("--bank-bins", nargs="+", default=None,
                            help="Ordered list of mass bin upper boundaries. "
                                 "An ordered list of type-boundary pairs, "
                                 "applied sequentially. Must provide a name "
                                 "(can be any unique string for tagging "
                                 "purposes), the parameter to bin "
                                 "on, and the membership condition via "
                                 "'lt' / 'gt' operators. "
                                 "Ex. name1:component:lt2 name2:total:lt15")
    bins_group.add_argument("--bank-file", default=None,
                            help="HDF format template bank file.")
    bins_group.add_argument("--f-lower", default=None,
                            help="Low frequency cutoff in Hz.")
    return bins_group

def bank_bins_from_cli(opts):
    """ Parses the CLI options related to binning templates in the bank.

    Parameters
    ----------
    opts : object
        Result of parsing the CLI with OptionParser.

    Results
    -------
    bins_idx : dict
        A dict with bin names as key and an array of their indices as value.
    bank : dict
        A dict of the datasets from the bank file.
    """
    bank = {}
    fp = h5py.File(opts.bank_file)
    for key in fp.keys():
        bank[key] = fp[key][:]
    bank["f_lower"] = float(opts.f_lower) if opts.f_lower else None
    if opts.bank_bins:
        bins_idx = coinc.background_bin_from_string(opts.bank_bins, bank)
    else:
        bins_idx = {"all" : numpy.arange(0, len(bank[fp.keys()[0]]))}
    fp.close()
    return bins_idx, bank

def insert_loudest_triggers_option_group(parser, coinc_options=True):
    """ Add options to the optparser object for selecting templates in bins.

    Parameters
    -----------
    parser : object
        OptionParser instance.
    """
    opt_group = insert_bank_bins_option_group(parser)
    opt_group.title = "Options for finding loudest triggers."
    if coinc_options:
        opt_group.add_argument("--statmap-file", default=None,
                               help="HDF format clustered coincident trigger "
                                    "result file.")
        opt_group.add_argument("--statmap-group", default="foreground",
                               help="Name of group in statmap file to "
                                    "get triggers.")
    opt_group.add_argument("--sngl-trigger-files", nargs="+", default=None,
                           action=types.MultiDetOptionAction,
                           help="HDF format merged single detector "
                                "trigger files.")
    opt_group.add_argument("--veto-file", default=None,
                           help="XML file with segment_definer and "
                                "segment table.")
    opt_group.add_argument("--veto-segment-name", default=None,
                           help="Name of segment to use as veto in "
                                 "XML file's segment_definer table.")
    opt_group.add_argument("--search-n-loudest", type=int, default=None,
                           help="Number of triggers to search over.")
    opt_group.add_argument("--n-loudest", type=int, default=None,
                           help="Number of triggers to return in results.")
    return opt_group

def loudest_triggers_from_cli(opts, coinc_parameters=None,
                              sngl_parameters=None, bank_parameters=None):
    """ Parses the CLI options related to find the loudest coincident or
    single detector triggers.

    Parameters
    ----------
    opts : object
        Result of parsing the CLI with OptionParser.
    coinc_parameters : list
        List of datasets in statmap file to retrieve.
    sngl_parameters : list
        List of datasets in single-detector trigger files to retrieve.
    bank_parameters : list
        List of datasets in template bank file to retrieve.

    Results
    -------
    bin_names : dict
        A list of bin names.
    bin_results : dict
        A list of dict holding trigger data data.
    """

    # list to hold trigger data
    bin_results = []

    # list of IFOs
    ifos = opts.sngl_trigger_files.keys()

    # get indices of bins in template bank
    bins_idx, bank_data = bank_bins_from_cli(opts)
    bin_names = bins_idx.keys()

    # if taking triggers from statmap file
    if opts.statmap_file and opts.bank_file and opts.sngl_trigger_files:

        # loop over each bin
        for bin_name in bin_names:
            data = {}

            # get template has and detection statistic for coincident events
            statmap = hdf.ForegroundTriggers(
                                 opts.statmap_file, opts.bank_file,
                                 sngl_files=opts.sngl_trigger_files.values(),
                                 n_loudest=opts.search_n_loudest,
                                 group=opts.statmap_group)
            template_hash = statmap.get_bankfile_array("template_hash")
            stat = statmap.get_coincfile_array("stat")

            # get indices of triggers in bin
            bin_idx = numpy.in1d(template_hash,
                              bank_data["template_hash"][bins_idx[bin_name]])

            # get indices for sorted detection statistic in bin
            sorting = stat[bin_idx].argsort()[::-1]

            # get variables for n-th loudest triggers
            for p in coinc_parameters:
                arr = statmap.get_coincfile_array(p)
                data[p] = arr[bin_idx][sorting][:opts.n_loudest]
            for p in sngl_parameters:
                for ifo in ifos:
                    key = "/".join([ifo, p])
                    arr = statmap.get_snglfile_array_dict(p)[ifo]
                    data[key] = arr[bin_idx][sorting][:opts.n_loudest]
            for p in bank_parameters:
                arr = statmap.get_bankfile_array(p)
                data[p] = arr[bin_idx][sorting][:opts.n_loudest]

            # append results
            bin_results.append(data)

    # if taking triggers from single detector file
    elif opts.bank_file and opts.sngl_trigger_files:

        # loop over each bin
        for bin_name in bin_names:
            data = {}

            # only use one IFO
            if len(opts.sngl_trigger_files.keys()) == 1:
                ifo = opts.sngl_trigger_files.keys()[0]
            else:
                raise ValueError("Too many IFOs")

            # get newSNR as statistic from single detector files
            sngls =  hdf.SingleDetTriggers(opts.sngl_trigger_files[ifo],
                                           opts.bank_file, opts.veto_file,
                                           opts.veto_segment_name, None, ifo)

            # cluster
            n_loudest = opts.search_n_loudest \
                           if opts.search_n_loudest else len(sngls.template_id)
            sngls.mask_to_n_loudest_clustered_events(n_loudest=n_loudest)
            template_hash = \
                  sngls.bank["template_hash"][:][sngls.template_id]

            # get indices of triggers in bin
            bin_idx = numpy.in1d(template_hash,
                              bank_data["template_hash"][bins_idx[bin_name]])

            # sort by detection statistic
            stats = sngls.stat
            sorting = stats[bin_idx].argsort()[::-1]

            # get indices for sorted detection statistic in bin
            for p in sngl_parameters:
                key = "/".join([ifo, p])
                arr = sngls.get_column(p)
                data[key] = arr[bin_idx][sorting][:opts.n_loudest]
            for p in bank_parameters:
                arr = sngls.bank[p][:]
                data[p] = \
                      arr[sngls.template_id][bin_idx][sorting][:opts.n_loudest]

            # append results
            bin_results.append(data)

    # else did not supply enough command line options
    else:
        raise ValueError("Must have --bank-file and --sngl-trigger-files")

    return bin_names, bin_results
