# Copyright (C) 2014 Alex Nitz, Andrew Miller, Tito Dal Canton
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
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
import copy
import logging
from pycbc.psd.read import *
from pycbc.psd.analytical import *
from pycbc.psd.analytical_space import *
from pycbc.psd.estimate import *
from pycbc.psd.variation import *
from pycbc.types import float32,float64
from pycbc.types import MultiDetOptionAppendAction, MultiDetOptionAction
from pycbc.types import DictOptionAction, MultiDetDictOptionAction
from pycbc.types import copy_opts_for_single_ifo
from pycbc.types import required_opts, required_opts_multi_ifo
from pycbc.types import ensure_one_opt, ensure_one_opt_multi_ifo

def from_cli(opt, length, delta_f, low_frequency_cutoff,
             strain=None, dyn_range_factor=1, precision=None):
    """Parses the CLI options related to the noise PSD and returns a
    FrequencySeries with the corresponding PSD. If necessary, the PSD is
    linearly interpolated to achieve the resolution specified in the CLI.

    Parameters
    ----------
    opt : object
        Result of parsing the CLI with OptionParser, or any object with the
        required attributes (psd_model, psd_file, asd_file, psd_estimation,
        psd_segment_length, psd_segment_stride, psd_inverse_length,
        psd_output).
    length : int
        The length in samples of the output PSD.
    delta_f : float
        The frequency step of the output PSD in hertz.
    low_frequency_cutoff: float
        The low frequency cutoff to use when calculating the PSD, in hertz. If
        the option `psd-low-frequency-cutoff` is supplied in `opt`, that value
        is used for the cutoff instead.
    strain : {None, TimeSeries}
        Time series containing the data from which the PSD should be measured,
        when psd_estimation is in use.
    dyn_range_factor : {1, float}
        For PSDs taken from models or text files, if `dyn_range_factor` is
        not None, then the PSD is multiplied by `dyn_range_factor` ** 2.
    precision : str, choices (None,'single','double')
        If not specified, or specified as None, the precision of the returned
        PSD will match the precision of the data, if measuring a PSD, or will
        match the default precision of the model if using an analytical PSD.
        If 'single' the PSD will be converted to float32, if not already in
        that precision. If 'double' the PSD will be converted to float64, if
        not already in that precision.

    Returns
    -------
    psd : FrequencySeries
        The frequency series containing the PSD.
    """
    if opt.psd_low_frequency_cutoff is not None:
        f_low = opt.psd_low_frequency_cutoff
    else:
        f_low = low_frequency_cutoff
    sample_rate = (length -1) * 2 * delta_f

    try:
        psd_estimation = opt.psd_estimation is not None
    except AttributeError:
        psd_estimation = False

    exclusive_opts = [opt.psd_model, opt.psd_file, opt.asd_file,
                      psd_estimation]
    if sum(map(bool, exclusive_opts)) != 1:
        err_msg = "You must specify exactly one of '--psd-file', "
        err_msg += "'--psd-model', '--asd-file', '--psd-estimation'"
        raise ValueError(err_msg)

    if (opt.psd_model or opt.psd_file or opt.asd_file):
        # PSD from lalsimulation or file
        if opt.psd_model:
            psd = from_string(opt.psd_model, length, delta_f, f_low,
                              **opt.psd_extra_args)
        elif opt.psd_file or opt.asd_file:
            if opt.asd_file:
                psd_file_name = opt.asd_file
            else:
                psd_file_name = opt.psd_file
            if psd_file_name.endswith(('.dat', '.txt')):
                is_asd_file = bool(opt.asd_file)
                psd = from_txt(psd_file_name, length,
                               delta_f, f_low, is_asd_file=is_asd_file)
            elif opt.asd_file:
                raise ValueError(
                    "ASD files must be in ASCII format (extension .dat or "
                    f".txt). Got {psd_file_name} instead"
                )
            elif psd_file_name.endswith(('.xml', '.xml.gz')):
                psd = from_xml(psd_file_name, length, delta_f, f_low,
                               ifo_string=opt.psd_file_xml_ifo_string,
                               root_name=opt.psd_file_xml_root_name)
        # Set values < flow to the value at flow (if flow > 0)
        kmin = int(f_low / psd.delta_f)
        if kmin > 0:
            psd[0:kmin] = psd[kmin]

        psd *= dyn_range_factor ** 2

    elif psd_estimation:
        # estimate PSD from data
        psd = welch(strain, avg_method=opt.psd_estimation,
                    seg_len=int(opt.psd_segment_length * sample_rate + 0.5),
                    seg_stride=int(opt.psd_segment_stride * sample_rate + 0.5),
                    num_segments=opt.psd_num_segments,
                    require_exact_data_fit=False)

        if delta_f != psd.delta_f:
            psd = interpolate(psd, delta_f, length)

    else:
        # Shouldn't be possible to get here
        raise RuntimeError("Shouldn't be possible to raise this!")

    if opt.psd_inverse_length:
        try:
            which_spectrum = opt.invpsd_trunc_which_spectrum
        except AttributeError:
            which_spectrum = 'invasd'
        try:
            fill_value = opt.invpsd_trunc_low_freq_fill_value
        except AttributeError:
            fill_value = 0.
        psd = inverse_spectrum_truncation(psd, 
            int(opt.psd_inverse_length * sample_rate),
            which_spectrum=which_spectrum,
            low_frequency_cutoff=f_low,
            low_frequency_fill_value=fill_value,
            trunc_method=opt.invpsd_trunc_method)

    if hasattr(opt, 'psd_output') and opt.psd_output:
        (psd.astype(float64) / (dyn_range_factor ** 2)).save(opt.psd_output)

    if precision is None:
        return psd
    elif precision == 'single':
        return psd.astype(float32)
    elif precision == 'double':
        return psd.astype(float64)
    raise ValueError(
        "If provided, the precision kwarg must be either 'single' or "
        f"'double'. You provided {precision}"
    )

def from_cli_single_ifo(opt, length, delta_f, low_frequency_cutoff, ifo,
             **kwargs):
    """
    Get the PSD for a single ifo when using the multi-detector CLI
    """
    single_det_opt = copy_opts_for_single_ifo(opt, ifo)
    return from_cli(single_det_opt, length, delta_f, low_frequency_cutoff,
                    **kwargs)

def from_cli_multi_ifos(opt, length_dict, delta_f_dict,
                        low_frequency_cutoff_dict, ifos, strain_dict=None,
                        **kwargs):
    """
    Get the PSD for all ifos when using the multi-detector CLI
    """
    psd = {}
    for ifo in ifos:
        if strain_dict is not None:
            strain = strain_dict[ifo]
        else:
            strain = None
        psd[ifo] = from_cli_single_ifo(opt, length_dict[ifo], delta_f_dict[ifo],
                                       low_frequency_cutoff_dict[ifo], ifo,
                                       strain=strain, **kwargs)
    return psd

def insert_psd_option_group(parser, output=True, include_data_options=True):
    """
    Adds the options used to call the pycbc.psd.from_cli function to an
    optparser as an OptionGroup. This should be used if you
    want to use these options in your code.

    Parameters
    -----------
    parser : object
        OptionParser instance.
    """
    psd_options = parser.add_argument_group(
                          "Options to select the method of PSD generation",
                          "The options --psd-model, --psd-file, --asd-file, "
                          "and --psd-estimation are mutually exclusive.")
    psd_options.add_argument("--psd-model",
                             help="Get PSD from given analytical model. ",
                             choices=get_psd_model_list())
    psd_options.add_argument("--psd-extra-args",
                             nargs='+', action=DictOptionAction,
                             metavar='PARAM:VALUE', default={}, type=float,
                             help="(optional) Extra arguments passed to "
                             "the PSD models.")
    psd_options.add_argument("--psd-file",
                             help="Get PSD using given PSD ASCII file")
    psd_options.add_argument("--asd-file",
                             help="Get PSD using given ASD ASCII file")
    psd_options.add_argument("--psd-inverse-length", type=float,
                             help="(Optional) The maximum length of the "
                             "impulse response of the overwhitening "
                             "filter (s)")
    psd_options.add_argument("--psd-low-frequency-cutoff", type=float,
                             help="(Optional) The low frequency cutoff for the "
                                  "PSD. If not specified, the low frequency "
                                  "cutoff of the matched filter/likelihood "
                                  "integral will be used.")
    # Truncation options if specifying psd-inverse-length
    psd_options.add_argument("--invpsd-trunc-method", choices=["hann"],
                             help="(Optional) What truncation method to use "
                                  "when applying psd-inverse-length. If not "
                                  "provided, a hard truncation will be used.")
    psd_options.add_argument("--invpsd-trunc-which-spectrum", type=str,
                             default='invasd', choices=["invasd", "invpsd"],
                             help="(Optional) Which spectrum to use to perform "
                                  "truncation when applying psd-inverse-length. "
                                  "The default ('invasd') truncates the inverse "
                                  "ASD, while 'invpsd' truncates the inverse "
                                  "PSD.")
    psd_options.add_argument("--invpsd-trunc-low-freq-fill-value", default=0.,
                             help="(Optional) Value to set the inverse PSD to "
                                  "for frequencies below the low frequency "
                                  "cutoff when applying psd-inverse-length. "
                                  "Default 0; accepts any float. Also accepts "
                                  "'fmin'; this sets values below the cutoff "
                                  "to the inverse PSD value specified by "
                                  "`psd-low-frequency-cutoff` (or the matched "
                                  "filter/likelihood integral otherwise).")
    # Options specific to XML PSD files
    psd_options.add_argument("--psd-file-xml-ifo-string",
                             help="If using an XML PSD file, use the PSD in "
                                  "the file's PSD dictionary with this "
                                  "ifo string. If not given and only one "
                                  "PSD present in the file return that, if "
                                  "not given and multiple (or zero) PSDs "
                                  "present an exception will be raised.")
    psd_options.add_argument("--psd-file-xml-root-name", default='psd',
                             help="If given use this as the root name for "
                                  "the PSD XML file. If this means nothing "
                                  "to you, then it is probably safe to "
                                  "ignore this option.")
    # Options for PSD variation
    psd_options.add_argument("--psdvar-segment", type=float,
                             metavar="SECONDS", help="Length of segment "
                             "for mean square calculation of PSD variation.")
    psd_options.add_argument("--psdvar-short-segment", type=float,
                             metavar="SECONDS", help="Length of short segment "
                             "for outliers removal in PSD variability "
                             "calculation.")
    psd_options.add_argument("--psdvar-long-segment", type=float,
                             metavar="SECONDS", help="Length of long segment "
                             "when calculating the PSD variability.")
    psd_options.add_argument("--psdvar-psd-duration", type=float,
                             metavar="SECONDS", help="Duration of short "
                             "segments for PSD estimation.")
    psd_options.add_argument("--psdvar-psd-stride", type=float,
                             metavar="SECONDS", help="Separation between PSD "
                             "estimation segments.")
    psd_options.add_argument("--psdvar-low-freq", type=float, metavar="HERTZ",
                             help="Minimum frequency to consider in strain "
                             "bandpass.")
    psd_options.add_argument("--psdvar-high-freq", type=float, metavar="HERTZ",
                             help="Maximum frequency to consider in strain "
                             "bandpass.")

    if include_data_options :
        psd_options.add_argument("--psd-estimation",
                                 help="Measure PSD from the data, using "
                                 "given average method.",
                                 choices=["mean", "median", "median-mean"])
        psd_options.add_argument("--psd-segment-length", type=float,
                                 help="(Required for --psd-estimation) The "
                                 "segment length for PSD estimation (s)")
        psd_options.add_argument("--psd-segment-stride", type=float,
                                 help="(Required for --psd-estimation) "
                                 "The separation between consecutive "
                                 "segments (s)")
        psd_options.add_argument("--psd-num-segments", type=int,
                                 help="(Optional, used only with "
                                 "--psd-estimation). If given, PSDs will "
                                 "be estimated using only this number of "
                                 "segments. If more data is given than "
                                 "needed to make this number of segments "
                                 "then excess data will not be used in "
                                 "the PSD estimate. If not enough data "
                                 "is given, the code will fail.")
    if output:
        psd_options.add_argument("--psd-output",
                          help="(Optional) Write PSD to specified file")

    return psd_options

def insert_psd_option_group_multi_ifo(parser):
    """
    Adds the options used to call the pycbc.psd.from_cli function to an
    optparser as an OptionGroup. This should be used if you
    want to use these options in your code.

    Parameters
    -----------
    parser : object
        OptionParser instance.
    """
    psd_options = parser.add_argument_group(
                          "Options to select the method of PSD generation",
                          "The options --psd-model, --psd-file, --asd-file, "
                          "and --psd-estimation are mutually exclusive.")
    psd_options.add_argument("--psd-model", nargs="+",
                          action=MultiDetOptionAction, metavar='IFO:MODEL',
                          help="Get PSD from given analytical model. "
                          "Choose from %s" %(', '.join(get_psd_model_list()),))
    psd_options.add_argument("--psd-extra-args",
                             nargs='+', action=MultiDetDictOptionAction,
                             metavar='DETECTOR:PARAM:VALUE', default={},
                             type=float, help="(optional) Extra arguments "
                             "passed to the PSD models.")
    psd_options.add_argument("--psd-file", nargs="+",
                          action=MultiDetOptionAction, metavar='IFO:FILE',
                          help="Get PSD using given PSD ASCII file")
    psd_options.add_argument("--asd-file", nargs="+",
                          action=MultiDetOptionAction, metavar='IFO:FILE',
                          help="Get PSD using given ASD ASCII file")
    psd_options.add_argument("--psd-estimation", nargs="+",
                          action=MultiDetOptionAction, metavar='IFO:FILE',
                          help="Measure PSD from the data, using given "
                          "average method. Choose from "
                          "mean, median or median-mean.")
    psd_options.add_argument("--psd-low-frequency-cutoff", nargs="+", type=float,
                             action=MultiDetOptionAction, metavar='IFO:FREQ',
                             help="(Optional) The low frequency cutoff for the "
                                  "PSD. If not specified, the low frequency "
                                  "cutoff of the matched filter/likelihood "
                                  "integral will be used.")
    psd_options.add_argument("--psd-segment-length", type=float, nargs="+",
                          action=MultiDetOptionAction, metavar='IFO:LENGTH',
                          help="(Required for --psd-estimation) The segment "
                               "length for PSD estimation (s)")
    psd_options.add_argument("--psd-segment-stride", type=float, nargs="+",
                          action=MultiDetOptionAction, metavar='IFO:STRIDE',
                          help="(Required for --psd-estimation) The separation"
                               " between consecutive segments (s)")
    psd_options.add_argument("--psd-num-segments", type=int, nargs="+",
                          action=MultiDetOptionAction, metavar='IFO:NUM',
                          help="(Optional, used only with --psd-estimation). "
                               "If given PSDs will be estimated using only "
                               "this number of segments. If more data is "
                               "given than needed to make this number of "
                               "segments than excess data will not be used in "
                               "the PSD estimate. If not enough data is given "
                               "the code will fail.")
    psd_options.add_argument("--psd-inverse-length", type=float, nargs="+",
                          action=MultiDetOptionAction, metavar='IFO:LENGTH',
                          help="(Optional) The maximum length of the impulse"
                          " response of the overwhitening filter (s)")
    psd_options.add_argument("--invpsd-trunc-method", choices=["hann"],
                             help="(Optional) What truncation method to use "
                                  "when applying psd-inverse-length. If not "
                                  "provided, a hard truncation will be used.")
    psd_options.add_argument("--invpsd-trunc-which-spectrum", type=str,
                             default='invasd', choices=["invasd", "invpsd"],
                             help="(Optional) Which spectrum to use to perform "
                                  "truncation when applying psd-inverse-length. "
                                  "The default ('invasd') truncates the inverse "
                                  "ASD, while 'invpsd' truncates the inverse "
                                  "PSD.")
    psd_options.add_argument("--invpsd-trunc-low-freq-fill-value", default=0.,
                             action=MultiDetOptionAction, metavar='IFO:VALUE',
                             nargs="+",
                             help="(Optional) Value to set the inverse PSD to "
                                  "for frequencies below the low frequency "
                                  "cutoff when applying psd-inverse-length. "
                                  "Default 0; accepts any float. Also accepts "
                                  "'fmin'; this sets values below the cutoff "
                                  "to the inverse PSD value specified by "
                                  "`psd-low-frequency-cutoff` (or the matched "
                                  "filter/likelihood integral otherwise).")
    psd_options.add_argument("--psd-output", nargs="+",
                          action=MultiDetOptionAction, metavar='IFO:FILE',
                          help="(Optional) Write PSD to specified file")

    # Options for PSD variation
    psd_options.add_argument("--psdvar-segment", type=float,
                             metavar="SECONDS", help="Length of segment "
                             "when calculating the PSD variability.")
    psd_options.add_argument("--psdvar-short-segment", type=float,
                             metavar="SECONDS", help="Length of short segment "
                             "for outliers removal in PSD variability "
                             "calculation.")
    psd_options.add_argument("--psdvar-long-segment", type=float,
                             metavar="SECONDS", help="Length of long segment "
                             "when calculating the PSD variability.")
    psd_options.add_argument("--psdvar-psd-duration", type=float,
                             metavar="SECONDS", help="Duration of short "
                             "segments for PSD estimation.")
    psd_options.add_argument("--psdvar-psd-stride", type=float,
                             metavar="SECONDS", help="Separation between PSD "
                             "estimation segments.")
    psd_options.add_argument("--psdvar-low-freq", type=float, metavar="HERTZ",
                             help="Minimum frequency to consider in strain "
                             "bandpass.")
    psd_options.add_argument("--psdvar-high-freq", type=float, metavar="HERTZ",
                             help="Maximum frequency to consider in strain "
                             "bandpass.")

    return psd_options

ensure_one_opt_groups = []
ensure_one_opt_groups.append(['--psd-file', '--psd-model',
                              '--psd-estimation', '--asd-file'])

def verify_psd_options(opt, parser):
    """Parses the CLI options and verifies that they are consistent and
    reasonable.

    Parameters
    ----------
    opt : object
        Result of parsing the CLI with OptionParser, or any object with the
        required attributes (psd_model, psd_file, asd_file, psd_estimation,
        psd_segment_length, psd_segment_stride, psd_inverse_length, psd_output).
    parser : object
        OptionParser instance.
    """
    try:
        psd_estimation = opt.psd_estimation is not None
    except AttributeError:
        psd_estimation = False

    for opt_group in ensure_one_opt_groups:
        ensure_one_opt(opt, parser, opt_group)

    if psd_estimation:
        required_opts(opt, parser,
                      ['--psd-segment-stride', '--psd-segment-length'],
                      required_by = "--psd-estimation")

def verify_psd_options_multi_ifo(opt, parser, ifos):
    """Parses the CLI options and verifies that they are consistent and
    reasonable.

    Parameters
    ----------
    opt : object
        Result of parsing the CLI with OptionParser, or any object with the
        required attributes (psd_model, psd_file, asd_file, psd_estimation,
        psd_segment_length, psd_segment_stride, psd_inverse_length, psd_output).
    parser : object
        OptionParser instance.
    """
    for ifo in ifos:
        for opt_group in ensure_one_opt_groups:
            ensure_one_opt_multi_ifo(opt, parser, ifo, opt_group)

        if opt.psd_estimation[ifo]:
            required_opts_multi_ifo(opt, parser, ifo,
                      ['--psd-segment-stride', '--psd-segment-length'],
                          required_by = "--psd-estimation")


def psd_estimation_data_length(opt, sample_rate, input_data_len):
    """Number of samples of data used for one Welch PSD estimation"""
    seg_stride = int(opt.psd_segment_stride * sample_rate)
    seg_len = int(opt.psd_segment_length * sample_rate)
    if opt.psd_num_segments is None:
        num_segments = int(input_data_len // seg_stride) - 1
    else:
        num_segments = int(opt.psd_num_segments)
    return (num_segments - 1) * seg_stride + seg_len


def generate_segment_psds(opt, gwstrain, analysis_segments, flen, delta_f, flow,
                          dyn_range_factor=1., precision=None):
    """Estimate a PSD for each analysis segment.

    For ``--psd-estimation`` the PSD for a segment is measured from a stretch of
    strain data placed, and slid wholly inside the available data, so that it
    overlaps as much as possible of the data the matched filter uses for that
    segment: centred on the analysed span when shorter than it, centred on the
    whole segment when at least as long as it, and in between covering the
    analysed span plus as much of the surrounding segment data as fits,
    preferring data before the analysed span over data after it.

    For a static PSD (``--psd-model`` / ``--psd-file`` / ``--asd-file``) the
    one PSD is returned for every segment.

    Parameters
    -----------
    opt : object
        Result of parsing the CLI, or any object with the attributes needed by
        ``from_cli`` (single-ifo).
    gwstrain : TimeSeries
        The timeseries of data on which to estimate PSDs.
    analysis_segments : list of (seg_start, seg_stop, ana_start, ana_stop)
        Sample indices into ``gwstrain`` (see ``StrainSegments.psd_segments``):
        ``seg_start``/``seg_stop`` bound the whole segment the matched filter
        uses; ``ana_start``/``ana_stop`` bound the portion analysed for
        triggers. Both are used to place each segment's PSD estimation stretch.
    flen : int
        The length in samples of the output PSDs.
    delta_f : float
        The frequency step of the output PSDs in hertz.
    flow: float
        The low frequency cutoff to use when calculating the PSD, in hertz.
    dyn_range_factor : {1, float}
        For PSDs taken from models or text files, the PSD is multiplied by
        ``dyn_range_factor ** 2`` (default 1).
    precision : str, choices (None,'single','double')
        Precision of the returned PSDs (see ``from_cli``).

    Returns
    --------
    psds_and_times : list of (start, end, PSD) tuples
        One entry per input segment, in the same order. ``start`` and ``end``
        are the ``gwstrain`` sample range used to estimate that segment's PSD
        (shared between segments that resolve to the same range).
    """
    if not opt.psd_estimation:
        psd = from_cli(opt, flen, delta_f, flow, strain=gwstrain,
                       dyn_range_factor=dyn_range_factor, precision=precision)
        return [(0, len(gwstrain), psd) for _ in analysis_segments]

    sample_rate = gwstrain.sample_rate
    input_data_len = len(gwstrain)
    psd_data_len = psd_estimation_data_length(opt, sample_rate, input_data_len)

    if input_data_len < psd_data_len:
        raise ValueError(
            "Input data (%.0fs) is shorter than the data needed to estimate a "
            "PSD ((psd-num-segments - 1) * psd-segment-stride + "
            "psd-segment-length = %.0fs)."
            % (input_data_len / sample_rate, psd_data_len / sample_rate))

    cache = {}
    psds_and_times = []
    uncovered = []
    for seg_start, seg_stop, ana_start, ana_stop in analysis_segments:
        # Data the matched filter actually uses for this segment, excluding any
        # zero-padding that falls outside the available strain.
        data_start = max(seg_start, 0)
        data_stop = min(seg_stop, input_data_len)
        ana_len = ana_stop - ana_start
        data_len = data_stop - data_start

        # Position the psd_data_len-sample estimation stretch:
        #  - shorter than the analysed span: centre it in the analysed span;
        #  - between the analysed span and the whole segment: cover the
        #    analysed span and spend the extra on segment data the filter uses,
        #    preferring data before the analysed span over data after it;
        #  - at least the whole segment: centre it on the segment.
        if psd_data_len <= ana_len:
            psd_start = (ana_start + ana_stop - psd_data_len) // 2
        elif psd_data_len < data_len:
            slack = psd_data_len - ana_len
            psd_start = ana_start - min(slack, ana_start - data_start)
        else:
            psd_start = (data_start + data_stop - psd_data_len) // 2
        psd_start = min(max(psd_start, 0), input_data_len - psd_data_len)
        psd_stop = psd_start + psd_data_len

        if psd_start > data_start or psd_stop < data_stop:
            uncovered.append((data_start, data_stop, psd_start, psd_stop))

        if (psd_start, psd_stop) not in cache:
            cache[(psd_start, psd_stop)] = from_cli(
                opt, flen, delta_f, flow, strain=gwstrain[psd_start:psd_stop],
                dyn_range_factor=dyn_range_factor, precision=precision)
        psds_and_times.append((psd_start, psd_stop, cache[(psd_start, psd_stop)]))

    if uncovered:
        epoch = float(gwstrain.start_time)
        a0, a1, p0, p1 = uncovered[0]
        logging.warning(
            "%d of %d analysis segment(s) are not fully covered by the data "
            "used to estimate their PSD, so the matched-filter SNR near the "
            "uncovered edge(s) will be slightly over-estimated. Increase "
            "--psd-num-segments (or analyse more data) so that "
            "(psd-num-segments - 1) * psd-segment-stride + psd-segment-length "
            ">= the segment length. First affected segment: data used spans "
            "[%.1f, %.1f] but its PSD was estimated from [%.1f, %.1f].",
            len(uncovered), len(analysis_segments),
            epoch + a0 / sample_rate, epoch + a1 / sample_rate,
            epoch + p0 / sample_rate, epoch + p1 / sample_rate)

    return psds_and_times

def associate_psds_to_segments(opt, fd_segments, gwstrain, flen, delta_f, flow,
                               dyn_range_factor=1., precision=None):
    """Measure a PSD for every analysis segment and store it on the segment.

    With ``--psd-model``, ``--psd-file`` or ``--asd-file`` a single PSD is used
    for every segment. With ``--psd-estimation`` each segment's PSD is measured
    from a stretch of ``gwstrain`` placed to overlap as much as possible of the
    data the matched filter uses for that segment (see ``generate_segment_psds``),
    so the data being analysed is always part of the data used to estimate its
    PSD.

    Parameters
    -----------
    opt : object
        Result of parsing the CLI with OptionParser, or any object with the
        required attributes (psd_model, psd_file, asd_file, psd_estimation,
        psd_segment_length, psd_segment_stride, psd_inverse_length, psd_output).
    fd_segments : StrainSegments.fourier_segments() object
        The fourier transforms of the various analysis segments. The psd
        attribute of each segment is updated to point to the appropriate PSD.
    gwstrain : Strain object
        The timeseries of raw data on which to estimate PSDs.
    flen : int
        The length in samples of the output PSDs.
    delta_f : float
        The frequency step of the output PSDs.
    flow: float
        The low frequency cutoff to use when calculating the PSD, in hertz.
    dyn_range_factor : {1, float}
        For PSDs taken from models or text files, if `dyn_range_factor` is
        not None, then the PSD is multiplied by `dyn_range_factor` ** 2.
    precision : str, choices (None,'single','double')
        If not specified, or specified as None, the precision of the returned
        PSD will match the precision of the data, if measuring a PSD, or will
        match the default precision of the model if using an analytical PSD.
        If 'single' the PSD will be converted to float32, if not already in
        that precision. If 'double' the PSD will be converted to float64, if
        not already in that precision.
    """
    analysis_segments = [fd_segment.psd_seg_bounds for fd_segment in fd_segments]
    psds_and_times = generate_segment_psds(
        opt, gwstrain, analysis_segments, flen, delta_f, flow,
        dyn_range_factor=dyn_range_factor, precision=precision)
    for fd_segment, (_, _, psd) in zip(fd_segments, psds_and_times):
        fd_segment.psd = psd

def associate_psds_to_single_ifo_segments(opt, fd_segments, gwstrain, flen,
                                          delta_f, flow, ifo,
                                          dyn_range_factor=1., precision=None):
    """
    Associate PSDs to segments for a single ifo when using the multi-detector
    CLI
    """
    single_det_opt = copy_opts_for_single_ifo(opt, ifo)
    associate_psds_to_segments(single_det_opt, fd_segments, gwstrain, flen,
                               delta_f, flow, dyn_range_factor=dyn_range_factor,
                               precision=precision)

def associate_psds_to_multi_ifo_segments(opt, fd_segments, gwstrain, flen,
                                         delta_f, flow, ifos,
                                         dyn_range_factor=1., precision=None):
    """
    Associate PSDs to segments for all ifos when using the multi-detector CLI
    """
    for ifo in ifos:
        if gwstrain is not None:
            strain = gwstrain[ifo]
        else:
            strain = None

        if fd_segments is not None:
            segments = fd_segments[ifo]
        else:
            segments = None

        associate_psds_to_single_ifo_segments(opt, segments, strain, flen,
                delta_f, flow, ifo, dyn_range_factor=dyn_range_factor,
                precision=precision)
