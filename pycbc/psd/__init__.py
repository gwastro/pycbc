#!/usr/bin/python
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
from glue import segments
from pycbc.psd.read import *
from pycbc.psd.analytical import *
from pycbc.psd.estimate import *
from pycbc.types import float32,float64
from pycbc.types import MultiDetOptionAppendAction, MultiDetOptionAction
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
        The frequency step of the output PSD.
    low_frequency_cutoff: float
        The low frequncy cutoff to use when calculating the PSD.
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
    f_low = low_frequency_cutoff
    sample_rate = int((length -1) * 2 * delta_f)

    if (opt.psd_model or opt.psd_file or opt.asd_file)\
                      and not opt.psd_estimation:
        # PSD from lalsimulation or file
        if opt.psd_model:
            psd = from_string(opt.psd_model, length, delta_f, f_low)
        elif opt.psd_file:
            psd = from_txt(opt.psd_file, length, 
                           delta_f, f_low, is_asd_file=False)
        elif opt.asd_file:
            psd = from_txt(opt.asd_file, length, 
                           delta_f, f_low, is_asd_file=True)
        # Set values < flow to the value at flow
        kmin = int(low_frequency_cutoff / psd.delta_f)
        psd[0:kmin] = psd[kmin]

        psd *= dyn_range_factor ** 2

    elif opt.psd_estimation and not (opt.psd_model or 
                                     opt.psd_file or opt.asd_file):
        # estimate PSD from data
        psd = welch(strain, avg_method=opt.psd_estimation,
                    seg_len=int(opt.psd_segment_length * sample_rate),
                    seg_stride=int(opt.psd_segment_stride * sample_rate),
                    num_segments=opt.psd_num_segments,
                    require_exact_data_fit=False)

        if delta_f != psd.delta_f:
            psd = interpolate(psd, delta_f)
    else:
        # no PSD options given
        return None

    if opt.psd_inverse_length:
        psd = inverse_spectrum_truncation(psd, 
            int(opt.psd_inverse_length * sample_rate),
            low_frequency_cutoff=f_low)

    if hasattr(opt, 'psd_output') and opt.psd_output:
        (psd.astype(float64) / (dyn_range_factor ** 2)).save(opt.psd_output)

    if precision is None:
        return psd
    elif precision == 'single':
        return psd.astype(float32)
    elif precision == 'double':
        return psd.astype(float64)
    else:
        err_msg = "If provided the precision kwarg must be either 'single' "
        err_msg += "or 'double'. You provided %s." %(precision)
        raise ValueError(err_msg)

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

def insert_psd_option_group(parser, output=True):
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
                          choices=get_lalsim_psd_list())
    psd_options.add_argument("--psd-file",
                          help="Get PSD using given PSD ASCII file")
    psd_options.add_argument("--asd-file",
                          help="Get PSD using given ASD ASCII file")
    psd_options.add_argument("--psd-estimation",
                          help="Measure PSD from the data, using given "
                          "average method.",
                          choices=["mean", "median", "median-mean"])
    psd_options.add_argument("--psd-segment-length", type=float, 
                          help="(Required for --psd-estimation) The segment "
                               "length for PSD estimation (s)")
    psd_options.add_argument("--psd-segment-stride", type=float, 
                          help="(Required for --psd-estimation) The separation"
                               " between consecutive segments (s)")
    psd_options.add_argument("--psd-num-segments", type=int, default=None,
                          help="(Optional, used only with --psd-estimation). "
                               "If given PSDs will be estimated using only "
                               "this number of segments. If more data is "
                               "given than needed to make this number of "
                               "segments than excess data will not be used in "
                               "the PSD estimate. If not enough data is given "
                               "the code will fail.")
    psd_options.add_argument("--psd-inverse-length", type=float, 
                          help="(Optional) The maximum length of the impulse"
                          " response of the overwhitening filter (s)")
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
                          "Choose from %s" %(', '.join(get_lalsim_psd_list()),))
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
    psd_options.add_argument("--psd-segment-length", type=float, nargs="+",
                          action=MultiDetOptionAction, metavar='IFO:LENGTH',
                          help="(Required for --psd-estimation) The segment "
                               "length for PSD estimation (s)")
    psd_options.add_argument("--psd-segment-stride", type=float, nargs="+",
                          action=MultiDetOptionAction, metavar='IFO:STRIDE',
                          help="(Required for --psd-estimation) The separation"
                               " between consecutive segments (s)")
    psd_options.add_argument("--psd-num-segments", type=int, nargs="+",
                          default=None,
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
    psd_options.add_argument("--psd-output", nargs="+",
                          action=MultiDetOptionAction, metavar='IFO:FILE',
                          help="(Optional) Write PSD to specified file")

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
    for opt_group in ensure_one_opt_groups:
        ensure_one_opt(opt, parser, opt_group)

    if opt.psd_estimation:
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

def generate_overlapping_psds(opt, gwstrain, flen, delta_f, flow,
                              dyn_range_factor=1., precision=None):
    """Generate a set of overlapping PSDs to cover a stretch of data. This
    allows one to analyse a long stretch of data with PSD measurements that
    change with time.

    Parameters
    -----------
    opt : object
        Result of parsing the CLI with OptionParser, or any object with the
        required attributes (psd_model, psd_file, asd_file, psd_estimation,
        psd_segment_length, psd_segment_stride, psd_inverse_length, psd_output).
    gwstrain : Strain object
        The timeseries of raw data on which to estimate PSDs.
    flen : int
        The length in samples of the output PSDs.
    delta_f : float
        The frequency step of the output PSDs.
    flow: float
        The low frequncy cutoff to use when calculating the PSD.
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
    --------
    psd_and_times : list of (start, end, psd) tuples
        This is a list of tuples containing one entry for each PSD. The first
        and second entries (start, end) in each tuple represent the index 
        range of the gwstrain data that was used to estimate that PSD. The
        third entry (psd) contains the PSD estimate between that interval.
    """
    if not opt.psd_estimation:
        psd = from_cli(opt, flen, delta_f, flow, strain=gwstrain,
                       dyn_range_factor=dyn_range_factor, precision=precision)
        psds_and_times = [ (0, len(gwstrain), psd) ]
        return psds_and_times

    # Figure out the data length used for PSD generation
    seg_stride = int(opt.psd_segment_stride * gwstrain.sample_rate)
    seg_len = int(opt.psd_segment_length * gwstrain.sample_rate)
    input_data_len = len(gwstrain)

    if opt.psd_num_segments is None:
        # FIXME: Should we make --psd-num-segments mandatory?
        #        err_msg = "You must supply --num-segments."
        #        raise ValueError(err_msg)
        num_segments = int(input_data_len // seg_stride) - 1
    else:
        num_segments = int(opt.psd_num_segments)

    psd_data_len = (num_segments - 1) * seg_stride + seg_len

    # How many unique PSD measurements is this?
    psds_and_times = []
    if input_data_len < psd_data_len:
        err_msg = "Input data length must be longer than data length needed "
        err_msg += "to estimate a PSD. You specified that a PSD should be "
        err_msg += "estimated with %d seconds. " %(psd_data_len)
        err_msg += "Input data length is %d seconds. " %(input_data_len)
        raise ValueError(err_msg)
    elif input_data_len == psd_data_len:
        num_psd_measurements = 1
        psd_stride = 0
    else:
        num_psd_measurements = int(2 * (input_data_len-1) / psd_data_len)
        psd_stride = int((input_data_len - psd_data_len) / num_psd_measurements)

    for idx in range(num_psd_measurements):
        if idx == (num_psd_measurements - 1):
            start_idx = input_data_len - psd_data_len
            end_idx = input_data_len
        else:
            start_idx = psd_stride * idx
            end_idx = psd_data_len + psd_stride * idx
        strain_part = gwstrain[start_idx:end_idx]
        psd = from_cli(opt, flen, delta_f, flow, strain=strain_part,
                       dyn_range_factor=dyn_range_factor, precision=precision)
        psds_and_times.append( (start_idx, end_idx, psd) )
    return psds_and_times

def associate_psds_to_segments(opt, fd_segments, gwstrain, flen, delta_f, flow,
                               dyn_range_factor=1., precision=None):
    """Generate a set of overlapping PSDs covering the data in GWstrain.
    Then associate these PSDs with the appropriate segment in strain_segments.

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
        The low frequncy cutoff to use when calculating the PSD.
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
    psds_and_times = generate_overlapping_psds(opt, gwstrain, flen, delta_f,
                                       flow, dyn_range_factor=dyn_range_factor,
                                       precision=precision)

    for fd_segment in fd_segments:
        best_psd = None
        psd_overlap = 0
        inp_seg = segments.segment(fd_segment.seg_slice.start,
                                   fd_segment.seg_slice.stop)
        for start_idx, end_idx, psd in psds_and_times:
            psd_seg = segments.segment(start_idx, end_idx)
            if psd_seg.intersects(inp_seg):
                curr_overlap = abs(inp_seg & psd_seg)
                if curr_overlap > psd_overlap:
                    psd_overlap = curr_overlap
                    best_psd = psd
        if best_psd is None:
            err_msg = "No PSDs found intersecting segment!"
            raise ValueError(err_msg)
        fd_segment.psd = best_psd
