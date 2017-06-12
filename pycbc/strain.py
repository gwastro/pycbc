#Copyright (C) 2013 Alex Nitz
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Generals
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""
This modules contains functions reading, generating, and segmenting strain data
"""
import copy
import logging, numpy, lal
import pycbc.noise
import pycbc.types
from pycbc.types import TimeSeries, zeros
from pycbc.types import Array, FrequencySeries, complex_same_precision_as
from pycbc.types import MultiDetOptionAppendAction, MultiDetOptionAction
from pycbc.types import MultiDetOptionActionSpecial
from pycbc.types import required_opts, required_opts_multi_ifo
from pycbc.types import ensure_one_opt, ensure_one_opt_multi_ifo
from pycbc.types import copy_opts_for_single_ifo
from pycbc.frame import read_frame, query_and_read_frame
from pycbc.inject import InjectionSet, SGBurstInjectionSet, RingdownInjectionSet
from pycbc.filter import resample_to_delta_t, highpass, make_frequency_series
from pycbc.filter.zpk import filter_zpk
import pycbc.psd
import pycbc.fft
import pycbc.events
import pycbc.frame
import pycbc.filter
from scipy.signal import kaiserord

def next_power_of_2(n):
    """Return the smallest integer power of 2 larger than the argument.

    Parameters
    ----------
    n : int
        A positive integer.

    Returns
    -------
    m : int
        Smallest integer power of 2 larger than n.
    """
    # TODO use 1 << n.bit_length() after Python 2.6 is gone
    return 1 << (len(bin(n)) - 2)

def detect_loud_glitches(strain, psd_duration=4., psd_stride=2.,
                         psd_avg_method='median', low_freq_cutoff=30.,
                         threshold=50., cluster_window=5., corrupt_time=4.,
                         high_freq_cutoff=None, output_intermediates=False):
    """Automatic identification of loud transients for gating purposes.

    This function first estimates the PSD of the input time series using the
    FindChirp Welch method. Then it whitens the time series using that
    estimate. Finally, it computes the magnitude of the whitened series,
    thresholds it and applies the FindChirp clustering over time to the
    surviving samples.

    Parameters
    ----------
    strain : TimeSeries
        Input strain time series to detect glitches over.
    psd_duration : {float, 4}
        Duration of the segments for PSD estimation in seconds.
    psd_stride : {float, 2}
        Separation between PSD estimation segments in seconds.
    psd_avg_method : {string, 'median'}
        Method for averaging PSD estimation segments.
    low_freq_cutoff : {float, 30}
        Minimum frequency to include in the whitened strain.
    threshold : {float, 50}
        Minimum magnitude of whitened strain for considering a transient to
        be present.
    cluster_window : {float, 5}
        Length of time window to cluster surviving samples over, in seconds.
    corrupt_time : {float, 4}
        Amount of time to be discarded at the beginning and end of the input
        time series.
    high_frequency_cutoff : {float, None}
        Maximum frequency to include in the whitened strain. If given, the
        input series is downsampled accordingly. If omitted, the Nyquist
        frequency is used.
    output_intermediates : {bool, False}
        Save intermediate time series for debugging.
    """
    logging.info('Autogating: tapering strain')
    taper_length = int(corrupt_time * strain.sample_rate)
    w = numpy.arange(taper_length) / float(taper_length)
    strain[0:taper_length] *= pycbc.types.Array(w, dtype=strain.dtype)
    strain[(len(strain)-taper_length):] *= pycbc.types.Array(w[::-1],
                                                             dtype=strain.dtype)

    # don't waste time trying to optimize a single FFT
    pycbc.fft.fftw.set_measure_level(0)

    if high_freq_cutoff:
        logging.info('Autogating: downsampling strain')
        strain = resample_to_delta_t(strain, 0.5 / high_freq_cutoff,
                                     method='ldas')
    if output_intermediates:
        strain.save_to_wav('strain_conditioned.wav')

    corrupt_length = int(corrupt_time * strain.sample_rate)

    # zero-pad strain to a power-of-2 length
    strain_pad_length = next_power_of_2(len(strain))
    pad_start = strain_pad_length/2 - len(strain)/2
    pad_end = pad_start + len(strain)
    strain_pad = pycbc.types.TimeSeries(
            pycbc.types.zeros(strain_pad_length, dtype=strain.dtype),
            delta_t=strain.delta_t, copy=False,
            epoch=strain.start_time-pad_start/strain.sample_rate)
    strain_pad[pad_start:pad_end] = strain[:]

    logging.info('Autogating: estimating PSD')
    psd = pycbc.psd.welch(strain[corrupt_length:(len(strain)-corrupt_length)],
                          seg_len=int(psd_duration * strain.sample_rate),
                          seg_stride=int(psd_stride * strain.sample_rate),
                          avg_method=psd_avg_method,
                          require_exact_data_fit=False)
    psd = pycbc.psd.interpolate(psd, 1./strain_pad.duration)
    psd = pycbc.psd.inverse_spectrum_truncation(
            psd, int(psd_duration * strain.sample_rate),
            low_frequency_cutoff=low_freq_cutoff,
            trunc_method='hann')
    kmin = int(low_freq_cutoff / psd.delta_f)
    psd[0:kmin] = numpy.inf
    if high_freq_cutoff:
        kmax = int(high_freq_cutoff / psd.delta_f)
        psd[kmax:] = numpy.inf

    logging.info('Autogating: time -> frequency')
    strain_tilde = pycbc.types.FrequencySeries(
            pycbc.types.zeros(len(strain_pad) / 2 + 1,
                              dtype=pycbc.types.complex_same_precision_as(strain)),
            delta_f=psd.delta_f, copy=False)
    pycbc.fft.fft(strain_pad, strain_tilde)

    logging.info('Autogating: whitening')
    if high_freq_cutoff:
        norm = high_freq_cutoff - low_freq_cutoff
    else:
        norm = strain.sample_rate/2. - low_freq_cutoff
    strain_tilde *= (psd * norm) ** (-0.5)

    logging.info('Autogating: frequency -> time')
    pycbc.fft.ifft(strain_tilde, strain_pad)

    pycbc.fft.fftw.set_measure_level(pycbc.fft.fftw._default_measurelvl)

    if output_intermediates:
        strain_pad[pad_start:pad_end].save_to_wav('strain_whitened.wav')

    logging.info('Autogating: computing magnitude')
    mag = abs(strain_pad[pad_start:pad_end])
    if output_intermediates:
        mag.save('strain_whitened_mag.npy')

    mag = mag.numpy()

    # remove strain corrupted by filters at the ends
    mag[0:corrupt_length] = 0
    mag[-1:-corrupt_length-1:-1] = 0

    logging.info('Autogating: finding loud peaks')
    indices = numpy.where(mag > threshold)[0]
    cluster_idx = pycbc.events.findchirp_cluster_over_window(
            indices, numpy.array(mag[indices]),
            int(cluster_window*strain.sample_rate))
    times = [idx * strain.delta_t + strain.start_time \
             for idx in indices[cluster_idx]]
    return times


def from_cli(opt, dyn_range_fac=1, precision='single',
             inj_filter_rejector=None):
    """Parses the CLI options related to strain data reading and conditioning.

    Parameters
    ----------
    opt : object
        Result of parsing the CLI with OptionParser, or any object with the
        required attributes  (gps-start-time, gps-end-time, strain-high-pass, 
        pad-data, sample-rate, (frame-cache or frame-files), channel-name, 
        fake-strain, fake-strain-seed, fake-strain-from-file, gating_file).
    dyn_range_fac: {float, 1}, optional
        A large constant to reduce the dynamic range of the strain.
    inj_filter_rejector: InjFilterRejector instance; optional, default=None
        If given send the InjFilterRejector instance to the inject module so
        that it can store a reduced representation of injections if
        necessary.

    Returns
    -------
    strain : TimeSeries
        The time series containing the conditioned strain data.
    """
    gating_info = {}

    if opt.frame_cache or opt.frame_files or opt.frame_type:
        if opt.frame_cache:
            frame_source = opt.frame_cache
        if opt.frame_files:
            frame_source = opt.frame_files

        logging.info("Reading Frames")
        
        if opt.frame_type:
            strain = query_and_read_frame(opt.frame_type, opt.channel_name,
                                          start_time=opt.gps_start_time-opt.pad_data,
                                          end_time=opt.gps_end_time+opt.pad_data)
        else:
            strain = read_frame(frame_source, opt.channel_name,
                            start_time=opt.gps_start_time-opt.pad_data,
                            end_time=opt.gps_end_time+opt.pad_data)

        if opt.zpk_z and opt.zpk_p and opt.zpk_k:
            logging.info("Highpass Filtering")
            strain = highpass(strain, frequency=opt.strain_high_pass)

            logging.info("Applying zpk filter")
            z = numpy.array(opt.zpk_z)
            p = numpy.array(opt.zpk_p)
            k = float(opt.zpk_k)
            strain = filter_zpk(strain.astype(numpy.float64), z, p, k)

        if opt.normalize_strain:
            logging.info("Dividing strain by constant")
            l = opt.normalize_strain
            strain = strain / l

        if opt.injection_file:
            logging.info("Applying injections")
            injector = InjectionSet(opt.injection_file)
            injections = \
                injector.apply(strain, opt.channel_name[0:2],
                               distance_scale=opt.injection_scale_factor,
                               inj_filter_rejector=inj_filter_rejector)

        if opt.sgburst_injection_file:
            logging.info("Applying sine-Gaussian burst injections")
            injector = SGBurstInjectionSet(opt.sgburst_injection_file)
            injector.apply(strain, opt.channel_name[0:2],
                             distance_scale=opt.injection_scale_factor)

        if opt.ringdown_injection_file:
            logging.info("Applying ringdown-only injection.")
            injector = RingdownInjectionSet(opt.ringdown_injection_file)
            injector.apply(strain, opt.channel_name[0:2])

        logging.info("Highpass Filtering")
        strain = highpass(strain, frequency=opt.strain_high_pass)

        if precision == 'single':
            logging.info("Converting to float32")
            strain = (strain * dyn_range_fac).astype(pycbc.types.float32)
        elif precision == "double":
            logging.info("Converting to float64")
            strain = (strain * dyn_range_fac).astype(pycbc.types.float64)
        else:
            raise ValueError("unrecognized precision {}".format(precision))

        if opt.gating_file is not None:
            logging.info("Gating glitches")
            gate_params = numpy.loadtxt(opt.gating_file)
            if len(gate_params.shape) == 1:
                gate_params = [gate_params]
            strain = gate_data(strain, gate_params)
            gating_info['file'] = \
                    [gp for gp in gate_params \
                     if (gp[0] + gp[1] + gp[2] >= strain.start_time) \
                     and (gp[0] - gp[1] - gp[2] <= strain.end_time)]

        if opt.autogating_threshold is not None:
            # the + 0 is for making a copy
            glitch_times = detect_loud_glitches(
                    strain + 0., threshold=opt.autogating_threshold,
                    cluster_window=opt.autogating_cluster,
                    low_freq_cutoff=opt.strain_high_pass,
                    high_freq_cutoff=opt.sample_rate/2,
                    corrupt_time=opt.pad_data+opt.autogating_pad)
            gate_params = [[gt, opt.autogating_width, opt.autogating_taper] \
                           for gt in glitch_times]
            if len(glitch_times) > 0:
                logging.info('Autogating at %s',
                             ', '.join(['%.3f' % gt for gt in glitch_times]))
            strain = gate_data(strain, gate_params)
            gating_info['auto'] = gate_params

        logging.info("Resampling data")
        strain = resample_to_delta_t(strain, 1.0/opt.sample_rate, method='ldas')

        logging.info("Highpass Filtering")
        strain = highpass(strain, frequency=opt.strain_high_pass)

        logging.info("Remove Padding")
        start = opt.pad_data*opt.sample_rate
        end = len(strain)-opt.sample_rate*opt.pad_data
        strain = strain[start:end]

    if opt.fake_strain or opt.fake_strain_from_file:
        logging.info("Generating Fake Strain")
        duration = opt.gps_end_time - opt.gps_start_time
        tlen = duration * opt.sample_rate
        pdf = 1.0/128
        plen = int(opt.sample_rate / pdf) / 2 + 1

        if opt.fake_strain_from_file:
            logging.info("Reading ASD from file")
            strain_psd = pycbc.psd.from_txt(opt.fake_strain_from_file, plen, pdf,
                                            opt.low_frequency_cutoff, is_asd_file=True)
        elif opt.fake_strain != 'zeroNoise':
            logging.info("Making PSD for strain")
            strain_psd = pycbc.psd.from_string(opt.fake_strain, plen, pdf,
                                               opt.low_frequency_cutoff)

        if opt.fake_strain == 'zeroNoise':
            logging.info("Making zero-noise time series")
            strain = TimeSeries(pycbc.types.zeros(tlen),
                                delta_t=1.0/opt.sample_rate)
        else:
            logging.info("Making colored noise")
            strain = pycbc.noise.noise_from_psd(tlen, 1.0/opt.sample_rate,
                                                strain_psd,
                                                seed=opt.fake_strain_seed)
        strain._epoch = lal.LIGOTimeGPS(opt.gps_start_time)

        if opt.injection_file:
            logging.info("Applying injections")
            injector = InjectionSet(opt.injection_file)
            injections = \
                injector.apply(strain, opt.channel_name[0:2],
                               distance_scale=opt.injection_scale_factor,
                               inj_filter_rejector=inj_filter_rejector)

        if opt.sgburst_injection_file:
            logging.info("Applying sine-Gaussian burst injections")
            injector =  SGBurstInjectionSet(opt.sgburst_injection_file)
            injector.apply(strain, opt.channel_name[0:2],
                             distance_scale=opt.injection_scale_factor)

        if opt.ringdown_injection_file:
            logging.info("Applying ringdown-only injection.")
            injector = RingdownInjectionSet(opt.ringdown_injection_file)
            injector.apply(strain, opt.channel_name[0:2])

        if precision == 'single':
            logging.info("Converting to float32")
            strain = (dyn_range_fac * strain).astype(pycbc.types.float32)

    if opt.taper_data:
        logging.info("Tapering data")
        # Use auto-gating stuff for this, a one-sided gate is a taper
        pd_taper_window = opt.taper_data
        gate_params = [(strain.start_time, 0., pd_taper_window)]
        gate_params.append( (strain.end_time, 0.,
                             pd_taper_window) )
        gate_data(strain, gate_params)

    if opt.injection_file:
        strain.injections = injections
    strain.gating_info = gating_info

    return strain  

def from_cli_single_ifo(opt, ifo, **kwargs):
    """
    Get the strain for a single ifo when using the multi-detector CLI
    """
    single_det_opt = copy_opts_for_single_ifo(opt, ifo)
    return from_cli(single_det_opt, **kwargs)

def from_cli_multi_ifos(opt, ifos, **kwargs):
    """
    Get the PSD for all ifos when using the multi-detector CLI
    """
    strain = {}
    for ifo in ifos:
        strain[ifo] = from_cli_single_ifo(opt, ifo, **kwargs)
    return strain


def insert_strain_option_group(parser, gps_times=True):
    """ Add strain-related options to the optparser object.

    Adds the options used to call the pycbc.strain.from_cli function to an
    optparser as an OptionGroup. This should be used if you
    want to use these options in your code.

    Parameters
    -----------
    parser : object
        OptionParser instance.
    """

    data_reading_group = parser.add_argument_group("Options for obtaining h(t)",
                  "These options are used for generating h(t) either by "
                  "reading from a file or by generating it. This is only "
                  "needed if the PSD is to be estimated from the data, ie. "
                  " if the --psd-estimation option is given.")

    # Required options
    
    if gps_times:
        data_reading_group.add_argument("--gps-start-time",
                                help="The gps start time of the data "
                                     "(integer seconds)", type=int)
        data_reading_group.add_argument("--gps-end-time",
                                help="The gps end time of the data "
                                     " (integer seconds)", type=int)
                                     
                                 
    data_reading_group.add_argument("--strain-high-pass", type=float,
                            help="High pass frequency")
    data_reading_group.add_argument("--pad-data",
              help="Extra padding to remove highpass corruption "
                   "(integer seconds)", type=int)
    data_reading_group.add_argument("--taper-data",
              help="Taper ends of data to zero using the supplied length as a "
                   "window (integer seconds)", type=int, default=0)
    data_reading_group.add_argument("--sample-rate", type=int,
                            help="The sample rate to use for h(t) generation (integer Hz).")
    data_reading_group.add_argument("--channel-name", type=str,
                   help="The channel containing the gravitational strain data")

    #Read from cache file
    data_reading_group.add_argument("--frame-cache", type=str, nargs="+",
                            help="Cache file containing the frame locations.")

    #Read from frame files
    data_reading_group.add_argument("--frame-files",
                            type=str, nargs="+",
                            help="list of frame files")

    #Use datafind to get frame files 
    data_reading_group.add_argument("--frame-type",
                            type=str,
                            help="(optional), replaces frame-files. Use datafind "
                                 "to get the needed frame file(s) of this type.")

    #Generate gaussian noise with given psd
    data_reading_group.add_argument("--fake-strain",
                help="Name of model PSD for generating fake gaussian noise.",
                     choices=pycbc.psd.get_lalsim_psd_list() + ['zeroNoise'])
    data_reading_group.add_argument("--fake-strain-seed", type=int, default=0,
                help="Seed value for the generation of fake colored"
                     " gaussian noise")
    data_reading_group.add_argument("--fake-strain-from-file",
                help="File containing ASD for generating fake noise from it.")

    #optional
    data_reading_group.add_argument("--injection-file", type=str,
                      help="(optional) Injection file used to add "
                           "waveforms into the strain")

    data_reading_group.add_argument("--sgburst-injection-file", type=str,
                      help="(optional) Injection file used to add "
                      "sine-Gaussian burst waveforms into the strain")

    data_reading_group.add_argument("--ringdown-injection-file", type=str,
                      help="(optional) Injection file used to add "
                           "ringdown-only waveforms into the strain.")

    data_reading_group.add_argument("--injection-scale-factor", type=float,
                    default=1, help="Divide injections by this factor "
                    "before injecting into the data.")

    data_reading_group.add_argument("--gating-file", type=str,
                    help="(optional) Text file of gating segments to apply."
                        " Format of each line is (all times in secs):"
                        "  gps_time zeros_half_width pad_half_width")

    data_reading_group.add_argument('--autogating-threshold', type=float,
                                    metavar='SIGMA',
                                    help='If given, find and gate glitches '
                                         'producing a deviation larger than '
                                         'SIGMA in the whitened strain time '
                                         'series.')
    data_reading_group.add_argument('--autogating-cluster', type=float,
                                    metavar='SECONDS', default=5.,
                                    help='Length of clustering window for '
                                         'detecting glitches for autogating.')
    data_reading_group.add_argument('--autogating-width', type=float,
                                    metavar='SECONDS', default=0.25,
                                    help='Half-width of the gating window.')
    data_reading_group.add_argument('--autogating-taper', type=float,
                                    metavar='SECONDS', default=0.25,
                                    help='Taper the strain before and after '
                                         'each gating window over a duration '
                                         'of SECONDS.')
    data_reading_group.add_argument('--autogating-pad', type=float,
                                    metavar='SECONDS', default=16,
                                    help='Ignore the given length of whitened '
                                         'strain at the ends of a segment, to '
                                         'avoid filters ringing.')

    data_reading_group.add_argument("--normalize-strain", type=float,
                    help="(optional) Divide frame data by constant.")

    data_reading_group.add_argument("--zpk-z", type=float, nargs="+",
                    help="(optional) Zero-pole-gain (zpk) filter strain. "
                        "A list of zeros for transfer function")

    data_reading_group.add_argument("--zpk-p", type=float, nargs="+",
                    help="(optional) Zero-pole-gain (zpk) filter strain. "
                        "A list of poles for transfer function")

    data_reading_group.add_argument("--zpk-k", type=float,
                    help="(optional) Zero-pole-gain (zpk) filter strain. "
                        "Transfer function gain")

    return data_reading_group

# FIXME: This repeats almost all of the options above. Any nice way of reducing
#        this?
def insert_strain_option_group_multi_ifo(parser):
    """
    Adds the options used to call the pycbc.strain.from_cli function to an
    optparser as an OptionGroup. This should be used if you
    want to use these options in your code.

    Parameters
    -----------
    parser : object
        OptionParser instance.
    """

    data_reading_group_multi = parser.add_argument_group("Options for obtaining"
                  " h(t)",
                  "These options are used for generating h(t) either by "
                  "reading from a file or by generating it. This is only "
                  "needed if the PSD is to be estimated from the data, ie. "
                  "if the --psd-estimation option is given. This group "
                  "supports reading from multiple ifos simultaneously.")

    # Required options
    data_reading_group_multi.add_argument("--gps-start-time", nargs='+',
                            action=MultiDetOptionAction, metavar='IFO:TIME',
                            help="The gps start time of the data "
                                 "(integer seconds)", type=int)
    data_reading_group_multi.add_argument("--gps-end-time", nargs='+', type=int,
                            action=MultiDetOptionAction, metavar='IFO:TIME',
                            help="The gps end time of the data "
                                 "(integer seconds)")
    data_reading_group_multi.add_argument("--strain-high-pass", nargs='+',
                            action=MultiDetOptionAction,
                            type=float, metavar='IFO:FREQUENCY',
                            help="High pass frequency")
    data_reading_group_multi.add_argument("--pad-data", nargs='+',
                            action=MultiDetOptionAction,
                            type=int, metavar='IFO:LENGTH',
                            help="Extra padding to remove highpass corruption "
                                "(integer seconds)")
    data_reading_group_multi.add_argument("--taper-data", nargs='+',
                            action=MultiDetOptionAction,
                            type=int, default=0, metavar='IFO:LENGTH',
                            help="Taper ends of data to zero using the "
                                "supplied length as a window (integer seconds)")
    data_reading_group_multi.add_argument("--sample-rate", type=int, nargs='+',
                            action=MultiDetOptionAction, metavar='IFO:RATE',
                            help="The sample rate to use for h(t) generation "
                                " (integer Hz).")
    data_reading_group_multi.add_argument("--channel-name", type=str, nargs='+',
                            action=MultiDetOptionActionSpecial,
                            metavar='IFO:CHANNEL',
                            help="The channel containing the gravitational "
                                "strain data")

    #Read from cache file
    data_reading_group_multi.add_argument("--frame-cache", type=str, nargs="+",
                            action=MultiDetOptionAppendAction,
                            metavar='IFO:FRAME_CACHE',
                            help="Cache file containing the frame locations.")

    #Read from frame files
    data_reading_group_multi.add_argument("--frame-files", type=str, nargs="+",
                            action=MultiDetOptionAppendAction,
                            metavar='IFO:FRAME_FILES',
                            help="list of frame files")

    # Use datafind to get frame files
    data_reading_group_multi.add_argument("--frame-type", type=str, nargs="+",
                                    action=MultiDetOptionAction,
                                    metavar='IFO:FRAME_TYPE',
                                    help="(optional) Replaces frame-files. "
                                         "Use datafind to get the needed frame "
                                         "file(s) of this type.")

    #Generate gaussian noise with given psd
    data_reading_group_multi.add_argument("--fake-strain", type=str, nargs="+",
                            action=MultiDetOptionAction, metavar='IFO:CHOICE',
                            help="Name of model PSD for generating fake "
                            "gaussian noise. Choose from %s or zeroNoise" \
                            %((', ').join(pycbc.psd.get_lalsim_psd_list()),) )
    data_reading_group_multi.add_argument("--fake-strain-seed", type=int,
                            default=0, nargs="+", action=MultiDetOptionAction,
                            metavar='IFO:SEED',
                            help="Seed value for the generation of fake "
                            "colored gaussian noise")
    data_reading_group_multi.add_argument("--fake-strain-from-file", nargs="+",
                            action=MultiDetOptionAction, metavar='IFO:FILE',
                            help="File containing ASD for generating fake "
                            "noise from it.")

    #optional
    data_reading_group_multi.add_argument("--injection-file", type=str,
                            nargs="+", action=MultiDetOptionAction,
                            metavar='IFO:FILE',
                            help="(optional) Injection file used to add "
                            "waveforms into the strain")

    data_reading_group_multi.add_argument("--sgburst-injection-file", type=str,
                      nargs="+", action=MultiDetOptionAction,
                      metavar='IFO:FILE',
                      help="(optional) Injection file used to add "
                      "sine-Gaussian burst waveforms into the strain")

    data_reading_group_multi.add_argument("--ringdown-injection-file", type=str,
                    nargs="+", action=MultiDetOptionAction, metavar='IFO:FILE',
                    help="(optional) Injection file used to add "
                           "ringdown-only waveforms into the strain.")

    data_reading_group_multi.add_argument("--injection-scale-factor",
                    type=float, nargs="+", action=MultiDetOptionAction,
                    metavar="IFO:VAL", default=1.,
                    help="Multiple injections by this factor "
                         "before injecting into the data.")

    data_reading_group_multi.add_argument("--gating-file", type=str,
                      nargs="+", action=MultiDetOptionAction,
                      metavar='IFO:FILE',
                      help="(optional) Text file of gating segments to apply."
                          " Format of each line is (all times in secs):"
                          "  gps_time zeros_half_width pad_half_width")

    data_reading_group_multi.add_argument('--autogating-threshold', type=float,
                                    nargs="+", action=MultiDetOptionAction,
                                    metavar='IFO:SIGMA',
                                    help='If given, find and gate glitches '
                                         'producing a deviation larger than '
                                         'SIGMA in the whitened strain time '
                                         'series.')
    data_reading_group_multi.add_argument('--autogating-cluster', type=float,
                                    nargs="+", action=MultiDetOptionAction,
                                    metavar='IFO:SECONDS', default=5.,
                                    help='Length of clustering window for '
                                         'detecting glitches for autogating.')
    data_reading_group_multi.add_argument('--autogating-width', type=float,
                                    nargs="+", action=MultiDetOptionAction,
                                    metavar='IFO:SECONDS', default=0.25,
                                    help='Half-width of the gating window.')
    data_reading_group_multi.add_argument('--autogating-taper', type=float,
                                    nargs="+", action=MultiDetOptionAction,
                                    metavar='IFO:SECONDS', default=0.25,
                                    help='Taper the strain before and after '
                                         'each gating window over a duration '
                                         'of SECONDS.')
    data_reading_group_multi.add_argument('--autogating-pad', type=float,
                                    nargs="+", action=MultiDetOptionAction,
                                    metavar='IFO:SECONDS', default=16,
                                    help='Ignore the given length of whitened '
                                         'strain at the ends of a segment, to '
                                         'avoid filters ringing.')

    data_reading_group_multi.add_argument("--normalize-strain", type=float,
                     nargs="+", action=MultiDetOptionAction,
                     metavar='IFO:VALUE',
                     help="(optional) Divide frame data by constant.")

    data_reading_group_multi.add_argument("--zpk-z", type=float,
                     nargs="+", action=MultiDetOptionAppendAction,
                     metavar='IFO:VALUE',
                     help="(optional) Zero-pole-gain (zpk) filter strain. "
                         "A list of zeros for transfer function")

    data_reading_group_multi.add_argument("--zpk-p", type=float,
                     nargs="+", action=MultiDetOptionAppendAction,
                     metavar='IFO:VALUE',
                     help="(optional) Zero-pole-gain (zpk) filter strain. "
                         "A list of poles for transfer function")

    data_reading_group_multi.add_argument("--zpk-k", type=float,
                     nargs="+", action=MultiDetOptionAppendAction,
                     metavar='IFO:VALUE',
                     help="(optional) Zero-pole-gain (zpk) filter strain. "
                         "Transfer function gain")

    return data_reading_group_multi


ensure_one_opt_groups = []
ensure_one_opt_groups.append(['--frame-cache','--fake-strain','--fake-strain-from-file','--frame-files', '--frame-type'])

required_opts_list = ['--gps-start-time', '--gps-end-time',
                      '--strain-high-pass', '--pad-data', '--sample-rate',
                      '--channel-name']

def verify_strain_options(opts, parser):
    """Sanity check provided strain arguments.

    Parses the strain data CLI options and verifies that they are consistent
    and reasonable.

    Parameters
    ----------
    opt : object
        Result of parsing the CLI with OptionParser, or any object with the
        required attributes (gps-start-time, gps-end-time, strain-high-pass, 
        pad-data, sample-rate, frame-cache, channel-name, fake-strain, 
        fake-strain-seed).
    parser : object
        OptionParser instance.
    """
    for opt_group in ensure_one_opt_groups:
        ensure_one_opt(opts, parser, opt_group)
    required_opts(opts, parser, required_opts_list)

def verify_strain_options_multi_ifo(opts, parser, ifos):
    """Sanity check provided strain arguments.

    Parses the strain data CLI options and verifies that they are consistent
    and reasonable.

    Parameters
    ----------
    opt : object
        Result of parsing the CLI with OptionParser, or any object with the
        required attributes (gps-start-time, gps-end-time, strain-high-pass, 
        pad-data, sample-rate, frame-cache, channel-name, fake-strain, 
        fake-strain-seed).
    parser : object
        OptionParser instance.
    ifos : list of strings
        List of ifos for which to verify options for
    """
    for ifo in ifos:
        for opt_group in ensure_one_opt_groups:
            ensure_one_opt_multi_ifo(opts, parser, ifo, opt_group)
        required_opts_multi_ifo(opts, parser, ifo, required_opts_list)


def gate_data(data, gate_params):
    """Apply a set of gating windows to a time series.
    
    Each gating window is
    defined by a central time, a given duration (centered on the given
    time) to zero out, and a given duration of smooth tapering on each side of
    the window. The window function used for tapering is a Tukey window.

    Parameters
    ----------
    data : TimeSeries
        The time series to be gated.
    gate_params : list
        List of parameters for the gating windows. Each element should be a
        list or tuple with 3 elements: the central time of the gating window,
        the half-duration of the portion to zero out, and the duration of the
        Tukey tapering on each side. All times in seconds. The total duration
        of the data affected by one gating window is thus twice the second
        parameter plus twice the third parameter.

    Returns
    -------
    data: TimeSeries
        The gated time series.
    """
    def inverted_tukey(M, n_pad):
        midlen = M - 2*n_pad
        if midlen < 0:
            raise ValueError("No zeros left after applying padding.")
        padarr = 0.5*(1.+numpy.cos(numpy.pi*numpy.arange(n_pad)/n_pad))
        return numpy.concatenate((padarr,numpy.zeros(midlen),padarr[::-1]))

    sample_rate = 1./data.delta_t
    temp = data.data

    for glitch_time, glitch_width, pad_width in gate_params:
        t_start = glitch_time - glitch_width - pad_width - data.start_time
        t_end = glitch_time + glitch_width + pad_width - data.start_time
        if t_start > data.duration or t_end < 0.:
            continue # Skip gate segments that don't overlap
        win_samples = int(2*sample_rate*(glitch_width+pad_width))
        pad_samples = int(sample_rate*pad_width)
        window = inverted_tukey(win_samples, pad_samples)
        offset = int(t_start * sample_rate)
        idx1 = max(0, -offset)
        idx2 = min(len(window), len(data)-offset)
        temp[idx1+offset:idx2+offset] *= window[idx1:idx2]

    return data

class StrainSegments(object):
    """ Class for managing manipulation of strain data for the purpose of
        matched filtering. This includes methods for segmenting and
        conditioning.
    """
    def __init__(self, strain, segment_length=None, segment_start_pad=0,
                 segment_end_pad=0, trigger_start=None, trigger_end=None,
                 filter_inj_only=False, injection_window=None,
                 allow_zero_padding=False):
        """ Determine how to chop up the strain data into smaller segments
            for analysis.
        """
        self._fourier_segments = None
        self.strain = strain

        self.delta_t = strain.delta_t
        self.sample_rate = strain.sample_rate

        if segment_length:
            seg_len = segment_length
        else:
            seg_len = strain.duration

        self.delta_f = 1.0 / seg_len
        self.time_len = seg_len * self.sample_rate
        self.freq_len = self.time_len / 2 + 1

        seg_end_pad = segment_end_pad
        seg_start_pad = segment_start_pad

        if not trigger_start:
            trigger_start = int(strain.start_time) + segment_start_pad
        else:
            if not allow_zero_padding:
                min_start_time = int(strain.start_time) + segment_start_pad
            else:
                min_start_time = int(strain.start_time)
            if trigger_start < min_start_time:
                err_msg = "Trigger start time must be within analysable "
                err_msg += "window. Asked to start from %d " %(trigger_start)
                err_msg += "but can only analyse from %d." %(min_start_time)
                raise ValueError(err_msg)

        if not trigger_end:
            trigger_end = int(strain.end_time) - segment_end_pad
        else:
            if not allow_zero_padding:
                max_end_time = int(strain.end_time) - segment_end_pad
            else:
                max_end_time = int(strain.end_time)
            if trigger_end > max_end_time:
                err_msg = "Trigger end time must be within analysable "
                err_msg += "window. Asked to end at %d " %(trigger_end)
                err_msg += "but can only analyse to %d." %(max_end_time)
                raise ValueError(err_msg)


        throwaway_size = seg_start_pad + seg_end_pad
        seg_width = seg_len - throwaway_size

        # The amount of time we can actually analyze given the
        # amount of padding that is needed
        analyzable = trigger_end - trigger_start
        data_start = (trigger_start - segment_start_pad) - \
                       int(strain.start_time)
        data_end = trigger_end + segment_end_pad - int(strain.start_time)
        data_dur = data_end - data_start
        data_start = data_start * strain.sample_rate
        data_end = data_end * strain.sample_rate

        #number of segments we need to analyze this data
        num_segs = int(numpy.ceil(float(analyzable) / float(seg_width)))

        # The offset we will use between segments
        seg_offset = int(numpy.ceil(analyzable / float(num_segs)))
        self.segment_slices = []
        self.analyze_slices = []

        # Determine how to chop up the strain into smaller segments
        for nseg in range(num_segs-1):
            # boundaries for time slices into the strain
            seg_start = int(data_start + (nseg*seg_offset) * strain.sample_rate)
            seg_end = int(seg_start + seg_len * strain.sample_rate)
            seg_slice = slice(seg_start, seg_end)
            self.segment_slices.append(seg_slice)

            # boundaries for the analyzable portion of the segment
            ana_start = int(seg_start_pad * strain.sample_rate)
            ana_end = int(ana_start + seg_offset * strain.sample_rate)
            ana_slice = slice(ana_start, ana_end)
            self.analyze_slices.append(ana_slice)

        # The last segment takes up any integer boundary slop
        seg_end = int(data_end)
        seg_start = int(seg_end - seg_len * strain.sample_rate)
        seg_slice = slice(seg_start, seg_end)
        self.segment_slices.append(seg_slice)

        remaining = (data_dur - ((num_segs - 1) * seg_offset + seg_start_pad))
        ana_start = int((seg_len - remaining) * strain.sample_rate)
        ana_end = int((seg_len - seg_end_pad) * strain.sample_rate)
        ana_slice = slice(ana_start, ana_end)
        self.analyze_slices.append(ana_slice)
        
        self.full_segment_slices = copy.deepcopy(self.segment_slices)

        #Remove segments that are outside trig start and end
        segment_slices_red = []
        analyze_slices_red = []
        trig_start_idx = (trigger_start - int(strain.start_time)) * strain.sample_rate
        trig_end_idx = (trigger_end - int(strain.start_time)) * strain.sample_rate

        if filter_inj_only and hasattr(strain, 'injections'):
            end_times = strain.injections.end_times()   
            end_times = [time for time in end_times if float(time) < trigger_end and float(time) > trigger_start]
            inj_idx = [(float(time) - float(strain.start_time)) * strain.sample_rate for time in end_times]

        for seg, ana in zip(self.segment_slices, self.analyze_slices):
            start = ana.start
            stop = ana.stop
            cum_start = start + seg.start
            cum_end = stop + seg.start

            # adjust first segment
            if trig_start_idx > cum_start:
                start += (trig_start_idx - cum_start)

            # adjust last segment
            if trig_end_idx < cum_end:
                stop -= (cum_end - trig_end_idx)

            if filter_inj_only and hasattr(strain, 'injections'):
                analyze_this = False
                inj_window = strain.sample_rate * 8
                for inj_id in inj_idx:
                    if inj_id < (cum_end + inj_window) and \
                            inj_id > (cum_start - inj_window):
                        analyze_this = True

                if not analyze_this:
                    continue

            if start < stop:
                segment_slices_red.append(seg)
                analyze_slices_red.append(slice(start, stop))

        self.segment_slices = segment_slices_red
        self.analyze_slices = analyze_slices_red

    def fourier_segments(self):
        """ Return a list of the FFT'd segments.
        Return the list of FrequencySeries. Additional properties are
        added that describe the strain segment. The property 'analyze'
        is a slice corresponding to the portion of the time domain equivelant
        of the segment to analyze for triggers. The value 'cumulative_index'
        indexes from the beginning of the original strain series.
        """
        if not self._fourier_segments:
            self._fourier_segments = []
            for seg_slice, ana in zip(self.segment_slices, self.analyze_slices):
                if seg_slice.start >= 0 and seg_slice.stop <= len(self.strain):
                    freq_seg = make_frequency_series(self.strain[seg_slice])
                # Assume that we cannot have a case where we both zero-pad on
                # both sides
                elif seg_slice.start < 0:
                    strain_chunk = self.strain[:seg_slice.stop]
                    strain_chunk.prepend_zeros(-seg_slice.start)
                    freq_seg = make_frequency_series(strain_chunk)
                elif seg_slice.stop > len(self.strain):
                    strain_chunk = self.strain[seg_slice.start:]
                    strain_chunk.append_zeros(seg_slice.stop - len(self.strain))
                    freq_seg = make_frequency_series(strain_chunk)
                freq_seg.analyze = ana
                freq_seg.cumulative_index = seg_slice.start + ana.start
                freq_seg.seg_slice = seg_slice
                self._fourier_segments.append(freq_seg)

        return self._fourier_segments

    @classmethod
    def from_cli(cls, opt, strain):
        """Calculate the segmentation of the strain data for analysis from
        the command line options.
        """
        return cls(strain, segment_length=opt.segment_length,
                   segment_start_pad=opt.segment_start_pad,
                   segment_end_pad=opt.segment_end_pad,
                   trigger_start=opt.trig_start_time,
                   trigger_end=opt.trig_end_time,
                   filter_inj_only=opt.filter_inj_only,
                   injection_window=opt.injection_window,
                   allow_zero_padding=opt.allow_zero_padding)

    @classmethod
    def insert_segment_option_group(cls, parser):
        segment_group = parser.add_argument_group(
                                   "Options for segmenting the strain",
                                  "These options are used to determine how to "
                                  "segment the strain into smaller chunks, "
                                  "and for determining the portion of each to "
                                  "analyze for triggers. ")
        segment_group.add_argument("--trig-start-time", type=int, default=0,
                    help="(optional) The gps time to start recording triggers")
        segment_group.add_argument("--trig-end-time", type=int, default=0,
                    help="(optional) The gps time to stop recording triggers")
        segment_group.add_argument("--segment-length", type=int,
                          help="The length of each strain segment in seconds.")
        segment_group.add_argument("--segment-start-pad", type=int,
                          help="The time in seconds to ignore of the "
                               "beginning of each segment in seconds. ")
        segment_group.add_argument("--segment-end-pad", type=int,
                          help="The time in seconds to ignore at the "
                               "end of each segment in seconds.")
        segment_group.add_argument("--allow-zero-padding", action='store_true',
                                   help="Allow for zero padding of data to "
                                        "analyze requested times, if needed.")
        # Injection optimization options
        segment_group.add_argument("--filter-inj-only", action='store_true',
                          help="Analyze only segments that contain an injection.")
        segment_group.add_argument("--injection-window", default=None,
                          type=float, help="""If using --filter-inj-only then
                          only search for injections within +/- injection
                          window of the injections's end time. This is useful
                          to speed up a coherent search or a search where we
                          initially filter at lower sample rate, and then
                          filter at full rate where needed. NOTE: Reverts to
                          full analysis if two injections are in the same
                          segment.""")


    @classmethod
    def from_cli_single_ifo(cls, opt, strain, ifo):
        """Calculate the segmentation of the strain data for analysis from
        the command line options.
        """
        return cls(strain, segment_length=opt.segment_length[ifo],
                   segment_start_pad=opt.segment_start_pad[ifo],
                   segment_end_pad=opt.segment_end_pad[ifo],
                   trigger_start=opt.trig_start_time[ifo],
                   trigger_end=opt.trig_end_time[ifo],
                   filter_inj_only=opt.filter_inj_only,
                   allow_zero_padding=opt.allow_zero_padding)

    @classmethod
    def from_cli_multi_ifos(cls, opt, strain_dict, ifos):
        """Calculate the segmentation of the strain data for analysis from
        the command line options.
        """
        strain_segments = {}
        for ifo in ifos:
            strain_segments[ifo] = cls.from_cli_single_ifo(opt,
                                                         strain_dict[ifo], ifo)
        return strain_segments

    @classmethod
    def insert_segment_option_group_multi_ifo(cls, parser):
        segment_group = parser.add_argument_group(
                                   "Options for segmenting the strain",
                                  "These options are used to determine how to "
                                  "segment the strain into smaller chunks, "
                                  "and for determining the portion of each to "
                                  "analyze for triggers. ")
        segment_group.add_argument("--trig-start-time", type=int, default=0,
                    nargs='+', action=MultiDetOptionAction, metavar='IFO:TIME',
                    help="(optional) The gps time to start recording triggers")
        segment_group.add_argument("--trig-end-time", type=int, default=0,
                    nargs='+', action=MultiDetOptionAction, metavar='IFO:TIME',
                    help="(optional) The gps time to stop recording triggers")
        segment_group.add_argument("--segment-length", type=int,
                    nargs='+', action=MultiDetOptionAction,
                    metavar='IFO:LENGTH',
                    help="The length of each strain segment in seconds.")
        segment_group.add_argument("--segment-start-pad", type=int,
                    nargs='+', action=MultiDetOptionAction, metavar='IFO:TIME',
                    help="The time in seconds to ignore of the "
                        "beginning of each segment in seconds. ")
        segment_group.add_argument("--segment-end-pad", type=int,
                    nargs='+', action=MultiDetOptionAction, metavar='IFO:TIME',
                    help="The time in seconds to ignore at the "
                         "end of each segment in seconds.")
        segment_group.add_argument("--allow-zero-padding", action='store_true',
                          help="Allow for zero padding of data to analyze "
                          "requested times, if needed.")
        segment_group.add_argument("--filter-inj-only", action='store_true',
                                   help="Analyze only segments that contain "
                                        "an injection.")

    required_opts_list = ['--segment-length',
                   '--segment-start-pad',
                   '--segment-end-pad',
                   ]

    @classmethod
    def verify_segment_options(cls, opt, parser):
        required_opts(opt, parser, cls.required_opts_list)

    @classmethod
    def verify_segment_options_multi_ifo(cls, opt, parser, ifos):
        for ifo in ifos:
            required_opts_multi_ifo(opt, parser, ifo, cls.required_opts_list)

class StrainBuffer(pycbc.frame.DataBuffer):
    def __init__(self, frame_src, channel_name, start_time,
                 max_buffer=512,
                 sample_rate=4096,
                 low_frequency_cutoff=20,
                 highpass_frequency=15.0,
                 highpass_reduction=200.0,
                 highpass_bandwidth=5.0,
                 psd_samples=30,
                 psd_segment_length=4,
                 psd_inverse_length=3.5,
                 trim_padding=0.25,
                 autogating_threshold=100,
                 autogating_cluster=0.25,
                 autogating_window=0.5,
                 autogating_pad=0.25,
                 state_channel=None,
                 data_quality_channel=None,
                 dyn_range_fac=pycbc.DYN_RANGE_FAC,
                 psd_abort_difference=None,
                 psd_recalculate_difference=None,
                 force_update_cache=True,
                 increment_update_cache=None,
                 analyze_flags=None,
                 data_quality_flags=None,
                 dq_padding=0):
        """ Class to produce overwhitened strain incrementally
        
        Parameters
        ----------
        frame_src: str of list of strings
            Strings that indicate where to read from files from. This can be a
            list of frame files, a glob, etc.
        channel_name: str
            Name of the channel to read from the frame files
        start_time: 
            Time to start reading from.
        max_buffer: {int, 512}, Optional
            Length of the buffer in seconds
        sample_rate: {int, 2048}, Optional
            Rate in Hz to sample the data.
        low_frequency_cutoff: {float, 20}, Optional
            The low frequency cutoff to use for inverse spectrum truncation
        highpass_frequency: {float, 15}, Optional
            The frequency to apply a highpass filter at before downsampling.
        highpass_reduction: {float, 200}, Optional
            The amount of reduction to apply to the low frequencies.
        highpass_bandwidth: {float, 5}, Optional
            The width of the transition region for the highpass filter.
        psd_samples: {int, 30}, Optional
            The number of sample to use for psd estimation
        psd_segment_length: {float, 4}, Optional
            The number of seconds in each psd sample.
        psd_inverse_length: {float, 3.5}, Optional
            The length in seconds for fourier transform of the inverse of the
            PSD to be truncated to.
        trim_padding: {float, 0.25}, Optional
            Amount of padding in seconds to give for truncated the overwhitened
            data stream.
        autogating_threshold: {float, 100}, Optional
            Sigma deviation required to cause gating of data
        autogating_cluster: {float, 0.25}, Optional
            Seconds to cluster possible gating locations
        autogating_window: {float, 0.5}, Optional
            Seconds to window out when gating a time
        autogating_pad: {float, 0.25}, Optional
            Seconds to pad either side of the gating window.
        state_channel: {str, None}, Optional
            Channel to use for state information about the strain
        data_quality_channel: {str, None}, Optional
            Channel to use for data quality information about the strain
        dyn_range_fac: {float, pycbc.DYN_RANGE_FAC}, Optional
            Scale factor to apply to strain
        psd_abort_difference: {float, None}, Optional
            The relative change in the inspiral range from the previous PSD
            estimate to trigger the data to be considered invalid.
        psd_recalculate_difference: {float, None}, Optional
            the relative change in the inspiral range from the previous PSD
            to trigger a re-estimatoin of the PSD.
        force_update_cache: {boolean, True}, Optional
            Re-check the filesystem for frame files on every attempt to 
            read more data.
        analyze_flags: list of strs
            The flags that must be on to mark the current data as valid for
            *any* use.
        data_quality_flags: list of strs
            The flags used to determine if to keep triggers.
        dq_padding: {float, 0}, optional
            Extra seconds to consider invalid before/after times with bad DQ.
        increment_update_cache: {str, None}, Optional
            Pattern to look for frame files in a GPS dependent directory. This
            is an alternate to the forced updated of the frame cache, and
            apptempts to predict the next frame file name without probing the
            filesystem.
        """ 
        super(StrainBuffer, self).__init__(frame_src, channel_name, start_time,
                                           max_buffer=max_buffer,
                                           force_update_cache=force_update_cache,
                                           increment_update_cache=increment_update_cache)

        self.low_frequency_cutoff = low_frequency_cutoff

        # Set up status buffers
        self.analyze_flags = analyze_flags
        self.data_quality_flags = data_quality_flags
        self.state = None
        self.dq = None
        self.dq_padding = dq_padding

        # State channel
        if state_channel is not None:
            valid_mask = 0
            for flag in self.analyze_flags:
                valid_mask = valid_mask | getattr(pycbc.frame, flag)
            self.state = pycbc.frame.StatusBuffer(
                frame_src,
                state_channel, start_time,
                max_buffer=max_buffer,
                valid_mask=valid_mask,
                force_update_cache=force_update_cache,
                increment_update_cache=increment_update_cache)

        # low latency dq channel
        if data_quality_channel is not None:
            valid_mask = 0
            for flag in self.data_quality_flags:
                valid_mask = valid_mask | getattr(pycbc.frame, flag)
            self.dq = pycbc.frame.StatusBuffer(
                frame_src,
                data_quality_channel, start_time,
                max_buffer=max_buffer,
                valid_mask=valid_mask,
                force_update_cache=force_update_cache,
                increment_update_cache=increment_update_cache)

        self.highpass_frequency = highpass_frequency
        self.highpass_reduction = highpass_reduction
        self.highpass_bandwidth = highpass_bandwidth

        self.autogating_threshold = autogating_threshold
        self.autogating_cluster = autogating_cluster
        self.autogating_pad = autogating_pad
        self.autogating_window = autogating_window

        self.sample_rate = sample_rate
        self.dyn_range_fac = dyn_range_fac

        self.psd_abort_difference = psd_abort_difference
        self.psd_recalculate_difference = psd_recalculate_difference
        self.psd_segment_length = psd_segment_length
        self.psd_samples = psd_samples
        self.psd_inverse_length = psd_inverse_length
        self.psd = None
        self.psds = {}

        strain_len = int(sample_rate * self.raw_buffer.delta_t * len(self.raw_buffer))
        self.strain = TimeSeries(zeros(strain_len, dtype=numpy.float32),
                                 delta_t=1.0/self.sample_rate, 
                                 epoch=start_time-max_buffer)

        # Determine the total number of corrupted samples for highpass
        # and PSD over whitening
        highpass_samples, self.beta = kaiserord(self.highpass_reduction,
          self.highpass_bandwidth / self.raw_buffer.sample_rate * 2 * numpy.pi)
        self.highpass_samples =  int(highpass_samples / 2)
        resample_corruption = 10 # If using the ldas method
        self.factor = int(1.0 / self.raw_buffer.delta_t / self.sample_rate)
        self.corruption = self.highpass_samples / self.factor + resample_corruption

        self.psd_corruption =  self.psd_inverse_length * self.sample_rate
        self.total_corruption = self.corruption + self.psd_corruption

        # Determine how much padding is needed after removing the parts
        # associated with PSD over whitening and highpass filtering
        self.trim_padding = int(trim_padding * self.sample_rate)
        if self.trim_padding > self.total_corruption:
            self.trim_padding = self.total_corruption

        self.psd_duration = (psd_samples - 1) / 2 * psd_segment_length

        self.reduced_pad = int(self.total_corruption - self.trim_padding)
        self.segments = {}
        
        # time to ignore output of frame (for initial buffering)
        self.add_hard_count()
        self.taper_immediate_strain = True

    @property
    def start_time(self):
        """ Return the start time of the current valid segment of data """
        return self.end_time - self.blocksize

    @property
    def end_time(self):
        """ Return the end time of the current valid segment of data """
        return float(self.strain.start_time + (len(self.strain) - self.total_corruption) / self.sample_rate)

    def add_hard_count(self):
        """ Reset the countdown timer, so that we don't analyze data long enough
        to generate a new PSD.
        """
        self.wait_duration = int(numpy.ceil(self.total_corruption / self.sample_rate +  self.psd_duration))
        self.invalidate_psd()

    def invalidate_psd(self):
        """ Make the current PSD invalid. A new one will be generated when
        it is next required """
        self.psd = None
        self.psds = {}

    def recalculate_psd(self):
        """ Recalculate the psd 
        """ 

        seg_len = self.sample_rate * self.psd_segment_length
        e = len(self.strain)
        s = e - ((self.psd_samples + 1) * self.psd_segment_length / 2) * self.sample_rate
        psd = pycbc.psd.welch(self.strain[s:e], seg_len=seg_len, seg_stride=seg_len / 2)

        from pycbc.waveform.spa_tmplt import spa_distance
        psd.dist = spa_distance(psd, 1.4, 1.4, self.low_frequency_cutoff) * pycbc.DYN_RANGE_FAC

        # If the new psd is similar to the old one, don't replace it
        if self.psd and self.psd_recalculate_difference:
            if abs(self.psd.dist - psd.dist) / self.psd.dist < self.psd_recalculate_difference:
                logging.info("Skipping Recalculation of PSD, %s-%s", self.psd.dist, psd.dist)
                return True
        
        # If the new psd is *really* different than the old one, return an error
        if self.psd and self.psd_abort_difference:
            if abs(self.psd.dist - psd.dist) / self.psd.dist > self.psd_abort_difference:
                logging.info("PSD is CRAZY, aborting!!!!, %s-%s", self.psd.dist, psd.dist)
                self.psd = psd
                self.psds = {}  
                return False

        # If the new estimate replaces the current one, invalide the ineterpolate PSDs
        self.psd = psd
        self.psds = {}     
        logging.info("Recalculating PSD, %s", psd.dist)  
        return True   

    def overwhitened_data(self, delta_f):
        """ Return overwhitened data

        Parameters
        ----------
        delta_f: float
            The sample step to generate overwhitened frequency domain data for
        
        Returns
        -------
        htilde: FrequencySeries
            Overwhited strain data
        """
        # we haven't alread computed htilde for this delta_f
        if delta_f not in self.segments:
            buffer_length = int(1.0 / delta_f)
            e = len(self.strain)
            s = int(e - buffer_length * self.sample_rate - self.reduced_pad * 2)
            fseries = make_frequency_series(self.strain[s:e])

            # we haven't calculate a resample psd for this delta_f
            if delta_f not in self.psds:
                psdt = pycbc.psd.interpolate(self.psd, fseries.delta_f)
                psdt = pycbc.psd.inverse_spectrum_truncation(psdt, 
                                       int(self.sample_rate * self.psd_inverse_length),
                                       low_frequency_cutoff=self.low_frequency_cutoff)
                psdt._delta_f = fseries.delta_f

                psd = pycbc.psd.interpolate(self.psd, delta_f)
                psd = pycbc.psd.inverse_spectrum_truncation(psd,
                                       int(self.sample_rate * self.psd_inverse_length),
                                       low_frequency_cutoff=self.low_frequency_cutoff)

                psd.psdt = psdt
                self.psds[delta_f] = psd

            psd = self.psds[delta_f]
            fseries /= psd.psdt

            # trim ends of strain
            if self.reduced_pad  != 0:
                overwhite = TimeSeries(zeros(e-s, dtype=self.strain.dtype),
                                             delta_t=self.strain.delta_t)
                pycbc.fft.ifft(fseries, overwhite)
                overwhite2 = overwhite[self.reduced_pad:len(overwhite)-self.reduced_pad]
                taper_window = self.trim_padding / 2.0 / overwhite.sample_rate
                gate_params = [(overwhite2.start_time, 0., taper_window),
                               (overwhite2.end_time, 0., taper_window)]
                gate_data(overwhite2, gate_params)
                fseries_trimmed = FrequencySeries(zeros(len(overwhite2) / 2 + 1,
                                                  dtype=fseries.dtype), delta_f=delta_f)
                pycbc.fft.fft(overwhite2, fseries_trimmed)
            else:
                fseries_trimmed = fseries

            fseries_trimmed.psd = psd
            self.segments[delta_f] = fseries_trimmed
            
        stilde = self.segments[delta_f]
        return stilde  

    def near_hwinj(self):
        """Check that the current set of triggers could be influenced by
        a hardware injection.

        Parameters
        ----------
        data_reader: dict of StrainBuffers
            A dict of StrainBuffer instances, indexed by ifos.
        """
        if not self.state:
            return False

        if not self.state.is_extent_valid(self.start_time, self.blocksize, pycbc.frame.NO_HWINJ):
            return True
        return False

    def null_advance_strain(self, blocksize):
        """ Advance and insert zeros

        Parameters
        ----------
        blocksize: int
            The number of seconds to attempt to read from the channel
        """
        sample_step = int(blocksize * self.sample_rate)
        csize = sample_step + self.corruption * 2
        self.strain.roll(-sample_step)

        # We should roll this off at some point too...
        self.strain[len(self.strain) - csize + self.corruption:] = 0
        self.strain.start_time += blocksize
        
        # The next time we need strain will need to be tapered
        self.taper_immediate_strain = True
       
    def advance(self, blocksize, timeout=10):
        """Advanced buffer blocksize seconds.

        Add blocksize seconds more to the buffer, push blocksize seconds
        from the beginning.

        Parameters
        ----------
        blocksize: int
            The number of seconds to attempt to read from the channel

        Returns
        -------
        status: boolean
            Returns True if this block is analyzable.         
        """
        ts = super(StrainBuffer, self).attempt_advance(blocksize, timeout=timeout)
        self.blocksize = blocksize

        # We have given up so there is no time series
        if ts is None:
            logging.info("Giving on up frame...")
            self.null_advance_strain(blocksize)
            if self.state:
                self.state.null_advance(blocksize)
            if self.dq:
                self.dq.null_advance(blocksize)
            return False

        # We collected some data so we are closer to being able to analyze data
        self.wait_duration -= blocksize

        # If the data we got was invalid, reset the counter on how much to collect
        # This behavior corresponds to how we handle CAT1 vetoes
        if self.state and self.state.advance(blocksize) is False:
            self.add_hard_count()
            self.null_advance_strain(blocksize)
            if self.dq:
                self.dq.null_advance(blocksize)
            logging.info("Time has invalid data, resetting buffer")
            return False

        # Also advance the dq vector in lockstep
        if self.dq:
            self.dq.advance(blocksize)            

        self.segments = {}

        # only condition with the needed raw data so we can continuously add
        # to the existing result

        # Precondition
        sample_step = int(blocksize * self.sample_rate)
        csize = sample_step + self.corruption * 2
        start = len(self.raw_buffer) - csize * self.factor
        strain = self.raw_buffer[start:]

        strain =  pycbc.filter.highpass_fir(strain, self.highpass_frequency,
                                       self.highpass_samples,
                                       beta=self.beta)
        strain = (strain * self.dyn_range_fac).astype(numpy.float32)
        
        strain = pycbc.filter.resample_to_delta_t(strain, 
                                           1.0/self.sample_rate, method='ldas')

        # remove corruption at beginning
        strain = strain[self.corruption:]
        
        # taper beginning if needed
        if self.taper_immediate_strain:
            logging.info("tapering start of strain block")
            strain = gate_data(strain, [(strain.start_time, 0., self.autogating_pad)])
            self.taper_immediate_strain = False

        # Stitch into continuous stream
        self.strain.roll(-sample_step)
        self.strain[len(self.strain) - csize + self.corruption:] = strain[:]
        self.strain.start_time += blocksize

        # apply gating if need be: NOT YET IMPLEMENTED
        if self.psd is None and self.wait_duration <=0:
            self.recalculate_psd()

        if self.wait_duration > 0:
            return False
        else:
            return True
