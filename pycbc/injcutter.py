# Copyright (C) 2013 Ian Harry
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
This module contains functions to filter injections with only useful templates.

This module implements a set of checks to test for each segment and template
combination whether injections contained within the segment are sufficiently
"similar" to the template to require a matched-filter. There are a few ways of
testing the "similarity" of templates and injections.
 * A chirp time threshold rejects templates if chirp time difference is large
 * A coarse match threshold rejects templates if a coarse overlap is small
"""

from pycbc import DYN_RANGE_FAC
from pycbc.pnutils import nearest_larger_binary_number
from pycbc.types import zeros

_injcutter_group_help = ("Options that, if injections are present in this "
                         "run, are responsible for performing pre-checks "
                         "between injections in the data being filtered and "
                         "the current search template "
                         "to determine if the template has any chance of "
                         "actually detecting the injection. The parameters "
                         "of this test are given by the various injcutter "
                         "options below. The --injcutter-chirp-time-threshold "
                         "and --injcutter-match-threshold options need to be "
                         "provided if those tests are desired. Other options "
                         "will take default values unless overriden. More "
                         "details on these options follow.")

_enable_cutter_help = ("If given 'injcutter' will be enabled.")

_injcutter_cthresh_help = ("If this value is not None and injcutter is "
                           "enabled then we will calculate the difference in "
                           "chirp time (tau_0) between the template and each "
                           "injection in the analysis segment. If the "
                           "difference is greate than this threshold for all "
                           "injections then filtering is not performed. By "
                           "default this will be None.")
_injcutter_mthresh_help = ("If this value is not None and injcutter is "
                           "enabled then we will calculate a 'coarse match' "
                           "between the template and each injection in the "
                           "analysis segment. If the match is less than this "
                           "threshold for all injections filtering is not "
                           "performed. Parameters for the 'coarse match' "
                           "follow. By default this value will be None.")
_injcutter_deltaf_help = ("If injcutter is enabled and "
                          "injcutter-match-threshold is not None, this option "
                          "specifies the frequency spacing that will be used "
                          "for injections, templates and PSD when computing "
                          "the 'coarse match'. Templates will be generated "
                          "directly with this spacing. The PSD and injections "
                          "will be resampled.")
_injcutter_fmax_help = ("If injcutter is enabled and "
                        "injcutter-match-threshold is not None, this option "
                        "specifies the maximum frequency that will be used "
                        "for injections, templates and PSD when computing "
                        "the 'coarse match'. Templates will be generated "
                        "directly with this max frequency. The PSD and "
                        "injections frequency series will be truncated.")
_injcutter_buffer_help = ("If injcutter is enabled, the injcutter tests "
                          "will determine if injections are 'in' the "
                          "specified analysis chunk by using the end times. "
                          "If this value is non-zero the analysis chunk is "
                          "extended on both sides by this amount before "
                          "determining if injections are within the given "
                          "window.")
_injcutter_flower_help = ("If injcutter is enabled, the injcutter tests "
                          "need a lower frequency for determine chirp times "
                          "or for doing matches. If this value is None the "
                          "lower frequency used for the full matched-filter "
                          "if used. Otherwise this value is used instead.")

def insert_injcutter_option_group(parser):
    injcutter_group = parser.add_argument_group(_injcutter_group_help)
    injcutter_group.add_argument("--enable-injcutter", action="store_true",
                                 default=False, help=_enable_cutter_help)
    injcutter_group.add_argument("--injcutter-chirp-time-threshold",
                                 type=float, default=None,
                                 help=_injcutter_cthresh_help)
    injcutter_group.add_argument("--injcutter-match-threshold", type=float,
                                 default=None, help=_injcutter_mthresh_help)
    injcutter_group.add_argument("--injcutter-coarsematch-deltaf",
                                 type=float, default=1.,
                                 help=_injcutter_deltaf_help)
    injcutter_group.add_argument("--injcutter-coarsematch-fmax", type=float,
                                 default=256., help=_injcutter_fmax_help)
    injcutter_group.add_argument("--injcutter-seg-buffer", type=int,
                                 default=10, help=_injcutter_buffer_help)
    injcutter_group.add_argument("--injcutter-f-lower", type=int,
                                 default=None, help=_injcutter_flower_help)


class InjCutter(object):
    """
    Class for holding parameters for using injection/template pre-filtering.
    """
    def __init__(self, injection_file, chirp_time_threshold,
                 match_threshold, f_lower, coarsematch_deltaf=1.,
                 seg_buffer=10):
        """ Initialise injcutter instance.
        """
        # Determine if injcutter is to be enabled
        if injection_file is None or injection_file == 'False' or\
            (chirp_time_threshold is None and match_threshold is None):
            self.enabled = False
            return
        else:
            self.enabled = True

        # Store parameters
        self.chirp_time_threshold = chirp_time_threshold
        self.match_threshold = match_threshold
        self.coarsematch_deltaf = coarsematch_deltaf
        self.coarsematch_fmax = coarsematch_fmax
        self.seg_buffer = seg_buffer
        self.f_lower = f_lower
        assert(self.f_lower is not None)

        # Variables for holding injcutter inputs (reduced injections, memory
        # for templates, reduced PSDs ...)
        self.short_injections = {}
        self._short_template_mem = None
        self._short_psd_storage = {}

    @classmethod
    def from_cli(self, opt):
        """
        Initialize an InjCutter instance from command-line options.
        """
        injection_file = opt.injection_file
        chirp_time_threshold = opt.injcutter_chirp_time_threshold
        match_threshold = opt.injcutter_match_threshold
        coarsematch_deltaf = opt.injcutter_coarsematch_deltaf
        coarsematch_fmax = opt.injcutter_coarsematch_fmax
        seg_buffer = opt.injcutter_seg_buffer
        if opt.injcutter_f_lower is not None:
            f_lower = opt.injcutter_f_lower
        else:
            # NOTE: Uses main low-frequency cutoff as default option. This may
            #       need some editing if using this in multi_inspiral, which I
            #       leave for future work, or if this is being used in another
            #       code which doesn't have --low-frequency-cutoff
            f_lower = opt.low_frequency_cutoff
        return cls(injection_file, chirp_time_threshold, match_threshold,
                   f_lower, coarsematch_deltaf=coarsematch_deltaf,
                   seg_buffer=seg_buffer)

    def generate_short_inj_from_inj(self, inj_waveform, simulation_id):
        """
        Generate and a store a truncated representation of inj_waveform.
        """
        if not self.enabled:
            # Do nothing!
            return
        if self.short_injections.has_key(simulation_id):
            err_msg = "An injection with simulation id "
            err_msg += str(simulation_id)
            err_msg += " has already been added to injcutter. This suggests "
            err_msg += "that your injection file contains injections with "
            err_msg += "duplicate simulation_ids. This is not allowed."
            raise ValueError(err_msg)
        curr_length = len(inj_waveform)
        new_length = int(nearest_larger_binary_number(inj_length))
        # Don't want length less than 1/delta_f
        while new_length * inj_waveform.delta_t < 1./self.coarsematch_deltaf:
            new_length = new_length * 2
        inj_waveform.resize(new_length)
        inj_tilde = inj_waveform.to_frequencyseries()
        # Dynamic range is important here!
        inj_tilde_np = inj_tilde.numpy() * DYN_RANGE_FAC
        delta_f = inj_tilde.get_delta_f()
        new_freq_len = int(self.coarsematch_deltaf / delta_f + 1)
        # This shouldn't be a problem if injections are generated at
        # 16384 Hz ... It is only a problem of injection sample rate
        # gives a lower Nyquist than the trunc_f_max. If this error is
        # ever raised one could consider zero-padding the injection.
        assert(new_freq_len <= len(inj_tilde))
        df_ratio = int(self.coarsematch_deltaf/delta_f)
        inj_tilde_np = inj_tilde_np[:new_freq_len:df_ratio]
        new_inj = FrequencySeries(inj_tilde_np, dtype=np.complex64,
                                  delta_f=self.coarsematch_deltaf)
        self.short_injections[simulation_id] = new_inj

    def template_segment_checker(self, bank, t_num, gwstrain, segment,
                                 start_time):
        """ Test if injections in segment are worth filtering with template.

        Using the current template, current segment, and injections within that
        segment. Test if the injections and sufficiently "similar" to any of
        the injections to justify actually performing a matched-filter call.
        Ther are two parts to this test: First we check if the chirp time of
        the template is within a provided window of any of the injections. If
        not then stop here, it is not worth filtering this template, segment
        combination for this injection set. If this check passes we compute a
        match between a coarse representation of the template and a coarse
        representation of each of the injections. If that match is above a
        user-provided value for any of the injections then filtering can
        proceed. This is currently only available if using frequency-domain
        templates.

        Parameters
        -----------
        FIXME

        Returns
        --------
        FIXME
        """
        if not self.enabled:
            # If disabled, always filter (ie. return True)
            return True

        # Get times covered by segment analyze
        seg_start_time = \
            segment.cumulative_index / float(gwstrain.sample_rate) + \
            start_time
        seg_end_time = \
            (segment.cumulative_index +
             (segment.analyze.stop - segment.analyze.start)) / \
            float(gwstrain.sample_rate) + start_time
        # And add buffer
        seg_start_time = seg_start_time - self.seg_buffer
        seg_end_time = seg_end_time + self.seg_buffer

        # Chirp time test
        if self.chirp_time_threshold is not None:
            m1 = bank.table[t_num]['mass1']
            m2 = bank.table[t_num]['mass2']
            tau0_temp, tau3_temp = \
                pycbc.pnutils.mass1_mass2_to_tau0_tau3(m1, m2, self.f_lower)
            for inj_idx, inj in enumerate(gwstrain.injections.table):
                end_time = inj.geocent_end_time + \
                    1E-9 * inj.geocent_end_time_ns
                if end_time > seg_end_time or end_time < seg_start_time:
                    continue
                tau0_inj, tau3_inj = \
                    pycbc.pnutils.mass1_mass2_to_tau0_tau3(inj.mass1,
                                                           inj.mass2,
                                                           self.f_lower)
                tau_diff = abs(tau0_temp - tau0_inj)
                if tau_diff <= ic_params['chirp_time_threshold']:
                    break
            else:
                # Get's here if all injections are outside chirp-time window
                return False

        # Coarse match test
        if self.match_threshold:
            if self._short_template_mem is None:
                # Set the memory for the short templates
                wav_len = 1 + int(self.coarsematch_fmax / \
                                  self.coarsematch_deltaf)
                self._short_template_mem = zeros(wav_len,
                                                 dtype=numpy.complex64)

            # Set the current short PSD to red_psd
            try:
                red_psd = self._short_psd_storage[id(segment.psd)]
            except KeyError:
                # PSD doesn't exist yet, so make it!
                curr_psd = segment.psd.numpy()
                step_size = int(self.coarsematch_deltaf / segment.psd.delta_f)
                max_idx = int(self.coarsematch_fmax / segment.psd.delta_f) + 1
                red_psd_data = curr_psd[:max_idx:step_size]
                red_psd = FrequencySeries(red_psd_data, copy=False,
                                          delta_f=self.coarsematch_deltaf)
                self._short_psd_storage[id(curr_psd)] = red_psd

            # Set htilde to be the current short template
            if not t_num == self._short_template_id:
                # Set the memory for the short templates if unset
                if self._short_template_mem is None:
                    wav_len = 1 + int(self.coarsematch_fmax / \
                                      self.coarsematch_deltaf)
                    self._short_template_mem = zeros(wav_len,
                                                     dtype=numpy.complex64)
                # Generate short waveform
                htilde = bank.generate_custom_size_waveform(
                    t_num, self.coarsematch_fmax, self.coarsematch_deltaf,
                    low_freq_cutoff=self.f_lower,
                    cached_mem=self._short_template_mem)
                self._short_template_id = t_num
                self._short_template_wav = htilde
            else:
                htilde = self._short_template_wav

            for inj_idx, inj in enumerate(gwstrain.injections.table):
                end_time = inj.geocent_end_time + \
                    1E-9 * inj.geocent_end_time_ns
                if end_time > seg_end_time or end_time < seg_start_time:
                    continue
                curr_inj = self.short_injections[inj.simulation_id]
                o, i = match(htilde, curr_inj, psd=red_psd,
                             low_frequency_cutoff=self.f_lower)
                if o > self.match_threshold:
                    break
            else:
                # Get's here if all injections are outside match threshold
                return False

        return True
