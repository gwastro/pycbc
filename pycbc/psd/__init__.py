#!/usr/bin/python
# Copyright (C) 2012 Alex Nitz, Andrew Miller, Tito Dal Canton
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


from pycbc.psd.read import *
from pycbc.psd.analytical import *
from pycbc.psd.estimate import *


def from_cli(opt, parser, strain=None, dyn_range_factor=None):
    """Parses the CLI options related to the noise PSD and returns a
    FrequencySeries with the corresponding PSD. If necessary, the PSD is
    linearly interpolated to achieve the resolution specified in the CLI.

    Parameters
    ----------
    opt : object
        Result of parsing the CLI with OptionParser, or any object with the
        required attributes (psd_model, psd_file, psd_estimation,
        psd_segment_length, psd_segment_stride, psd_inverse_length, sample_rate
        and low_frequency_cutoff).
    parser : object
        OptionParser instance.
    strain : {None, TimeSeries}
        Time series containing the data from which the PSD should be measured,
        when psd_estimation is in use.
    dyn_range_factor : {None, float}
        For PSDs taken from models or text files, if `dyn_range_factor` is
        not None, then the PSD is multiplied by `dyn_range_factor` ** 2.

    Returns
    -------
    psd : FrequencySeries
        The frequency series containing the PSD.
    """
    f_low = opt.low_frequency_cutoff
    delta_f = 1. / opt.segment_length

    if (opt.psd_model or opt.psd_file) and not opt.psd_estimation:
        # PSD from lalsimulation or file
        psd_size = int(opt.segment_length * opt.sample_rate) / 2 + 1
        if opt.psd_model:
            psd = from_string(opt.psd_model, psd_size, delta_f, f_low)
        elif opt.psd_file:
            psd = from_asd_txt(opt.psd_file, psd_size, delta_f, f_low)
        if dyn_range_factor:
            psd *= dyn_range_factor ** 2
    elif opt.psd_estimation and not (opt.psd_model or opt.psd_file):
        # estimate PSD from data
        if not opt.psd_segment_length or not opt.psd_segment_stride:
            parser.error('PSD estimation requires --psd-segment-length '
                         'and --psd-segment-stride')
        psd = welch(strain, avg_method=opt.psd_estimation,
                    seg_len=int(opt.psd_segment_length * opt.sample_rate),
                    seg_stride=int(opt.psd_segment_stride * opt.sample_rate))
        if delta_f != psd.delta_f:
            psd = interpolate(psd, delta_f)
        if opt.psd_inverse_length:
            psd = inverse_spectrum_truncation(
                    psd, int(opt.sample_rate * opt.psd_inverse_length),
                    low_frequency_cutoff=f_low)
    else:
        parser.error('Please specify either --psd-model, --psd-file or '
                     '--psd-estimation')

    return psd
