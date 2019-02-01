# Copyright (C) 2018 Alex Nitz
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

"""This module provides model classes that assume the noise is Gaussian.
"""

import numpy
import scipy.special

from pycbc import filter as pyfilter
from pycbc.waveform import get_fd_waveform
from pycbc.detector import Detector

from .base_data import BaseDataModel

# In this model we only calculate terms up to a constant.
# We are primarily interested in the posterior result


class SingleTemplate(BaseDataModel):
    r"""Model that assumes we know all the intrinsic parameters.

    This model assumes we know all the intrinsic parameters, and are only
    maximizing over the extrinsic ones. We also assume a dominant mode waveform
    approximant only and non-precessing.
    """
    name = 'single_template'

    def __init__(self, data, psds,
                 low_frequency_cutoff=None,
                 high_frequency_cutoff=None,
                 sample_rate=32768,
                 **kwargs):

        super(SingleTemplate, self).__init__(data=data, **kwargs)

        if low_frequency_cutoff is not None:
            low_frequency_cutoff = float(low_frequency_cutoff)

        if high_frequency_cutoff is not None:
            high_frequency_cutoff = float(high_frequency_cutoff)

        # Generate template waveforms
        df = data[data.keys()[0]].delta_f
        p = self.static_params.copy()
        if 'distance' in p:
            p.pop('distance')
        if 'inclination' in p:
            p.pop('inclination')

        hp, _ = get_fd_waveform(delta_f=df, distance=1, inclination=0, **p)

        if high_frequency_cutoff is None:
            high_frequency_cutoff = len(data[data.keys()[0]]-1) * df

        # Extend data and template to high sample rate
        flen = int(sample_rate / df) / 2 + 1
        hp.resize(flen)
        for ifo in data:
            data[ifo].resize(flen)

        # Calculate high sample rate SNR time series
        self.sh = {}
        self.hh = {}
        self.det = {}
        for ifo in data:
            self.det[ifo] = Detector(ifo)
            snr, _, _ = pyfilter.matched_filter_core(
                hp, data[ifo],
                psd=psds[ifo],
                low_frequency_cutoff=low_frequency_cutoff,
                high_frequency_cutoff=high_frequency_cutoff)

            self.sh[ifo] = 4 * df * snr
            self.hh[ifo] = -0.5 * pyfilter.sigmasq(
                hp, psd=psds[ifo],
                low_frequency_cutoff=low_frequency_cutoff,
                high_frequency_cutoff=high_frequency_cutoff)
        self.time = None

    def _loglikelihood(self):
        return self.loglr

    def _lognl(self):
        return 0

    def _loglr(self):
        r"""Computes the log likelihood ratio

        Returns
        -------
        float
            The value of the log likelihood ratio.
        """
        # calculate <d-h|d-h> = <h|h> - 2<h|d> + <d|d> up to a constant
        p = self.current_params.copy()
        p.update(self.static_params)

        if self.time is None:
            self.time = p['tc']

        shloglr = hhloglr = 0
        for ifo in self.sh:
            fp, fc = self.det[ifo].antenna_pattern(p['ra'], p['dec'],
                                                   p['polarization'],
                                                   self.time)
            dt = self.det[ifo].time_delay_from_earth_center(p['ra'],
                                                            p['dec'],
                                                            self.time)
            ip = numpy.cos(p['inclination'])
            ic = 0.5 * (1.0 + ip * ip)
            htf = (fp * ip + 1.0j * fc * ic) / p['distance']

            sh = self.sh[ifo].at_time(p['tc'] + dt) * htf
            shloglr += sh
            hhloglr += self.hh[ifo] * abs(htf) ** 2.0

        vloglr = numpy.log(scipy.special.i0e(abs(shloglr)))
        vloglr += abs(shloglr) + hhloglr

        return float(vloglr)
