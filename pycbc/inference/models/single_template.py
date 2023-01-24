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

import logging
import numpy
import itertools

from pycbc import filter as pyfilter
from pycbc.waveform import get_fd_waveform
from pycbc.detector import Detector

from .gaussian_noise import BaseGaussianNoise
from .tools import DistMarg


class SingleTemplate(DistMarg, BaseGaussianNoise):
    r"""Model that assumes we know all the intrinsic parameters.

    This model assumes we know all the intrinsic parameters, and are only
    maximizing over the extrinsic ones. We also assume a dominant mode waveform
    approximant only and non-precessing.


    Parameters
    ----------
    variable_params : (tuple of) string(s)
        A tuple of parameter names that will be varied.
    data : dict
        A dictionary of data, in which the keys are the detector names and the
        values are the data (assumed to be unwhitened). All data must have the
        same frequency resolution.
    low_frequency_cutoff : dict
        A dictionary of starting frequencies, in which the keys are the
        detector names and the values are the starting frequencies for the
        respective detectors to be used for computing inner products.
    sample_rate : int, optional
        The sample rate to use. Default is 32768.
    polarization_samples: int, optional
        Parameter to specify how finely to marginalize over polarization angle.
        If None, then polarization must be a parameter.
    \**kwargs :
        All other keyword arguments are passed to
        :py:class:`BaseGaussianNoise`; see that class for details.
    """
    name = 'single_template'

    def __init__(self, variable_params, data, low_frequency_cutoff,
                 sample_rate=32768,
                 marginalize_phase=True,
                 **kwargs):
        variable_params, kwargs = self.setup_marginalization(
                                   variable_params,
                                   marginalize_phase=marginalize_phase,
                                   **kwargs)
        super(SingleTemplate, self).__init__(
            variable_params, data, low_frequency_cutoff, **kwargs)

        sample_rate = float(sample_rate)

        # Generate template waveforms
        df = data[self.detectors[0]].delta_f
        self.df = df
        p = self.static_params.copy()
        for k in self.static_params:
            if p[k] == 'REPLACE':
                p.pop(k)
        if 'distance' in p:
            _ = p.pop('distance')
        if 'inclination' in p:
            _ = p.pop('inclination')

        hp, _ = get_fd_waveform(delta_f=df, distance=1, inclination=0, **p)

        # Extend template to high sample rate
        flen = int(round(sample_rate / df) / 2 + 1)
        hp.resize(flen)

        # Calculate high sample rate SNR time series
        self.sh = {}
        self.hh = {}
        self.snr = {}
        self.det = {}
        for ifo in self.data:
            flow = self.kmin[ifo] * df
            fhigh = self.kmax[ifo] * df
            # Extend data to high sample rate
            self.data[ifo].resize(flen)
            self.det[ifo] = Detector(ifo)
            snr, _, norm = pyfilter.matched_filter_core(
                hp, self.data[ifo],
                psd=self.psds[ifo],
                low_frequency_cutoff=flow,
                high_frequency_cutoff=fhigh)

            self.sh[ifo] = 4 * df * snr
            self.snr[ifo] = snr * norm

            self.hh[ifo] = pyfilter.sigmasq(
                hp, psd=self.psds[ifo],
                low_frequency_cutoff=flow,
                high_frequency_cutoff=fhigh)

        self.waveform = hp
        self.htfs = {}  # Waveform phase / distance transformation factors
        self.dts = {}

        # Retrict to analyzing around peaks if chosen and choose what
        # ifos to draw from
        self.setup_peak_lock(snrs=self.snr,
                             sample_rate=sample_rate,
                             **kwargs)
        self.draw_ifos(self.snr)

    @property
    def multi_signal_support(self):
        """ The list of classes that this model supports in a multi-signal
        likelihood
        """
        # Check if this model *can* be included in a multi-signal model.
        # All marginalizations must currently be disabled to work!
        if (self.marginalize_vector_params or
            self.marginalize_distance or
            self.marginalize_phase):
            logging.info("Cannot use single template model inside of"
                         "multi_signal if marginalizations are enabled")
        return [type(self)]

    def calculate_hihjs(self, models):
        """ Pre-calculate the hihj inner products on a grid
        """
        self.hihj = {}
        for m1, m2 in itertools.combinations(models, 2):
            self.hihj[(m1, m2)] = {}
            h1 = m1.waveform
            h2 = m2.waveform
            for ifo in self.data:
                flow = self.kmin[ifo] * self.df
                fhigh = self.kmax[ifo] * self.df
                h1h2, _, _ = pyfilter.matched_filter_core(
                h1, h2,
                psd=self.psds[ifo],
                low_frequency_cutoff=flow,
                high_frequency_cutoff=fhigh)
                self.hihj[(m1, m2)][ifo] = 4 * self.df * h1h2

    def multi_loglikelihood(self, models):
        """ Calculate a multi-model (signal) likelihood
        """
        models = [self] + models
        loglr = 0
        # handle sum[<d|h_i> - 0.5 <h_i|h_i>]
        for m in models:
            loglr += m.loglr

        if not hasattr(self, 'hihj'):
            self.calculate_hihjs(models)

        # finally add in the lognl term from this model
        for m1, m2 in itertools.combinations(models, 2):
            for det in self.data:
                hihj_vec = self.hihj[(m1, m2)][det]
                dt = m1.dts[det] - m2.dts[det] + hihj_vec.start_time
                if dt < hihj_vec.start_time:
                    dt += hihj_vec.duration

                h1h2 = hihj_vec.at_time(dt, nearest_sample=True)
                h1h2 *= m1.htfs[det] * m2.htfs[det].conj()
                loglr += - h1h2.real # This is -0.5 * re(<h1|h2> + <h2|h1>)
        return loglr + self.lognl

    def _loglr(self):
        r"""Computes the log likelihood ratio

        Returns
        -------
        float
            The value of the log likelihood ratio.
        """
        # calculate <d-h|d-h> = <h|h> - 2<h|d> + <d|d> up to a constant
        p = self.current_params

        phase = 1
        if 'coa_phase' in p:
            phase = numpy.exp(-1.0j * 2 * p['coa_phase'])

        sh_total = hh_total = 0

        ic = numpy.cos(p['inclination'])
        ip = 0.5 * (1.0 + ic * ic)
        pol_phase = numpy.exp(-2.0j * p['polarization'])

        self.snr_draw(self.snr)

        for ifo in self.sh:
            dt = self.det[ifo].time_delay_from_earth_center(p['ra'], p['dec'],
                                                            p['tc'])
            self.dts[ifo] = p['tc'] + dt

            fp, fc = self.det[ifo].antenna_pattern(p['ra'], p['dec'],
                                                   0, p['tc'])
            f = (fp + 1.0j * fc) * pol_phase

            # Note, this includes complex conjugation already
            # as our stored inner products were hp* x data
            htf = (f.real * ip + 1.0j * f.imag * ic) / p['distance'] * phase
            self.htfs[ifo] = htf
            sh = self.sh[ifo].at_time(self.dts[ifo], interpolate='quadratic')
            sh_total += sh * htf
            hh_total += self.hh[ifo] * abs(htf) ** 2.0

        loglr = self.marginalize_loglr(sh_total, hh_total)
        return loglr
