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
                 polarization_samples=None,
                 marginalize_phase=True,
                 **kwargs):

        #polarization array to marginalize over if num_samples given
        self.pflag = 0
        self.polarization = None
        if polarization_samples is not None:
            self.polarization = numpy.linspace(0, 2*numpy.pi,
                                               int(polarization_samples))
            self.pflag = 1
        marg_vector = 'polarization' if self.pflag else False

        variable_params, kwargs = self.setup_distance_marginalization(
                                   variable_params,
                                   marginalize_phase=marginalize_phase,
                                   marginalize_vector=marg_vector,
                                   marginalize_vector_params=self.polarization,
                                   **kwargs)
        super(SingleTemplate, self).__init__(
            variable_params, data, low_frequency_cutoff, **kwargs)

        # Generate template waveforms
        df = data[self.detectors[0]].delta_f
        p = self.static_params.copy()
        if 'distance' in p:
            _ = p.pop('distance')
        if 'inclination' in p:
            _ = p.pop('inclination')
        hp, _ = get_fd_waveform(delta_f=df, distance=1, inclination=0, **p)

        # Extend template to high sample rate
        flen = int(int(sample_rate) / df) / 2 + 1
        hp.resize(flen)

        # Calculate high sample rate SNR time series
        self.sh = {}
        self.hh = {}
        self.det = {}
        for ifo in self.data:
            flow = self.kmin[ifo] * df
            fhigh = self.kmax[ifo] * df
            # Extend data to high sample rate
            self.data[ifo].resize(flen)
            self.det[ifo] = Detector(ifo)
            snr, _, _ = pyfilter.matched_filter_core(
                hp, self.data[ifo],
                psd=self.psds[ifo],
                low_frequency_cutoff=flow,
                high_frequency_cutoff=fhigh)

            self.sh[ifo] = 4 * df * snr
            self.hh[ifo] = pyfilter.sigmasq(
                hp, psd=self.psds[ifo],
                low_frequency_cutoff=flow,
                high_frequency_cutoff=fhigh)
        self.time = None

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
        if self.pflag == 0:
            polarization = p['polarization']
        elif self.pflag == 1:
            polarization = self.polarization

        if self.time is None:
            self.time = p['tc']

        phase = 1
        if 'coa_phase' in p:
            phase = numpy.exp(1.0j * 2 * p['coa_phase'])

        sh_total = hh_total = 0
        for ifo in self.sh:
            fp, fc = self.det[ifo].antenna_pattern(p['ra'], p['dec'],
                                                   polarization, self.time)
            dt = self.det[ifo].time_delay_from_earth_center(p['ra'], p['dec'],
                                                            self.time)
            ic = numpy.cos(p['inclination'])
            ip = 0.5 * (1.0 + ic * ic)
            htf = (fp * ip + 1.0j * fc * ic) / p['distance'] * phase
            sh = self.sh[ifo].at_time(p['tc'] + dt) * htf
            sh_total += sh
            hh_total += self.hh[ifo] * abs(htf) ** 2.0

        return self.marginalize_loglr(sh_total, hh_total)
