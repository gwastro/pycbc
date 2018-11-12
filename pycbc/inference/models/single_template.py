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
from pycbc.waveform import NoWaveformError, get_fd_waveform
from pycbc.types import Array
from pycbc.detector import Detector

from .base import BaseModel

# In this model we only calculate terms up to a constant.
# We are primarily interested in the posterior result

class SingleTemplate(BaseModel):
    r"""Model that assumes we know all the intrinsic parameters.

    This model assumes we know all the intrinsic parameters, and are only
    maximizing over the extrinsic ones. We also assume a dominant mode waveform
    approximant only and non-precessing.
    """
    name = 'single_template'

    def __init__(self, data, psds, f_lower, f_upper=None, **kwargs):
        # set up the boiler-plate attributes; note: we'll compute the
        # log evidence later
        super(SingleTemplate, self).__init__(**kwargs)          
        sample_rate = 32768
               
        # Generate template waveforms
        df = data[data.keys()[0]].delta_f
        hp, _ = get_fd_waveform(delta_f=df, **self.static_params)
 
        if f_upper is None:
            f_upper = len(data[data.keys()[0]]-1) * df
 
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
            snr, _, norm = pyfilter.matched_filter_core(hp, data[ifo], 
                                            psd=psds[ifo],                   
                                            low_frequency_cutoff=f_lower, 
                                            high_frequency_cutoff=f_upper)
                                            
            self.sh[ifo] = 4 * df * snr
            self.hh[ifo] = - 2.0 * df / norm

    def _lognl(self):
        """Computes the log likelihood assuming the data is noise.

        Since this is a constant for Gaussian noise, this is only computed once
        then stored.
        """
        # FIXME we won't bother calculating this for now since it is a 
        # constant in this model.
        return 0

    def loglr(self):
        r"""Computes the log likelihood ratio,

        Returns
        -------
        float
            The value of the log likelihood ratio.
        """
        # calculate <d-h|d-h> = <h|h> - 2<h|d> + <d|d> up to a constant
        p = self.current_params.copy()
        p.update(self.static_params)

        vloglr = 0
        for ifo in self.sh:
            fp, fc = self.det[ifo].antenna_pattern(p['ra'], p['dec'],
                                         p['polarization'], p['tc'])   
            dt = self.det[ifo].time_delay_from_earth_center(p['ra'],
                                         p['dec'], p['tc'])                          
            ip = numpy.cos(p['inclination'])
            ic = 0.5 * (1.0 + ip * ip)                     
            htf = (fp * ip + 1.0j * fc * ic) / p['distance']
            htf *= numpy.exp(1.0j * p['coa_phase'])
            
            vloglr += (self.sh[ifo].at_time(p['tc'] + dt) * htf).real
            vloglr += self.hh[ifo] * abs(htf) ** 2.0     
        return float(vloglr)

    def _loglikelihood(self):
        r"""Computes the log likelihood of the paramaters,

        """
        return self.loglr()
