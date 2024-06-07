# Copyright (C) 2018  Collin Capano
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

"""
This modules provides models that have analytic solutions for the
log likelihood.
"""

import copy
import logging

import pycbc.types

from .base import BaseModel

import pycbc.psd

from pycbc.waveform.pre_merger_waveform import (
    pre_process_data_lisa_pre_merger,
    generate_waveform_lisa_pre_merger,
)
from pycbc.psd.lisa_pre_merger import generate_pre_merger_psds
from pycbc.waveform.waveform import parse_mode_array
from pycbc.waveform.utils import apply_fseries_time_shift
from .tools import marginalize_likelihood


class LISAPreMergerModel(BaseModel):
    r"""Model for pre-merger inference in LISA.

    Parameters
    ----------
    variable_params : (tuple of) string(s)
        A tuple of parameter names that will be varied.
    static_params: Dict[str: Any]
        Dictionary of static parameters used for waveform generation.
    psd_file : str
        Path to the PSD file. Uses the same PSD file for LISA_A and LISA_E
        channels.
    **kwargs :
        All other keyword arguments are passed to ``BaseModel``.
    """
    name = "lisa_pre_merger"

    def __init__(
        self,
        variable_params,
        static_params=None,
        psd_file=None,
        **kwargs
    ):
        # Pop relevant values from kwargs
        cutoff_time = int(kwargs.pop('cutoff_time'))
        kernel_length = int(kwargs.pop('kernel_length'))
        psd_kernel_length = int(kwargs.pop('psd_kernel_length'))
        window_length = int(kwargs.pop('window_length'))
        extra_forward_zeroes = int(kwargs.pop('extra_forward_zeroes'))
        tlen = int(kwargs.pop('tlen'))
        sample_rate = float(kwargs.pop('sample_rate'))
        data_file = kwargs.pop('data_file')
        
        # set up base likelihood parameters
        super().__init__(variable_params, **kwargs)

        self.static_params = parse_mode_array(static_params)

        if psd_file is None:
            raise ValueError("Must specify a PSD file!")

        # Zero phase PSDs for whitening
        # Only store the frequency-domain PSDs
        logging.info("Generating pre-merger PSDs")
        self.whitening_psds = {}
        self.whitening_psds['LISA_A'] = generate_pre_merger_psds(
            psd_file,
            sample_rate=sample_rate,
            duration=tlen,
            kernel_length=psd_kernel_length,
        )["FD"]
        self.whitening_psds['LISA_E'] = generate_pre_merger_psds(
            psd_file,
            sample_rate=sample_rate,
            duration=tlen,
            kernel_length=psd_kernel_length,
        )["FD"]

        # Store data for doing likelihoods.
        self.kernel_length = kernel_length
        self.window_length = window_length
        self.sample_rate = sample_rate
        self.cutoff_time = cutoff_time
        self.extra_forward_zeroes = extra_forward_zeroes

        # Load the data from the file
        data = {}
        for channel in ["LISA_A", "LISA_E"]:
            data[channel] = pycbc.types.timeseries.load_timeseries(
                data_file,
                group=f"/{channel}",
            )

        # Pre-process the pre-merger data
        # Returns time-domain data
        # Uses UIDs: 4235(0), 4236(0)
        logging.info("Pre-processing pre-merger data")
        pre_merger_data = pre_process_data_lisa_pre_merger(
            data,
            sample_rate=sample_rate,
            psds_for_whitening=self.whitening_psds,
            window_length=self.window_length, 
            cutoff_time=self.cutoff_time,
            extra_forward_zeroes=self.extra_forward_zeroes,
        )

        self.lisa_a_strain = pre_merger_data["LISA_A"]
        self.lisa_e_strain = pre_merger_data["LISA_E"]

        # Frequency-domain data for computing log-likelihood
        self.lisa_a_strain_fd = pycbc.strain.strain.execute_cached_fft(
            self.lisa_a_strain,
            copy_output=True,
            uid=3223965
        )
        self.lisa_e_strain_fd = pycbc.strain.strain.execute_cached_fft(
            self.lisa_e_strain,
            copy_output=True,
            uid=3223967
        )
        # Data epoch
        self._epoch = self.lisa_a_strain_fd._epoch

    def get_waveforms(self, params):
        """Generate the waveforms given the parameters.
        
        Note: `params` should already include the static parameters.
        """
        # Generate the pre-merger waveform
        # These waveforms are whitened
        # Uses UIDs: 1234(0), 1235(0), 1236(0), 1237(0)
        ws = generate_waveform_lisa_pre_merger(
            params,
            psds_for_whitening=self.whitening_psds,
            window_length=self.window_length,
            sample_rate=self.sample_rate,
            cutoff_time=self.cutoff_time,
            extra_forward_zeroes=self.extra_forward_zeroes,
        )

        wf = {}
        # Adjust epoch to match data and shift merger to the
        # correct time.
        # Can safely set copy=False since ws won't be used again.
        dt = params["tc"] - self._epoch
        for channel in ws.keys():
            wf[channel] = apply_fseries_time_shift(
                ws[channel], dt, copy=False,
            )
            wf[channel]._epoch = self._epoch
        return wf

    def _loglikelihood(self):
        """Compute the pre-merger log-likelihood."""
        cparams = copy.deepcopy(self.static_params)
        cparams.update(self.current_params)

        # Generate the waveforms
        wforms = self.get_waveforms(cparams)

        # Compute <h|d> for each channel
        snr_A = pycbc.filter.overlap_cplx(
            wforms["LISA_A"],
            self.lisa_a_strain_fd,
            normalized=False,
        )
        snr_E = pycbc.filter.overlap_cplx(
            wforms["LISA_E"],
            self.lisa_e_strain_fd,
            normalized=False,
        )
        # Compute <h|h> for each channel
        a_norm = pycbc.filter.sigmasq(wforms["LISA_A"])
        e_norm = pycbc.filter.sigmasq(wforms["LISA_E"])

        hs = snr_A + snr_E
        hh = (a_norm + e_norm)

        return marginalize_likelihood(hs, hh, phase=False)
