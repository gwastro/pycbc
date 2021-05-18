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

"""This module provides model classes that assume the noise is Gaussian and
introduces a gate to remove given times from the data, using the inpainting
method to fill the removed part such that it does not enter the likelihood.
"""

import numpy
from .gaussian_noise import (BaseGaussianNoise, create_waveform_generator)
from pycbc.waveform import (NoWaveformError, FailedWaveformError)
from pycbc.types import FrequencySeries
from .base_data import BaseDataModel
from pycbc.detector import Detector
from pycbc.pnutils import hybrid_meco_frequency
from pycbc.waveform.utils import time_from_frequencyseries


class GatedGaussianNoise(BaseGaussianNoise):
    r"""Model that applies a time domain gate, assuming stationary Gaussian
    noise.
    """
    name = 'gated_gaussian_noise'

    def __init__(self, variable_params, data, low_frequency_cutoff, psds=None,
                 high_frequency_cutoff=None, normalize=False,
                 static_params=None, **kwargs):
        # we'll want the time-domain data, so store that
        self._td_data = {}
        # set up the boiler-plate attributes
        super(GatedGaussianNoise, self).__init__(
            variable_params, data, low_frequency_cutoff, psds=psds,
            high_frequency_cutoff=high_frequency_cutoff, normalize=normalize,
            static_params=static_params, **kwargs)
        # caches for debuging
        self.current_wfs = {}
        self.current_gated_wfs = {}
        self.current_gated_data = {}
        # create the waveform generator
        self.waveform_generator = create_waveform_generator(
            self.variable_params, self.data,
            waveform_transforms=self.waveform_transforms,
            recalibration=self.recalibration,
            gates=self.gates, **self.static_params)

    @BaseDataModel.data.setter
    def data(self, data):
        """Store a copy of the FD and TD data."""
        BaseDataModel.data.fset(self, data)
        # store the td version
        self._td_data = {det: d.to_timeseries() for det, d in data.items()}

    @property
    def td_data(self):
        return self._td_data

    @property
    def psds(self):
        """Dictionary of detectors -> PSD frequency series.
        If no PSD was provided for a detector, this will just be a frequency
        series of ones.
        """
        return self._psds

    @BaseGaussianNoise.psds.setter
    def psds(self, psds):
        """Sets the psds, and calculates the weight and norm from them.
        The data and the low and high frequency cutoffs must be set first.
        """
        # check that the data has been set
        if self._data is None:
            raise ValueError("No data set")
        if self._f_lower is None:
            raise ValueError("low frequency cutoff not set")
        if self._f_upper is None:
            raise ValueError("high frequency cutoff not set")
        # make sure the relevant caches are cleared
        self._psds.clear()
        self._invpsds.clear()
        self._weight.clear()
        self._whitened_data.clear()
        for det, d in self._data.items():
            if psds is None:
                # No psd means assume white PSD
                p = FrequencySeries(numpy.ones(int(self._N/2+1)),
                                    delta_f=d.delta_f)
            else:
                # copy for storage
                p = psds[det].copy()
            self._psds[det] = p
            # we'll store the weight to apply to the inner product
            invp = 1./p
            self._invpsds[det] = invp
            self._weight[det] = numpy.sqrt(4 * invp.delta_f * invp)
            self._whitened_data[det] = d.copy()
            self._whitened_data[det] *= self._weight[det]

    def det_lognorm(self, det):
        # FIXME: just returning 0 for now, but should be the determinant
        # of the truncated covariance matrix
        return 0.

    @property
    def normalize(self):
        """Determines if the loglikelihood includes the normalization term.
        """
        return self._normalize

    @normalize.setter
    def normalize(self, normalize):
        """Clears the current stats if the normalization state is changed.
        """
        self._normalize = normalize

    @property
    def _extra_stats(self):
        """Adds ``lognl``, plus ```optimal_snrsq`` in each
        detector."""
        return []

    def _nowaveform_logl(self):
        """Convenience function to set logl values if no waveform generated.
        """
        return -numpy.inf

    def _loglr(self):
        r"""Computes the log likelihood ratio.
        Returns
        -------
        float
            The value of the log likelihood ratio evaluated at the given point.
        """
        return self.loglikelihood - self.lognl

    def get_gate_times(self):
        """Gets the time to apply a gate based on the current sky position.
        Returns
        -------
        dict :
            Dictionary of detector names -> (gate start, gate width)
        """
        params = self.current_params
        try:
            gatefunc = self.current_params['gatefunc']
        except KeyError:
            gatefunc = None
        if gatefunc == 'hmeco':
            return self.get_gate_times_hmeco()
        # gate input for ringdown analysis which consideres a start time
        # and an end time
        gatestart = params['t_gate_start']
        gateend = params['t_gate_end']
        dgate = gateend-gatestart
        # we'll need the sky location for determining time shifts
        ra = self.current_params['ra']
        dec = self.current_params['dec']
        gatetimes = {}
        for det, invpsd in self._invpsds.items():
            thisdet = Detector(det)
            # account for the time delay between the waveforms of the
            # different detectors
            gatestartdelay = gatestart + thisdet.time_delay_from_earth_center(
                ra, dec, gatestart)
            gateenddelay = gateend + thisdet.time_delay_from_earth_center(
                ra, dec, gateend)
            dgatedelay = gateenddelay - gatestartdelay
            gatetimes[det] = (gatestartdelay, dgatedelay)
        return gatetimes

    def get_gate_times_hmeco(self):
        """Gets the time to apply a gate based on the current sky position.
        Returns
        -------
        dict :
            Dictionary of detector names -> (gate start, gate width)
        """
        params = self.current_params
        # gate input for ringdown analysis which consideres a start time
        # and an end time
        dgate = params['gate_window']
        gatestart = params['t_gate_start']

        # generate the template waveform
        try:
            wfs = self.waveform_generator.generate(**params)
        except NoWaveformError:
            return self._nowaveform_logl()
        except FailedWaveformError as e:
            if self.ignore_failed_waveforms:
                return self._nowaveform_logl()
            else:
                raise e

        spin1 = params['spin1z']
        spin2 = params['spin2z']
        meco_f = hybrid_meco_frequency(params['mass1'], params['mass2'],
                                       spin1, spin2, qm1=None, qm2=None)
        gatetimes = {}
        for det, h in wfs.items():
            invpsd = self._invpsds[det]
            h.resize(len(invpsd))
            ht = h.to_timeseries()
            f_low = int((self._f_lower[det]+1)/h.delta_f)
            sample_freqs = h.sample_frequencies[f_low:].numpy()
            f_idx = numpy.where(sample_freqs <= meco_f)[0][-1]
            # find time corresponding to meco frequency
            t_from_freq = time_from_frequencyseries(
                h[f_low:], sample_frequencies=sample_freqs)
            if t_from_freq[f_idx]>0:
                Gatestartdelay = t_from_freq[f_idx] + float(t_from_freq.epoch)
            else:
                Gatestartdelay = t_from_freq[f_idx] + ht.sample_times[-1]
            gatestartdelay = min(Gatestartdelay, params['t_gate_start'])
            gatetimes[det] = (gatestartdelay, dgate)
        return gatetimes


    def _lognl(self):
        """Calculates the log of the noise likelihood.
        """
        params = self.current_params
        # clear variables
        lognl = 0.
        self._det_lognls.clear()
        # get the times of the gates
        gate_times = self.get_gate_times()
        self.current_nproj = {}
        for det, invpsd in self._invpsds.items():
            norm = self.det_lognorm(det)
            gatestartdelay, dgatedelay = gate_times[det]
            # we always filter the entire segment starting from kmin, since the
            # gated series may have high frequency components
            slc = slice(self._kmin[det], self._kmax[det])
            # gate the data
            data = self.td_data[det]
            gated_dt = data.gate(gatestartdelay + dgatedelay/2,
                                 window=dgatedelay/2, copy=True,
                                 invpsd=invpsd, method='paint')
            self.current_nproj[det] = (gated_dt.proj, gated_dt.projslc)
            # convert to the frequency series
            gated_d = gated_dt.to_frequencyseries()
            # overwhiten
            gated_d *= invpsd
            d = self.data[det]
            # inner product
            ip = 4 * invpsd.delta_f * d[slc].inner(gated_d[slc]).real # <d, d>
            dd = norm - 0.5*ip
            # store
            self._det_lognls[det] = dd
            lognl += dd
        return float(lognl)

    def det_lognl(self, det):
        # make sure lognl has been called
        _ = self._trytoget('lognl', self._lognl)
        # the det_lognls dict should have been updated, so can call it now
        return self._det_lognls[det]

    def _loglikelihood(self):
        r"""Computes the log likelihood after removing the power within the
        given time window,
        .. math::
            \log p(d|\Theta) = -\frac{1}{2} \sum_i
             \left< d_i - h_i(\Theta) | d_i - h_i(\Theta) \right>,
        at the current parameter values :math:`\Theta`.
        Returns
        -------
        float
            The value of the log likelihood.
        """
        params = self.current_params
        # generate the template waveform
        try:
            wfs = self.waveform_generator.generate(**params)
        except NoWaveformError:
            return self._nowaveform_logl()
        except FailedWaveformError as e:
            if self.ignore_failed_waveforms:
                return self._nowaveform_logl()
            else:
                raise e
        # get the times of the gates
        gate_times = self.get_gate_times()
        # clear variables
        logl = 0.
        self.current_wfs = wfs
        self.current_gated_wfs.clear()
        self.current_gated_data.clear()
        self.current_proj = {}
        for det, h in wfs.items():
            invpsd = self._invpsds[det]
            norm = self.det_lognorm(det)
            gatestartdelay, dgatedelay = gate_times[det]
            # we always filter the entire segment starting from kmin, since the
            # gated series may have high frequency components
            slc = slice(self._kmin[det], self._kmax[det])
            # calculate the residual
            data = self.td_data[det]
            h.resize(len(self.data[det]))
            ht = h.to_timeseries()
            res = data - ht
            rtilde = res.to_frequencyseries()
            gated_res = res.gate(gatestartdelay + dgatedelay/2,
                                 window=dgatedelay/2, copy=True,
                                 invpsd=invpsd, method='paint')
            self.current_proj[det] = (gated_res.proj, gated_res.projslc)
            gated_rtilde = gated_res.to_frequencyseries()
            # overwhiten
            gated_rtilde *= invpsd
            rr = 4 * invpsd.delta_f * rtilde[slc].inner(gated_rtilde[slc]).real
            logl += norm - 0.5*rr
        return float(logl)

    def get_waveforms(self, whiten=False):
        params = self.current_params
        wfs = self.waveform_generator.generate(**params)
        if whiten:
            for det, h in wfs.items():
                invpsd = self._invpsds[det]
                if whiten == 2:
                    h *= invpsd
                else:
                    h *= invpsd**0.5
                wfs[det] = h
        return wfs

    def get_gated_waveforms(self, whiten=False):
        params = self.current_params
        wfs = self.waveform_generator.generate(**params)
        gate_times = self.get_gate_times()
        for det, h in wfs.items():
            invpsd = self._invpsds[det]
            gatestartdelay, dgatedelay = gate_times[det]
            data = self.td_data[det]
            ht = h.to_timeseries()
            ht = ht.gate(gatestartdelay + dgatedelay/2,
                         window=dgatedelay/2, copy=False,
                         invpsd=invpsd, method='paint')
            h = ht.to_frequencyseries()
            if whiten == 2:
                h *= invpsd
            elif whiten:
                h *= invpsd**0.5
            wfs[det] = h
        return wfs

    def get_residual(self, whiten=False):
        params = self.current_params
        wfs = self.waveform_generator.generate(**params)
        out = {}
        for det, h in wfs.items():
            invpsd = self._invpsds[det]
            d = self.data[det]
            res = d - h
            if whiten == 2:
                res *= invpsd
            elif whiten:
                res *= invpsd**0.5
            out[det] = res
        return out

    def get_gated_residual(self, whiten=False):
        params = self.current_params
        wfs = self.waveform_generator.generate(**params)
        gate_times = self.get_gate_times()
        out = {}
        for det, h in wfs.items():
            invpsd = self._invpsds[det]
            gatestartdelay, dgatedelay = gate_times[det]
            data = self.td_data[det]
            ht = h.to_timeseries()
            res = data - ht
            res = res.gate(gatestartdelay + dgatedelay/2,
                           window=dgatedelay/2, copy=True,
                           invpsd=invpsd, method='paint')
            res = res.to_frequencyseries()
            if whiten == 2:
                res *= invpsd
            elif whiten:
                res *= invpsd**0.5
            out[det] = res
        return out

    def get_data(self, whiten=False):
        data = {det: d.copy() for det, d in self.data.items()}
        if whiten:
            for det, dtilde in data.items():
                invpsd = self._invpsds[det]
                if whiten == 2:
                    dtilde *= invpsd
                else:
                    dtilde *= invpsd**0.5
                data[det] = dtilde
        return data

    def get_gated_data(self, whiten=False):
        gate_times = self.get_gate_times()
        out = {}
        for det, d in self.td_data.items():
            invpsd = self._invpsds[det]
            gatestartdelay, dgatedelay = gate_times[det]
            d = d.gate(gatestartdelay + dgatedelay/2,
                       window=dgatedelay/2, copy=True,
                       invpsd=invpsd, method='paint')
            dtilde = d.to_frequencyseries()
            if whiten == 2:
                dtilde *= invpsd
            elif whiten:
                dtilde *= invpsd**0.5
            out[det] = dtilde
        return out

    def write_metadata(self, fp):
        """Adds writing the psds.
        The analyzed detectors, their analysis segments, and the segments
        used for psd estimation are written to the file's ``attrs``, as
        ``analyzed_detectors``, ``{{detector}}_analysis_segment``, and
        ``{{detector}}_psd_segment``, respectively.
        Parameters
        ----------
        fp : pycbc.inference.io.BaseInferenceFile instance
            The inference file to write to.
        """
        super(BaseGaussianNoise, self).write_metadata(fp)
        # write the analyzed detectors and times
        fp.attrs['analyzed_detectors'] = self.detectors
        for det, data in self.data.items():
            key = '{}_analysis_segment'.format(det)
            fp.attrs[key] = [float(data.start_time), float(data.end_time)]
        if self._psds is not None:
            fp.write_psd(self._psds)
        # write the times used for psd estimation (if they were provided)
        for det in self.psd_segments:
            key = '{}_psd_segment'.format(det)
            fp.attrs[key] = list(map(float, self.psd_segments[det]))
        try:
            attrs = fp[fp.samples_group].attrs
        except KeyError:
            # group doesn't exist, create it
            fp.create_group(fp.samples_group)
            attrs = fp[fp.samples_group].attrs
        for det in self.detectors:
            # Save each IFO's low frequency cutoff used in the likelihood
            # computation as an attribute
            fp.attrs['{}_likelihood_low_freq'.format(det)] = self._f_lower[det]
            # Save the IFO's high frequency cutoff used in the likelihood
            # computation as an attribute if one was provided the user
            if self._f_upper[det] is not None:
                fp.attrs['{}_likelihood_high_freq'.format(det)] = \
                                                        self._f_upper[det]
