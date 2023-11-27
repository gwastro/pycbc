# Copyright (C) 2020  Collin Capano and Shilpa Kastha
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

from abc import abstractmethod
import logging
import numpy
import scipy
from scipy import special

from pycbc.waveform import (NoWaveformError, FailedWaveformError)
from pycbc.types import FrequencySeries
from pycbc.detector import Detector
from pycbc.pnutils import hybrid_meco_frequency
from pycbc.psd import interpolate
from pycbc import types
from pycbc.waveform.utils import time_from_frequencyseries
from pycbc.waveform import generator
from pycbc.filter import highpass
from .gaussian_noise import (BaseGaussianNoise, create_waveform_generator)
from .base_data import BaseDataModel
from .data_utils import fd_data_from_strain_dict


class BaseGatedGaussian(BaseGaussianNoise):
    r"""Base model for gated gaussian.

    Provides additional routines for applying a time-domain gate to data.
    See :py:class:`GatedGaussianNoise` for more details.
    """
    def __init__(self, variable_params, data, low_frequency_cutoff, psds=None,
                 high_frequency_cutoff=None, normalize=False,
                 static_params=None, highpass_waveforms=False,
                 zero_before_gate=False, zero_after_gate=False,
                 **kwargs):
        # we'll want the time-domain data, so store that
        self._td_data = {}
        # cache the overwhitened data
        self._overwhitened_data = {}
        # cache the current gated data
        self._gated_data = {}
        # cache terms related to normalization and gating 
        self._invasds = {}
        self._Rss = {}
        self._lognorm = {}
        # cache samples and linear regression for determinant extrapolation
        self._cov_samples = {}
        self._cov_regressions = {}
        # highpass waveforms with the given frequency
        self.highpass_waveforms = highpass_waveforms
        if self.highpass_waveforms:
            logging.info("Will highpass waveforms at %f Hz",
                         highpass_waveforms)
        self.zero_before_gate = zero_before_gate
        self.zero_after_gate = zero_after_gate
        # set up the boiler-plate attributes
        super().__init__(
            variable_params, data, low_frequency_cutoff, psds=psds,
            high_frequency_cutoff=high_frequency_cutoff, normalize=normalize,
            static_params=static_params, **kwargs)

    @classmethod
    def from_config(cls, cp, data_section='data', data=None, psds=None,
                    **kwargs):
        """Adds addiotional keyword arguments based on config file.
        
        Additional keyword arguments are:

           * ``highpass_waveforms`` : waveforms will be highpassed.
           * ``zero_before_gate``: data/waveforms will be zeroed out at all
             times before the gate.
           * ``zero_after_gate``: data/waveforms will zeroed out at all times
             after the gate.
        """
        if cp.has_option('model', 'zero-before-gate'):
            kwargs['zero_before_gate'] = True
        if cp.has_option('model', 'zero-after-gate'):
            kwargs['zero_after_gate'] = True
        if cp.has_option(data_section, 'strain-high-pass') and \
            'highpass_waveforms' not in kwargs:
            kwargs['highpass_waveforms'] = float(cp.get(data_section,
                                                        'strain-high-pass'))
        return super().from_config(cp, data_section=data_section,
                                   data=data, psds=psds,
                                   **kwargs)

    @BaseDataModel.data.setter
    def data(self, data):
        """Store a copy of the FD and TD data."""
        BaseDataModel.data.fset(self, data)
        # store the td version
        self._td_data = {det: d.to_timeseries() for det, d in data.items()}

    @property
    def td_data(self):
        """The data in the time domain."""
        return self._td_data

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
        self._invasds.clear()
        self._gated_data.clear()
        # store the psds
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
            self._invasds[det] = invp**0.5
            # store the autocorrelation function and covariance matrix for each detector
            Rss = p.astype(types.complex_same_precision_as(p)).to_timeseries()
            self._Rss[det] = Rss
            # calculate and store the linear regressions to extrapolate determinant values
            if self.normalize:
                cov = scipy.linalg.toeplitz(Rss/2) # full covariance matrix
                samples, fit = self.logdet_fit(cov, p)
                self._cov_samples[det] = samples
                self._cov_regressions[det] = fit
        self._overwhitened_data = self.whiten(self.data, 2, inplace=False)

    def logdet_fit(self, cov, p):
        """Construct a linear regression from a sample of truncated covariance matrices.
        
        Returns the sample points used for linear fit generation as well as the linear fit parameters.
        """
        # initialize lists for matrix sizes and determinants
        sample_sizes = []
        sample_dets = []
        # set sizes of sample matrices; ensure exact calculations are only on small matrices
        s = cov.shape[0]
        max_size = 8192
        if s > max_size:
            sample_sizes = [s, max_size, max_size//2, max_size//4]
        else:
            sample_sizes = [s, s//2, s//4, s//8]
        for i in sample_sizes:
            # calculate logdet of the full matrix using circulant eigenvalue approximation
            if i == s:
                ld = 2*numpy.log(p/(2*p.delta_t)).sum()
                sample_dets.append(ld)
            # generate three more sample matrices using exact calculations
            else:
                gate_size = s - i
                start = (s - gate_size)//2
                end = start + gate_size
                tc = numpy.delete(numpy.delete(cov, slice(start, end), 0), slice(start, end), 1)
                ld = numpy.linalg.slogdet(tc)[1]
                sample_dets.append(ld)
        # generate a linear regression using the four points (size, logdet)
        x = numpy.vstack([sample_sizes, numpy.ones(len(sample_sizes))]).T
        m, b = numpy.linalg.lstsq(x, sample_dets, rcond=None)[0]
        return (sample_sizes, sample_dets), (m, b)
            
    def gate_indices(self, det):
        """Calculate the indices corresponding to start and end of gate.
        """
        # get time series start and delta_t
        ts = self._Rss[det]
        start_time_gc = float(self.td_data[det].start_time) # need float for conversion input
        delta_t = ts.delta_t
        # get gate start and length from get_gate_times
        gate_start, gate_length = self.get_gate_times()[det]
        # convert to indices
        lindex = float(gate_start - start_time_gc) // delta_t
        rindex = lindex + (gate_length // delta_t)
        lindex = lindex if lindex >= 0 else 0
        rindex = rindex if rindex <= len(ts) else len(ts)
        return lindex, rindex
    
    def det_lognorm(self, det, start_index=None, end_index=None):
        """Calculate the normalization term from the truncated covariance matrix.
        Determinant is estimated using a linear fit to logdet vs truncated matrix size
        """
        if not self.normalize:
            return 0
        try:
            # check if the key already exists; if so, return its value
            lognorm = self._lognorm[(det, start_index, end_index)]
        except KeyError:
            # if not, extrapolate the normalization term
            start_index, end_index = self.gate_indices(det)
            # get the size of the matrix
            n = len(self._Rss[det])
            trunc_size = n - (end_index - start_index)
            # call the linear regression
            m, b = self._cov_regressions[det]
            # extrapolate from linear fit
            ld = m*trunc_size + b
            lognorm = -0.5*(numpy.log(2*numpy.pi)*trunc_size + ld) # full normalization term
            # cache the result
            self._lognorm[(det, start_index, end_index)] = lognorm
        return lognorm

    @property
    def normalize(self):
        """Determines if the loglikelihood includes the normalization term.
        """
        return self._normalize

    ### This is called before psds, so doing direct calls to self._psds, etc. throws an error for now ###
    
    @normalize.setter
    def normalize(self, normalize):
        """Clears the current stats if the normalization state is changed.
        """
        self._normalize = normalize

    @staticmethod
    def _nowaveform_logl():
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
        return self._loglikelihood() - self._lognl()

    def whiten(self, data, whiten, inplace=False):
        """Whitens the given data.

        Parameters
        ----------
        data : dict
            Dictionary of detector names -> FrequencySeries.
        whiten : {0, 1, 2}
            Integer indicating what level of whitening to apply. Levels are:
            0: no whitening; 1: whiten; 2: overwhiten.
        inplace : bool, optional
            If True, modify the data in place. Otherwise, a copy will be
            created for whitening.


        Returns
        -------
        dict :
            Dictionary of FrequencySeries after the requested whitening has
            been applied.
        """
        if not inplace:
            data = {det: d.copy() for det, d in data.items()}
        if whiten:
            for det, dtilde in data.items():
                if whiten == 1:
                    dtilde *= self._invasds[det]
                elif whiten == 2:
                    dtilde *= self._invpsds[det]
                else:
                    raise ValueError("whiten must be either 0, 1, or 2")
        return data

    def get_waveforms(self):
        """The waveforms generated using the current parameters.

        If the waveforms haven't been generated yet, they will be generated,
        resized to the same length as the data, and cached. If the
        ``highpass_waveforms`` attribute is set, a highpass filter will
        also be applied to the waveforms.

        Returns
        -------
        dict :
            Dictionary of detector names -> FrequencySeries.
        """
        if self._current_wfs is None:
            params = self.current_params
            wfs = self.waveform_generator.generate(**params)
            for det, h in wfs.items():
                # make the same length as the data
                h.resize(len(self.data[det]))
                # apply high pass
                if self.highpass_waveforms:
                    h = highpass(
                        h.to_timeseries(),
                        frequency=self.highpass_waveforms).to_frequencyseries()
                wfs[det] = h
            self._current_wfs = wfs
        return self._current_wfs

    @abstractmethod
    def get_gated_waveforms(self, white=False):
        """Generates and gates waveforms using the current parameters.

        Parameters
        ----------
        white : bool, optional
            Gate and in-paint the waveforms using the inverse ASD instead of
            the inverse PSD. This will cause the whitened waveform to be
            identically zero in the gated region rather than the overwhitened
            waveform. Default is False.

        Returns
        -------
        dict :
            Dictionary of detector names -> FrequencySeries.
        """
        pass

    def get_residuals(self):
        """Generates the residuals ``d-h`` using the current parameters.

        Returns
        -------
        dict :
            Dictionary of detector names -> FrequencySeries.
        """
        wfs = self.get_waveforms()
        
        out = {}
        for det, h in wfs.items():
            d = self.data[det]
            out[det] = d - h
        return out

    def get_data(self):
        """Return a copy of the data.

        Returns
        -------
        dict :
            Dictionary of detector names -> FrequencySeries.
        """
        return {det: d.copy() for det, d in self.data.items()}

    def get_gated_data(self):
        """Return a copy of the gated data.

        The gated data will be cached for faster retrieval.

        Returns
        -------
        dict :
            Dictionary of detector names -> FrequencySeries.
        """
        gate_times = self.get_gate_times()
        out = {}
        for det, d in self.td_data.items():
            # make sure the cache at least has the detectors in it
            try:
                cache = self._gated_data[det]
            except KeyError:
                cache = self._gated_data[det] = {}
            invasd = self._invasds[det]
            gatestartdelay, dgatedelay = gate_times[det]
            try:
                dtilde = cache[gatestartdelay, dgatedelay]
            except KeyError:
                # doesn't exist yet, or the gate times changed
                cache.clear()
                d = d.gate(gatestartdelay + dgatedelay/2,
                           window=dgatedelay/2, copy=True,
                           kernel=invasd, method='paint',
                           zero_before_gate=self.zero_before_gate,
                           zero_after_gate=self.zero_after_gate)
                dtilde = d.to_frequencyseries()
                # save for next time
                cache[gatestartdelay, dgatedelay] = dtilde
            out[det] = dtilde
        return out

    def get_gate_times(self):
        """Gets the time to apply a gate based on the current sky position.

        If the parameter ``gatefunc`` is set to ``'hmeco'``, the gate times
        will be calculated based on the hybrid MECO of the given set of
        parameters; see ``get_gate_times_hmeco`` for details. Otherwise, the
        gate times will just be retrieved from the ``t_gate_start`` and
        ``t_gate_end`` parameters.

        .. warning::
            Since the normalization of the likelihood is currently not
            being calculated, it is recommended that you do not use
            ``gatefunc``, instead using fixed gate times.

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
        # we'll need the sky location for determining time shifts
        ra = self.current_params['ra']
        dec = self.current_params['dec']
        gatetimes = {}
        for det in self.data:
            thisdet = Detector(det)
            # account for the time delay between the waveforms of the
            # different detectors
            gatestartdelay = gatestart + thisdet.time_delay_from_earth_center(ra, dec, gatestart)
            gateenddelay = gateend + thisdet.time_delay_from_earth_center(ra, dec, gateend)
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
        # generate the template waveform
        try:
            wfs = self.get_waveforms()
        except NoWaveformError:
            return self._nowaveform_logl()
        except FailedWaveformError as e:
            if self.ignore_failed_waveforms:
                return self._nowaveform_logl()
            raise e
        # get waveform parameters
        params = self.current_params
        spin1 = params['spin1z']
        spin2 = params['spin2z']
        # gate input for ringdown analysis which consideres a start time
        # and an end time
        dgate = params['gate_window']
        meco_f = hybrid_meco_frequency(params['mass1'], params['mass2'], spin1, spin2)
        # figure out the gate times
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
            if t_from_freq[f_idx] > 0:
                gatestartdelay = t_from_freq[f_idx] + float(t_from_freq.epoch)
            else:
                gatestartdelay = t_from_freq[f_idx] + ht.sample_times[-1]
            gatestartdelay = min(gatestartdelay, params['t_gate_start'])
            gatetimes[det] = (gatestartdelay, dgate)
        return gatetimes

    def _lognl(self):
        """Calculates the log of the noise likelihood.
        """
        # clear variables
        lognl = 0.
        self._det_lognls.clear()
        # get the times of the gates
        gate_times = self.get_gate_times()
        for det, invasd in self._invasds.items():
            start_index, end_index = self.gate_indices(det)
            norm = self.det_lognorm(det, start_index, end_index) # linear estimation
            gatestartdelay, dgatedelay = gate_times[det]
            # we always filter the entire segment starting from kmin, since the
            # gated series may have high frequency components
            slc = slice(self._kmin[det], self._kmax[det])
            # gate the data
            data = self.td_data[det]
            gated_dt = data.gate(gatestartdelay + dgatedelay/2,
                                 window=dgatedelay/2, copy=True,
                                 kernel=invasd, method='paint',
                                 zero_before_gate=self.zero_before_gate,
                                 zero_after_gate=self.zero_after_gate)
            # convert to the frequency series
            gated_d = gated_dt.to_frequencyseries()
            # whiten
            gated_d *= invasd
            # inner product
            ip = 4 * invasd.delta_f * gated_d[slc].inner(gated_d[slc]).real
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

    @staticmethod
    def _fd_data_from_strain_dict(opts, strain_dict, psd_strain_dict):
        """Wrapper around :py:func:`data_utils.fd_data_from_strain_dict`.

        Ensures that if the PSD is estimated from data, the inverse spectrum
        truncation uses a Hann window, and that the low frequency cutoff is
        zero.
        """
        if opts.psd_inverse_length and opts.invpsd_trunc_method is None:
            # make sure invpsd truncation is set to hanning
            logging.info("Using Hann window to truncate inverse PSD")
            opts.invpsd_trunc_method = 'hann'
        lfs = None
        if opts.psd_estimation:
            # make sure low frequency cutoff is zero
            logging.info("Setting low frequency cutoff of PSD to 0")
            lfs = opts.low_frequency_cutoff.copy()
            opts.low_frequency_cutoff = {d: 0. for d in lfs}
        out = fd_data_from_strain_dict(opts, strain_dict, psd_strain_dict)
        # set back
        if lfs is not None:
            opts.low_frequency_cutoff = lfs
        return out

    def write_metadata(self, fp, group=None):
        """Adds writing the psds, and analyzed detectors.

        The analyzed detectors, their analysis segments, and the segments
        used for psd estimation are written as
        ``analyzed_detectors``, ``{{detector}}_analysis_segment``, and
        ``{{detector}}_psd_segment``, respectively. These are either written
        to the specified ``group``'s attrs, or to the top level attrs if
        ``group`` is None.

        Parameters
        ----------
        fp : pycbc.inference.io.BaseInferenceFile instance
            The inference file to write to.
        group : str, optional
            If provided, the metadata will be written to the attrs specified
            by group, i.e., to ``fp[group].attrs``. Otherwise, metadata is
            written to the top-level attrs (``fp.attrs``).
        """
        BaseDataModel.write_metadata(self, fp, group=group)
        attrs = fp.getattrs(group=group)
        # write the analyzed detectors and times
        attrs['analyzed_detectors'] = self.detectors # store fitting values here
        for det, data in self.data.items():
            key = '{}_analysis_segment'.format(det)
            attrs[key] = [float(data.start_time), float(data.end_time)]
            # store covariance determinant extrapolation information (checkpoint)
            if self.normalize:
                attrs['{}_cov_sample'.format(det)] = self._cov_samples[det]
                attrs['{}_cov_regression'.format(det)] = self._cov_regressions[det]
        if self._psds is not None and not self.no_save_data:
            fp.write_psd(self._psds, group=group)
        # write the times used for psd estimation (if they were provided)
        for det in self.psd_segments:
            key = '{}_psd_segment'.format(det)
            attrs[key] = list(map(float, self.psd_segments[det]))
        # save the frequency cutoffs
        for det in self.detectors:
            attrs['{}_likelihood_low_freq'.format(det)] = self._f_lower[det]
            if self._f_upper[det] is not None:
                attrs['{}_likelihood_high_freq'.format(det)] = \
                    self._f_upper[det]


class GatedGaussianNoise(BaseGatedGaussian):
    r"""Model that applies a time domain gate, assuming stationary Gaussian
    noise.

    The gate start and end times are set by providing ``t_gate_start`` and
    ``t_gate_end`` parameters, respectively. This will cause the gated times
    to be excised from the analysis. For more details on the likelihood
    function and its derivation, see
    `arXiv:2105.05238 <https://arxiv.org/abs/2105.05238>`_.

    .. warning::
        The normalization of the likelihood depends on the gate times. However,
        at the moment, the normalization is not calculated, as it depends on
        the determinant of the truncated covariance matrix (see Eq. 4 of
        arXiv:2105.05238). For this reason it is recommended that you only
        use this model for fixed gate times.

    """
    name = 'gated_gaussian_noise'

    def __init__(self, variable_params, data, low_frequency_cutoff, psds=None,
                 high_frequency_cutoff=None, normalize=False,
                 static_params=None, **kwargs):
        # set up the boiler-plate attributes
        super().__init__(
            variable_params, data, low_frequency_cutoff, psds=psds,
            high_frequency_cutoff=high_frequency_cutoff, normalize=normalize,
            static_params=static_params, **kwargs)
        # create the waveform generator
        self.waveform_generator = create_waveform_generator(
            self.variable_params, self.data,
            waveform_transforms=self.waveform_transforms,
            recalibration=self.recalibration,
            gates=self.gates, **self.static_params)

    @property
    def _extra_stats(self):
        """No extra stats are stored."""
        return []

    ### This needs to be called after running model.update() for changes to take effect ###
    
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
        # generate the template waveform
        try:
            wfs = self.get_waveforms()
        except NoWaveformError:
            return self._nowaveform_logl()
        except FailedWaveformError as e:
            if self.ignore_failed_waveforms:
                return self._nowaveform_logl()
            raise e
        # get the times of the gates
        gate_times = self.get_gate_times()
        logl = 0.
        for det, h in wfs.items():
            invasd = self._invasds[det]
            start_index, end_index = self.gate_indices(det)
            norm = self.det_lognorm(det, start_index, end_index)
            gatestartdelay, dgatedelay = gate_times[det]
            # we always filter the entire segment starting from kmin, since the
            # gated series may have high frequency components
            slc = slice(self._kmin[det], self._kmax[det])
            # calculate the residual
            data = self.data[det]
            res = (data - h).to_timeseries()
            res = res.gate(gatestartdelay + dgatedelay/2,
                           window=dgatedelay/2, copy=False,
                           invasd=invasd, method='paint',
                           zero_before_gate=self.zero_before_gate,
                           zero_after_gate=self.zero_after_gate)
            rtilde = res.to_frequencyseries()
            # whiten
            rtilde *= invasd
            rr = 4 * invasd.delta_f * rtilde[slc].inner(rtilde[slc]).real
            logl += norm - 0.5*rr
        return float(logl)
    
    @property
    def multi_signal_support(self):
        """ The list of classes that this model supports in a multi-signal
        likelihood
        """
        return [type(self)]
    
    def multi_loglikelihood(self, models):
        """ Calculate a multi-model (signal) likelihood
        """
        # Generate the waveforms for each submodel
        wfs = []
        for m in models + [self]:
            wf = m.get_waveforms()
            wfs.append(wf)
        # combine into a single waveform
        combine = {}
        for det in self.data:
            # get max waveform length
            mlen = max([len(x[det]) for x in wfs])
            [x[det].resize(mlen) for x in wfs]
            combine[det] = sum([x[det] for x in wfs])
        self._current_wfs = combine
        return self._loglikelihood()

    def get_gated_waveforms(self):
        wfs = self.get_waveforms()
        gate_times = self.get_gate_times()
        out = {}
        # apply the gate
        for det, h in wfs.items():
            invasd = self._invasds[det]
            gatestartdelay, dgatedelay = gate_times[det]
            ht = h.to_timeseries()
            ht = ht.gate(gatestartdelay + dgatedelay/2,
                         window=dgatedelay/2, copy=False,
                         kernel=invasd, method='paint',
                         zero_before_gate=self.zero_before_gate,
                         zero_after_gate=self.zero_after_gate)
            h = ht.to_frequencyseries()
            out[det] = h
        return out

    def get_gated_residuals(self):
        """Generates the gated residuals ``d-h`` using the current parameters.

        Returns
        -------
        dict :
            Dictionary of detector names -> FrequencySeries.
        """
        params = self.current_params
        wfs = self.waveform_generator.generate(**params)
        gate_times = self.get_gate_times()
        out = {}
        for det, h in wfs.items():
            invasd = self._invasds[det]
            gatestartdelay, dgatedelay = gate_times[det]
            data = self.data[det]
            res = (data - h).to_timeseries()
            res = res.gate(gatestartdelay + dgatedelay/2,
                           window=dgatedelay/2, copy=False,
                           invasd=invasd, method='paint',
                           zero_before_gate=self.zero_before_gate,
                           zero_after_gate=self.zero_after_gate)
            rtilde = res.to_frequencyseries()
            out[det] = rtilde
        return out


class GatedGaussianMargPol(BaseGatedGaussian):
    r"""Gated gaussian model with numerical marginalization over polarization.

    This implements the GatedGaussian likelihood with an explicit numerical
    marginalization over polarization angle. This is accomplished using
    a fixed set of integration points distribution uniformation between
    0 and 2pi. By default, 1000 integration points are used.
    The 'polarization_samples' argument can be passed to set an alternate
    number of integration points.
    """
    name = 'gated_gaussian_margpol'

    def __init__(self, variable_params, data, low_frequency_cutoff, psds=None,
                 high_frequency_cutoff=None, normalize=False,
                 static_params=None,
                 polarization_samples=1000, **kwargs):
        # set up the boiler-plate attributes
        super().__init__(
            variable_params, data, low_frequency_cutoff, psds=psds,
            high_frequency_cutoff=high_frequency_cutoff, normalize=normalize,
            static_params=static_params, **kwargs)
        # the polarization parameters
        self.polarization_samples = polarization_samples
        self.pol = numpy.linspace(0, 2*numpy.pi, self.polarization_samples)
        self.dets = {}
        # create the waveform generator
        self.waveform_generator = create_waveform_generator(
            self.variable_params, self.data,
            waveform_transforms=self.waveform_transforms,
            recalibration=self.recalibration,
            generator_class=generator.FDomainDetFrameTwoPolGenerator,
            **self.static_params)

    def get_waveforms(self):
        if self._current_wfs is not None:
            return self._current_wfs
        params = self.current_params
        wfs = self.waveform_generator.generate(**params)
        for det, (hp, hc) in wfs.items():
            # make the same length as the data
            hp.resize(len(self.data[det]))
            hc.resize(len(self.data[det]))
            # apply high pass
            if self.highpass_waveforms:
                hp = highpass(
                    hp.to_timeseries(),
                    frequency=self.highpass_waveforms).to_frequencyseries()
                hc = highpass(
                    hc.to_timeseries(),
                    frequency=self.highpass_waveforms).to_frequencyseries()
            wfs[det] = (hp, hc)
        self._current_wfs = wfs
        return self._current_wfs

    def get_gated_waveforms(self):
        wfs = self.get_waveforms()
        gate_times = self.get_gate_times()
        out = {}
        for det in wfs:
            invasd = self._invasds[det]
            gatestartdelay, dgatedelay = gate_times[det]
            # the waveforms are a dictionary of (hp, hc)
            pols = []
            for h in wfs[det]:
                ht = h.to_timeseries() 
                try:
                    ht = ht.gate(gatestartdelay + dgatedelay/2,
                             window=dgatedelay/2, copy=False,
                             kernel=invasd, method='paint',
                             zero_before_gate=self.zero_before_gate,
                             zero_after_gate=self.zero_after_gate)
                    h = ht.to_frequencyseries()
                except ValueError as e:
                    numpy.save('fail_params.out', self.current_params,
                               allow_pickle=True)
                    ht.save('fail_wf.hdf')
                    raise e
                pols.append(h)
            out[det] = tuple(pols)
        return out
    
    def get_gate_times_hmeco(self):
        """Gets the time to apply a gate based on the current sky position.
        Returns
        -------
        dict :
            Dictionary of detector names -> (gate start, gate width)
        """
        # generate the template waveform
        try:
            wfs = self.get_waveforms()
        except NoWaveformError:
            return self._nowaveform_logl()
        except FailedWaveformError as e:
            if self.ignore_failed_waveforms:
                return self._nowaveform_logl()
            raise e
        # get waveform parameters
        params = self.current_params
        spin1 = params['spin1z']
        spin2 = params['spin2z']
        # gate input for ringdown analysis which consideres a start time
        # and an end time
        dgate = params['gate_window']
        meco_f = hybrid_meco_frequency(params['mass1'], params['mass2'], spin1, spin2)
        # figure out the gate times
        gatetimes = {}
        # for now only calculating time from plus polarization; should be all that's necessary
        for det, (hp, hc) in wfs.items():
            invpsd = self._invpsds[det]
            hp.resize(len(invpsd))
            ht = hp.to_timeseries()
            f_low = int((self._f_lower[det]+1)/hp.delta_f)
            sample_freqs = hp.sample_frequencies[f_low:].numpy()
            f_idx = numpy.where(sample_freqs <= meco_f)[0][-1]
            # find time corresponding to meco frequency
            t_from_freq = time_from_frequencyseries(
                hp[f_low:], sample_frequencies=sample_freqs)
            if t_from_freq[f_idx] > 0:
                gatestartdelay = t_from_freq[f_idx] + float(t_from_freq.epoch)
            else:
                gatestartdelay = t_from_freq[f_idx] + ht.sample_times[-1]
            gatestartdelay = min(gatestartdelay, params['t_gate_start'])
            gatetimes[det] = (gatestartdelay, dgate)
        return gatetimes

    @property
    def _extra_stats(self):
        """Adds the maxL polarization and corresponding likelihood."""
        return ['maxl_polarization', 'maxl_logl']

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
        # get the gated waveforms
        try:
            wfs = self.get_gated_waveforms()
        except NoWaveformError:
            return self._nowaveform_logl()
        except FailedWaveformError as e:
            if self.ignore_failed_waveforms:
                return self._nowaveform_logl()
            raise e
        # get the gated data
        data = self.get_gated_data()
        # cycle over
        loglr = 0.
        lognl = 0.
        for det, (hp, hc) in wfs.items():
            # get the antenna patterns
            if det not in self.dets:
                self.dets[det] = Detector(det)
            fp, fc = self.dets[det].antenna_pattern(self.current_params['ra'],
                                                    self.current_params['dec'],
                                                    self.pol,
                                                    self.current_params['tc'])
            start_index, end_index = self.gate_indices(det)
            norm = self.det_lognorm(det, start_index, end_index)
            # we always filter the entire segment starting from kmin, since the
            # gated series may have high frequency components
            slc = slice(self._kmin[det], self._kmax[det])
            # whiten the data and waveforms
            invasd = self._invasds[det]
            hp = hp[slc] * invasd[slc]
            hc = hc[slc] * invasd[slc]
            d = data[det][slc] * invasd[slc]
            # get the various gated inner products
            hpd = hp.inner(d).real  # <hp, d>
            hcd = hc.inner(d).real  # <hc, d>
            dhp = d.inner(hp).real  # <d, hp>
            dhc = d.inner(hc).real  # <d, hc>
            hphp = hp.inner(hp).real  # <hp, hp>
            hchc = hc.inner(hc).real  # <hc, hc>
            hphc = hp.inner(hc).real  # <hp, hc>
            hchp = hc.inner(hp).real  # <hc, hp>
            dd = d.inner(d).real  # <d, d>
            # since the antenna patterns are real,
            # <h, d>/2 + <d, h>/2 = fp*(<hp, d>/2 + <d, hp>/2)
            #                     + fc*(<hc, d>/2 + <d, hc>/2)
            hd = fp*(hpd + dhp) + fc*(hcd + dhc)
            # <h, h>/2 = <fp*hp + fc*hc, fp*hp + fc*hc>/2
            #          = fp*fp*<hp, hp>/2 + fc*fc*<hc, hc>/2
            #            + fp*fc*<hp, hc>/2 + fc*fp*<hc, hp>/2
            hh = fp*fp*hphp + fc*fc*hchc + fp*fc*(hphc + hchp)
            # sum up; note that the factor is 2df instead of 4df to account
            # for the factor of 1/2
            loglr += norm + 2*invasd.delta_f*(hd - hh)
            lognl += -2 * invasd.delta_f * dd
        # store the maxl polarization
        idx = loglr.argmax()
        setattr(self._current_stats, 'maxl_polarization', self.pol[idx])
        setattr(self._current_stats, 'maxl_logl', loglr[idx] + lognl)
        # compute the marginalized log likelihood
        marglogl = special.logsumexp(loglr) + lognl - numpy.log(len(self.pol))
        return float(marglogl)
    
    @property
    def multi_signal_support(self):
        """ The list of classes that this model supports in a multi-signal
        likelihood
        """
        return [type(self)]
    
    def multi_loglikelihood(self, models):
        """ Calculate a multi-model (signal) likelihood
        """
        # Generate the waveforms for each submodel
        wfs = []
        for m in models + [self]:
            wf = m.get_waveforms()
            wfs.append(wf)

        # combine into a single waveform
        combine = {}
        for det in self.data:
            # get max waveform length
            mlenp = max([len(x[det][0]) for x in wfs])
            mlenc = max([len(x[det][1]) for x in wfs])
            mlen = max([mlenp, mlenc])
            # resize plus and cross
            [x[det][0].resize(mlen) for x in wfs]
            [x[det][1].resize(mlen) for x in wfs]
            # combine waveforms
            combine[det] = (sum([x[det][0] for x in wfs]), sum([x[det][1] for x in wfs]))

        self._current_wfs = combine
        return self._loglikelihood()
