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

"""This module provides model classes that assume the noise is Gaussian.
"""

import logging
from abc import ABCMeta
from six import add_metaclass
import numpy

from pycbc import filter as pyfilter
from pycbc.waveform import NoWaveformError
from pycbc.waveform import generator
from pycbc.types import Array, FrequencySeries
from pycbc.strain import gates_from_cli
from pycbc.strain.calibration import Recalibrate
from pycbc.inject import InjectionSet
from pycbc.types.optparse import MultiDetOptionAction

from .base import ModelStats
from .base_data import BaseDataModel
from .data_utils import (data_opts_from_config, data_from_cli,
                         fd_data_from_strain_dict, gate_overwhitened_data)


@add_metaclass(ABCMeta)
class BaseGaussianNoise(BaseDataModel):
    r"""Model for analyzing GW data with assuming a wide-sense stationary
    Gaussian noise model.

    This model will load gravitational wave data and calculate the log noise
    likelihood ``_lognl`` and normalization. It also implements the
    ``_loglikelihood`` function as the sum of the log likelihood ratio and the
    ``lognl``. It does not implement a log likelihood ratio function
    ``_loglr``, however, since that can differ depending on the signal model.
    Models that analyze GW data assuming it is stationary Gaussian should
    therefore inherit from this class and implement their own ``_loglr``
    function.

    For more details on the inner product used, the log likelihood of the
    noise, and the normalization factor, see :py:class:`GaussianNoise`.

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
    psds : dict, optional
        A dictionary of FrequencySeries keyed by the detector names. The
        dictionary must have a psd for each detector specified in the data
        dictionary. If provided, the inner products in each detector will be
        weighted by 1/psd of that detector.
    high_frequency_cutoff : dict, optional
        A dictionary of ending frequencies, in which the keys are the
        detector names and the values are the ending frequencies for the
        respective detectors to be used for computing inner products. If not
        provided, the minimum of the largest frequency stored in the data
        and a given waveform will be used.
    normalize : bool, optional
        If True, the normalization factor :math:`alpha` will be included in the
        log likelihood. See :py:class:`GaussianNoise` for details. Default is
        to not include it.
    static_params : dict, optional
        A dictionary of parameter names -> values to keep fixed.
    \**kwargs :
        All other keyword arguments are passed to ``BaseDataModel``.

    Attributes
    ----------
    data : dict
        Dictionary of detectors -> frequency-domain data.
    low_frequency_cutoff
    high_frequency_cutoff
    kmin
    kmax
    psds
    psd_segments
    weight
    whitened_data
    normalize
    lognorm
    """
    def __init__(self, variable_params, data, low_frequency_cutoff, psds=None,
                 high_frequency_cutoff=None, normalize=False,
                 static_params=None, **kwargs):
        # set up the boiler-plate attributes
        super(BaseGaussianNoise, self).__init__(variable_params, data,
                                                static_params=static_params,
                                                **kwargs)
        # check if low frequency cutoff has been provided for every IFO with
        # data
        for ifo in self.data.keys():
            if low_frequency_cutoff[ifo] is None:
                raise ValueError(
                    "A low-frequency-cutoff must be provided for every "
                    "detector for which data has been provided. If "
                    "loading the model settings from "
                    "a config file, please provide "
                    "`{DETECTOR}:low-frequency-cutoff` options for "
                    "every detector in the `[model]` section, where "
                    "`{DETECTOR} is the name of the detector,"
                    "or provide a single low-frequency-cutoff option"
                    "which will be used for all detectors")

        # check that the data sets all have the same delta fs and delta ts
        dts = numpy.array([d.delta_t for d in self.data.values()])
        dfs = numpy.array([d.delta_f for d in self.data.values()])
        if not all(dts == dts[0]):
            raise ValueError("all data must have the same sample rate")
        if not all(dfs == dfs[0]):
            raise ValueError("all data must have the same segment length")
        # store the number of samples in the time domain
        self._N = int(1./(dts[0]*dfs[0]))
        # Set low frequency cutoff
        self.low_frequency_cutoff = self._f_lower = low_frequency_cutoff
        # set upper frequency cutoff
        self._f_upper = None
        self.high_frequency_cutoff = high_frequency_cutoff
        # Set the cutoff indices
        self._kmin = {}
        self._kmax = {}

        for (det, d) in self._data.items():
            kmin, kmax = pyfilter.get_cutoff_indices(self._f_lower[det],
                                                     self._f_upper[det],
                                                     d.delta_f, self._N)
            self._kmin[det] = kmin
            self._kmax[det] = kmax
        # store the psd segments
        self._psd_segments = {}
        if psds is not None:
            self.set_psd_segments(psds)
        # store the psds and calculate the inner product weight
        self._psds = {}
        self._weight = {}
        self._lognorm = {}
        self._det_lognls = {}
        self._whitened_data = {}
        # set the normalization state
        self._normalize = False
        self.normalize = normalize
        # store the psds and whiten the data
        self.psds = psds

    @property
    def high_frequency_cutoff(self):
        """The high frequency cutoff of the inner product.

        If a high frequency cutoff was not provided for a detector, it will
        be ``None``.
        """
        return self._f_upper

    @high_frequency_cutoff.setter
    def high_frequency_cutoff(self, high_frequency_cutoff):
        """Sets the high frequency cutoff.

        Parameters
        ----------
        high_frequency_cutoff : dict
            Dictionary mapping detector names to frequencies. If a high
            frequency cutoff is not provided for one or more detectors, the
            Nyquist frequency will be used for those detectors.
        """
        self._f_upper = {}
        if high_frequency_cutoff is not None and bool(high_frequency_cutoff):
            for det in self._data:
                if det in high_frequency_cutoff:
                    self._f_upper[det] = high_frequency_cutoff[det]
                else:
                    self._f_upper[det] = None
        else:
            for det in self._data.keys():
                self._f_upper[det] = None

    @property
    def kmin(self):
        """Dictionary of starting indices for the inner product.

        This is determined from the lower frequency cutoff and the ``delta_f``
        of the data using
        :py:func:`pycbc.filter.matchedfilter.get_cutoff_indices`.
        """
        return self._kmin

    @property
    def kmax(self):
        """Dictionary of ending indices for the inner product.

        This is determined from the high frequency cutoff and the ``delta_f``
        of the data using
        :py:func:`pycbc.filter.matchedfilter.get_cutoff_indices`. If no high
        frequency cutoff was provided, this will be the indice corresponding to
        the Nyquist frequency.
        """
        return self._kmax

    @property
    def psds(self):
        """Dictionary of detectors -> PSD frequency series.

        If no PSD was provided for a detector, this will just be a frequency
        series of ones.
        """
        return self._psds

    @psds.setter
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
        self._weight.clear()
        self._lognorm.clear()
        self._det_lognls.clear()
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
            w = Array(numpy.zeros(len(p)))
            # only set weight in band we will analyze
            kmin = self._kmin[det]
            kmax = self._kmax[det]
            w[kmin:kmax] = numpy.sqrt(4.*p.delta_f/p[kmin:kmax])
            self._weight[det] = w
            self._whitened_data[det] = d.copy()
            self._whitened_data[det][kmin:kmax] *= w[kmin:kmax]
        # set the lognl and lognorm; we'll get this by just calling lognl
        _ = self.lognl

    @property
    def psd_segments(self):
        """Dictionary giving times used for PSD estimation for each detector.

        If a detector's PSD was not estimated from data, or the segment wasn't
        provided, that detector will not be in the dictionary.
        """
        return self._psd_segments

    def set_psd_segments(self, psds):
        """Sets the PSD segments from a dictionary of PSDs.

        This attempts to get the PSD segment from a ``psd_segment`` attribute
        of each detector's PSD frequency series. If that attribute isn't set,
        then that detector is not added to the dictionary of PSD segments.

        Parameters
        ----------
        psds : dict
            Dictionary of detector name -> PSD frequency series. The segment
            used for each PSD will try to be retrieved from the PSD's
            ``.psd_segment`` attribute.
        """
        for det, p in psds.items():
            try:
                self._psd_segments[det] = p.psd_segment
            except AttributeError:
                continue

    @property
    def weight(self):
        r"""Dictionary of detectors -> frequency series of inner-product
        weights.

        The weights are :math:`\sqrt{4 \Delta f / S_n(f)}`. This is set when
        the PSDs are set.
        """
        return self._weight

    @property
    def whitened_data(self):
        r"""Dictionary of detectors -> whitened data frequency series.

        The whitened data is the data multiplied by the inner-product weight.
        Note that this includes the :math:`\sqrt{4 \Delta f}` factor. This
        is set when the PSDs are set.
        """
        return self._whitened_data

    def det_lognorm(self, det):
        """The log of the likelihood normalization in the given detector.

        If ``self.normalize`` is False, will just return 0.
        """
        if not self.normalize:
            return 0.
        try:
            return self._lognorm[det]
        except KeyError:
            # hasn't been calculated yet
            p = self._psds[det]
            dt = self._whitened_data[det].delta_t
            kmin = self._kmin[det]
            kmax = self._kmax[det]
            lognorm = -float(self._N*numpy.log(numpy.pi*self._N*dt)/2.
                             + numpy.log(p[kmin:kmax]).sum())
            self._lognorm[det] = lognorm
            return self._lognorm[det]

    @property
    def normalize(self):
        """Determines if the loglikelihood includes the normalization term.
        """
        return self._normalize

    @normalize.setter
    def normalize(self, normalize):
        """Clears the current stats if the normalization state is changed.
        """
        if normalize != self._normalize:
            self._current_stats = ModelStats()
            self._lognorm.clear()
            self._det_lognls.clear()
        self._normalize = normalize

    @property
    def lognorm(self):
        """The log of the normalization of the log likelihood."""
        return sum(self.det_lognorm(det) for det in self._data)

    def det_lognl(self, det):
        r"""Returns the log likelihood of the noise in the given detector:

        .. math::

            \log p(d_i|n_i) = \log \alpha_i -
                \frac{1}{2} \left<d_i | d_i\right>.


        Parameters
        ----------
        det : str
            The name of the detector.

        Returns
        -------
        float :
            The log likelihood of the noise in the requested detector.
        """
        try:
            return self._det_lognls[det]
        except KeyError:
            # hasn't been calculated yet; calculate & store
            kmin = self._kmin[det]
            kmax = self._kmax[det]
            d = self._whitened_data[det]
            lognorm = self.det_lognorm(det)
            lognl = lognorm - 0.5 * d[kmin:kmax].inner(d[kmin:kmax]).real
            self._det_lognls[det] = lognl
            return self._det_lognls[det]

    def _lognl(self):
        """Computes the log likelihood assuming the data is noise.

        Since this is a constant for Gaussian noise, this is only computed once
        then stored.
        """
        return sum(self.det_lognl(det) for det in self._data)

    def _loglikelihood(self):
        r"""Computes the log likelihood of the paramaters,

        .. math::

            \log p(d|\Theta, h) = \log \alpha -\frac{1}{2}\sum_i
                \left<d_i - h_i(\Theta) | d_i - h_i(\Theta)\right>,

        at the current parameter values :math:`\Theta`.

        Returns
        -------
        float
            The value of the log likelihood evaluated at the given point.
        """
        # since the loglr has fewer terms, we'll call that, then just add
        # back the noise term that canceled in the log likelihood ratio
        return self.loglr + self.lognl

    def write_metadata(self, fp):
        """Adds writing the psds and lognl, since it's a constant.

        The lognl is written to the sample group's ``attrs``.

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
        attrs['lognl'] = self.lognl
        for det in self.detectors:
            # Save lognl for each IFO as attributes in the samples group
            attrs['{}_lognl'.format(det)] = self.det_lognl(det)
            # Save each IFO's low frequency cutoff used in the likelihood
            # computation as an attribute
            fp.attrs['{}_likelihood_low_freq'.format(det)] = self._f_lower[det]
            # Save the IFO's high frequency cutoff used in the likelihood
            # computation as an attribute if one was provided the user
            if self._f_upper[det] is not None:
                fp.attrs['{}_likelihood_high_freq'.format(det)] = \
                                                        self._f_upper[det]

    @classmethod
    def from_config(cls, cp, data_section='data', **kwargs):
        r"""Initializes an instance of this class from the given config file.

        In addition to ``[model]``, a ``data_section`` (default ``[data]``)
        must be in the configuration file. The data section specifies settings
        for loading data and estimating PSDs. See the `online documentation
        <http://pycbc.org/pycbc/latest/html/inference.html#setting-data>`_ for
        more details.

        The following options are read from the ``[model]`` section, in
        addition to ``name`` (which must be set):

        * ``{{DET}}-low-frequency-cutoff = FLOAT`` :
          The low frequency cutoff to use for each detector {{DET}}. A cutoff
          must be provided for every detector that may be analyzed (any
          additional detectors are ignored).
        * ``{{DET}}-high-frequency-cutoff = FLOAT`` :
          (Optional) A high frequency cutoff for each detector. If not
          provided, the Nyquist frequency is used.
        * ``check-for-valid-times =`` :
          (Optional) If provided, will check that there are no data quality
          flags on during the analysis segment and the segment used for PSD
          estimation in each detector. To check for flags,
          :py:func:`pycbc.dq.query_flag` is used, with settings pulled from the
          ``dq-*`` options in the ``[data]`` section. If a detector has bad
          data quality during either the analysis segment or PSD segment, it
          will be removed from the analysis.
        * ``shift-psd-times-to-valid =`` :
          (Optional) If provided, the segment used for PSD estimation will
          automatically be shifted left or right until a continous block of
          data with no data quality issues can be found. If no block can be
          found with a maximum shift of +/- the requested psd segment length,
          the detector will not be analyzed.
        * ``err-on-missing-detectors =`` :
          Raises an error if any detector is removed from the analysis because
          a valid time could not be found. Otherwise, a warning is printed
          to screen and the detector is removed from the analysis.
        * ``normalize =``:
          (Optional) Turn on the normalization factor.

        Parameters
        ----------
        cp : WorkflowConfigParser
            Config file parser to read.
        data_section : str, optional
            The name of the section to load data options from.
        \**kwargs :
            All additional keyword arguments are passed to the class. Any
            provided keyword will over ride what is in the config file.
        """
        args = cls._init_args_from_config(cp)
        # check if normalize is set
        if cp.has_option('model', 'normalize'):
            args['normalize'] = True
        # get any other keyword arguments provided in the model section
        ignore_args = ['name', 'normalize']
        for option in cp.options("model"):
            if option in ("low-frequency-cutoff", "high-frequency-cutoff"):
                ignore_args.append(option)
                name = option.replace('-', '_')
                args[name] = cp.get_cli_option('model', name,
                                               nargs='+', type=float,
                                               action=MultiDetOptionAction)
        if 'low_frequency_cutoff' not in args:
            raise ValueError("low-frequency-cutoff must be provided in the"
                             " model section, but is not found!")

        # data args
        bool_args = ['check-for-valid-times', 'shift-psd-times-to-valid',
                     'err-on-missing-detectors']
        data_args = {arg.replace('-', '_'): True for arg in bool_args
                     if cp.has_option('model', arg)}
        ignore_args += bool_args
        # load the data
        opts = data_opts_from_config(cp, data_section,
                                     args['low_frequency_cutoff'])
        strain_dict, psd_strain_dict = data_from_cli(opts, **data_args)
        # convert to frequency domain and get psds
        stilde_dict, psds = fd_data_from_strain_dict(opts, strain_dict,
                                                     psd_strain_dict)
        # save the psd data segments if the psd was estimated from data
        if opts.psd_estimation is not None:
            _tdict = psd_strain_dict or strain_dict
            for det in psds:
                psds[det].psd_segment = (_tdict[det].start_time,
                                         _tdict[det].end_time)
        # gate overwhitened if desired
        if opts.gate_overwhitened and opts.gate is not None:
            stilde_dict = gate_overwhitened_data(stilde_dict, psds, opts.gate)
        args.update({'data': stilde_dict, 'psds': psds})
        # any extra args
        args.update(cls.extra_args_from_config(cp, "model",
                                               skip_args=ignore_args))
        # get the injection file
        # Note: PyCBC's multi-ifo parser uses key:ifo for
        # the injection file, even though we will use the same
        # injection file for all detectors. This
        # should be fixed in a future version of PyCBC. Once it is,
        # update this. Until then, just use the first file.
        if opts.injection_file:
            injection_file = tuple(opts.injection_file.values())[0]
            # None if not set
        else:
            injection_file = None
        args['injection_file'] = injection_file
        # update any static params that are set to FROM_INJECTION
        replace_params = get_static_params_from_injection(
            args['static_params'], injection_file)
        args['static_params'].update(replace_params)
        # get ifo-specific instances of calibration model
        if cp.has_section('calibration'):
            logging.info("Initializing calibration model")
            recalib = {
                ifo: Recalibrate.from_config(cp, ifo, section='calibration')
                for ifo in opts.instruments}
            args['recalibration'] = recalib
        # get gates for templates
        gates = gates_from_cli(opts)
        if gates:
            args['gates'] = gates
        return cls(**args)


class GaussianNoise(BaseGaussianNoise):
    r"""Model that assumes data is stationary Gaussian noise.

    With Gaussian noise the log likelihood functions for signal
    :math:`\log p(d|\Theta, h)` and for noise :math:`\log p(d|n)` are given by:

    .. math::

        \log p(d|\Theta, h) &=  \log\alpha -\frac{1}{2} \sum_i
            \left< d_i - h_i(\Theta) | d_i - h_i(\Theta) \right> \\
        \log p(d|n) &= \log\alpha -\frac{1}{2} \sum_i \left<d_i | d_i\right>

    where the sum is over the number of detectors, :math:`d_i` is the data in
    each detector, and :math:`h_i(\Theta)` is the model signal in each
    detector. The (discrete) inner product is given by:

    .. math::

        \left<a_i | b_i\right> = 4\Re \Delta f
            \sum_{k=k_{\mathrm{min}}}^{k_{\mathrm{max}}}
            \frac{\tilde{a}_i^{*}[k] \tilde{b}_i[k]}{S^{(i)}_n[k]},

    where :math:`\Delta f` is the frequency resolution (given by 1 / the
    observation time :math:`T`), :math:`k` is an index over the discretely
    sampled frequencies :math:`f = k \Delta_f`, and :math:`S^{(i)}_n[k]` is the
    PSD in the given detector. The upper cutoff on the inner product
    :math:`k_{\max}` is by default the Nyquist frequency
    :math:`k_{\max} = N/2+1`, where :math:`N = \lfloor T/\Delta t \rfloor`
    is the number of samples in the time domain, but this can be set manually
    to a smaller value.

    The normalization factor :math:`\alpha` is:

    .. math::

        \alpha = \prod_{i} \frac{1}{\left(\pi T\right)^{N/2}
            \prod_{k=k_\mathrm{min}}^{k_{\mathrm{max}}} S^{(i)}_n[k]},

    where the product is over the number of detectors. By default, the
    normalization constant is not included in the log likelihood, but it can
    be turned on using the ``normalize`` keyword argument.

    Note that the log likelihood ratio has fewer terms than the log likelihood,
    since the normalization and :math:`\left<d_i|d_i\right>` terms cancel:

    .. math::

        \log \mathcal{L}(\Theta) = \sum_i \left[
            \left<h_i(\Theta)|d_i\right> -
            \frac{1}{2} \left<h_i(\Theta)|h_i(\Theta)\right> \right]

    Upon initialization, the data is whitened using the given PSDs. If no PSDs
    are given the data and waveforms returned by the waveform generator are
    assumed to be whitened.

    For more details on initialization parameters and definition of terms, see
    :py:class:`models.BaseDataModel`.

    Parameters
    ----------
    variable_params : (tuple of) string(s)
        A tuple of parameter names that will be varied.
    data : dict
        A dictionary of data, in which the keys are the detector names and the
        values are the data (assumed to be unwhitened). The list of keys must
        match the waveform generator's detectors keys, and the epoch of every
        data set must be the same as the waveform generator's epoch.
    low_frequency_cutoff : dict
        A dictionary of starting frequencies, in which the keys are the
        detector names and the values are the starting frequencies for the
        respective detectors to be used for computing inner products.
    psds : dict, optional
        A dictionary of FrequencySeries keyed by the detector names. The
        dictionary must have a psd for each detector specified in the data
        dictionary. If provided, the inner products in each detector will be
        weighted by 1/psd of that detector.
    high_frequency_cutoff : dict, optional
        A dictionary of ending frequencies, in which the keys are the
        detector names and the values are the ending frequencies for the
        respective detectors to be used for computing inner products. If not
        provided, the minimum of the largest frequency stored in the data
        and a given waveform will be used.
    normalize : bool, optional
        If True, the normalization factor :math:`alpha` will be included in the
        log likelihood. Default is to not include it.
    static_params : dict, optional
        A dictionary of parameter names -> values to keep fixed.
    \**kwargs :
        All other keyword arguments are passed to ``BaseDataModel``.

    Examples
    --------
    Create a signal, and set up the model using that signal:

    >>> from pycbc import psd as pypsd
    >>> from pycbc.inference.models import GaussianNoise
    >>> from pycbc.waveform.generator import (FDomainDetFrameGenerator,
    ...                                       FDomainCBCGenerator)
    >>> seglen = 4
    >>> sample_rate = 2048
    >>> N = seglen*sample_rate/2+1
    >>> fmin = 30.
    >>> static_params = {'approximant': 'IMRPhenomD', 'f_lower': fmin,
    ...                  'mass1': 38.6, 'mass2': 29.3,
    ...                  'spin1z': 0., 'spin2z': 0., 'ra': 1.37, 'dec': -1.26,
    ...                  'polarization': 2.76, 'distance': 3*500.}
    >>> variable_params = ['tc']
    >>> tsig = 3.1
    >>> generator = FDomainDetFrameGenerator(
    ...     FDomainCBCGenerator, 0., detectors=['H1', 'L1'],
    ...     variable_args=variable_params,
    ...     delta_f=1./seglen, **static_params)
    >>> signal = generator.generate(tc=tsig)
    >>> psd = pypsd.aLIGOZeroDetHighPower(N, 1./seglen, 20.)
    >>> psds = {'H1': psd, 'L1': psd}
    >>> low_frequency_cutoff = {'H1': fmin, 'L1': fmin}
    >>> model = GaussianNoise(variable_params, signal, low_frequency_cutoff,
                              psds=psds, static_params=static_params)

    Set the current position to the coalescence time of the signal:

    >>> model.update(tc=tsig)

    Now compute the log likelihood ratio and prior-weighted likelihood ratio;
    since we have not provided a prior, these should be equal to each other:

    >>> print('{:.2f}'.format(model.loglr))
    282.43
    >>> print('{:.2f}'.format(model.logplr))
    282.43

    Print all of the default_stats:

    >>> print(',\n'.join(['{}: {:.2f}'.format(s, v)
    ...                   for (s, v) in sorted(model.current_stats.items())]))
    H1_cplx_loglr: 177.76+0.00j,
    H1_optimal_snrsq: 355.52,
    L1_cplx_loglr: 104.67+0.00j,
    L1_optimal_snrsq: 209.35,
    logjacobian: 0.00,
    loglikelihood: 0.00,
    loglr: 282.43,
    logprior: 0.00

    Compute the SNR; for this system and PSD, this should be approximately 24:

    >>> from pycbc.conversions import snr_from_loglr
    >>> x = snr_from_loglr(model.loglr)
    >>> print('{:.2f}'.format(x))
    23.77

    Since there is no noise, the SNR should be the same as the quadrature sum
    of the optimal SNRs in each detector:

    >>> x = (model.det_optimal_snrsq('H1') +
    ...      model.det_optimal_snrsq('L1'))**0.5
    >>> print('{:.2f}'.format(x))
    23.77

    Toggle on the normalization constant:

    >>> model.normalize = True
    >>> model.loglikelihood
    835397.8757405131

    Using the same model, evaluate the log likelihood ratio at several points
    in time and check that the max is at tsig:

    >>> import numpy
    >>> times = numpy.linspace(tsig-1, tsig+1, num=101)
    >>> loglrs = numpy.zeros(len(times))
    >>> for (ii, t) in enumerate(times):
    ...     model.update(tc=t)
    ...     loglrs[ii] = model.loglr
    >>> print('tsig: {:.2f}, time of max loglr: {:.2f}'.format(
    ...     tsig, times[loglrs.argmax()]))
    tsig: 3.10, time of max loglr: 3.10

    Create a prior and use it (see distributions module for more details):

    >>> from pycbc import distributions
    >>> uniform_prior = distributions.Uniform(tc=(tsig-0.2,tsig+0.2))
    >>> prior = distributions.JointDistribution(variable_params, uniform_prior)
    >>> model = GaussianNoise(variable_params,
    ...     signal, low_frequency_cutoff, psds=psds, prior=prior,
    ...     static_params=static_params)
    >>> model.update(tc=tsig)
    >>> print('{:.2f}'.format(model.logplr))
    283.35
    >>> print(',\n'.join(['{}: {:.2f}'.format(s, v)
    ...                   for (s, v) in sorted(model.current_stats.items())]))
    H1_cplx_loglr: 177.76+0.00j,
    H1_optimal_snrsq: 355.52,
    L1_cplx_loglr: 104.67+0.00j,
    L1_optimal_snrsq: 209.35,
    logjacobian: 0.00,
    loglikelihood: 0.00,
    loglr: 282.43,
    logprior: 0.92

    """
    name = 'gaussian_noise'

    def __init__(self, variable_params, data, low_frequency_cutoff, psds=None,
                 high_frequency_cutoff=None, normalize=False,
                 static_params=None, **kwargs):
        # set up the boiler-plate attributes
        super(GaussianNoise, self).__init__(
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
        """Adds ``loglr``, plus ``cplx_loglr`` and ``optimal_snrsq`` in each
        detector."""
        return ['loglr'] + \
               ['{}_cplx_loglr'.format(det) for det in self._data] + \
               ['{}_optimal_snrsq'.format(det) for det in self._data]

    def _nowaveform_loglr(self):
        """Convenience function to set loglr values if no waveform generated.
        """
        for det in self._data:
            setattr(self._current_stats, 'loglikelihood', -numpy.inf)
            setattr(self._current_stats, '{}_cplx_loglr'.format(det),
                    -numpy.inf)
            # snr can't be < 0 by definition, so return 0
            setattr(self._current_stats, '{}_optimal_snrsq'.format(det), 0.)
        return -numpy.inf

    def _loglr(self):
        r"""Computes the log likelihood ratio,

        .. math::

            \log \mathcal{L}(\Theta) = \sum_i
                \left<h_i(\Theta)|d_i\right> -
                \frac{1}{2}\left<h_i(\Theta)|h_i(\Theta)\right>,

        at the current parameter values :math:`\Theta`.

        Returns
        -------
        float
            The value of the log likelihood ratio.
        """
        params = self.current_params
        try:
            wfs = self.waveform_generator.generate(**params)
        except NoWaveformError:
            return self._nowaveform_loglr()
        lr = 0.
        for det, h in wfs.items():
            # the kmax of the waveforms may be different than internal kmax
            kmax = min(len(h), self._kmax[det])
            if self._kmin[det] >= kmax:
                # if the waveform terminates before the filtering low frequency
                # cutoff, then the loglr is just 0 for this detector
                cplx_hd = 0j
                hh = 0.
            else:
                slc = slice(self._kmin[det], kmax)
                # whiten the waveform
                h[self._kmin[det]:kmax] *= self._weight[det][slc]
                # the inner products
                cplx_hd = self._whitened_data[det][slc].inner(h[slc])  # <h, d>
                hh = h[slc].inner(h[slc]).real  # < h, h>
            cplx_loglr = cplx_hd - 0.5*hh
            # store
            setattr(self._current_stats, '{}_optimal_snrsq'.format(det), hh)
            setattr(self._current_stats, '{}_cplx_loglr'.format(det),
                    cplx_loglr)
            lr += cplx_loglr.real
        # also store the loglikelihood, to ensure it is populated in the
        # current stats even if loglikelihood is never called
        self._current_stats.loglikelihood = lr + self.lognl
        return float(lr)

    def det_cplx_loglr(self, det):
        """Returns the complex log likelihood ratio in the given detector.

        Parameters
        ----------
        det : str
            The name of the detector.

        Returns
        -------
        complex float :
            The complex log likelihood ratio.
        """
        # try to get it from current stats
        try:
            return getattr(self._current_stats, '{}_cplx_loglr'.format(det))
        except AttributeError:
            # hasn't been calculated yet; call loglr to do so
            self._loglr()
            # now try returning again
            return getattr(self._current_stats, '{}_cplx_loglr'.format(det))

    def det_optimal_snrsq(self, det):
        """Returns the opitmal SNR squared in the given detector.

        Parameters
        ----------
        det : str
            The name of the detector.

        Returns
        -------
        float :
            The opimtal SNR squared.
        """
        # try to get it from current stats
        try:
            return getattr(self._current_stats, '{}_optimal_snrsq'.format(det))
        except AttributeError:
            # hasn't been calculated yet; call loglr to do so
            self._loglr()
            # now try returning again
            return getattr(self._current_stats, '{}_optimal_snrsq'.format(det))


#
# =============================================================================
#
#                               Support functions
#
# =============================================================================
#
def get_static_params_from_injection(static_params, injection_file):
    """Gets FROM_INJECTION static params from injection.

    Parameters
    ----------
    static_params : dict
        Dictionary of static params.
    injection_file : str
        Name of the injection file to use.

    Returns
    -------
    dict :
        A dictionary mapping parameter names to values retrieved from the
        injection file. The dictionary will only contain parameters that were
        set to ``"FROM_INJECTION"`` in the ``static_params``.
        Parameter names -> values The injection parameters.
    """
    # pull out the parameters that need replacing
    replace_params = {p: None for (p, val) in static_params.items()
                      if val == 'FROM_INJECTION'}
    if replace_params != {}:
        # make sure there's actually an injection file
        if injection_file is None:
            raise ValueError("one or more static params are set to "
                             "FROM_INJECTION, but no injection file "
                             "provided in data section")
        inj = InjectionSet(injection_file)
        # make sure there's only one injection provided
        if inj.table.size > 1:
            raise ValueError("Some static params set to FROM_INJECTION, but "
                             "more than one injection exists in the injection "
                             "file.")
        for param in replace_params:
            try:
                injval = inj.table[param][0]
            except NameError:
                # means the parameter doesn't exist
                raise ValueError("Static param {} with placeholder "
                                 "FROM_INJECTION has no counterpart in "
                                 "injection file.".format(param))
            replace_params[param] = injval
    return replace_params


def create_waveform_generator(variable_params, data, waveform_transforms=None,
                              recalibration=None, gates=None,
                              **static_params):
    """Creates a waveform generator for use with a model.

    Parameters
    ----------
    variable_params : list of str
        The names of the parameters varied.
    data : dict
        Dictionary mapping detector names to either a
        :py:class:`<pycbc.types.TimeSeries TimeSeries>` or
        :py:class:`<pycbc.types.FrequencySeries FrequencySeries>`.
    waveform_transforms : list, optional
        The list of transforms applied to convert variable parameters into
        parameters that will be understood by the waveform generator.
    recalibration : dict, optional
        Dictionary mapping detector names to
        :py:class:`<pycbc.calibration.Recalibrate>` instances for
        recalibrating data.
    gates : dict of tuples, optional
        Dictionary of detectors -> tuples of specifying gate times. The
        sort of thing returned by :py:func:`pycbc.gate.gates_from_cli`.

    Returns
    -------
    pycbc.waveform.FDomainDetFrameGenerator
        A waveform generator for frequency domain generation.
    """
    # the waveform generator will get the variable_params + the output
    # of the waveform transforms, so we'll add them to the list of
    # parameters
    if waveform_transforms is not None:
        wfoutputs = set.union(*[t.outputs
                                for t in waveform_transforms])
    else:
        wfoutputs = set()
    variable_params = list(variable_params) + list(wfoutputs)
    # figure out what generator to use based on the approximant
    try:
        approximant = static_params['approximant']
    except KeyError:
        raise ValueError("no approximant provided in the static args")
    generator_function = generator.select_waveform_generator(approximant)
    # get data parameters; we'll just use one of the data to get the
    # values, then check that all the others are the same
    delta_f = None
    for d in data.values():
        if delta_f is None:
            delta_f = d.delta_f
            delta_t = d.delta_t
            start_time = d.start_time
        else:
            if not all([d.delta_f == delta_f, d.delta_t == delta_t,
                        d.start_time == start_time]):
                raise ValueError("data must all have the same delta_t, "
                                 "delta_f, and start_time")
    waveform_generator = generator.FDomainDetFrameGenerator(
        generator_function, epoch=start_time,
        variable_args=variable_params, detectors=list(data.keys()),
        delta_f=delta_f, delta_t=delta_t,
        recalib=recalibration, gates=gates,
        **static_params)
    return waveform_generator
