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
import numpy
from pycbc import filter as pyfilter
from pycbc.waveform import NoWaveformError
from pycbc.waveform import generator
from pycbc.types import Array

from .base_data import BaseDataModel


class GaussianNoise(BaseDataModel):
    r"""Model that assumes data is stationary Gaussian noise.

    With Gaussian noise the log likelihood functions for signal
    :math:`\log p(d|\Theta)` and for noise :math:`\log p(d|n)` are given by:

    .. math::

        \log p(d|\Theta) &=  -\frac{1}{2} \sum_i
            \left< h_i(\Theta) - d_i | h_i(\Theta) - d_i \right> \\
        \log p(d|n) &= -\frac{1}{2} \sum_i \left<d_i | d_i\right>

    where the sum is over the number of detectors, :math:`d_i` is the data in
    each detector, and :math:`h_i(\Theta)` is the model signal in each
    detector. The inner product is given by:

    .. math::

        \left<a | b\right> = 4\Re \int_{0}^{\infty}
            \frac{\tilde{a}(f) \tilde{b}(f)}{S_n(f)} \mathrm{d}f,

    where :math:`S_n(f)` is the PSD in the given detector.

    Note that the log prior-weighted likelihood ratio has one fewer term
    than the log posterior, since the :math:`\left<d_i|d_i\right>` term cancels
    in the likelihood ratio:

    .. math::

        \log \hat{\mathcal{L}} = \log p(\Theta) + \sum_i \left[
            \left<h_i(\Theta)|d_i\right> -
            \frac{1}{2} \left<h_i(\Theta)|h_i(\Theta)\right> \right]

    Upon initialization, the data is whitened using the given PSDs. If no PSDs
    are given the data and waveforms returned by the waveform generator are
    assumed to be whitened. The likelihood function of the noise,

    .. math::

        p(d|n) = \frac{1}{2} \sum_i \left<d_i|d_i\right>,

    is computed on initialization and stored as the `lognl` attribute.

    By default, the data is assumed to be equally sampled in frequency, but
    unequally sampled data can be supported by passing the appropriate
    normalization using the ``norm`` keyword argument.

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
    psds : {None, dict}
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
    norm : {None, float or array}
        An extra normalization weight to apply to the inner products. Can be
        either a float or an array. If ``None``, ``4*data.values()[0].delta_f``
        will be used.
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
    282.29
    >>> print('{:.2f}'.format(model.logplr))
    282.29

    Print all of the default_stats:

    >>> print(',\n'.join(['{}: {:.2f}'.format(s, v)
    ...                   for (s, v) in sorted(model.current_stats.items())]))
    H1_cplx_loglr: 177.66+0.00j,
    H1_optimal_snrsq: 355.33,
    L1_cplx_loglr: 104.63+0.00j,
    L1_optimal_snrsq: 209.26,
    logjacobian: 0.00,
    loglikelihood: 0.00,
    loglr: 282.29,
    logprior: 0.00

    Compute the SNR; for this system and PSD, this should be approximately 24:

    >>> from pycbc.conversions import snr_from_loglr
    >>> x = snr_from_loglr(model.loglr)
    >>> print('{:.2f}'.format(x))
    23.76

    Since there is no noise, the SNR should be the same as the quadrature sum
    of the optimal SNRs in each detector:

    >>> x = (model.det_optimal_snrsq('H1') +
    ...      model.det_optimal_snrsq('L1'))**0.5
    >>> print('{:.2f}'.format(x))
    23.76

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
    >>> model = pycbc.inference.models.GaussianNoise(variable_params,
    ...     signal, generator, low_frequency_cutoff, psds=psds, prior=prior)
    >>> model.update(tc=tsig)
    >>> print('{:.2f}'.format(model.logplr))
    283.21
    >>> print(',\n'.join(['{}: {:.2f}'.format(s, v)
    ...                   for (s, v) in sorted(model.current_stats.items())]))
    H1_cplx_loglr: 177.66+0.00j,
    H1_optimal_snrsq: 355.33,
    L1_cplx_loglr: 104.63+0.00j,
    L1_optimal_snrsq: 209.26,
    logjacobian: 0.00,
    loglikelihood: 0.00,
    loglr: 282.29,
    logprior: 0.92

    """
    name = 'gaussian_noise'

    def __init__(self, variable_params, data, low_frequency_cutoff, psds=None,
                 high_frequency_cutoff=None, norm=None, static_params=None,
                 **kwargs):
        # set up the boiler-plate attributes
        super(GaussianNoise, self).__init__(variable_params, data,
                                            static_params=static_params,
                                            **kwargs)
        # check if low frequency cutoff has been provided for every IFO.
        if set(low_frequency_cutoff.keys()) != set(self.data.keys()):
            raise KeyError("IFOs for which data has been provided should "
                           "match IFOs for which low-frequency-cutoff has "
                           "been provided for likelihood inner product "
                           "calculation. If loading the model settings from "
                           "a config file, please provide an "
                           "`IFO-low-frequency-cutoff` input for each "
                           "detector in the `[model]` section, where IFO is "
                           "the name of the detector.")
        # create the waveform generator
        self._waveform_generator = create_waveform_generator(
            self.variable_params, self.data, recalibration=self.recalibration,
            gates=self.gates, **self.static_params)
        # check that the data sets all have the same lengths
        dlens = numpy.array([len(d) for d in self.data.values()])
        if not all(dlens == dlens[0]):
            raise ValueError("all data must be of the same length")
        # we'll use the first data set for setting values
        d = list(self.data.values())[0]
        N = len(d)
        # Set low frequency cutoff
        self._f_lower = low_frequency_cutoff
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
        if norm is None:
            norm = 4*d.delta_f
        # we'll store the weight to apply to the inner product
        if psds is None:
            self._psds = None
            w = Array(numpy.sqrt(norm)*numpy.ones(N))
            self._weight = {det: w for det in data}
        else:
            # store a copy of the psds
            self._psds = {ifo: d.copy() for (ifo, d) in psds.items()}
            # temporarily suppress numpy divide by 0 warning
            numpysettings = numpy.seterr(divide='ignore')
            self._weight = {det: Array(numpy.sqrt(norm/psds[det]))
                            for det in data}
            numpy.seterr(**numpysettings)
        self._kmin = {}
        self._kmax = {}
        for det in self._data:
            # whiten the data
            kmin, kmax = pyfilter.get_cutoff_indices(self._f_lower[det],
                                                     self._f_upper[det],
                                                     d.delta_f, (N-1)*2)
            self._data[det][kmin:kmax] *= self._weight[det][kmin:kmax]
            self._kmin[det] = kmin
            self._kmax[det] = kmax

    @property
    def waveform_generator(self):
        """The waveform generator used."""
        return self._waveform_generator

    @property
    def low_frequency_cutoff(self):
        """The low frequency cutoff of the inner product."""
        return self._f_lower

    @property
    def high_frequency_cutoff(self):
        """The high frequency cutoff of the inner product."""
        return self._f_upper

    @property
    def _extra_stats(self):
        """Adds ``loglr``, plus ``cplx_loglr`` and ``optimal_snrsq`` in each
        detector."""
        return ['loglr'] + \
               ['{}_cplx_loglr'.format(det) for det in self._data] + \
               ['{}_optimal_snrsq'.format(det) for det in self._data]

    def _lognl(self):
        """Computes the log likelihood assuming the data is noise.

        Since this is a constant for Gaussian noise, this is only computed once
        then stored.
        """
        try:
            return self.__lognl
        except AttributeError:
            det_lognls = {}
            for (det, d) in self._data.items():
                kmin = self._kmin[det]
                kmax = self._kmax[det]
                det_lognls[det] = -0.5 * d[kmin:kmax].inner(d[kmin:kmax]).real
            self.__det_lognls = det_lognls
            self.__lognl = sum(det_lognls.values())
            return self.__lognl

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
            wfs = self._waveform_generator.generate(**params)
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
                cplx_hd = self.data[det][slc].inner(h[slc])  # <h, d>
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

    def _loglikelihood(self):
        r"""Computes the log likelihood of the paramaters,

        .. math::

            p(d|\Theta) = -\frac{1}{2}\sum_i
                \left<h_i(\Theta) - d_i | h_i(\Theta) - d_i\right>,

        at the current parameter values :math:`\Theta`.

        Returns
        -------
        float
            The value of the log likelihood evaluated at the given point.
        """
        # since the loglr has fewer terms, we'll call that, then just add
        # back the noise term that canceled in the log likelihood ratio
        return self.loglr + self.lognl

    def det_lognl(self, det):
        """Returns the log likelihood of the noise in the given detector.

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
            return self.__det_lognls[det]
        except AttributeError:
            # hasn't been calculated yet, call lognl to calculate & store
            self._lognl()
            # now try returning
            return self.__det_lognls[det]

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

    def write_metadata(self, fp):
        """Adds writing the psds and lognl, since it's a constant.

        The lognl is written to the sample group's ``attrs``.

        Parameters
        ----------
        fp : pycbc.inference.io.BaseInferenceFile instance
            The inference file to write to.
        """
        super(GaussianNoise, self).write_metadata(fp)
        if self._psds is not None:
            fp.write_psd(self._psds)
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
    def from_config(cls, cp, **kwargs):
        r"""Initializes an instance of this class from the given config file.

        Parameters
        ----------
        cp : WorkflowConfigParser
            Config file parser to read.
        \**kwargs :
            All additional keyword arguments are passed to the class. Any
            provided keyword will over ride what is in the config file.
        """
        args = cls._init_args_from_config(cp)
        args['low_frequency_cutoff'] = low_frequency_cutoff_from_config(cp)
        args['high_frequency_cutoff'] = high_frequency_cutoff_from_config(cp)
        # get any other keyword arguments provided in the model section
        ignore_args = ['name']
        for option in cp.options("model"):
            if any([option.endswith("-low-frequency-cutoff"),
                    option.endswith("-high-frequency-cutoff")]):
                ignore_args.append(option)
        args.update(cls.extra_args_from_config(cp, "model",
                                               skip_args=ignore_args))
        args.update(kwargs)
        return cls(**args)


#
# =============================================================================
#
#                               Support functions
#
# =============================================================================
#


def create_waveform_generator(variable_params, data,
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


def low_frequency_cutoff_from_config(cp):
    """Gets the low frequency cutoff for all detectors to be used from the
    given config file.

    The low-frequency-cutoff for each detector should be provided using an
    option ``IFO-low-frequency-cutoff`` in the ``[model]`` section, where IFO
    is the name of the detector. The low frequency cutoff value is then casted
    to float. If the casting to float fails, an error is raised.

    Parameters
    ----------
    cp : WorkflowConfigParser
        Config file parser to read.

    Returns
    -------
    dict :
        Dictionary with the IFO names as the keys and the respective low
        frequency cutoffs to be used in the likelihood calculation as the
        values.
    """
    low_frequency_cutoff = {}
    for option in cp.options("model"):
        if option.endswith("-low-frequency-cutoff"):
            ifo = option.rsplit("-low-frequency-cutoff")[0].upper()
            try:
                low_frequency_cutoff[ifo] = float(cp.get("model", option))
            except Exception as e:
                logging.warning("Low frequency cutoff of %s could not be "
                                "converted to float", ifo)
                raise e
    return low_frequency_cutoff


def high_frequency_cutoff_from_config(cp):
    """Gets the high frequency cutoff, if one is provided for a detector from
    the given config file for likelihood computation.

    This looks for options ``IFO-high-frequency-cutoff`` in the ``[model]``
    section, where IFO is the name of the detector, and casts it to float.

    Parameters
    ----------
    cp : WorkflowConfigParser
        Config file parser to read.

    Returns
    -------
    dict or None :
        Dictionary with the IFO names as the keys and the respective high
        frequency cutoffs to be used in the likelihood calculation as the
        values. If high frequency cutoffs for no detectors are provided,
        returns an empty dictionary.
    """
    high_frequency_cutoff = {}
    for option in cp.options("model"):
        if option.endswith("-high-frequency-cutoff"):
            ifo = option.rsplit("-high-frequency-cutoff")[0].upper()
            try:
                high_frequency_cutoff[ifo] = float(cp.get("model", option))
            except Exception as e:
                logging.warning("High frequency cutoff of %s could not be "
                                "converted to float", ifo)
                raise e
    return high_frequency_cutoff
