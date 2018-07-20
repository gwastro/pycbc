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

import numpy
from scipy import special

from pycbc import filter
from pycbc.waveform import NoWaveformError
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
    ``BaseModel``.

    Parameters
    ----------
    variable_params : (tuple of) string(s)
        A tuple of parameter names that will be varied.
    waveform_generator : generator class
        A generator class that creates waveforms. This must have a ``generate``
        function which takes parameter values as keyword arguments, a
        detectors attribute which is a dictionary of detectors keyed by their
        names, and an epoch which specifies the start time of the generated
        waveform.
    data : dict
        A dictionary of data, in which the keys are the detector names and the
        values are the data (assumed to be unwhitened). The list of keys must
        match the waveform generator's detectors keys, and the epoch of every
        data set must be the same as the waveform generator's epoch.
    f_lower : float
        The starting frequency to use for computing inner products.
    psds : {None, dict}
        A dictionary of FrequencySeries keyed by the detector names. The
        dictionary must have a psd for each detector specified in the data
        dictionary. If provided, the inner products in each detector will be
        weighted by 1/psd of that detector.
    f_upper : {None, float}
        The ending frequency to use for computing inner products. If not
        provided, the minimum of the largest frequency stored in the data
        and a given waveform will be used.
    norm : {None, float or array}
        An extra normalization weight to apply to the inner products. Can be
        either a float or an array. If ``None``, ``4*data.values()[0].delta_f``
        will be used.
    **kwargs :
        All other keyword arguments are passed to ``BaseDataModel``.

    Examples
    --------
    Create a signal, and set up the model on that signal:

    >>> from pycbc import psd as pypsd
    >>> from pycbc.waveform.generator import (FDomainDetFrameGenerator,
                                              FDomainCBCGenerator)
    >>> import gwin
    >>> seglen = 4
    >>> sample_rate = 2048
    >>> N = seglen*sample_rate/2+1
    >>> fmin = 30.
    >>> m1, m2, s1z, s2z, tsig, ra, dec, pol, dist = (
            38.6, 29.3, 0., 0., 3.1, 1.37, -1.26, 2.76, 3*500.)
    >>> variable_params = ['tc']
    >>> generator = FDomainDetFrameGenerator(
            FDomainCBCGenerator, 0.,
            variable_args=variable_params, detectors=['H1', 'L1'],
            delta_f=1./seglen, f_lower=fmin,
            approximant='SEOBNRv2_ROM_DoubleSpin',
            mass1=m1, mass2=m2, spin1z=s1z, spin2z=s2z,
            ra=ra, dec=dec, polarization=pol, distance=dist)
    >>> signal = generator.generate(tc=tsig)
    >>> psd = pypsd.aLIGOZeroDetHighPower(N, 1./seglen, 20.)
    >>> psds = {'H1': psd, 'L1': psd}
    >>> model = gwin.models.GaussianNoise(
            variable_params, signal, generator, fmin, psds=psds)

    Set the current position to the coalescence time of the signal:

    >>> model.update(tc=tsig)

    Now compute the log likelihood ratio and prior-weighted likelihood ratio;
    since we have not provided a prior, these should be equal to each other:

    >>> model.loglr
    278.9612860719217
    >>> model.logplr
    278.9612860719217

    Print all of the default_stats:

    >>> model.current_stats
    {'H1_cplx_loglr': (175.56552899471038+0j),
     'H1_optimal_snrsq': 351.13105798942075,
     'L1_cplx_loglr': (103.39575707721129+0j),
     'L1_optimal_snrsq': 206.79151415442257,
     'logjacobian': 0.0,
     'loglr': 278.9612860719217,
     'logprior': 0.0}

    Compute the SNR; for this system and PSD, this should be approximately 24:

    >>> from pycbc.conversions import snr_from_loglr
    >>> snr_from_loglr(model.loglr)
    23.62038467391764

    Since there is no noise, the SNR should be the same as the quadrature sum
    of the optimal SNRs in each detector:

    >>> (model.det_optimal_snrsq('H1') + model.det_optimal_snrsq('L1'))**0.5
    23.62038467391764

    Using the same model, evaluate the log likelihood ratio at several points
    in time and check that the max is at tsig:

    >>> import numpy
    >>> times = numpy.arange(seglen*sample_rate)/float(sample_rate)
    >>> loglrs = numpy.zeros(len(times))
    >>> for (ii, t) in enumerate(times):
            model.update(tc=t)
            loglrs[ii] = model.loglr
    >>> print('tsig: {}, time of max loglr: {}'.format(
            tsig, times[loglrs.argmax()]))
    tsig: 3.1, time of max loglr: 3.10009765625

    Create a prior and use it (see distributions module for more details):

    >>> from pycbc import distributions
    >>> uniform_prior = distributions.Uniform(tc=(tsig-0.2,tsig+0.2))
    >>> prior = distributions.JointDistribution(variable_params, uniform_prior)
    >>> model = gwin.models.GaussianNoise(variable_params,
            signal, generator, 20., psds=psds, prior=prior)
    >>> model.update(tc=tsig)
    >>> model.logplr
    279.8775768037958
    >>> model.current_stats
    {'H1_cplx_loglr': (175.56552899471038+0j),
     'H1_optimal_snrsq': 351.13105798942075,
     'L1_cplx_loglr': (103.39575707721127+0j),
     'L1_optimal_snrsq': 206.79151415442254,
     'logjacobian': 0.0,
     'loglr': 278.9612860719217,
     'logprior': 0.9162907318741542}

    """
    name = 'gaussian_noise'

    def __init__(self, variable_params, data, waveform_generator,
                 f_lower, psds=None, f_upper=None, norm=None,
                 **kwargs):
        # set up the boiler-plate attributes; note: we'll compute the
        # log evidence later
        super(GaussianNoise, self).__init__(variable_params, data,
                                            waveform_generator, **kwargs)
        # check that the data and waveform generator have the same detectors
        if (sorted(waveform_generator.detectors.keys()) !=
                sorted(self._data.keys())):
            raise ValueError(
                "waveform generator's detectors ({0}) does not "
                "match data ({1})".format(
                    ','.join(sorted(waveform_generator.detector_names)),
                    ','.join(sorted(self._data.keys()))))
        # check that the data and waveform generator have the same epoch
        if any(waveform_generator.epoch != d.epoch
               for d in self._data.values()):
            raise ValueError("waveform generator does not have the same epoch "
                             "as all of the data sets.")
        # check that the data sets all have the same lengths
        dlens = numpy.array([len(d) for d in data.values()])
        if not all(dlens == dlens[0]):
            raise ValueError("all data must be of the same length")
        # we'll use the first data set for setting values
        d = data.values()[0]
        N = len(d)
        # figure out the kmin, kmax to use
        kmin, kmax = filter.get_cutoff_indices(f_lower, f_upper, d.delta_f,
                                               (N-1)*2)
        self._kmin = kmin
        self._kmax = kmax
        if norm is None:
            norm = 4*d.delta_f
        # we'll store the weight to apply to the inner product
        if psds is None:
            w = Array(numpy.sqrt(norm)*numpy.ones(N))
            self._weight = {det: w for det in data}
        else:
            # temporarily suppress numpy divide by 0 warning
            numpysettings = numpy.seterr(divide='ignore')
            self._weight = {det: Array(numpy.sqrt(norm/psds[det]))
                            for det in data}
            numpy.seterr(**numpysettings)
        # whiten the data
        for det in self._data:
            self._data[det][kmin:kmax] *= self._weight[det][kmin:kmax]

    @property
    def default_stats(self):
        """The stats that ``get_current_stats`` returns by default."""
        return ['logjacobian', 'logprior', 'loglr'] + \
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
                kmin = self._kmin
                kmax = self._kmax
                det_lognls[det] = -0.5 * d[kmin:kmax].inner(d[kmin:kmax]).real
            self.__det_lognls = det_lognls
            self.__lognl = sum(det_lognls.values())
            return self.__lognl

    def _nowaveform_loglr(self):
        """Convenience function to set loglr values if no waveform generated.
        """
        for det in self._data:
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
            kmax = min(len(h), self._kmax)
            if self._kmin >= kmax:
                # if the waveform terminates before the filtering low frequency
                # cutoff, then the loglr is just 0 for this detector
                cplx_hd = 0j
                hh = 0.
            else:
                slc = slice(self._kmin, kmax)
                # whiten the waveform
                h[self._kmin:kmax] *= self._weight[det][slc]
                # the inner products
                cplx_hd = self.data[det][slc].inner(h[slc])  # <h, d>
                hh = h[slc].inner(h[slc]).real  # < h, h>
            cplx_loglr = cplx_hd - 0.5*hh
            # store
            setattr(self._current_stats, '{}_optimal_snrsq'.format(det), hh)
            setattr(self._current_stats, '{}_cplx_loglr'.format(det),
                    cplx_loglr)
            lr += cplx_loglr.real
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
            self.loglr
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
            self.loglr
            # now try returning again
            return getattr(self._current_stats, '{}_optimal_snrsq'.format(det))


class MarginalizedPhaseGaussianNoise(GaussianNoise):
    r"""The likelihood is analytically marginalized over phase.

    This class can be used with signal models that can be written as:

    .. math::

        \tilde{h}(f; \Theta, \phi) = A(f; \Theta)e^{i\Psi(f; \Theta) + i \phi},

    where :math:`\phi` is an arbitrary phase constant. This phase constant
    can be analytically marginalized over with a uniform prior as follows:
    assuming the noise is stationary and Gaussian (see `GaussianNoise`
    for details), the posterior is:

    .. math::

        p(\Theta,\phi|d)
            &\propto p(\Theta)p(\phi)p(d|\Theta,\phi) \\
            &\propto p(\Theta)\frac{1}{2\pi}\exp\left[
                -\frac{1}{2}\sum_{i}^{N_D} \left<
                    h_i(\Theta,\phi) - d_i, h_i(\Theta,\phi) - d_i
                \right>\right].

    Here, the sum is over the number of detectors :math:`N_D`, :math:`d_i`
    and :math:`h_i` are the data and signal in detector :math:`i`,
    respectively, and we have assumed a uniform prior on :math:`phi \in [0,
    2\pi)`. With the form of the signal model given above, the inner product
    in the exponent can be written as:

    .. math::

        -\frac{1}{2}\left<h_i - d_i, h_i- d_i\right>
            &= \left<h_i, d_i\right> -
               \frac{1}{2}\left<h_i, h_i\right> -
               \frac{1}{2}\left<d_i, d_i\right> \\
            &= \Re\left\{O(h^0_i, d_i)e^{-i\phi}\right\} -
               \frac{1}{2}\left<h^0_i, h^0_i\right> -
               \frac{1}{2}\left<d_i, d_i\right>,

    where:

    .. math::

        h_i^0 &\equiv \tilde{h}_i(f; \Theta, \phi=0); \\
        O(h^0_i, d_i) &\equiv 4 \int_0^\infty
            \frac{\tilde{h}_i^*(f; \Theta,0)\tilde{d}_i(f)}{S_n(f)}\mathrm{d}f.

    Gathering all of the terms that are not dependent on :math:`\phi` together:

    .. math::

        \alpha(\Theta, d) \equiv \exp\left[-\frac{1}{2}\sum_i
            \left<h^0_i, h^0_i\right> + \left<d_i, d_i\right>\right],

    we can marginalize the posterior over :math:`\phi`:

    .. math::

        p(\Theta|d)
            &\propto p(\Theta)\alpha(\Theta,d)\frac{1}{2\pi}
                     \int_{0}^{2\pi}\exp\left[\Re \left\{
                         e^{-i\phi} \sum_i O(h^0_i, d_i)
                     \right\}\right]\mathrm{d}\phi \\
            &\propto p(\Theta)\alpha(\Theta, d)\frac{1}{2\pi}
                     \int_{0}^{2\pi}\exp\left[
                         x(\Theta,d)\cos(\phi) + y(\Theta, d)\sin(\phi)
                     \right]\mathrm{d}\phi.

    The integral in the last line is equal to :math:`2\pi I_0(\sqrt{x^2+y^2})`,
    where :math:`I_0` is the modified Bessel function of the first kind. Thus
    the marginalized log posterior is:

    .. math::

        \log p(\Theta|d) \propto \log p(\Theta) +
            I_0\left(\left|\sum_i O(h^0_i, d_i)\right|\right) -
            \frac{1}{2}\sum_i\left[ \left<h^0_i, h^0_i\right> -
                                    \left<d_i, d_i\right> \right]

    This class computes the above expression for the log likelihood.
    """
    name = 'marginalized_phase'

    @property
    def default_stats(self):
        """The stats that ``get_current_stats`` returns by default."""
        return ['logjacobian', 'logprior', 'loglr'] + \
               ['{}_optimal_snrsq'.format(det) for det in self._data]

    def _loglr(self):
        r"""Computes the log likelihood ratio,

        .. math::

            \log \mathcal{L}(\Theta) =
                I_0 \left(\left|\sum_i O(h^0_i, d_i)\right|\right) -
                \frac{1}{2}\left<h^0_i, h^0_i\right>,

        at the current point in parameter space :math:`\Theta`.

        Returns
        -------
        float
            The value of the log likelihood ratio evaluated at the given point.
        """
        params = self.current_params
        try:
            wfs = self._waveform_generator.generate(**params)
        except NoWaveformError:
            return self._nowaveform_loglr()
        hh = 0.
        hd = 0j
        for det, h in wfs.items():
            # the kmax of the waveforms may be different than internal kmax
            kmax = min(len(h), self._kmax)
            if self._kmin >= kmax:
                # if the waveform terminates before the filtering low frequency
                # cutoff, then the loglr is just 0 for this detector
                hh_i = 0.
                hd_i = 0j
            else:
                # whiten the waveform
                h[self._kmin:kmax] *= self._weight[det][self._kmin:kmax]
                # calculate inner products
                hh_i = h[self._kmin:kmax].inner(h[self._kmin:kmax]).real
                hd_i = self.data[det][self._kmin:kmax].inner(
                    h[self._kmin:kmax])
            # store
            setattr(self._current_stats, '{}_optimal_snrsq'.format(det), hh_i)
            # TODO: is it possible to store a something like cplx_loglr for
            # each detector?
            hh += hh_i
            hd += hd_i
        hd = abs(hd)
        return numpy.log(special.i0e(hd)) + hd - 0.5*hh
