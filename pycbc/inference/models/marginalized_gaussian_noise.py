# Copyright (C) 2018  Charlie Hoy, Collin Capano
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
allows for the likelihood to be marginalized over phase and/or time and/or
distance.
"""

import itertools
import numpy
from scipy import special

from pycbc.waveform import generator
from pycbc.waveform import (NoWaveformError, FailedWaveformError)
from pycbc.detector import Detector
from .gaussian_noise import (BaseGaussianNoise,
                             create_waveform_generator,
                             GaussianNoise)
from .tools import marginalize_likelihood, DistMarg


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
    respectively, and we have assumed a uniform prior on :math:`\phi \in [0,
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
    the marginalized posterior is:

    .. math::

        p(\Theta|d) \propto
            I_0\left(\left|\sum_i O(h^0_i, d_i)\right|\right)
             p(\Theta)\exp\left[\frac{1}{2}\sum_i\left( \left<h^0_i, h^0_i\right> -
                                    \left<d_i, d_i\right> \right)\right]
    """
    name = 'marginalized_phase'

    def __init__(self, variable_params, data, low_frequency_cutoff, psds=None,
                 high_frequency_cutoff=None, normalize=False,
                 static_params=None, **kwargs):
        # set up the boiler-plate attributes
        super(MarginalizedPhaseGaussianNoise, self).__init__(
            variable_params, data, low_frequency_cutoff, psds=psds,
            high_frequency_cutoff=high_frequency_cutoff, normalize=normalize,
            static_params=static_params, **kwargs)

    @property
    def _extra_stats(self):
        """Adds ``loglr``, plus ``cplx_loglr`` and ``optimal_snrsq`` in each
        detector."""
        return ['loglr', 'maxl_phase'] + \
               ['{}_optimal_snrsq'.format(det) for det in self._data]

    def _nowaveform_loglr(self):
        """Convenience function to set loglr values if no waveform generated.
        """
        setattr(self._current_stats, 'loglikelihood', -numpy.inf)
        # maxl phase doesn't exist, so set it to nan
        setattr(self._current_stats, 'maxl_phase', numpy.nan)
        for det in self._data:
            # snr can't be < 0 by definition, so return 0
            setattr(self._current_stats, '{}_optimal_snrsq'.format(det), 0.)
        return -numpy.inf

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
            if self.all_ifodata_same_rate_length:
                wfs = self.waveform_generator.generate(**params)
            else:
                wfs = {}
                for det in self.data:
                    wfs.update(self.waveform_generator[det].generate(**params))

        except NoWaveformError:
            return self._nowaveform_loglr()
        except FailedWaveformError as e:
            if self.ignore_failed_waveforms:
                return self._nowaveform_loglr()
            else:
                raise e
        hh = 0.
        hd = 0j
        for det, h in wfs.items():
            # the kmax of the waveforms may be different than internal kmax
            kmax = min(len(h), self._kmax[det])
            if self._kmin[det] >= kmax:
                # if the waveform terminates before the filtering low frequency
                # cutoff, then the loglr is just 0 for this detector
                hh_i = 0.
                hd_i = 0j
            else:
                # whiten the waveform
                h[self._kmin[det]:kmax] *= \
                    self._weight[det][self._kmin[det]:kmax]
                # calculate inner products
                hh_i = h[self._kmin[det]:kmax].inner(
                    h[self._kmin[det]:kmax]).real
                hd_i = h[self._kmin[det]:kmax].inner(
                    self._whitened_data[det][self._kmin[det]:kmax])
            # store
            setattr(self._current_stats, '{}_optimal_snrsq'.format(det), hh_i)
            hh += hh_i
            hd += hd_i
        self._current_stats.maxl_phase = numpy.angle(hd)
        return marginalize_likelihood(hd, hh, phase=True)


class MarginalizedTime(DistMarg, BaseGaussianNoise):
    r""" This likelihood numerically marginalizes over time

    This likelihood is optimized for marginalizing over time, but can also
    handle marginalization over polarization, phase (where appropriate),
    and sky location. The time series is interpolated using a
    quadratic apparoximation for sub-sample times.
    """
    name = 'marginalized_time'

    def __init__(self, variable_params,
                 data, low_frequency_cutoff, psds=None,
                 high_frequency_cutoff=None, normalize=False,
                 **kwargs):

        self.kwargs = kwargs
        variable_params, kwargs = self.setup_marginalization(
                               variable_params,
                               **kwargs)

        # set up the boiler-plate attributes
        super(MarginalizedTime, self).__init__(
            variable_params, data, low_frequency_cutoff, psds=psds,
            high_frequency_cutoff=high_frequency_cutoff, normalize=normalize,
            **kwargs)
        # Determine if all data have the same sampling rate and segment length
        if self.all_ifodata_same_rate_length:
            # create a waveform generator for all ifos
            self.waveform_generator = create_waveform_generator(
                self.variable_params, self.data,
                waveform_transforms=self.waveform_transforms,
                recalibration=self.recalibration,
                generator_class=generator.FDomainDetFrameTwoPolNoRespGenerator,
                gates=self.gates, **kwargs['static_params'])
        else:
            # create a waveform generator for each ifo respestively
            self.waveform_generator = {}
            for det in self.data:
                self.waveform_generator[det] = create_waveform_generator(
                    self.variable_params, {det: self.data[det]},
                    waveform_transforms=self.waveform_transforms,
                    recalibration=self.recalibration,
                    generator_class=generator.FDomainDetFrameTwoPolNoRespGenerator,
                    gates=self.gates, **kwargs['static_params'])

        self.dets = {}

    def _nowaveform_loglr(self):
        """Convenience function to set loglr values if no waveform generated.
        """
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
        from pycbc.filter import matched_filter_core

        params = self.current_params
        try:
            if self.all_ifodata_same_rate_length:
                wfs = self.waveform_generator.generate(**params)
            else:
                wfs = {}
                for det in self.data:
                    wfs.update(self.waveform_generator[det].generate(**params))
        except NoWaveformError:
            return self._nowaveform_loglr()
        except FailedWaveformError as e:
            if self.ignore_failed_waveforms:
                return self._nowaveform_loglr()
            else:
                raise e

        sh_total = hh_total = 0.
        snr_estimate = {}
        cplx_hpd = {}
        cplx_hcd = {}
        hphp = {}
        hchc = {}
        hphc = {}
        for det, (hp, hc) in wfs.items():
            # the kmax of the waveforms may be different than internal kmax
            kmax = min(max(len(hp), len(hc)), self._kmax[det])
            slc = slice(self._kmin[det], kmax)

            # whiten both polarizations
            hp[self._kmin[det]:kmax] *= self._weight[det][slc]
            hc[self._kmin[det]:kmax] *= self._weight[det][slc]

            hp.resize(len(self._whitened_data[det]))
            hc.resize(len(self._whitened_data[det]))
            cplx_hpd[det], _, _ = matched_filter_core(
                                 hp,
                                 self._whitened_data[det],
                                 low_frequency_cutoff=self._f_lower[det],
                                 high_frequency_cutoff=self._f_upper[det],
                                 h_norm=1)
            cplx_hcd[det], _, _ = matched_filter_core(
                                 hc,
                                 self._whitened_data[det],
                                 low_frequency_cutoff=self._f_lower[det],
                                 high_frequency_cutoff=self._f_upper[det],
                                 h_norm=1)

            hphp[det] = hp[slc].inner(hp[slc]).real
            hchc[det] = hc[slc].inner(hc[slc]).real
            hphc[det] = hp[slc].inner(hc[slc]).real

            snr_proxy = ((cplx_hpd[det] / hphp[det] ** 0.5).squared_norm() +
                         (cplx_hcd[det] / hchc[det] ** 0.5).squared_norm())
            snr_estimate[det] = (0.5 * snr_proxy) ** 0.5

        self.draw_ifos(snr_estimate, log=False, **self.kwargs)
        self.snr_draw(snr_estimate)

        for det in wfs:
            if det not in self.dets:
                self.dets[det] = Detector(det)
            fp, fc = self.dets[det].antenna_pattern(
                                    params['ra'],
                                    params['dec'],
                                    params['polarization'],
                                    params['tc'])
            dt = self.dets[det].time_delay_from_earth_center(params['ra'],
                                                             params['dec'],
                                                             params['tc'])
            dtc = params['tc'] + dt
            cplx_hd = fp * cplx_hpd[det].at_time(dtc,
                                                 interpolate='quadratic')
            cplx_hd += fc * cplx_hcd[det].at_time(dtc,
                                                  interpolate='quadratic')
            hh = (fp * fp * hphp[det] +
                  fc * fc * hchc[det] +
                  2.0 * fp * fc * hphc[det])

            sh_total += cplx_hd
            hh_total += hh

        loglr = self.marginalize_loglr(sh_total, hh_total)
        return loglr


class MarginalizedPolarization(DistMarg, BaseGaussianNoise):
    r""" This likelihood numerically marginalizes over polarization angle

    This class implements the Gaussian likelihood with an explicit numerical
    marginalization over polarization angle. This is accomplished using
    a fixed set of integration points distribution uniformation between
    0 and 2pi. By default, 1000 integration points are used.
    The 'polarization_samples' argument can be passed to set an alternate
    number of integration points.
    """
    name = 'marginalized_polarization'

    def __init__(self, variable_params, data, low_frequency_cutoff, psds=None,
                 high_frequency_cutoff=None, normalize=False,
                 polarization_samples=1000,
                 **kwargs):

        variable_params, kwargs = self.setup_marginalization(
                               variable_params,
                               polarization_samples=polarization_samples,
                               **kwargs)

        # set up the boiler-plate attributes
        super(MarginalizedPolarization, self).__init__(
            variable_params, data, low_frequency_cutoff, psds=psds,
            high_frequency_cutoff=high_frequency_cutoff, normalize=normalize,
            **kwargs)
        # Determine if all data have the same sampling rate and segment length
        if self.all_ifodata_same_rate_length:
            # create a waveform generator for all ifos
            self.waveform_generator = create_waveform_generator(
                self.variable_params, self.data,
                waveform_transforms=self.waveform_transforms,
                recalibration=self.recalibration,
                generator_class=generator.FDomainDetFrameTwoPolGenerator,
                gates=self.gates, **kwargs['static_params'])
        else:
            # create a waveform generator for each ifo respestively
            self.waveform_generator = {}
            for det in self.data:
                self.waveform_generator[det] = create_waveform_generator(
                    self.variable_params, {det: self.data[det]},
                    waveform_transforms=self.waveform_transforms,
                    recalibration=self.recalibration,
                    generator_class=generator.FDomainDetFrameTwoPolGenerator,
                    gates=self.gates, **kwargs['static_params'])

        self.dets = {}

    @property
    def _extra_stats(self):
        """Adds ``loglr``, ``maxl_polarization``, and the ``optimal_snrsq`` in
        each detector.
        """
        return ['loglr', 'maxl_polarization', 'maxl_loglr'] + \
               ['{}_optimal_snrsq'.format(det) for det in self._data]

    def _nowaveform_loglr(self):
        """Convenience function to set loglr values if no waveform generated.
        """
        setattr(self._current_stats, 'loglr', -numpy.inf)
        # maxl phase doesn't exist, so set it to nan
        setattr(self._current_stats, 'maxl_polarization', numpy.nan)
        for det in self._data:
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
            if self.all_ifodata_same_rate_length:
                wfs = self.waveform_generator.generate(**params)
            else:
                wfs = {}
                for det in self.data:
                    wfs.update(self.waveform_generator[det].generate(**params))
        except NoWaveformError:
            return self._nowaveform_loglr()
        except FailedWaveformError as e:
            if self.ignore_failed_waveforms:
                return self._nowaveform_loglr()
            else:
                raise e

        lr = sh_total = hh_total = 0.
        for det, (hp, hc) in wfs.items():
            if det not in self.dets:
                self.dets[det] = Detector(det)
            fp, fc = self.dets[det].antenna_pattern(
                                    params['ra'],
                                    params['dec'],
                                    params['polarization'],
                                    params['tc'])

            # the kmax of the waveforms may be different than internal kmax
            kmax = min(max(len(hp), len(hc)), self._kmax[det])
            slc = slice(self._kmin[det], kmax)

            # whiten both polarizations
            hp[self._kmin[det]:kmax] *= self._weight[det][slc]
            hc[self._kmin[det]:kmax] *= self._weight[det][slc]

            # h = fp * hp + hc * hc
            # <h, d> = fp * <hp,d> + fc * <hc,d>
            # the inner products
            cplx_hpd = hp[slc].inner(self._whitened_data[det][slc])  # <hp, d>
            cplx_hcd = hc[slc].inner(self._whitened_data[det][slc])  # <hc, d>

            cplx_hd = fp * cplx_hpd + fc * cplx_hcd

            # <h, h> = <fp * hp + fc * hc, fp * hp + fc * hc>
            # = Real(fpfp * <hp,hp> + fcfc * <hc,hc> + \
            #  fphc * (<hp, hc> + <hc, hp>))
            hphp = hp[slc].inner(hp[slc]).real  # < hp, hp>
            hchc = hc[slc].inner(hc[slc]).real  # <hc, hc>

            # Below could be combined, but too tired to figure out
            # if there should be a sign applied if so
            hphc = hp[slc].inner(hc[slc]).real  # <hp, hc>
            hchp = hc[slc].inner(hp[slc]).real  # <hc, hp>

            hh = fp * fp * hphp + fc * fc * hchc + fp * fc * (hphc + hchp)
            # store
            setattr(self._current_stats, '{}_optimal_snrsq'.format(det), hh)
            sh_total += cplx_hd
            hh_total += hh

        lr, idx, maxl = self.marginalize_loglr(sh_total, hh_total,
                  return_peak=True)

        # store the maxl polarization
        setattr(self._current_stats,
                'maxl_polarization',
                params['polarization'])
        setattr(self._current_stats, 'maxl_loglr', maxl)

        # just store the maxl optimal snrsq
        for det in wfs:
            p = '{}_optimal_snrsq'.format(det)
            setattr(self._current_stats, p,
                    getattr(self._current_stats, p)[idx])

        return lr


class MarginalizedHMPolPhase(BaseGaussianNoise):
    r"""Numerically marginalizes waveforms with higher modes over polarization
    `and` phase.

    This class implements the Gaussian likelihood with an explicit numerical
    marginalization over polarization angle and orbital phase. This is
    accomplished using a fixed set of integration points distributed uniformly
    between 0 and 2:math:`\pi` for both the polarization and phase. By default,
    100 integration points are used for each parameter, giving :math:`10^4`
    evaluation points in total. This can be modified using the
    ``polarization_samples`` and ``coa_phase_samples`` arguments.

    This only works with waveforms that return separate spherical harmonic
    modes for each waveform. For a list of currently supported approximants,
    see :py:func:`pycbc.waveform.waveform_modes.fd_waveform_mode_approximants`
    and :py:func:`pycbc.waveform.waveform_modes.td_waveform_mode_approximants`.

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
    polarization_samples : int, optional
        How many points to use in polarization. Default is 100.
    coa_phase_samples : int, optional
        How many points to use in phase. Defaults is 100.
    \**kwargs :
        All other keyword arguments are passed to
        :py:class:`BaseGaussianNoise
        <pycbc.inference.models.gaussian_noise.BaseGaussianNoise>`.

    """
    name = 'marginalized_hmpolphase'

    def __init__(self, variable_params, data, low_frequency_cutoff, psds=None,
                 high_frequency_cutoff=None, normalize=False,
                 polarization_samples=100,
                 coa_phase_samples=100,
                 static_params=None, **kwargs):
        # set up the boiler-plate attributes
        super(MarginalizedHMPolPhase, self).__init__(
            variable_params, data, low_frequency_cutoff, psds=psds,
            high_frequency_cutoff=high_frequency_cutoff, normalize=normalize,
            static_params=static_params, **kwargs)
        # create the waveform generator
        self.waveform_generator = create_waveform_generator(
            self.variable_params, self.data,
            waveform_transforms=self.waveform_transforms,
            recalibration=self.recalibration,
            generator_class=generator.FDomainDetFrameModesGenerator,
            gates=self.gates, **self.static_params)
        pol = numpy.linspace(0, 2*numpy.pi, polarization_samples)
        phase = numpy.linspace(0, 2*numpy.pi, coa_phase_samples)
        # remap to every combination of the parameters
        # this gets every combination by mappin them to an NxM grid
        # one needs to be transposed so that they run allong opposite
        # dimensions
        n = coa_phase_samples * polarization_samples
        self.nsamples = n
        self.pol = numpy.resize(pol, n)
        phase = numpy.resize(phase, n)
        phase = phase.reshape(coa_phase_samples, polarization_samples)
        self.phase = phase.T.flatten()
        self._phase_fac = {}
        self.dets = {}

    def phase_fac(self, m):
        r"""The phase :math:`\exp[i m \phi]`."""
        try:
            return self._phase_fac[m]
        except KeyError:
            # hasn't been computed yet, calculate it
            self._phase_fac[m] = numpy.exp(1.0j * m * self.phase)
            return self._phase_fac[m]

    @property
    def _extra_stats(self):
        """Adds ``maxl_polarization`` and the ``maxl_phase``
        """
        return ['maxl_polarization', 'maxl_phase', ]

    def _nowaveform_loglr(self):
        """Convenience function to set loglr values if no waveform generated.
        """
        # maxl phase doesn't exist, so set it to nan
        setattr(self._current_stats, 'maxl_polarization', numpy.nan)
        setattr(self._current_stats, 'maxl_phase', numpy.nan)
        return -numpy.inf

    def _loglr(self, return_unmarginalized=False):
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
        except FailedWaveformError as e:
            if self.ignore_failed_waveforms:
                return self._nowaveform_loglr()
            else:
                raise e

        # ---------------------------------------------------------------------
        # Some optimizations not yet taken:
        # * higher m calculations could have a lot of redundancy
        # * fp/fc need not be calculated except where polarization is different
        # * may be possible to simplify this by making smarter use of real/imag
        # ---------------------------------------------------------------------
        lr = 0.
        hds = {}
        hhs = {}
        for det, modes in wfs.items():
            if det not in self.dets:
                self.dets[det] = Detector(det)

            fp, fc = self.dets[det].antenna_pattern(params['ra'],
                                                    params['dec'],
                                                    self.pol,
                                                    params['tc'])

            # loop over modes and prepare the waveform modes
            # we will sum up zetalm = glm <ulm, d> + i glm <vlm, d>
            # over all common m so that we can apply the phase once
            zetas = {}
            rlms = {}
            slms = {}
            for mode in modes:
                l, m = mode
                ulm, vlm = modes[mode]

                # whiten the waveforms
                # the kmax of the waveforms may be different than internal kmax
                kmax = min(max(len(ulm), len(vlm)), self._kmax[det])
                slc = slice(self._kmin[det], kmax)
                ulm[self._kmin[det]:kmax] *= self._weight[det][slc]
                vlm[self._kmin[det]:kmax] *= self._weight[det][slc]

                # the inner products
                # <ulm, d>
                ulmd = ulm[slc].inner(self._whitened_data[det][slc]).real
                # <vlm, d>
                vlmd = vlm[slc].inner(self._whitened_data[det][slc]).real

                # add inclination, and pack into a complex number
                import lal
                glm = lal.SpinWeightedSphericalHarmonic(
                    params['inclination'], 0, -2, l, m).real

                if m not in zetas:
                    zetas[m] = 0j
                zetas[m] += glm * (ulmd + 1j*vlmd)

                # Get condense set of the parts of the waveform that only diff
                # by m, this is used next to help calculate <h, h>
                r = glm * ulm
                s = glm * vlm

                if m not in rlms:
                    rlms[m] = r
                    slms[m] = s
                else:
                    rlms[m] += r
                    slms[m] += s

            # now compute all possible <hlm, hlm>
            rr_m = {}
            ss_m = {}
            rs_m = {}
            sr_m = {}
            combos = itertools.combinations_with_replacement(rlms.keys(), 2)
            for m, mprime in combos:
                r = rlms[m]
                s = slms[m]
                rprime = rlms[mprime]
                sprime = slms[mprime]
                rr_m[mprime, m] = r[slc].inner(rprime[slc]).real
                ss_m[mprime, m] = s[slc].inner(sprime[slc]).real
                rs_m[mprime, m] = s[slc].inner(rprime[slc]).real
                sr_m[mprime, m] = r[slc].inner(sprime[slc]).real
                # store the conjugate for easy retrieval later
                rr_m[m, mprime] = rr_m[mprime, m]
                ss_m[m, mprime] = ss_m[mprime, m]
                rs_m[m, mprime] = sr_m[mprime, m]
                sr_m[m, mprime] = rs_m[mprime, m]
            # now apply the phase to all the common ms
            hpd = 0.
            hcd = 0.
            hphp = 0.
            hchc = 0.
            hphc = 0.
            for m, zeta in zetas.items():
                phase_coeff = self.phase_fac(m)

                # <h+, d> = (exp[i m phi] * zeta).real()
                # <hx, d> = -(exp[i m phi] * zeta).imag()
                z = phase_coeff * zeta
                hpd += z.real
                hcd -= z.imag

                # now calculate the contribution to <h, h>
                cosm = phase_coeff.real
                sinm = phase_coeff.imag

                for mprime in zetas:
                    pcprime = self.phase_fac(mprime)

                    cosmprime = pcprime.real
                    sinmprime = pcprime.imag
                    # needed components
                    rr = rr_m[m, mprime]
                    ss = ss_m[m, mprime]
                    rs = rs_m[m, mprime]
                    sr = sr_m[m, mprime]
                    # <hp, hp>
                    hphp += rr * cosm * cosmprime \
                        + ss * sinm * sinmprime \
                        - rs * cosm * sinmprime \
                        - sr * sinm * cosmprime
                    # <hc, hc>
                    hchc += rr * sinm * sinmprime \
                        + ss * cosm * cosmprime \
                        + rs * sinm * cosmprime \
                        + sr * cosm * sinmprime
                    # <hp, hc>
                    hphc += -rr * cosm * sinmprime \
                        + ss * sinm * cosmprime \
                        + sr * sinm * sinmprime \
                        - rs * cosm * cosmprime

            # Now apply the polarizations and calculate the loglr
            # We have h = Fp * hp + Fc * hc
            # loglr = <h, d> - <h, h>/2
            #       = Fp*<hp, d> + Fc*<hc, d>
            #          - (1/2)*(Fp*Fp*<hp, hp> + Fc*Fc*<hc, hc>
            #                   + 2*Fp*Fc<hp, hc>)
            # (in the last line we have made use of the time series being
            #  real, so that <a, b> = <b, a>).
            hd = fp * hpd + fc * hcd
            hh = fp * fp * hphp + fc * fc * hchc + 2 * fp * fc * hphc
            hds[det] = hd
            hhs[det] = hh
            lr += hd - 0.5 * hh

        if return_unmarginalized:
            return self.pol, self.phase, lr, hds, hhs

        lr_total = special.logsumexp(lr) - numpy.log(self.nsamples)

        # store the maxl values
        idx = lr.argmax()
        setattr(self._current_stats, 'maxl_polarization', self.pol[idx])
        setattr(self._current_stats, 'maxl_phase', self.phase[idx])
        return float(lr_total)
