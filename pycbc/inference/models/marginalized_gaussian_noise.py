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

import numpy
from scipy import special

from pycbc.waveform import generator
from pycbc.waveform import (NoWaveformError, FailedWaveformError)
from pycbc.detector import Detector
from .gaussian_noise import (BaseGaussianNoise, create_waveform_generator)


class MarginalizedPhaseGaussianNoise(BaseGaussianNoise):
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
            wfs = self.waveform_generator.generate(**params)
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
                hd_i = self._whitened_data[det][self._kmin[det]:kmax].inner(
                    h[self._kmin[det]:kmax])
            # store
            setattr(self._current_stats, '{}_optimal_snrsq'.format(det), hh_i)
            hh += hh_i
            hd += hd_i
        hd = abs(hd)
        self._current_stats.maxl_phase = numpy.angle(hd)
        return numpy.log(special.i0e(hd)) + hd - 0.5*hh


class MarginalizedPolarization(BaseGaussianNoise):
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
                 static_params=None, **kwargs):
        # set up the boiler-plate attributes
        super(MarginalizedPolarization, self).__init__(
            variable_params, data, low_frequency_cutoff, psds=psds,
            high_frequency_cutoff=high_frequency_cutoff, normalize=normalize,
            static_params=static_params, **kwargs)
        # create the waveform generator
        self.waveform_generator = create_waveform_generator(
            self.variable_params, self.data,
            waveform_transforms=self.waveform_transforms,
            recalibration=self.recalibration,
            generator_class=generator.FDomainDetFrameTwoPolGenerator,
            gates=self.gates, **self.static_params)

        self.polarization_samples = polarization_samples
        self.pol = numpy.linspace(0, 2*numpy.pi, self.polarization_samples)
        self.dets = {}

    @property
    def _extra_stats(self):
        """Adds ``loglr``, ``maxl_polarization``, and the ``optimal_snrsq`` in
        each detector.
        """
        return ['loglr', 'maxl_polarization'] + \
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
            wfs = self.waveform_generator.generate(**params)
        except NoWaveformError:
            return self._nowaveform_loglr()
        except FailedWaveformError as e:
            if self.ignore_failed_waveforms:
                return self._nowaveform_loglr()
            else:
                raise e

        lr = 0.
        for det, (hp, hc) in wfs.items():
            if det not in self.dets:
                self.dets[det] = Detector(det)
            fp, fc = self.dets[det].antenna_pattern(self.current_params['ra'],
                                                    self.current_params['dec'],
                                                    self.pol,
                                                    self.current_params['tc'])

            # the kmax of the waveforms may be different than internal kmax
            kmax = min(max(len(hp), len(hc)), self._kmax[det])
            slc = slice(self._kmin[det], kmax)

            # whiten both polarizations
            hp[self._kmin[det]:kmax] *= self._weight[det][slc]
            hc[self._kmin[det]:kmax] *= self._weight[det][slc]

            # h = fp * hp + hc * hc
            # <h, d> = fp * <hp,d> + fc * <hc,d>
            # the inner products
            cplx_hpd = self._whitened_data[det][slc].inner(hp[slc])  # <hp, d>
            cplx_hcd = self._whitened_data[det][slc].inner(hc[slc])  # <hc, d>

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
            lr += cplx_hd.real - 0.5 * hh

        lr_total = special.logsumexp(lr) - numpy.log(len(self.pol))

        # store the maxl polarization
        idx = lr.argmax()
        setattr(self._current_stats, 'maxl_polarization', self.pol[idx])

        # just store the maxl optimal snrsq
        for det in wfs:
            p = '{}_optimal_snrsq'.format(det)
            setattr(self._current_stats, p,
                    getattr(self._current_stats, p)[idx])

        return float(lr_total)
