# Copyright (C) 2016  Collin Capano
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


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
"""
This modules provides classes and functions for evaluating the log likelihood
for parameter estimation.
"""

from pycbc import filter
from pycbc.types import Array
import prior as pyprior
import numpy

class _BaseLikelihoodEvaluator:
    """Base container class for generating waveforms, storing the data, and
    computing log likelihoods. The likelihood function defined here does
    nothing. The likelihood function either needs to be monkey patched, or
    other classes inherit from this class and define their own
    likelihood functions.

    Parameters
    ----------
    waveform_generator : generator class
        A generator class that creates waveforms. This must have a generate
        function which takes a set of parameter values as arguments, and a
        detectors attribute which is a dictionary of detectors keyed by their
        names.
    data : dict
        A dictionary of data, in which the keys are the detector names and the
        values are the data (assumed to be unwhitened). The list of keys must
        match the waveform generator's detectors keys.
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
        either a float or an array. If None, 4*data.values()[0].delta_f will be
        used.
    prior : callable
        A callable class or function that computes the prior.
    """

    def __init__(self, waveform_generator, data, f_lower, psds=None,
            f_upper=None, norm=None, prior=None):
        self._waveform_generator = waveform_generator
        self._data = data
        # check that the data and waveform generator have the same detectors
        if sorted(waveform_generator.detectors.keys()) != \
                sorted(self._data.keys()):
            raise ValueError("waveform generator's detectors (%s) " %(
                ','.join(sorted(waveform_generator.detector_names))) +
                "does not match data (%s)" %(
                ','.join(sorted(self._data.keys()))))
        # check that the data sets all have the same lengths
        dlens = numpy.array([len(d) for d in data.values()])
        if not all(dlens == dlens[0]):
            raise ValueError("all data must be of the same length")
        N = dlens[0]
        # we'll use the first data set for setting values
        d = data.values()[0]
        # figure out the kmin, kmax to use
        kmin, kmax = filter.get_cutoff_indices(f_lower, f_upper, d.delta_f,
            (N-1)*2)
        self._kmin = kmin
        self._kmax = kmax
        if norm is None:
            norm = 4*d.delta_f
        # we'll store the weight to apply to the inner product
        if psds is None:
            w = norm*Array(numpy.ones(N))
            # FIXME: use the following when we've switched to 2.7
            #self._weight = {det: w for det in data} 
            self._weight = dict([(det, w) for det in data])
        else:
            # temporarily suppress numpy divide by 0 warning
            numpy.seterr(divide='ignore')
            # FIXME: use the following when we've switched to 2.7
            #self._weight = {det: norm/psds[det] for det in data}
            self._weight = dict([(det, norm/psds[det]) for det in data])
            numpy.seterr(divide='warn')
        # compute <d, d>
        # FIXME: use the following when we've switched to 2.7
        #self._dd = {det:
        #    d[kmin:kmax].inner(d[kmin:kmax]*self._weight[det][kmin:kmax]).real/2.
         #   for det,d in self._data.items()}
        self._dd = dict([(det,
            d[kmin:kmax].inner(d[kmin:kmax]*self._weight[det][kmin:kmax]).real/2.)
            for det,d in self._data.items()])
        # store prior
        if prior:
            self._prior = prior
        else:
            self._prior = pyprior.no_prior

    @property
    def waveform_generator(self):
        """Returns the waveform generator that was set."""
        return self._waveform_generator

    @property
    def data(self):
        """Returns the data that was set."""
        return self._data

    def prior(self, params):
        """This function should return the prior of the given params.
        """
        return self._prior(params)

    def loglikelihood(self, params):
        """This function should return the log likelihood of the given params.
        """
        raise ValueError("Likelihood function not set.")

    def __call__(self, params):
        return self.loglikelihood(params)



class GaussianLikelihood(_BaseLikelihoodEvaluator):
    r"""Computes the log likelihood for the given parameters using:

    .. math:: \log \mathcal{L} \propto -\frac{1}{2}\sum_{i} \left<\mathbf{h}_i(\mathbf{\theta}) - \mathbf{s}_i | \mathbf{h}_i(\mathbf{\theta}) - \mathbf{s}_i \right>

    where :math:`\mathbf{h}_i(\mathbf{\theta})` and :math:`\mathbf{s}_i`
    are the model waveform vector with parameters :math:`\mathbf{\theta}`
    and the data, respectively, in the :math:`i`th detector (the sum is
    over the detectors), and the inner product is given by:

    .. math:: \left<a | b\right> = 4\Re \int \tilde{a}(f) \tilde{b}(f) \mathrm{d}f.

    No analytic marginalization is done, and the data is assumed to not be
    whitened.

    For details on initialization parameters, see _BaseLikelihoodEvaluator.

    Examples
    --------
    Create a signal, and compute the log likelihood for a template with the same
    parameters at the same time (the log likelihood should be zero in this case):
    >>> from pycbc import psd as pypsd, inference, waveform
    >>> seglen = 4
    >>> m1, m2, s1z, s2z, tsig, ra, dec, pol = 38.6, 29.3, 0.33, -0.94, 3.1, 1.37, -1.26, 2.76
    >>> generator = waveform.FDomainDetFrameGenerator(waveform.FDomainCBCGenerator, variable_args=['tc'], detectors=['H1', 'L1'], delta_f=1./seglen, f_lower=20., approximant='SEOBNRv2_ROM_DoubleSpin', mass1=m1, mass2=m2, spin1z=s1z, spin2z=s2z, ra=ra, dec=dec, polarization=pol)
    >>> signal = generator.generate(tsig)
    >>> psd = pypsd.aLIGOZeroDetHighPower(seglen*2048/2+1, 1./seglen, 20.)
    >>> psds = {'H1': psd, 'L1': psd}
    >>> likelihood_eval = inference.GaussianLikelihood(generator, signal, 20., psds=psds)
    >>> likelihood_eval.loglikelihood([tsig])
    ArrayWithAligned(0.0)

    Using the same likelihood evaluator, evaluate the log likelihood at several
    points in time, check that the max is at tsig, and plot:
    >>> times = numpy.arange(seglen*2048)/2048.
    >>> lls = numpy.array([likelihood_eval.loglikelihood([t]) for t in times])
    >>> times[lls.argmax()]
    3.10009765625
    >>> fig = pyplot.figure(); ax = fig.add_subplot(111)
    >>> ax.plot(times, lls)
    [<matplotlib.lines.Line2D at 0x12780ff90>]
    >>> fig.show()
    """

    def loglikelihood(self, params):
        """Computes the log-likelihood at the given point in parameter space.

        Parameters
        ----------
        params: array-like
            An array of numerical values to pass to the waveform generator.

        Returns
        -------
        float
            The value of the log-likelhood evaluated at the given point in
            parameter space.
        """
        # get prior
        prior = self.prior(params)
        # prior will return -numpy.inf if params are invalid
        if prior == -numpy.inf:
            return -numpy.inf
        hs = self._waveform_generator.generate(*params)
        # the kmax of the waveforms may be different than internal kmax
        kmax = min(len(hs.values()[0]), self._kmax)
        return prior + sum([
            # <h, d>
            self.data[det][self._kmin:kmax].inner(
                h[self._kmin:kmax]*self._weight[det][self._kmin:kmax]).real
            # - <h, h>/2.
            - h[self._kmin:kmax].inner(
                h[self._kmin:kmax]*self._weight[det][self._kmin:kmax]).real/2.
            # - <d, d>/2.
            - self._dd[det]
            for det,h in hs.items()])


likelihood_evaluators = {'gaussian': GaussianLikelihood}

__all__ = ['GaussianLikelihood', 'likelihood_evaluators']
