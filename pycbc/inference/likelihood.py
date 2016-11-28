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
import numpy

class _NoPrior(object):
    """Dummy class to just return 0 if no prior is provided in a
    likelihood generator.
    """
    @staticmethod
    def apply_boundary_conditions(params):
        return params

    def __call__(self, params):
        return 0.


def snr_from_loglr(loglr):
    """Returns SNR computed from the given log likelihood ratio(s). This is
    defined as `sqrt(2*loglr)`.If the log likelihood ratio is < 0, returns 0.

    Parameters
    ----------
    loglr : array or float
        The log likelihood ratio(s) to evaluate.

    Returns
    -------
    array or float
        The SNRs computed from the log likelihood ratios.
    """
    singleval = isinstance(loglr, float)
    if singleval:
        loglr = numpy.array([loglr])
    # temporarily quiet sqrt(-1) warnings
    numpysettings = numpy.seterr(invalid='ignore')
    snrs = numpy.sqrt(2*loglr)
    numpy.seterr(**numpysettings)
    snrs[numpy.isnan(snrs)] = 0.
    if singleval:
        snrs = snrs[0]
    return snrs

class _BaseLikelihoodEvaluator(object):
    r"""Base container class for generating waveforms, storing the data, and
    computing posteriors.

    The nomenclature used by this class and those that inherit from it is as
    follows: Given some model parameters :math:`\Theta` and some data
    :math:`d` with noise model :math:`n`, we define:

     * the *likelihood function*: :math:`p(d|\Theta)`

     * the *noise likelihood*: :math:`p(d|n)`

     * the *likelihood ratio*: :math:`\mathcal{L}(\Theta) = \frac{p(d|\Theta)}{p(d|n)}`

     * the *prior*: :math:`p(\Theta)`

     * the *posterior*: :math:`p(\Theta|d) \propto p(d|\Theta)p(\Theta)`

     * the *prior-weighted likelihood ratio*: :math:`\hat{\mathcal{L}}(\Theta) = \frac{p(d|\Theta)p(\Theta)}{p(d|n)}
   
     * the *SNR*: :math:`\rho(\Theta) = \sqrt{2\log\mathcal{L}(\Theta)}`; for
       two detectors, this is approximately the same quantity as the coincident
       SNR used in the CBC search.
   
    .. note::

        Although the posterior probability is only proportional to
        :math:`p(d|\Theta)p(\Theta)`, here we refer to this quantity as the
        posterior. Also note that for a given noise model, the prior-weighted
        likelihood ratio is proportional to the posterior, and so the two can
        usually be swapped for each other.

    When performing parameter estimation we work with the log of these values
    since we are mostly concerned with their values around the maxima. If
    we have multiple detectors, each with data :math:`d_i`, then these values
    simply sum over the detectors. For example, the log likelihood ratio is:

    .. math::
        \log \mathcal{L}(\Theta) = \sum_i \left[\log p(\Theta|d_i) - \log p(n|d_i)\right]
   
    This class provides boiler-plate methods and attributes for evaluating the
    log likelihood ratio, log prior, and log likelihood. This class
    makes no assumption about the detectors' noise model :math:`n`. As such,
    the methods for computing these values raise `NotImplementedError`s. These
    functions need to be monkey patched, or other classes that inherit from
    this class need to define their own functions.

    Instances of this class can be called like a function. The default is for
    this class to call its `logposterior` function, but this can be changed by
    with the `set_callfunc` method.

    Parameters
    ----------
    waveform_generator : generator class
        A generator class that creates waveforms. This must have a generate
        function which takes a set of parameter values as arguments, a
        detectors attribute which is a dictionary of detectors keyed by their
        names, and an epoch which specifies the start time of the generated
        waveform.
    data : dict
        A dictionary of data, in which the keys are the detector names and the
        values are the data (assumed to be unwhitened). The list of keys must
        match the waveform generator's detectors keys, and the epoch of every
        data set must be the same as the waveform generator's epoch.
    prior : callable
        A callable class or function that computes the log of the prior. If
        None provided, will use `_noprior`, which returns 0 for all parameter
        values.

    Attributes
    ----------
    waveform_generator : dict
        The waveform generator that the class was initialized with.
    data : dict
        The data that the class was initialized with.
    lognl : {None, float}
        The log of the noise likelihood summed over the number of detectors.
    return_meta : {True, bool}
        If True, `logposterior` and `logplr` will return the value of the
        prior and the loglikelihood ratio, along with the posterior/plr.

    Methods
    -------
    prior :
        A function that returns the log of the prior given a list of
        parameters.
    loglikelihood :
        A function that returns the log of the likelihood function of a given
        list of parameters.
    logposterior :
        A function that returns the log of the posterior of a given list of
        parameters.
    loglr :
        A function that returns the log of the likelihood ratio of a given list
        of parameters.
    logplr :
        A function that returns the log of the prior-weighted likelihood ratio
        of a given list of parameters.
    snr :
        A function that returns the square root of twice the log likelihood
        ratio. If the log likelihood ratio is < 0, will return 0.
    set_callfunc :
        Set the function to use when the class is called as a function.
    """
    name = None

    def __init__(self, waveform_generator, data, prior=None,
                 return_meta=True):
        self._waveform_generator = waveform_generator
        # we'll store a copy of the data which we'll later whiten in place
        self._data = dict([[ifo, 1*data[ifo]] for ifo in data])
        # check that the data and waveform generator have the same detectors
        if sorted(waveform_generator.detectors.keys()) != \
                sorted(self._data.keys()):
            raise ValueError("waveform generator's detectors (%s) " %(
                ','.join(sorted(waveform_generator.detector_names))) +
                "does not match data (%s)" %(
                ','.join(sorted(self._data.keys()))))
        # check that the data and waveform generator have the same epoch
        if any(waveform_generator.epoch != d.epoch for d in self._data.values()):
            raise ValueError("waveform generator does not have the same epoch "
                "as all of the data sets.")
        # check that the data sets all have the same lengths
        dlens = numpy.array([len(d) for d in data.values()])
        if not all(dlens == dlens[0]):
            raise ValueError("all data must be of the same length")
        # store prior
        if prior is None:
            self._prior = _NoPrior()
        else:
            # check that the variable args of the prior evaluator is the same
            # as the waveform generator
            if prior.variable_args != self._waveform_generator.variable_args:
                raise ValueError("variable args of prior and waveform "
                    "generator do not match")
            self._prior = prior
        # initialize the log nl to 0
        self._lognl = None
        self.return_meta = return_meta

    @property
    def waveform_generator(self):
        """Returns the waveform generator that was set."""
        return self._waveform_generator

    @property
    def data(self):
        """Returns the data that was set."""
        return self._data

    @property
    def lognl(self):
        """Returns the log of the noise likelihood."""
        return self._lognl

    def set_lognl(self, lognl):
        """Set the value of the log noise likelihood."""
        self._lognl = lognl

    def prior(self, params):
        """This function should return the prior of the given params.
        """
        return self._prior(params)

    def loglikelihood(self, params):
        """Returns the natural log of the likelihood function.
        """
        raise NotImplementedError("Likelihood function not set.")

    def loglr(self, params):
        """Returns the natural log of the likelihood ratio.
        """
        raise NotImplementedError("Likelihood ratio function not set.")


    # the names and order of data returned by _formatreturn when
    # return_metadata is True
    metadata_fields = ["prior", "loglr"]

    def _formatreturn(self, val, prior=None, loglr=None):
        """Adds the prior to the return value if return_meta is True.
        Otherwise, just returns the value.

        Parameters
        ----------
        val : float
            The value to return.
        prior : {None, float}
            The value of the prior.
        loglr : {None, float}
            The value of the log likelihood-ratio.

        Returns
        -------
        val : float
            The given value to return.
        *If return_meta is True:*
        metadata : (prior, loglr)
            A tuple of the prior and log likelihood ratio.
        """
        if self.return_meta:
            return val, (prior, loglr)
        else:
            return val

    def logplr(self, params):
        """Returns the log of the prior-weighted likelihood ratio.
        """
        # if the prior returns -inf, just return
        logp = self._prior(params)
        if logp == -numpy.inf:
            return self._formatreturn(logp, prior=logp)
        llr = self.loglr(params)
        return self._formatreturn(llr + logp, prior=logp, loglr=llr)

    def logposterior(self, params):
        """Returns the log of the posterior of the given params.
        """
        # if the prior returns -inf, just return
        logp = self._prior(params)
        if logp == -numpy.inf:
            return self._formatreturn(logp, prior=logp)
        ll = self.loglikelihood(params)
        return self._formatreturn(ll + logp, prior=logp, loglr=ll-self._lognl)

    def snr(self, params):
        """Returns the "SNR" of the given params. This will return
        imaginary values if the log likelihood ratio is < 0.
        """
        return snr_from_loglr(self.loglr(params))

    _callfunc = logposterior

    @classmethod
    def set_callfunc(cls, funcname):
        """Sets the function used when the class is called as a function.

        Parameters
        ----------
        funcname : str
            The name of the function to use; must be the name of an instance
            method.
        """
        cls._callfunc = getattr(cls, funcname)

    def __call__(self, params):
        # apply any boundary conditions to the parameters before
        # generating/evaluating
        return self._callfunc(self._prior.apply_boundary_conditions(params))



class GaussianLikelihood(_BaseLikelihoodEvaluator):
    r"""Computes log likelihoods assuming the detectors' noise is Gaussian.

    With Gaussian noise the log likelihood functions for signal
    :math:`\log p(d|\Theta)` and for noise :math:`log p(d|n)` are given by:

    .. math::

        \log p(d|\Theta) = -\frac{1}{2} \sum_i \left<h_i(\Theta) - d_i | h_i(\Theta - d_i\right>

        \log p(d|n) = -\frac{1}{2} \sum_i \left<d_i | d_i\right>

    where the sum is over the number of detectors, :math:`d_i` is the data in
    each detector, and :math:`h_i(\Theta)` is the model signal in each
    detector. The inner product is given by:

    .. math::

        \left<a | b\right> = 4\Re \int_{0}^{\infty} \frac{\tilde{a}(f) \tilde{b}(f)}{S_n(f)} \mathrm{d}f,

    where :math:`S_n(f)` is the PSD in the given detector.
    
    Note that the log prior-weighted likelihood ratio has one less term
    than the log posterior, since the :math:`\left<d_i|d_i\right>` term cancels
    in the likelihood ratio:

    .. math::

        \log \hat{\mathcal{L}} = \log p(\Theta) + \sum_i \left\[\left<h_i(\Theta)|d_i\right> - \frac{1}{2} \left<h_i(\Theta)|h_i(\Theta)\right>\right]

    For this reason, by default this class returns `logplr` when called as a
    function instead of `logposterior`. This can be changed via the
    `set_callfunc` method.

    Upon initialization, the data is whitened using the given PSDs. If no PSDs
    are given the data and waveforms returned by the waveform generator are
    assumed to be whitened. The likelihood function of the noise,
    
    .. math::
    
        p(d|n) = \frac{1}{2} \sum_i \left<d_i|d_i\right>,

    is computed on initialization and stored as the `lognl` attribute.
    
    By default, the data is assumed to be equally sampled in frequency, but
    unequally sampled data can be supported by passing the appropriate
    normalization using the `norm` keyword argument.

    For more details on initialization parameters and definition of terms, see
    `_BaseLikelihoodEvaluator`.

    Parameters
    ----------
    waveform_generator : generator class
        A generator class that creates waveforms. This must have a generate
        function which takes a set of parameter values as arguments, a
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
        either a float or an array. If `None`, `4*data.values()[0].delta_f
        will be used.
    prior : callable
        A callable class or function that computes the prior.
    return_meta : {True, bool}
        If True, `logposterior` and `logplr` will return the value of the
        prior and the loglikelihood ratio, along with the posterior/plr.

    Examples
    --------
    Create a signal, and set up the likelihood evaluator on that signal:
    >>> seglen = 4
    >>> sample_rate = 2048
    >>> N = seglen*sample_rate/2+1
    >>> fmin = 30.
    >>> m1, m2, s1z, s2z, tsig, ra, dec, pol, dist = 38.6, 29.3, 0., 0., 3.1, 1.37, -1.26, 2.76, 3*500.
    >>> generator = waveform.FDomainDetFrameGenerator(waveform.FDomainCBCGenerator, 0., variable_args=['tc'], detectors=['H1', 'L1'], delta_f=1./seglen, f_lower=fmin, approximant='SEOBNRv2_ROM_DoubleSpin', mass1=m1, mass2=m2, spin1z=s1z, spin2z=s2z, ra=ra, dec=dec, polarization=pol, distance=dist)
    >>> signal = generator.generate(tsig)
    >>> psd = pypsd.aLIGOZeroDetHighPower(N, 1./seglen, 20.)
    >>> psds = {'H1': psd, 'L1': psd}
    >>> likelihood_eval = inference.GaussianLikelihood(generator, signal, fmin, psds=psds, return_meta=False)

    Now compute the log likelihood ratio and prior-weighted likelihood ratio;
    since we have not provided a prior, these should be equal to each other:
    >>> likelihood_eval.loglr([tsig]), likelihood_eval.logplr([tsig])
    (ArrayWithAligned(277.92945279883855), ArrayWithAligned(277.92945279883855))

    Compute the log likelihood and log posterior; since we have not
    provided a prior, these should both be equal to zero:
    >>> likelihood_eval.loglikelihood([tsig]), likelihood_eval.logposterior([tsig])
    (ArrayWithAligned(0.0), ArrayWithAligned(0.0))

    Compute the SNR; for this system and PSD, this should be approximately 24:
    >>> likelihood_eval.snr([tsig])
    ArrayWithAligned(23.576660187517593)

    Using the same likelihood evaluator, evaluate the log prior-weighted
    likelihood ratio at several points in time, check that the max is at tsig,
    and plot (note that we use the class as a function here, which defaults
    to calling `logplr`):
    >>> from matplotlib import pyplot
    >>> times = numpy.arange(seglen*sample_rate)/float(sample_rate)
    >>> lls = numpy.array([likelihood_eval([t]) for t in times])
    >>> times[lls.argmax()]
    3.10009765625
    >>> fig = pyplot.figure(); ax = fig.add_subplot(111)
    >>> ax.plot(times, lls)
    [<matplotlib.lines.Line2D at 0x1274b5c50>]
    >>> fig.show()

    Create a prior and use it (see prior module for more details):
    >>> from pycbc.inference import prior
    >>> uniform_prior = prior.Uniform(tc=(tsig-0.2,tsig+0.2))
    >>> prior_eval = prior.PriorEvaluator(['tc'], uniform_prior)
    >>> likelihood_eval = inference.GaussianLikelihood(generator, signal, 20., psds=psds, prior=prior_eval, return_meta=False)
    >>> likelihood_eval.logplr([tsig]), likelihood_eval.logposterior([tsig])
    (ArrayWithAligned(278.84574353071264), ArrayWithAligned(0.9162907318741418))
    """
    name = 'gaussian'

    def __init__(self, waveform_generator, data, f_lower, psds=None,
            f_upper=None, norm=None, prior=None, return_meta=True):
        # set up the boiler-plate attributes; note: we'll compute the
        # log evidence later
        super(GaussianLikelihood, self).__init__(waveform_generator, data,
            prior=prior, return_meta=return_meta)
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
            # FIXME: use the following when we've switched to 2.7
            #self._weight = {det: w for det in data} 
            self._weight = dict([(det, w) for det in data])
        else:
            # temporarily suppress numpy divide by 0 warning
            numpysettings = numpy.seterr(divide='ignore')
            # FIXME: use the following when we've switched to 2.7
            #self._weight = {det: Array(numpy.sqrt(norm/psds[det])) for det in data}
            self._weight = dict([(det, Array(numpy.sqrt(norm/psds[det]))) for det in data])
            numpy.seterr(**numpysettings)
        # whiten the data
        for det in self._data:
            self._data[det][kmin:kmax] *= self._weight[det][kmin:kmax]
        # compute the log likelihood function of the noise and save it
        self.set_lognl(-0.5*sum([
            d[kmin:kmax].inner(d[kmin:kmax]).real
            for d in self._data.values()]))
        # set default call function to logplor
        self.set_callfunc('logplr')

    @property
    def lognl(self):
        return self._lognl

    def loglr(self, params):
        r"""Computes the log likelihood ratio,
        
        .. math::
            
            \log \mathcal{L}(\Theta) = \sum_i \left<h_i(\Theta)|d_i\right> - \frac{1}{2}\left<h_i(\Theta)|h_i(\Theta)\right>,

        at the given point in parameter space :math:`\Theta`.

        Parameters
        ----------
        params: array-like
            An array of numerical values to pass to the waveform generator.

        Returns
        -------
        float
            The value of the log likelihood ratio evaluated at the given point.
        """
        lr = 0
        for det,h in self._waveform_generator.generate(*params).items():
            # the kmax of the waveforms may be different than internal kmax
            kmax = min(len(h), self._kmax)
            # whiten the waveform
            h[self._kmin:kmax] *= self._weight[det][self._kmin:kmax]
            lr += (
                # <h, d>
                self.data[det][self._kmin:kmax].inner(h[self._kmin:kmax]).real
                # - <h, h>/2.
                - 0.5*h[self._kmin:kmax].inner(h[self._kmin:kmax]).real
                )
        return lr

    def loglikelihood(self, params):
        r"""Computes the log likelihood of the paramaters,
        
        .. math::
        
        p(d|\Theta) = -\frac{1}{2}\sum_i \left<h_i(\Theta) - d_i | h_i(\Theta) - d_i\right>

        Parameters
        ----------
        params: array-like
            An array of numerical values to pass to the waveform generator.

        Returns
        -------
        float
            The value of the log likelihood evaluated at the given point.
        """
        # since the loglr has fewer terms, we'll call that, then just add
        # back the noise term that canceled in the log likelihood ratio
        return self.loglr(params) + self._lognl


    def logposterior(self, params):
        """Computes the log-posterior probability at the given point in
        parameter space.

        Parameters
        ----------
        params: array-like
            An array of numerical values to pass to the waveform generator.

        Returns
        -------
        float
            The value of the log-posterior evaluated at the given point in
            parameter space.
        metadata : tuple
            If `return_meta`, the prior and likelihood ratio as a tuple.
            Otherwise, just returns the log-posterior.
        """
        # since the logplr has fewer terms, we'll call that, then just add
        # back the noise term that canceled in the log likelihood ratio
        logplr = self.logplr(params)
        if self.return_meta:
            logplr, (pr, lr) = logplr
        else:
            pr = lr = None
        return self._formatreturn(logplr + self._lognl, prior=pr, loglr=lr)

likelihood_evaluators = {GaussianLikelihood.name: GaussianLikelihood}

__all__ = ['_BaseLikelihoodEvaluator', 'GaussianLikelihood',
           'likelihood_evaluators']
