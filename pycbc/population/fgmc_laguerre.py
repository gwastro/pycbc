# Copyright (C) 2016 Jolien Creighton
#           (C) 2021 Jolien Creighton & Thomas Dent
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.

"""
Based ultimately on code used for O1 rate calculations, see
https://git.ligo.org/RatesAndPopulations/lvc-rates-and-pop/-/blob/master/bin/O1_scripts/lvc_rates_calc_posterior
and technical documentation at https://dcc.ligo.org/LIGO-T1700029/public
"""

import numpy
import scipy.stats as sst
import scipy.special as ssp
import scipy.integrate as sig
import scipy.optimize as sop


class augmented_rv_continuous(sst.rv_continuous):

    def __init__(self, unit='dimensionless', texunit=r'\mbox{dimensionless}',
                 texsymb=r'x', **kwargs):
        '''
        Parameters
        ----------
        unit : string, optional
            units of independent variable
        texunit : string, optional
            units of independent variable, in tex format
        texsymb : string, optional
            symbol of independent variable, in tex format
        '''

        super(augmented_rv_continuous, self).__init__(**kwargs)
        self._hpd_interval_vec = numpy.vectorize(self._hpd_interval_scalar)
        self.unit = unit
        self.texunit = texunit
        self.texsymb = texsymb

    def _hpd_interval_scalar(self, alpha):

        def width(a):
            return self.ppf(alpha + self.cdf(a)) - a

        a = self.ppf(1e-6)  # a is displaced slightly from 0
        b = self.ppf(alpha)
        if self.pdf(a) >= self.pdf(b):  # upper limit
            return self.a, b
        a = sop.fminbound(width, a, self.ppf(1.0 - alpha))
        b = self.ppf(alpha + self.cdf(a))
        return a, b

    def hpd_interval(self, alpha):
        '''
        Confidence interval of highest probability density.

        Parameters
        ----------
        alpha : array_like of float
            Probability that an rv will be drawn from the returned range.
            Each value should be in the range [0, 1].

        Returns
        -------
        a, b : ndarray of float
            end-points of range that contain ``100 * alpha %`` of the rv's
            possible values.
        '''
        if isinstance(alpha, (float, numpy.number)):
            a, b = self._hpd_interval_scalar(alpha)
        else:
            a, b = self._hpd_interval_vec(alpha)
        return a, b


class count_posterior(augmented_rv_continuous):
    '''
    Count posterior distribution.
    '''

    def __init__(self, logbf, laguerre_n, Lambda0, prior=-0.5,
                 name='count posterior', unit='signals/experiment',
                 texunit=r'\mathrm{signals}/\mathrm{experiment}',
                 texsymb=r'\Lambda_1'):
        '''
        Parameters
        ----------
        logbf : array_like
            logs of normalized foreground over background pdf ratios of events
        laguerre_n: int
            degree of generalized Laguerre polynomial for quadrature formula
        Lambda0 : float
            background rate (default=len(bayesfac))
        prior : float or count_posterior, optional
            prior distribution power law of improper prior if float
            or count posterior distribution if count_posterior
            (default=-0.5: Jeffreys prior)
        '''
        super(count_posterior, self).__init__(a=0.0, b=numpy.inf, name=name,
                                              unit=unit, texunit=texunit,
                                              texsymb=texsymb)
        self.Lambda0 = Lambda0
        # weighted Bayes factor
        self.k = numpy.exp(numpy.array(logbf)) / self.Lambda0

        # power-law priors
        self.alpha = prior
        if prior == 0:
            self.prior = lambda x: 1.0
        elif prior > 0:
            self.prior = lambda x: x ** prior
        else:
            # regularize at x = 0
            self.prior = lambda x: (x + self.xtol) ** prior

        # pre-compute Gaussian-Generalized-Laguerre quadrature
        # abscissas and weights, along with pdf at these abscissas
        self.x, w = ssp.la_roots(laguerre_n, self.alpha)
        self.p = numpy.array([ww * numpy.prod(1.0 + self.k * xx)
                              for xx, ww in zip(self.x, w)])
        self.norm = 1.0 / sum(self.p)
        self.p *= self.norm

    def _pdf(self, x):
        # discourage underflows by evaluating ln L and using ln(1+x) function
        logL = -x + numpy.sum(numpy.log1p(self.k * x))
        P = numpy.exp(logL) * self.prior(x)
        return self.norm * P

    def _cdf(self, x):
        return sig.quad(self._pdf, 0.0, x)

    def expect(self, func):
        '''
        Calculate expected value of a function with respect to the
            distribution.

        The expected value of a function ``f(x)`` with respect to a
        distribution ``dist`` is defined as::

            E[x] = Integral(f(x) * dist.pdf(x))

        Parameters
        ----------
        func : callable
            Function for which integral is calculated. Takes only one argument.

        Returns
        -------
        expect : float
            The calculated expected value.
        '''
        # FIXME: not as feature rich as the expect method this overrides
        return sum(pp * func(xx) for xx, pp in zip(self.x, self.p))

    def _munp(self, n):
        return self.expect(lambda x: x**n)

    def p_bg(self, logbf):
        '''
        Calculate the false alarm probabilities of the events.

        Parameters
        ----------
        logbf : array_like
            Logs of foreground over background probability ratios of events.
        '''
        # get weighted bayes factor
        k = numpy.exp(numpy.asarray(logbf)) / self.Lambda0

        P0 = numpy.dot(1./(1. + numpy.outer(k, self.x)), self.p)
        if isinstance(k, (float, int, numpy.number)):
            return P0.item()
        if isinstance(k, numpy.ndarray) and k.ndim == 0:
            return P0.item()
        # except in special cases above, return array of values
        return P0

__all__ = ['augmented_rv_continuous', 'count_posterior']
