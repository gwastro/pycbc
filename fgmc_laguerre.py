"""
Based ultimately on code used for O1 rate calculations, see
https://git.ligo.org/RatesAndPopulations/lvc-rates-and-pop/-/blob/master/bin/O1_scripts/lvc_rates_calc_posterior
and technical documentation at
https://dcc.ligo.org/LIGO-T1700029/public
"""

import numpy
import scipy.stats as sst
import scipy.special as ssp
import scipy.integrate as sig
import scipy.optimize as sop
import json

from matplotlib import figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

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

    def _munp(self, n, *args):
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
        elif isinstance(k, numpy.ndarray) and k.ndim == 0:
            return P0.item()
        else:
            return P0


def plotodds(rankstats, p_b):

    fig = figure.Figure()
    FigureCanvas(fig)
    ax = fig.gca()
    ax.loglog()
    ax.grid(True)
    ax.plot(rankstats, (1.0 - p_b) / p_b, 'k.')
    ax.plot([rankstats.min(), rankstats.max()], [1.0, 1.0], 'c--')
    ax.set_title(r'Foreground/Background Odds')
    ax.set_xlabel(r'ranking statistic')
    ax.set_ylabel(r'$P_1/P_0$')
    ax.set_xlim(0.99 * rankstats.min(), 1.2 * rankstats.max())
    return fig


def plotpbg(rankstats, p_b):

    fig = figure.Figure()
    FigureCanvas(fig)
    ax = fig.gca()
    ax.loglog()
    ax.grid(True)
    ax.plot(rankstats, p_b, 'k.')
    ax.set_title(r'Probability of background origin')
    ax.set_xlabel(r'ranking statistic')
    ax.set_ylabel(r'$P_0$')
    ax.set_xlim(0.99 * rankstats.min(), 1.2 * rankstats.max())
    return fig


def plotoddsifar(ifar, p_b):

    fig = figure.Figure()
    FigureCanvas(fig)
    ax = fig.gca()
    ax.loglog()
    ax.grid(True)
    ax.plot(ifar, (1.0 - p_b) / p_b, 'k.')
    ax.plot([ifar.min(), ifar.max()], [1.0, 1.0], 'c--')
    ax.set_title(r'Foreground/Background Odds')
    ax.set_xlabel(r'IFAR')
    ax.set_ylabel(r'$P_1/P_0$')
    ax.set_xlim(0.9 * ifar.min(), 1.1 * ifar.max())
    return fig


def plotfdr(p_b, ntop):

    # get smallest N p_terr values
    p_b = numpy.sort(p_b)[:ntop]
    # cumulative probable noise/signal counts
    cum_false = p_b.cumsum()
    cum_true = (1. - p_b).cumsum()
    fig = figure.Figure()
    FigureCanvas(fig)
    ax = fig.gca()
    ax.semilogy()
    ax.grid(True)
    ax.plot(p_b, cum_false / cum_true, 'b+')
    ax.plot(p_b, 1. / (numpy.arange(len(p_b)) + 1), 'c--', label=r'1 noise event')
    leg = ax.legend()
    ax.set_xlabel(r'$p_{\rm terr}$')
    ax.set_ylabel(r'Cumulative $p_{\rm terr}$ / Cumulative $p_{\rm astro}$')
    ax.set_xlim(0., 1.05 * p_b.max())
    return fig


def finalize_plot(fig, args, extensions, name, pltype, tag):
    # Helper function
    for extn in extensions:
        filename = args.pldir + '_'.join(name.split()) + '_' + pltype + tag + extn
        if args.verbose:
            print('writing %s ...' % filename)
        fig.savefig(filename)


def odds_summary(args, rankstats, ifars, p_b, ntop, times=None, mchirps=None,
                 name='events', plot_extensions=None):

    print('\nSummary of Top %i %s' % (ntop, name.title()))

    # do sort in reverse order
    statsort = numpy.argsort(1. / numpy.array(rankstats))
    topn = statsort[:ntop]  # indices giving top n
    topgps = []
    topstat = []
    topifar = []
    toppastro = []
    for n, i in enumerate(topn):
        gps = times[i] if times is not None else ''
        ifar = ifars[i]
        stat = rankstats[i]
        mchirpstring = 'mchirp %.3F' % mchirps[i] if mchirps is not None else ''
        topgps.append(gps)
        topstat.append(stat)
        topifar.append(ifar)
        print('#%d event:' % (n + 1), str(gps), mchirpstring)
        print('    rankstat = %-8.3f' % stat)
        print('    IFAR     = %.2f' % ifar)
        print('    odds     = %g' % ((1. - p_b[i]) / p_b[i]))
        toppastro.append(1. - p_b[i])

    if args.p_astro_txt is not None:
        numpy.savetxt(args.p_astro_txt,
                      numpy.column_stack((topgps, topstat, topifar, toppastro)),
                      fmt=['%.3F', '%.2F', '%.2F', '%.5F'],
                      delimiter=',',
                      header='GPS seconds, stat, IFAR/yr, p_astro')

    if hasattr(args, 'json_tag') and args.json_tag is not None:
        # save to catalog-style files
        def dump_json(gps, p_a, p_b):
            jfile = args.plot_dir + 'H1L1V1-PYCBC_%s-%s-1.json' % \
                      (args.json_tag, str(int(gps)))  # truncate to integer GPS
            with open(jfile, 'w') as jf:
                json.dump({'Astro': p_a, 'Terrestrial': p_b}, jf)
        if hasattr(args, 'json_min_ifar') and args.json_min_ifar is not None:
            for g, ifar, pt in zip(times, ifars, p_b):
                if ifar < args.json_min_ifar:
                    continue
                dump_json(g, 1. - pt, pt)
        else:
            for g, pa in zip(topgps, toppastro):
                dump_json(g, pa, 1. - pa)

    if plot_extensions is not None:
        plottag = args.plot_tag or ''
        if plottag is not '':
            plottag = '_' + plottag
        fig = plotodds(rankstats, p_b)
        finalize_plot(fig, args, plot_extensions, name, 'odds', plottag)
        fig = plotpbg(rankstats, p_b)
        finalize_plot(fig, args, plot_extensions, name, 'pbg', plottag)
        fig = plotoddsifar(ifars, p_b)
        finalize_plot(fig, args, plot_extensions, name, 'ifarodds', plottag)
        fig = plotfdr(p_b, ntop)
        finalize_plot(fig, args, plot_extensions, name, 'fdr', plottag)        


def plotdist(rv, plot_lim=None, middle=None, credible_intervals=None, style='linear'):

    fig = figure.Figure()
    FigureCanvas(fig)
    ax = fig.gca()

    name = rv.name if hasattr(rv, 'name') else None
    symb = rv.texsymb if hasattr(rv, 'texsymb') else r'x'
    unit = rv.texunit if hasattr(rv, 'texunit') else None

    xlabel = r'$' + symb + '$'
    if unit is not None:
        xlabel += r' ($' + unit + r'$)'

    a, b = rv.interval(0.9999)

    if style == 'loglog':
        ax.loglog()
        ylabel = r'$p(' + symb + r')$'
        space = lambda a, b: numpy.logspace(numpy.log10(a), numpy.log10(b), 100)
        func = numpy.vectorize(lambda x: rv.pdf(x))
        xmin = a
        ymin = rv.pdf(b)

    elif style == 'semilogx':
        ax.semilogx()
        ylabel = r'$' + symb + r'\,p(' + symb + r')$'
        space = lambda a, b: numpy.logspace(numpy.log10(a), numpy.log10(b), 100)
        func = numpy.vectorize(lambda x: x * rv.pdf(x))
        xmin = a
        ymin = 0.0

    else:  # linear
        ax.yaxis.set_ticklabels([])
        ylabel = r'$p(' + symb + r')$'
        space = lambda a, b: numpy.linspace(a, b, 100)
        func = numpy.vectorize(lambda x: rv.pdf(x))
        xmin = 0.0
        ymin = 0.0

    x = space(a, b)
    y = func(x)

    ax.plot(x, y, color='k', linestyle='-')
    if plot_lim is not None:
        xmin, xmax = plot_lim
        ax.set_xlim(xmin=xmin, xmax=xmax)
    else:
        ax.set_xlim(xmin=xmin)
    ax.set_ylim(ymin=ymin)

    if y.max() < 2 and y.max() > 1:
        ax.set_ylim(ymax = 2.)

    if name is not None:
        ax.set_title(name.title())
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if middle is not None:
        ax.plot([middle, middle], [ymin, func(middle)], 'k--')

    if credible_intervals is not None:
        # alpha : density of fill shading
        alpha = 1.0 / (1.0 + len(credible_intervals))
        for lo, hi in credible_intervals.values():
            lo = max(a, lo)
            hi = min(b, hi)
            x = space(lo, hi)
            y = func(x)
            ax.fill_between(x, y, ymin, color='k', alpha=alpha)

    return fig


def dist_summary(args, rv, plot_styles=['linear', 'loglog', 'semilogx'],
                 plot_extensions=None, middle=None, credible_intervals=None):

    name = rv.name if hasattr(rv, 'name') else 'posterior'
    unit = rv.unit if hasattr(rv, 'unit') else ''
    median = rv.median()
    mode = rv.mode() if hasattr(rv, 'mode') else None

    print('Summary of ' + name.title())
    print('mean   =', rv.mean(), unit)
    print('median =', median, unit)
    if mode is not None:
        print('mode   =', mode, unit)
    print('stddev =', rv.std(), unit)

    if credible_intervals is not None and len(credible_intervals) > 0:
        print('equal-tailed credible intervals:')
        equal_tailed_credible_intervals = {}
        for cred in credible_intervals:
            lo, hi = rv.interval(cred)
            equal_tailed_credible_intervals[cred] = (lo, hi)
            print('%g%%' % (cred * 100), 'credible interval =', '[%g, %g]' %
                                         (lo, hi), unit)

        if hasattr(rv, 'hpd_interval'):
            print('highest probability density credible intervals:')
            hpd_credible_intervals = {}
            for cred in credible_intervals:
                hpdlo, hpdhi = rv.hpd_interval(cred)
                hpd_credible_intervals[cred] = (hpdlo, hpdhi)
                print('%g%%' % (cred * 100), 'credible interval =', '[%g, %g]' %
                                             (hpdlo, hpdhi), unit)
        else:
            hpd_credible_intervals = None

    if len(credible_intervals) == 0:
        credible_intervals = None
        intervals = None

    if middle == 'mode' and mode is not None:
        middle = mode
        if credible_intervals is not None:
            # use hpd intervals with mode
            intervals = hpd_credible_intervals
    else:
        middle = median
        if credible_intervals is not None:
            # use equal tailed intervals with median
            intervals = equal_tailed_credible_intervals

    # plot distributions
    if plot_extensions is not None:
        plottag = args.plot_tag or ''
        if plottag is not '':
            plottag = '_' + plottag
        for style in plot_styles:
            fig = plotdist(rv, plot_lim=args.plot_limits, middle=middle,
                           credible_intervals=intervals, style=style)
            finalize_plot(fig, args, plot_extensions, name, style, plottag)

    if credible_intervals is not None and len(credible_intervals) == 1:
        return median, lo - median, hi - median

