# Copyright (C) 2021 Thomas Dent
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.

"""
A set of helper functions for evaluating event rates, densities etc.

See https://dcc.ligo.org/LIGO-T2100060/public for technical explanations
"""

from os.path import basename
import h5py
import bisect
from itertools import chain as it_chain, combinations as it_comb
import numpy as np

from pycbc import conversions as conv
from pycbc import events
from pycbc.events.coinc import mean_if_greater_than_zero as coinc_meanigz
from pycbc.events import triggers


def filter_bin_lo_hi(values, lo, hi):
    in_bin = np.sign((values - lo) * (hi - values))
    if np.any(in_bin == 0):
        raise RuntimeError('Edge case! Bin edges', lo, hi,
                           'value(s)', values[in_bin == 0])
    return in_bin == 1


def filter_tmplt_mchirp(bankf, lo_mchirp, hi_mchirp):
    with h5py.File(bankf) as bank:
        mchirp = conv.mchirp_from_mass1_mass2(bank['mass1'][:], bank['mass2'][:])
    # Boolean over template id
    return filter_bin_lo_hi(mchirp, lo_mchirp, hi_mchirp)


def read_full_data(fullf, rhomin, tmplt_filter=None):
    """Read the zero- and time-lagged triggers identified by a specific
       set of templates.

       Parameters
       ----------
       fullf:
              File that stores zerolag and slide triggers
       bankf:
              File with template mass/spin information
       rhomin: float
              Ranking statistic threshold
       tmplt_filter: array of Booleans
              Filter over the array of templates stored in bankf

       Returns
       -------
       dictionary
              containing foreground triggers and background information
    """
    with h5py.File(fullf, 'r') as full_data:
        # apply template filter
        tid_bkg = full_data['background_exc/template_id'][:]
        tid_fg = full_data['foreground/template_id'][:]
        bkg_inbin = tmplt_filter[tid_bkg]  # Boolean over bg events
        fg_inbin = tmplt_filter[tid_fg]  # Boolean over fg events
        zerolagstat = full_data['foreground/stat'][:][fg_inbin]
        zerolagifar = full_data['foreground/ifar'][:][fg_inbin]
        # arbitrarily choose time from one of the ifos
        zerolagtime = full_data['foreground/time1'][:][fg_inbin]

        cstat_back_exc = full_data['background_exc/stat'][:][bkg_inbin]
        dec_factors = full_data['background_exc/decimation_factor'][:][bkg_inbin]

    # filter on stat value
    above = zerolagstat > rhomin
    back_above = cstat_back_exc > rhomin
    return {'zerolagstat': zerolagstat[above],
            'zerolagifar': zerolagifar[above],
            'zerolagtime': zerolagtime[above],
            'dec_factors': dec_factors[back_above],
            'cstat_back_exc': cstat_back_exc[back_above],
            'file_name': fullf}


def read_full_data_mchirp(fullf, bankf, rhomin, mc_lo, mc_hi):
    tmp_filter = filter_tmplt_mchirp(bankf, mc_lo, mc_hi)
    return read_full_data(fullf, rhomin, tmp_filter)


def log_rho_bg(trigs, counts, bins):
    """
    trigs: zerolag event statistic values
    counts: background histogram
    bins: bin edges of the background histogram

    Returns:
    log of background PDF at the zerolag statistic values,
    fractional uncertainty due to Poisson count (set to 100% for empty bins)
    """
    trigs = np.atleast_1d(trigs)
    if len(trigs) == 0:  # corner case
        return np.array([]), np.array([])

    assert np.all(trigs >= np.min(bins)), "can't have triggers below bin lower limit"

    N = sum(counts)
    log_rhos = []
    fracerr = []

    # If any zerolag triggers that are louder than the max bin, put one
    # fictitious count in a bin that extends from the limits of the slide triggers
    # out to the loudest trigger.
    if np.any(trigs >= np.max(bins)):
        N = N + 1

    for t in trigs:
        if t >= np.max(bins):
            # For a trigger louder than the max bin, put one fictitious count in
            # a bin that extends from the limits of the slide triggers out to the
            # loudest trigger.  Fractional error is 100%
            log_rhos.append(-np.log(N) - np.log(np.max(trigs) - bins[-1]))
            fracerr.append(1.)
        else:
            i = bisect.bisect(bins, t) - 1
            # If there are no counts for a foreground trigger put a fictitious
            # count in the background bin
            if counts[i] == 0:
                counts[i] = 1
            log_rhos.append(np.log(counts[i]) - np.log(bins[i+1] - bins[i])
                            - np.log(N))
            fracerr.append(counts[i] ** -0.5)
    return np.array(log_rhos), np.array(fracerr)


def log_rho_fg_analytic(trigs, rhomin):
    # PDF of a rho^-4 distribution defined above the threshold rhomin
    return np.log(3.) + 3. * np.log(rhomin) - 4 * np.log(trigs)


def log_rho_fg(trigs, injstats, bins):
    """
    trigs: zerolag event statistic values
    injstats: injection event statistic values
    bins: histogram bins

    Returns:
    log of signal PDF at the zerolag statistic values,
    fractional uncertainty from Poisson count
    """
    trigs = np.atleast_1d(trigs)
    if len(trigs) == 0:  # corner case
        return np.array([])

    assert np.min(trigs) >= np.min(bins)
    # allow 'very loud' triggers
    bmax = np.max(bins)
    if np.max(trigs) >= bmax:
        print('Replacing stat values lying above highest bin')
        print(str(bmax))
        trigs = np.where(trigs >= bmax, bmax - 1e-6, trigs)
    assert np.max(trigs) < np.max(bins)  # check it worked

    counts, bins = np.histogram(injstats, bins)
    N = sum(counts)
    dens = counts / np.diff(bins) / float(N)
    fracerr = counts ** -0.5

    tinds = np.searchsorted(bins, trigs) - 1
    return np.log(dens[tinds]), fracerr[tinds]


def get_start_dur(path):
    fname = basename(path)  # remove directory path
    # file name is IFOS-DESCRIPTION-START-DURATION.type
    pieces = fname.split('.')[0].split('-')
    return pieces[2], pieces[3]


def in_coinc_time_incl(f, cstring, test_times):
    """ filter to all times where coincs of type given by cstring exist
    """
    in_time = np.zeros(len(test_times))
    starts = np.array(f['segments/%s/start' % cstring][:])
    ends = np.array(f['segments/%s/end' % cstring][:])
    idx_within_segment = events.indices_within_times(test_times, starts, ends)
    in_time[idx_within_segment] = np.ones_like(idx_within_segment)
    return in_time


# what to change for more/fewer ifos
_ifoset = ('H1', 'L1', 'V1')


# all combinations of ifos with length mincount or more
# each returned as a tuple in same order as ifos
def alltimes(ifos, mincount=1):
    assert mincount <= len(ifos)
    assert len(set(ifos)) == len(ifos)  # can't work with duplicates
    return it_chain.from_iterable(it_comb(ifos, r) for r in
                                  np.arange(mincount, len(ifos) + 1))


_alltimes = frozenset(alltimes(_ifoset, mincount=1))
_alltimestring = frozenset([''.join(t) for t in _alltimes])
_allctimes = frozenset(alltimes(_ifoset, mincount=2))


def ifos_from_combo(ct):
    # extract ifos in alphabetical order from a coinc time string
    return sorted([ct[i:i + 2] for i in range(0, len(ct), 2)])


def type_in_time(ct, cty):
    # returns True if the given coinc type can exist in the coinc time ct
    return all(i in ct for i in cty)


class EventRate(object):
    def __init__(self, args, coinc_times, coinc_types=None, bin_param='mchirp',
                 bin_lo=None, bin_hi=None):
        """
        coinc_times: iterable of strings indicating combinations of ifos operating
        coinc_types: list of strings indicating coinc event types to be considered

        """
        # allow for single-ifo time although not supported in pipeline yet
        if hasattr(args, 'min_ifos'):
            self.mincount = args.min_ifos
        else:
            self.mincount = 2
        if hasattr(args, 'network') and sorted(args.network) != list(_ifoset):
            self.ifos = sorted(args.network)
        else:
            self.ifos = _ifoset
        self.allctimes = frozenset(alltimes(self.ifos, mincount=self.mincount))
        self.allctimestring = frozenset([''.join(t) for t in self.allctimes])
        for ct in coinc_times:
            assert ct in list(self.allctimestring)
        self.ctimes = coinc_times

        if coinc_types is None:
            # all types possible during given times
            self.coinc_types = self.allctimestring
        else:
            # any coinc type must also be a time (?)
            for ty in coinc_types:
                assert ty in list(self.allctimes)
            self.coinc_types = frozenset([''.join(t) for t in coinc_types])
        if args.verbose:
            print('Using', self.coinc_types, 'coincs in',
                  self.allctimestring, 'times')
        self.args = args
        self.thr = self.args.stat_threshold
        self.bin_param = bin_param
        self.lo = bin_lo
        self.hi = bin_hi
        self.bank = None
        self.massspins = None
        self.tpars = None
        self.tmplt_filter = None
        self.in_bin = None
        self.incl_livetimes = {}
        self.livetimes = {}

    def add_bank(self, bank_file):
        self.bank = bank_file
        with h5py.File(self.bank, 'r') as b:
            tids = np.arange(len(b['mass1'][:]))
            # tuples of m1, m2, s1z, s2z in template id order
            self.massspins = triggers.get_mass_spin(b, tids)

    def filter_templates(self):
        """
        calculate array of Booleans in template id order to filter events
        """
        assert self.massspins is not None
        assert self.lo is not None
        assert self.hi is not None
        if self.args.verbose:
            print('Cutting on %s between %f - %f' %
                  (self.bin_param, self.lo, self.hi))
        self.tpars = triggers.get_param(self.bin_param, None, *self.massspins)
        self.in_bin = filter_bin_lo_hi(self.tpars, self.lo, self.hi)

    def make_bins(self, maxval, choice='bg'):
        # allow options to be strings describing bin formulae as well as floats?
        try:
            linbw = getattr(self.args, choice + '_bin_width')
            logbw = getattr(self.args, choice + '_log_bin_width')
        except AttributeError:
            pass
        if linbw is not None:
            n_bins = int((maxval - self.thr) / float(linbw))
            bins = np.linspace(self.thr - 0.0001, maxval, n_bins + 1)
        elif logbw is not None:
            n_bins = int(np.log(maxval / self.thr) / float(logbw))
            bins = np.logspace(np.log10(self.thr) - 0.0001, np.log10(maxval),
                               n_bins + 1)
        else:
            raise RuntimeError("Can't make bins without a %s bin width option!"
                               % choice)
        if self.args.verbose:
            print(str(n_bins) + ' ' + choice + ' stat bins')
        return bins

    def get_ctypes(self, tdict):
        # tdict is a ifo -> trigger time dictionary
        ifotimes = zip(*[tdict[i] for i in self.ifos])
        cty = []
        for ts in ifotimes:
            # if an ifo doesn't participate, time is sentinel value -1
            cty.append(''.join([i for i, t in zip(self.ifos, ts) if t > 0]))
        # return is array of coinc types strings
        return np.array(cty)

    def moreifotimes(self, ctstring):
        # get list of coinc times with more ifos than ctstring
        allctime_moreifos = [ct for ct in self.allctimestring if
                             len(ct) > len(ctstring)]
        # only return those when at least the same ifos are operating
        ret = []
        ifos = ifos_from_combo(ctstring)
        for act in allctime_moreifos:
            if all(i in act for i in ifos):
                ret.append(act)
        return ret

    def in_coinc_time_excl(self, f, cstring, test_times):
        """ filter to all times where exactly the ifos in cstring are observing
        """
        if len(cstring) == max(len(s) for s in self.allctimestring):
            # ctime string already uniquely specifies time
            return in_coinc_time_incl(f, cstring, test_times)
        in_time = in_coinc_time_incl(f, cstring, test_times)
        # if 'more-ifo' coinc times exist, exclude them
        for combo in self.moreifotimes(cstring):
            in_moreifo_time = in_coinc_time_incl(f, combo, test_times)
            # subtract one if in more-ifo time
            in_time -= in_moreifo_time
        # if subtraction yields anything other than 1, set to 0
        np.putmask(in_time, in_time != 1, 0)
        return in_time

    def get_livetimes(self, fi):
        with h5py.File(fi, 'r') as f:
            for ct in self.ctimes:
                # 'inclusive' time when at least the ifos specified by ct are on
                fgt = conv.sec_to_year(f[ct].attrs['foreground_time'])
                # index dict on chunk start time / coinc type
                self.incl_livetimes[(get_start_dur(fi)[0], ct)] = fgt
                # subtract times during which 1 more ifo was on,
                # ie subtract H1L1* time from H1L1; subtract H1* time from H1; etc
                for combo in self.moreifotimes(ct):
                    if len(combo) == len(ct) + 2:
                        fgt -= conv.sec_to_year(f[combo].attrs['foreground_time'])
                # index dict on chunk start time / coinc time
                self.livetimes[(get_start_dur(fi)[0], ct)] = fgt


class ForegroundEvents(EventRate):
    def __init__(self, args, coinc_times, coinc_types=None, bin_param='mchirp',
                 bin_lo=None, bin_hi=None):
        EventRate.__init__(self, args, coinc_times, coinc_types=coinc_types,
                           bin_param=bin_param, bin_lo=bin_lo, bin_hi=bin_hi)
        self.thr = self.args.stat_threshold
        # set of arrays in parallel containing zerolag event properties
        self.starttimes = []
        self.gpst = np.array([])
        self.stat = np.array([])
        self.ifar = np.array([])
        self.masspars = np.array([])
        self.start = np.array([])
        self.ctime = np.array([], dtype=object)  # allow unequal length strings
        self.ctype = np.array([], dtype=object)
        self.bg_pdf = np.array([])
        self.sg_pdf = np.array([])

    def add_zerolag(self, full_file):
        start = get_start_dur(full_file)[0]
        self.starttimes.append(start)
        with h5py.File(full_file, 'r') as f:
            # get stat values & threshold
            _stats = f['foreground/stat'][:]
            _keepstat = _stats > self.thr

            # get templates & apply filter
            _tids = f['foreground/template_id'][:]
            # we need the template filter to have already been made
            assert self.in_bin is not None
            _keep = np.logical_and(_keepstat, self.in_bin[_tids])
            massp = self.tpars[_tids][_keep]  # filtered template params

            # assign times and coinc types
            _times = {}
            for i in self.ifos:
                _times[i] = f['foreground/' + i + '/time'][:][_keep]
            # if an ifo doesn't participate, time is sentinel value -1
            # event time is mean of remaining positive GPS times
            meantimes = np.array([coinc_meanigz(ts)[0]
                                  for ts in zip(*_times.values())])
            _ctype = self.get_ctypes(_times)
            if len(_ctype) == 0:
                if self.args.verbose:
                    print('No events in ' + start)
                return
            # filter events
            in_ctypes = np.array([cty in self.coinc_types for cty in _ctype])
            meantimes = meantimes[in_ctypes]
            # get coinc time as strings
            # (strings may have different lengths)
            _ctime = np.repeat(np.array([''], dtype=object), len(meantimes))
            for ct in self.allctimestring:
                intime = self.in_coinc_time_excl(f, ct, meantimes)
                _ctime[intime == 1] = ct
                if self.args.verbose:
                    print('Got %i events in %s time' % (len(_ctime[intime == 1]), ct))
            # store
            self.stat = np.append(self.stat, _stats[_keep][in_ctypes])
            try:  # injection analyses only have 'ifar_exc', not 'ifar'
                self.ifar = np.append(self.ifar,
                                     f['foreground/ifar'][:][_keep][in_ctypes])
            except KeyError:
                self.ifar = np.append(self.ifar,
                                 f['foreground/ifar_exc'][:][_keep][in_ctypes])
            self.gpst = np.append(self.gpst, meantimes)
            self.masspars = np.append(self.masspars, massp)
            self.start = np.append(self.start, int(start) *
                                   np.ones_like(meantimes))
            self.ctime = np.append(self.ctime, _ctime)
            self.ctype = np.append(self.ctype, _ctype[in_ctypes])

    def get_bg_pdf(self, bg_rate):
        assert isinstance(bg_rate, BackgroundEventRate)
        self.bg_pdf = np.zeros_like(self.stat)  # initialize

        # do the calculation by chunk / coinc time / coinc type
        for st in self.starttimes:
            for ct in self.allctimestring:
                for cty in self.coinc_types:
                    if not type_in_time(ct, cty):
                        continue
                    _idx = np.logical_and((self.ctime == ct), (self.ctype == cty))
                    _idx = np.logical_and(_idx, (self.start == int(st)))
                    _vals = self.stat[_idx]
                    if len(_vals) == 0:
                        continue
                    # evaluate bg pdf for specific chunk, coinc time & type
                    _pdf = bg_rate.eval_pdf(st, ct, cty, _vals)
                    # store
                    self.bg_pdf[_idx] = _pdf
                    if self.args.verbose:
                        print('Found bg PDFs for ' + cty + ' coincs from ' + st)

    def get_sg_pdf(self, sg_rate):
        assert isinstance(sg_rate, SignalEventRate)
        self.sg_pdf = np.zeros_like(self.stat)

        for st in self.starttimes:
            for ct in self.allctimestring:
                for cty in self.coinc_types:
                    if not type_in_time(ct, cty):
                        continue
                    _idx = np.logical_and((self.ctime == ct), (self.ctype == cty))
                    _idx = np.logical_and(_idx, (self.start == int(st)))
                    _vals = self.stat[_idx]
                    if len(_vals) == 0:
                        continue
                    # norm of PDF is chunk-dependent so need the chunk start time
                    _pdf = sg_rate.eval_pdf(st, ct, cty, _vals)
                    # store
                    self.sg_pdf[_idx] = _pdf
                    if self.args.verbose:
                        print('Found sg PDFs for %s coincs in %s time from %s' %
                              (cty, ct, st))


class BackgroundEventRate(EventRate):
    def __init__(self, args, coinc_times, coinc_types=None, bin_param='mchirp',
                 bin_lo=None, bin_hi=None):
        EventRate.__init__(self, args, coinc_times, coinc_types=coinc_types,
                           bin_param=bin_param, bin_lo=bin_lo, bin_hi=bin_hi)
        self.thr = self.args.stat_threshold
        # BG values in dict indexed on tuple (chunk start, coinc type)
        self.bg_vals = {}
        self.bg_dec = {}
        # BG livetimes
        self.bg_livetimes = {}
        # BG hist stored as bin heights / edges
        self.bg_hist = {}
        # Expected counts of BG events
        self.exp_bg = {}
        # Total expected BG count
        self.norm = 0

    def add_background(self, full_file):
        start = get_start_dur(full_file)[0]
        self.get_livetimes(full_file)

        with h5py.File(full_file, 'r') as ff:
            # get stat values and threshold
            _bgstat = ff['background_exc/stat'][:]
            _keepstat = _bgstat > self.thr

            # get template ids and filter
            _bgtid = ff['background_exc/template_id'][:]
            # need the template filter to have already been made
            assert self.in_bin is not None
            _keep = np.logical_and(_keepstat, self.in_bin[_bgtid])
            _bgstat = _bgstat[_keep]
            _bgdec = ff['background_exc/decimation_factor'][:][_keep]

            # assign coinc types
            _times = {}
            for i in self.ifos:
                # NB times are time-shifted between ifos
                _times[i] = ff['background_exc/' + i + '/time'][:][_keep]
            _ctype = self.get_ctypes(_times)
            for cty in self.coinc_types:
                self.bg_vals[(start, cty)] = _bgstat[_ctype == cty]
                self.bg_dec[(start, cty)] = _bgdec[_ctype == cty]
                # get bg livetime for noise rate estimate
                # - convert to years
                self.bg_livetimes[(start, cty)] = conv.sec_to_year(
                                          ff[cty].attrs['background_time_exc'])

                # make histogram
                bins = self.make_bins(np.max(_bgstat[_ctype == cty]), 'bg')
                # hack to make larger bins for H1L1V1
                if cty == 'H1L1V1':
                    if self.args.verbose:
                        print('Halving bg bins for triple bg hist')
                    bins = bins[::2].copy()  # take every 2nd bin edge
                self.bg_hist[(start, cty)] = \
                      np.histogram(_bgstat[_ctype == cty],
                                   weights=_bgdec[_ctype == cty], bins=bins)
                # get expected number of bg events for this chunk and coinc type
                self.exp_bg[(start, cty)] = _bgdec[_ctype == cty].sum() * \
                      self.incl_livetimes[(start, cty)] / \
                      self.bg_livetimes[(start, cty)]

    def plot_bg(self):
        from matplotlib import pyplot as plt
        for chunk_type, hist in self.bg_hist.items():
            print('Plotting', chunk_type, 'background PDF ...')
            xplot = np.linspace(self.thr, self.args.plot_max_stat, 500)
            heights, bins = hist[0], hist[1]
            logpdf, _ = log_rho_bg(xplot, heights, bins)
            plt.plot(xplot, np.exp(logpdf))
            # plot error bars at bin centres
            lpdf, fracerr = log_rho_bg(0.5 * (bins[:-1] + bins[1:]), heights, bins)
            plt.errorbar(0.5 * (bins[:-1] + bins[1:]), np.exp(lpdf),
                         yerr=np.exp(lpdf) * fracerr, fmt='none')
            plt.semilogy()
            plt.grid(True)
            plt.xlim(xmax=self.args.plot_max_stat + 0.5)
            plt.ylim(ymin=0.7 * np.exp(logpdf.min()))
            plt.xlabel('Ranking statistic')
            plt.ylabel('Background PDF')
            plt.savefig(self.args.plot_dir + '%s-bg_pdf-%s' %
                        (chunk_type[1], chunk_type[0]) + '.png')
            plt.close()

    def get_norms(self):
        for count in self.exp_bg.values():
            self.norm += count

    def eval_pdf(self, chunk, ctime, ctype, statvals):
        # given statistic values all in the same data chunk and coinc type,
        # evaluate the background pdf normalized over all chunks & types
        assert self.norm > 0
        chunk_type = (chunk, ctype)
        # fraction of expected noise events in given chunk & coinc type
        frac_chunk_type = self.exp_bg[chunk_type] / self.norm
        # fraction of inj in specified chunk, coinc type *and* time
        frac_in_time = self.livetimes[(chunk, ctime)] /\
                                                self.incl_livetimes[chunk_type]
        # unpack heights / bins from bg hist object
        local_pdfs, _ = log_rho_bg(statvals, *self.bg_hist[chunk_type])
        return local_pdfs + np.log(frac_chunk_type * frac_in_time)


class SignalEventRate(EventRate):
    def __init__(self, args, coinc_times, coinc_types=None, bin_param='mchirp',
                 bin_lo=None, bin_hi=None):
        EventRate.__init__(self, args, coinc_times, coinc_types=coinc_types,
                           bin_param=bin_param, bin_lo=bin_lo, bin_hi=bin_hi)
        self.thr = self.args.stat_threshold
        self.starts = []  # bookkeeping
        # for the moment roll all inj chunks together
        # but sort both by coinc time and coinc type
        self.inj_vals = {}  # dict indexed on tuple (coinc time, coinc type)
        self.fg_bins = {}
        self.norm = 0

    def add_injections(self, inj_file, fg_file):
        # fg_file only needed for coinc time info :/
        self.starts.append(get_start_dur(inj_file)[0])
        self.get_livetimes(inj_file)

        with h5py.File(inj_file, 'r') as jf:
            # get stat values and threshold
            _injstat = jf['found_after_vetoes/stat'][:]
            _keepstat = _injstat > self.thr

            # get template ids and filter
            _injtid = jf['found_after_vetoes/template_id'][:]
            assert self.in_bin is not None
            _keep = np.logical_and(_keepstat, self.in_bin[_injtid])
            _injstat = _injstat[_keep]

            # assign coinc types
            _times = {}
            for i in self.ifos:
                _times[i] = jf['found_after_vetoes/' + i + '/time'][:][_keep]
            meantimes = np.array([coinc_meanigz(ts)[0]
                                  for ts in zip(*_times.values())])
            _ctype = self.get_ctypes(_times)
            # get coinc time as strings
            # (strings may have different lengths)
            _ctime = np.repeat(np.array([''], dtype=object), len(meantimes))
            for ct in self.allctimestring:
                # get coinc time info from segments in fg file
                intime = self.in_coinc_time_excl(
                                        h5py.File(fg_file, 'r'), ct, meantimes)
                _ctime[intime == 1] = ct  # do we need this?
                if self.args.verbose:
                    print('Got %i ' % (intime == 1).sum() + 'inj in %s time' % ct)
                # filter by coinc type and add to array
                for cty in self.coinc_types:
                    if not type_in_time(ct, cty):
                        continue
                    my_vals = _injstat[np.logical_and(_ctype == cty, intime == 1)]
                    if self.args.verbose:
                        print('%d ' % len(my_vals) + 'are %s coincs' % cty)
                    if (ct, cty) not in self.inj_vals:  # initialize
                        self.inj_vals[(ct, cty)] = np.array([])
                    if len(my_vals) > 0:
                        self.inj_vals[(ct, cty)] = \
                                  np.append(self.inj_vals[(ct, cty)], my_vals)
                del intime, my_vals

    def make_all_bins(self):
        for ct in self.allctimestring:
            for cty in self.coinc_types:
                if not type_in_time(ct, cty):
                    continue
                vals = self.inj_vals[(ct, cty)]
                # get norm of fg histogram by taking bins out to max injection stat
                binmax = vals.max() * 1.01
                self.fg_bins[(ct, cty)] = self.make_bins(binmax, 'inj')

    def plot_inj(self):
        from matplotlib import pyplot as plt
        for ct in self.allctimestring:
            for cty in self.coinc_types:
                if not type_in_time(ct, cty):
                    continue
                print('Plotting ' + cty + ' signal PDF in ' + ct + ' time ...')
                samples = self.inj_vals[(ct, cty)]
                bins = self.fg_bins[(ct, cty)]
                xplot = np.logspace(np.log10(self.thr),
                                    np.log10(samples.max()), 500)
                logpdf, _ = log_rho_fg(xplot, samples, bins)
                plt.plot(xplot, np.exp(logpdf))
                # plot error bars at bin centres
                lpdf, fracerr = log_rho_fg(0.5 * (bins[:-1] + bins[1:]),
                                           samples, bins)
                plt.errorbar(0.5 * (bins[:-1] + bins[1:]), np.exp(lpdf),
                             yerr=np.exp(lpdf) * fracerr, fmt='none')
                plt.semilogy()
                plt.grid(True)
                # zoom in on the 'interesting' range
                plt.xlim(xmin=self.thr, xmax=2. * self.args.plot_max_stat)
                plt.ylim(ymin=0.7 * np.exp(logpdf.min()))
                plt.title(r'%i injs plotted, \# of bins %i' %
                                                (len(samples), len(bins) - 1))
                plt.xlabel('Ranking statistic')
                plt.ylabel('Signal PDF')
                plt.savefig(self.args.plot_dir + '%s-fg_pdf-%s' % (ct, cty)
                            + '.png')
                plt.close()

    def get_norms(self):
        for vals in self.inj_vals.values():
            # injections don't have weights/decimation
            self.norm += float(len(vals))

    def eval_pdf(self, chunk, ctime, ctype, statvals):
        # given statistic values in the same chunk, coinc time and coinc type,
        # evaluate the signal pdf normalized over all chunks, times and types
        assert self.norm > 0
        time_type = (ctime, ctype)
        # fraction of inj in specified coinc time and type
        frac_time_type = float(len(self.inj_vals[time_type])) / self.norm
        # total livetime for specified coinc time
        total_coinc_time = sum([self.livetimes[(ch, ctime)] for ch in self.starts])
        # fraction of inj in specified chunk *and* coinc time/type
        this_norm = frac_time_type * self.livetimes[(chunk, ctime)] / \
                                   total_coinc_time
        local_pdfs, _ = log_rho_fg(statvals, self.inj_vals[time_type],
                                   self.fg_bins[time_type])
        return local_pdfs + np.log(this_norm)

__all__ = ['filter_bin_lo_hi', 'filter_tmplt_mchirp', 'read_full_data',
           'read_full_data_mchirp', 'log_rho_bg', 'log_rho_fg_analytic',
           'log_rho_fg', 'get_start_dur', 'in_coinc_time_incl', 'alltimes',
           'ifos_from_combo', 'type_in_time', 'EventRate', 'ForegroundEvents',
           'BackgroundEventRate', 'SignalEventRate']
