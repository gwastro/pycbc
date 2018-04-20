"""
A set of helper functions for evaluating rates.
"""

from scipy import integrate, optimize
import numpy as np, h5py
from numpy import log
import scipy.stats as ss

import bisect
from pycbc.conversions import mchirp_from_mass1_mass2

def process_full_data(fname, rhomin, mass1, mass2, lo_tlt_mass, hi_tlt_mass):
    """Read the zero and time-lagged triggers identified by BBH templates.

       Parameters
       ----------
       hdfile:
              File that stores all the triggers
       rhomin: float
              Minimum value of SNR threhold (will need including ifar)
       mass1: array
              First mass of the waveform in the template bank
       mass2: array
              Second mass of the waveform in the template bank
       min_bh_mass: float
              Minimimum mass for the black hole

       Returns
       -------
       dictionary
              containing foreground triggers and background information
    """
    with h5py.File(fname, 'r') as bulk:

        id_bkg = bulk['background_exc/template_id'][:]
        id_fg = bulk['foreground/template_id'][:]

        bound = np.sign((mass1[id_bkg] - lo_tlt_mass) * (hi_tlt_mass - mass1[id_bkg]))
        bound += np.sign((mass2[id_bkg] - lo_tlt_mass) * (hi_tlt_mass - mass2[id_bkg]))
        idx_bkg = np.where(bound == 2)
        bound = np.sign((mass1[id_fg] - lo_tlt_mass) * (hi_tlt_mass - mass1[id_fg]))
        bound += np.sign((mass2[id_fg] - lo_tlt_mass) * (hi_tlt_mass - mass2[id_fg]))
        idx_fg = np.where(bound == 2)

        zerolagstat = bulk['foreground/stat'][:][idx_fg]
        cstat_back_exc = bulk['background_exc/stat'][:][idx_bkg]
        dec_factors = bulk['background_exc/decimation_factor'][:][idx_bkg]

    return {'zerolagstat': zerolagstat[zerolagstat > rhomin],
           'dec_factors': dec_factors[cstat_back_exc > rhomin],
           'cstat_back_exc': cstat_back_exc[cstat_back_exc > rhomin]}

def save_bkg_falloff(fname_statmap, fname_bank, path, rhomin, lo_tlt_mass, hi_tlt_mass):
    ''' Read the STATMAP files to derive snr falloff for the background events.
        Save the output to a txt file
        Bank file is also provided to restrict triggers to BBH templates.

        Parameters
        ----------
        fname_statmap: string
               STATMAP file containing trigger information
        fname_bank: string
               File name of the template bank
        path: string
               Destination where txt file is saved
        rhomin: float
               Minimum value of SNR threhold (will need including ifar)
        lo_tlt_mass: float
               Low mass for template for trigger to be considered
        hi_tlt_mass: float
               High mass for template for trigger to be considered
    '''

    with h5py.File(fname_bank, 'r') as bulk:
        mass1_bank = bulk['mass1'][:]
        mass2_bank = bulk['mass2'][:]
        full_data = process_full_data(fname_statmap, rhomin,
                           mass1_bank, mass2_bank, lo_tlt_mass, hi_tlt_mass)

    max_bg_stat = np.max(full_data['cstat_back_exc'])
    bg_bins = np.linspace(rhomin, max_bg_stat, 76)
    bg_counts = np.histogram(full_data['cstat_back_exc'],
                         weights=full_data['dec_factors'], bins=bg_bins)[0]

    zerolagstat = full_data['zerolagstat']
    coincs = zerolagstat[zerolagstat >= rhomin]

    np.savetxt(path+"/background_bins.txt",
               np.column_stack([bg_bins[:-1],
               bg_bins[1:], bg_counts]), fmt='%.4e',
               header="bin min, bin max, count")

    np.savetxt(path+"/coincs.txt", coincs,
               fmt='%.4e', header="coincs above threshold %.2f" % rhomin)

def log_rho_bg(trigs, bins, counts):
    ''' Calculate the log of background fall-off

        Parameters
        ----------
        trigs: array
               SNR values of all the triggers
        bins: string
               bins for histogrammed triggers
        path: string
               counts for histogrammed triggers

        Returns
        -------
        array
    '''

    trigs = np.atleast_1d(trigs)

    N = sum(counts)

    assert np.all(trigs >= np.min(bins)), 'triggers can not besmaller than bin lower limit'

    # If there are any triggers that are louder than the max bin, put one
    # fictituous count in a bin that extends from the limits of the slide
    # triggers out to the loudest trigger.

    # If there is no counts for a foreground trigger put a fictious count
    # in the background bin
    if np.any(trigs >= np.max(bins)):
        N = N + 1
        #log_plimit = -np.log(N) - np.log(np.max(trigs) - bins[-1]) CHECK IT

    log_rhos = []
    for t in trigs:
        if t >= np.max(bins):
            log_rhos.append(-log(N)-log(np.max(trigs) - bins[-1]))
        else:
            i = bisect.bisect(bins, t) - 1

            if counts[i] == 0:
                counts[i] = 1
            log_rhos.append(log(counts[i]) - log(bins[i+1] - bins[i]) - log(N))
    return np.array(log_rhos)

def log_rho_fg_mc(t, injstats, bins):
    counts, bins = np.histogram(injstats, bins)

    N = sum(counts)
    dens = counts / np.diff(bins) / float(N)

    assert np.min(t) >= np.min(bins)
    assert np.max(t) < np.max(bins)

    tinds = np.searchsorted(bins, t) - 1

    return log(dens[tinds])

def fg_mc(log_fg_ratios, mu_log_vt, sigma_log_vt, Rf, maxfg):
    '''
    Function to fit the likelihood Fixme
    '''

    Lb = np.random.uniform(0., maxfg, len(Rf))
    pquit = 0

    while pquit < 0.1:
        # quit when the posterior on Lf is very close to its prior

        nsamp = len(Lb)
        Rf_sel = np.random.choice(Rf, nsamp)
        vt = np.random.lognormal(mu_log_vt, sigma_log_vt, len(Rf_sel))

        Lf = Rf_sel * vt

        log_Lf, log_Lb = log(Lf), log(Lb)

        plR = 0
        for lfr in log_fg_ratios:

            plR += np.logaddexp(lfr + log_Lf, log_Lb)

        plR -= (Lf + Lb)

        plRn = plR - max(plR)
        idx = np.exp(plRn) > np.random.random(len(plRn))

        pquit = ss.stats.ks_2samp(Lb, Lb[idx])[1]

        Lb = Lb[idx]

    return Rf_sel[idx], Lf[idx], Lb

def _optm(x, alpha, mu, sigma):
    '''Return probability density of skew-lognormal
       See scipy.optimize.curve_fit
    '''
    return ss.skewnorm.pdf(x, alpha, mu, sigma)

def fit(R):
    ''' Fit skew - lognormal to the rate samples achived from a prior analysis
        Parameters
        ----------
        R: array
           Rate samples
        Returns
        -------
        ff[0]: float
            The skewness
        ff[1]: float
            The mean
        ff[2]: float
            The standard deviation
    '''

    R = np.log(R)
    mu_norm, sigma_norm = np.mean(R), np.std(R)

    xs = np.linspace(min(R), max(R), 200)
    kde = ss.gaussian_kde(R)
    pxs = kde(xs)

    # Initial guess has been taken as the mean and std-dev of the data
    # And a guess assuming small skewness
    ff = optimize.curve_fit(_optm, xs, pxs, p0 = [0.1, mu_norm, sigma_norm])[0]
    return ff[0], ff[1], ff[2]

def skew_lognormal_samples(alpha, mu, sigma, minrp, maxrp):
    ''' Returns a large number of Skew lognormal samples
        Parameters
        ----------
        alpha: float
           Skewness of the distribution
        mu: float
           Mean of the distribution
        sigma: float
           Scale of the distribution
        minrp: float
           Minimum value for the samples
        maxrp: float
           Maximum value for the samples
        Returns
        -------
        Rfs: array
            Large number of samples (may need fixing)
    '''

    nsamp = 100000000
    lRu = np.random.uniform(minrp, maxrp, nsamp)
    plRu = ss.skewnorm.pdf(lRu, alpha, mu, sigma)
    rndn = np.random.random(nsamp)
    maxp = max(plRu)
    idx = np.where(plRu/maxp > rndn)
    log_Rf = lRu[idx]
    Rfs = np.exp(log_Rf)

    return Rfs

# The flat in log and power-law mass distribution models  #

# PDF for the two caninical models
def prob_lnm(m1, m2, s1z, s2z, **kwargs):
    ''' Return probability density for uniform in log
        Parameters
        ----------
        m1: array
            Component masses 1
        m2: array
            Component masses 2
        s1z: array
            Aligned spin 1(Not in use currently)
        s2z:
            Aligned spin 2(Not in use currently)
        **kwargs: string
            Keyword arguments as model parameters
        Returns
        -------
        p_m1_m2: array
            The probability density for m1, m2 pair
    '''

    min_mass = kwargs.get('min_mass', 5.)
    max_mass = kwargs.get('max_mass', 95.)
    max_mtotal = min_mass + max_mass
    m1, m2 = np.array(m1), np.array(m2)

    C_lnm = integrate.quad(lambda x: (log(max_mtotal - x) - log(min_mass))/x, min_mass, max_mass)[0]

    xx = np.minimum(m1, m2)
    m1 = np.maximum(m1, m2)
    m2 = xx

    bound = np.sign(max_mtotal - m1 - m2)
    bound += np.sign(max_mass - m1) * np.sign(m2 - min_mass)
    idx = np.where(bound != 2)

    p_m1_m2 = (1/C_lnm)*(1./m1)*(1./m2)
    p_m1_m2[idx] = 0

    return p_m1_m2

def prob_imf(m1, m2, s1z, s2z, **kwargs):
    ''' Return probability density for power-law
        Parameters
        ----------
        m1: array
            Component masses 1
        m2: array
            Component masses 2
        s1z: array
            Aligned spin 1(Not in use currently)
        s2z:
            Aligned spin 2(Not in use currently)
        **kwargs: string
            Keyword arguments as model parameters

        Returns
        -------
        p_m1_m2: array
           the probability density for m1, m2 pair
    '''

    min_mass = kwargs.get('min_mass', 5.)
    max_mass = kwargs.get('max_mass', 95.)
    alpha = kwargs.get('alpha', -2.35)
    max_mtotal = min_mass + max_mass
    m1, m2 = np.array(m1), np.array(m2)

    C_imf = max_mass**(alpha + 1)/(alpha + 1)
    C_imf -= min_mass**(alpha + 1)/(alpha + 1)

    xx = np.minimum(m1, m2)
    m1 = np.maximum(m1, m2)
    m2 = xx

    bound = np.sign(max_mtotal - m1 - m2)
    bound += np.sign(max_mass - m1) * np.sign(m2 - min_mass)
    idx = np.where(bound != 2)

    p_m1_m2 = np.zeros_like(m1)
    idx = np.where(m1 <= max_mtotal/2.)
    p_m1_m2[idx] = (1./C_imf) * m1[idx]**alpha /(m1[idx] - min_mass)
    idx = np.where(m1 > max_mtotal/2.)
    p_m1_m2[idx] = (1./C_imf) * m1[idx]**alpha /(max_mass - m1[idx])
    p_m1_m2[idx] = 0

    return p_m1_m2/2.

# Functions to generate samples for the two canonical models
def draw_imf_samples(**kwargs):
    ''' Draw samples for power-law model

        Parameters
        ----------
        **kwargs: string
           Keyword arguments as model parameters and number of samples

        Returns
        -------
        array
           The first mass
        array
           The second mass
    '''

    alpha_salpeter = kwargs.get('alpha', -2.35)
    nsamples = kwargs.get('nsamples', 1)
    min_mass = kwargs.get('min_mass', 5.)
    max_mass = kwargs.get('max_mass', 95.)
    max_mtotal = min_mass + max_mass

    a = (max_mass/min_mass)**(alpha_salpeter + 1.0) - 1.0
    beta = 1.0 / (alpha_salpeter + 1.0)

    k = nsamples * int(1.5 + log(1 + 100./nsamples))
    aa = min_mass * (1.0 + a * np.random.random(k))**beta
    bb = np.random.uniform(min_mass, aa, k)

    idx = np.where(aa + bb < max_mtotal)
    m1, m2 = (np.maximum(aa, bb))[idx], (np.minimum(aa, bb))[idx]

    return np.resize(m1, nsamples), np.resize(m2, nsamples)

def draw_lnm_samples(**kwargs):
    ''' Draw samples for uniform-in-log model

        Parameters
        ----------
        **kwargs: string
           Keyword arguments as model parameters and number of samples

        Returns
        -------
        array
           The first mass
        array
           The second mass
    '''

    #PDF doesnt match with sampler
    nsamples = kwargs.get('nsamples', 1)
    min_mass = kwargs.get('min_mass', 5.)
    max_mass = kwargs.get('max_mass', 95.)
    max_mtotal = min_mass + max_mass
    lnmmin = log(min_mass)
    lnmmax = log(max_mass)

    k = nsamples * int(1.5 + log(1 + 100./nsamples))
    aa = np.exp(np.random.uniform(lnmmin, lnmmax, k))
    bb = np.exp(np.random.uniform(lnmmin, lnmmax, k))

    idx = np.where(aa + bb < max_mtotal)
    m1, m2 = (np.maximum(aa, bb))[idx], (np.minimum(aa, bb))[idx]

    return np.resize(m1, nsamples), np.resize(m2, nsamples)

# Functions to generate chirp mass samples for the two canonical models
def mchirp_sampler_lnm(**kwargs):
    ''' Draw chirp mass samples for uniform-in-log model

        Parameters
        ----------
        **kwargs: string
           Keyword arguments as model parameters and number of samples

        Returns
        -------
        mchirp-astro: array
           The chirp mass samples for the population
    '''
    m1, m2 = draw_lnm_samples(**kwargs)
    mchirp_astro = mchirp_from_mass1_mass2(m1, m2)

    return mchirp_astro

def mchirp_sampler_imf(**kwargs):
    ''' Draw chirp mass samples for power-law model

        Parameters
        ----------
        **kwargs: string
           Keyword arguments as model parameters and number of samples

        Returns
        -------
        mchirp-astro: array
           The chirp mass samples for the population
    '''
    m1, m2 = draw_imf_samples(**kwargs)
    mchirp_astro = mchirp_from_mass1_mass2(m1, m2)

    return mchirp_astro
