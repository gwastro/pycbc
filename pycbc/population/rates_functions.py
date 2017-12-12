"""
A set of helper functions for evaluating rates.
"""

import glob
from scipy import integrate, optimize
from scipy.special import erf
import numpy as np, h5py
import scipy.stats as ss

import bisect

def process_full_data(hdffile, rhomin, mass1, mass2, min_bh_mass):
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
    with h5py.File(hdffile, 'r') as bulk:
        
        tmplt_bkg = bulk['background_exc/template_id'][:]
        tmplt_fg = bulk['foreground/template_id'][:]

        idx_bkg = np.where( (mass1[tmplt_bkg] > 5.) & (mass2[tmplt_bkg] > min_bh_mass) )
        idx_fg = np.where( (mass1[tmplt_fg] > 5.) & (mass2[tmplt_fg] > min_bh_mass) )
                                                  
        zerolagstat = bulk['foreground/stat'][:][idx_fg]
        cstat_back_exc = bulk['background_exc/stat'][:][idx_bkg]
        dec_factors = bulk['background_exc/decimation_factor'][:][idx_bkg]

    return {'zerolagstat': zerolagstat[zerolagstat > rhomin],
           'dec_factors': dec_factors[cstat_back_exc > rhomin],
           'cstat_back_exc': cstat_back_exc[cstat_back_exc > rhomin]}


def merge_full_data(all_bkg):
    """Merge triggers over all chunks

       Parameters
       ----------
       all_bkg: dictionary

       Returns
       -------
       dictionary
              merged dictionaries
    """
    
    merged_bkg = {}
    for data in ['zerolagstat', 'dec_factors', 'cstat_back_exc']:
        merged_bkg[data] = np.concatenate([all_bkg[c][data] for c in all_bkg.keys()])
        
    return merged_bkg

def save_bkg_falloff(folder_name_statmap, folder_name_bank, path, rhomin, min_bh_mass = 5.):
    ''' Read the STATMAP files to derive snr falloff for the background events. 
        Save the output to a txt file
        Bank file is also provided to restrict triggers to BBH templates.

        Parameters
        ----------
        folder_name_statmap: string
               Folder name where STATMAP files are located
        folder_name_bank: string
               Folder name where the bank files are located
        path: string
               Destination where txt file is saved
        rhomin: float
               Minimum value of SNR threhold (will need including ifar)
        min_bh_mass: float
              Minimimum mass for the black hole
    '''

    full_data = {}
    statmap_files = glob.glob(folder_name_statmap)
    bank_files = glob.glob(folder_name_bank)
 
    i = 0
    for sfile in statmap_files:
        for bfile in bank_files:
            
            if  sfile[-22:] == bfile[-22:]: # check if the encoded GPS match (better FIXME ?)
                with h5py.File(bfile, 'r') as bulk:
                    
                    print "Loading zero-lag results from file: %s" % sfile
                    print "Bank file is: %s" % bfile
        
                    mass1_bank = bulk['mass1'][:]
                    mass2_bank = bulk['mass2'][:]
                    full_data[str(i)] = process_full_data(sfile, rhomin, mass1_bank, mass2_bank, min_bh_mass)
                    i =+ 1
                break
                       
    if len(statmap_files) > len(full_data.keys()):
        raise Exception('Missing bank file!')
     
        
    full_data = merge_full_data(full_data)

    max_bg_stat = np.max(full_data['cstat_back_exc'])
    bg_bins = np.linspace(rhomin, max_bg_stat, 76)
    bg_counts, junk = np.histogram(full_data['cstat_back_exc'], 
                                weights=full_data['dec_factors'],
                                bins=bg_bins)    

    zerolagstat = full_data['zerolagstat']
    coincs = zerolagstat[zerolagstat >= rhomin]

    np.savetxt(path+"/background_bins.txt", np.column_stack([bg_bins[:-1], bg_bins[1:], bg_counts]), 
            fmt='%.4e', header="bin min, bin max, count")
    np.savetxt(path+"/coincs.txt", coincs, fmt='%.4e', header="coincs above threshold %.2f" % rhomin)

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
    
    assert np.all(trigs >= np.min(bins)), 'cannot have triggers smaller than bin lower limit'
    
    # If there are any triggers that are louder than the max bin, put one 
    # fictituous count in a bin that extends from the limits of the slide triggers
    # out to the loudest trigger.
    
    # If there is no counts for a foreground trigger put a fictious count in the background bin
    if np.any(trigs >= np.max(bins)):
        N = N + 1
        #log_plimit = -np.log(N) - np.log(np.max(trigs) - bins[-1])  NEEDS CHECKING
    
    log_rhos = []
    for t in trigs:
        if t >= np.max(bins):
            log_rhos.append(-np.log(N)-np.log(np.max(trigs) - bins[-1])) 
        else:
            i = bisect.bisect(bins, t) - 1
            #print t, counts[i], np.log(bins[i+1] - bins[i])
            if counts[i] == 0:
                counts[i] = 1
            log_rhos.append(np.log(counts[i]) - np.log(bins[i+1] - bins[i]) - np.log(N))
    return np.array(log_rhos)


#def log_rho_fg_analytic(t):
#    return np.log(3.0) + 3.0*np.log(rhomin) - 4.0*np.log(t)


def log_rho_fg_mc(t, injstats, bins):
    counts, bins = np.histogram(injstats, bins)
    
    N = sum(counts)
    dens = counts / np.diff(bins) / float(N)

    assert np.min(t) >= np.min(bins)
    assert np.max(t) < np.max(bins)
    
    tinds = np.searchsorted(bins, t) - 1
    
    return np.log(dens[tinds])
    
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
    
    def pdf(x):
        return 1/np.sqrt(2*np.pi) * np.exp(-x**2/2)

    def cdf(x):
        return (1 + erf(x/np.sqrt(2))) / 2

    def skew(x, alpha, mu, sigma):
        t = (x - mu) / sigma
        return 2 / sigma * pdf(t) * cdf(alpha*t)

    def optm(l, y, x):
            alpha, mu, sigma = l
            return kde(x)- skew(x, alpha, mu, sigma)    
    
    R = np.log(R)
    xs = np.linspace(min(R), max(R), 200)
    
    mu_norm, sigma_norm = np.mean(R), np.std(R)
    
    kde = ss.gaussian_kde(R)
    p_xs = kde(xs)
    
    ff, flag = optimize.leastsq(optm, [0.1, mu_norm, sigma_norm], args=(p_xs, xs))
    return ff[0], ff[1], ff[2]

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
    
    C_lnm = integrate.quad(lambda x: (np.log(max_mtotal - x) - np.log(min_mass))/x, min_mass, max_mass)[0]    
    
    xx = np.minimum(m1, m2)
    m1 = np.maximum(m1, m2)
    m2 = xx    
    
    bound =np.sign(max_mtotal - m1 - m2) + np.sign(max_mass - m1) * np.sign(m2 - min_mass)
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
    
    C_imf = max_mass**(alpha + 1)/(alpha + 1) - min_mass**(alpha + 1)/(alpha + 1)    
    
    xx = np.minimum(m1, m2)
    m1 = np.maximum(m1, m2)
    m2 = xx    
    
    bound = np.sign(max_mtotal - m1 - m2) + np.sign(max_mass - m1) * np.sign(m2 - min_mass)
    idx = np.where(bound != 2)
    
    p_m1_m2 = np.zeros_like(m1)
    p_m1_m2[m1 <= max_mtotal/2.] = (1./C_imf) * m1[m1 <= max_mtotal/2.]**alpha /(m1[m1 <= max_mtotal/2.] - min_mass)
    p_m1_m2[m1 > max_mtotal/2.] = (1./C_imf) * m1[m1 > max_mtotal/2.]**alpha /(max_mass - m1[m1 > max_mtotal/2.])
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
    
    k = nsamples * int(1.5 + np.log(1 + 100./nsamples))
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
    lnmmin = np.log(min_mass)
    lnmmax = np.log(max_mass)
    
    k = nsamples * int(1.5 + np.log(1 + 100./nsamples))
    aa = np.exp(np.random.uniform(lnmmin, lnmmax, k))
    bb = np.exp(np.random.uniform(lnmmin, lnmmax, k))
  
    idx = np.where(aa + bb < max_mtotal)
    m1, m2 = (np.maximum(aa, bb))[idx], (np.minimum(aa, bb))[idx]
    
    return np.resize(m1, nsamples), np.resize(m2, nsamples)

def m1m2_to_mcheta(m1, m2):
    ''' Get chirp mass and eta for m1 and m2

        Parameters
        ----------
        m1: array
            Component masses 1
        m2: array
            Component masses 2

        Returns
        -------
        array
           Chirp mass
        array
           Symmetric mass ratio 
    '''
    m1, m2 = np.array(m1), np.array(m2)
    return (m1*m2)**.6/(m1+m2)**.2, m1*m2/(m1+m2)**2

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
    mchirp_astro, junk = m1m2_to_mcheta(m1, m2)
    
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
    mchirp_astro, junk = m1m2_to_mcheta(m1, m2)
    
    return mchirp_astro
