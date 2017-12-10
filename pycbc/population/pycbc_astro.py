#!/usr/bin/env python

# Copyright (C) 2017 Vaibhav Tiwari

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

"""
Given an xml(.gz) table with CBC injections, this script separates them into:
(1) (potentially) found injections
(2) injections that we expect to miss
The two sets are stored into two separate output files. 
"""

__author__ = "Vaibhav Tiwari"
__email__ = "vaibhav.tiwari@ligo.org"
__version__ = "0.0"
__date__ = "31.10.2017"

import logging
import argparse
import pycbc.version

import os
import h5py
import pystan
import numpy as np
import scipy.stats as ss
import scale_injections as si
from rates_functions import *
import matplotlib.pyplot as plt

import seaborn as sns
import matplotlib as mpl

# Parse command line
parser = argparse.ArgumentParser(description='Estimate rates for flat in log and power-law mass distribution.')
#parser.add_argument("--version", action="version",
#                    version=pycbc.version.git_verbose_msg)
#group = parser.add_mutually_exclusive_group(required=True)

parser.add_argument('--sim-folder', dest='sim_folder', required=True, 
                  help='Folder containing simulation files - delimiters can be used.')
parser.add_argument('--bank-folder', dest='bank_folder', required=True, 
                  help='Folder containing simulation files - delimiters can be used.')
parser.add_argument('--statmap-folder', dest='statmap_folder', required=True,
                  help="Folder containing the STATMAP files (files that save background and foregound trigger information.")
parser.add_argument('--stan-file', dest='stan_file', required=True,
                  help="Stan file to perform the fit over simulation and background events.")
parser.add_argument('--prior-samples', dest='prior_samples', required=True,
                  help="File storing samples of prior for the analysis - posterior from the previous analysis")
parser.add_argument('--output-folder', dest='output_folder', required=True,
                  help="Output folder where files are to be saved.")
parser.add_argument('--min-snr', dest='min_snr', required=False, default = 8.0, 
                  help="Minimum SNR value at which rates are to be calculated")
    
opts = parser.parse_args()
path = opts.output_folder + '/astro_output'
if path is not None and not os.path.exists(path):
    os.makedirs(path)

# Read the simulation files
injections = si.read_injections(opts.sim_folder)

# Read the chirp-mass samples -- Imported from rates_function
mchirp_sampler, prob, label = {}, {}, {}
distrs = ['lnm', 'imf']
mchirp_sampler['lnm'] = mchirp_sampler_lnm
mchirp_sampler['imf'] = mchirp_sampler_imf

prob['lnm'] = prob_lnm
prob['imf'] = prob_imf
label['lnm'] = 'Flat in log'
label['imf'] = 'Power-law'

# Estimate the rates and make supporting plots
vt, vol_time, sig_vt, snr_falloff = {}, {}, {}, {}
for dist in distrs:
    vt[dist] = si.estimate_vt(injections, mchirp_sampler[dist], prob[dist])
    
    vol_time[dist], sig_vt[dist] = si.get_summed_vt(vt[dist])
    snr_falloff[dist] = si.get_accumulated_falloff(vt[dist])
    
#for dist in distrs:
#    hist4, bins= np.histogram(snr_falloff[dist], bins = 30, density=True) 
#    width, center = (bins[1] - bins[0]), (bins[:-1] + bins[1:]) / 2
#    plt.bar(center, hist4, align='center', width=width, alpha = 0.8, label = label[dist])
#    plt.legend(loc = 'best')
#plt.title('SNR Fall Off') 


#Sabe background data and coincidences
save_bkg_falloff(opts.statmap_folder, opts.bank_folder, path, opts.min_snr, 5.0)

#Load background data and coincidences/ make some plots
coincs = np.loadtxt(path + "/coincs.txt")
bg_l, bg_h, bg_counts = np.loadtxt(path + "/background_bins.txt", unpack=True)
bg_bins = np.append(bg_l, bg_h[-1])

fg_stats = np.concatenate([snr_falloff[dist] for dist in distrs])
fg_bins = np.logspace(np.log10(opts.min_snr), np.log10(np.max(fg_stats)), 101)
fg_stats = fg_stats[fg_stats > opts.min_snr]

ts = np.linspace(opts.min_snr, 11.0, 1000)
plt.figure()
plt.plot(ts, np.exp(log_rho_bg(ts, bg_bins, bg_counts)), label=r'$p_0$')
plt.plot(ts, np.exp(log_rho_fg_mc(ts, fg_stats, fg_bins)), label=r'$p_1$')
plt.yscale('log')
plt.xlabel(r"$x'$")
plt.ylabel(r"$p(x')$")
plt.legend(loc='lower left')
plt.axis(ymin=1e-3)
plt.tight_layout()
plt.savefig(path+'/fg-bg.png')

#Load prior samples and fit a skew-log-normal to it
with h5py.File(opts.prior_samples, "r") as f:
    Rfl = np.array(f['flat/Rf'])
    Rpl = np.array(f['power-law/Rf'])
    
alpha, mu, sigma = {}, {}, {}
alpha['lnm'], mu['lnm'], sigma['lnm'] = fit(Rfl)
alpha['imf'], mu['imf'], sigma['imf'] = fit(Rpl)

#Estimate rates
rate_data = {}
rate_samples = {}
stan_file = '/home/vaibhav/software/my_svn/rates_related/scale_injections/log_dist/bkg_and_rates_prior/rate-single.stan'
log_fg_ratios = log_rho_fg_mc(coincs, fg_stats, fg_bins) - log_rho_bg(coincs, bg_bins, bg_counts)
for dist in distrs:
    rate_data[dist] = {'Ncoincs': len(log_fg_ratios),
           'log_fg_ratios': log_fg_ratios,
           'mu_vt': vol_time[dist][0]/1e9,
           'sigma_vt': sig_vt[dist][0]/1e9,
           'alpha': alpha[dist],
           'mu': mu[dist],
           'sigma': sigma[dist]}
    rate_fit = pystan.stan(file=stan_file, data=rate_data[dist], iter=20000)
    rate_samples[dist] = rate_fit.extract(permuted=True)
    rate_fit

r50, r95, r05 = {}, {}, {}    
for dist in distrs:
    r50[dist], r95[dist], r05[dist] = np.percentile(rate_samples[dist]['Rf'], [50, 95, 5])    
    
#Save rate posteriors
with h5py.File(path+'/rate-posterior.hdf5', 'w') as out:

    pl = out.create_group('power-law')
    pl.create_dataset('Lf', data=rate_samples['imf']['Lf'], compression='gzip')
    pl.create_dataset('Rf', data=rate_samples['imf']['Rf'], compression='gzip')
    
    flat = out.create_group('flat')
    flat.create_dataset('Lf', data=rate_samples['lnm']['Lf'], compression='gzip')
    flat.create_dataset('Rf', data=rate_samples['lnm']['Rf'], compression='gzip')
    
    d = out.create_group('data')
    d.create_dataset('log_fg_bg_ratio', data=log_fg_ratios, compression='gzip')
    d.create_dataset('newsnr', data=coincs, compression='gzip')

# Make prior/posterior plot
post_alpha, post_mu, post_sigma = {}, {}, {}
post_alpha['lnm'], post_mu['lnm'], post_sigma['lnm'] = fit(rate_samples['lnm']['Rf'])
post_alpha['imf'], post_mu['imf'], post_sigma['imf'] = fit(rate_samples['imf']['Rf'])

plt.figure()
log_R = np.log(Rfl)
xs = np.linspace(min(log_R), max(log_R), 200)
plt.plot(xs, ss.skewnorm.pdf(xs, alpha['lnm'], mu['lnm'], sigma['lnm']), '--', label='Flat Prior', color=sns.color_palette()[0])

log_R = np.log(Rpl)
xs = np.linspace(min(log_R), max(log_R), 200)
plt.plot(xs, ss.skewnorm.pdf(xs, alpha['imf'], mu['imf'], sigma['imf']), '--', label='Power Prior', color=sns.color_palette()[1])

log_R = np.log(Rfl)
xs = np.linspace(min(log_R), max(log_R), 200)
plt.plot(xs, ss.skewnorm.pdf(xs, post_alpha['lnm'], mu['lnm'], post_sigma['lnm']), label='Flat Posterior', color=sns.color_palette()[0])

log_R = np.log(Rpl)
xs = np.linspace(min(log_R), max(log_R), 200)
plt.plot(xs, ss.skewnorm.pdf(xs, post_alpha['imf'], post_mu['imf'], post_sigma['imf']), label='Power Posterior', color=sns.color_palette()[1])

plt.xlabel(r'$log(R)$ ($\mathrm{Gpc}^{-3} \, \mathrm{yr}^{-1}$)')
plt.ylabel(r'$R p(R)$')

plt.legend(loc='best')
plt.savefig(path+'/rate_prior_posterior.png')

#Estimate p_astro for top 5 events

from numpy import logaddexp, log, newaxis, expm1
log_pastros_imf = logaddexp.reduce(log(rate_samples['imf']['Lf'][:, newaxis]) + log_fg_ratios[newaxis,:] - logaddexp(log(rate_samples['imf']['Lf'][:, newaxis]) + log_fg_ratios[newaxis,:], log(rate_samples['imf']['Lb'][:, newaxis])), axis=0) - log(rate_samples['imf']['Lf'].shape[0])
log_pastros_lnm = logaddexp.reduce(log(rate_samples['lnm']['Lf'][:, newaxis]) + log_fg_ratios[newaxis,:] - logaddexp(log(rate_samples['lnm']['Lf'][:, newaxis]) + log_fg_ratios[newaxis,:], log(rate_samples['lnm']['Lb'][:, newaxis])), axis=0) - log(rate_samples['lnm']['Lf'].shape[0])
plt.figure()
plt.plot(log_fg_ratios, -expm1(log_pastros_imf), '.', label='IMF')
plt.plot(log_fg_ratios, -expm1(log_pastros_lnm), '.', label='LNM')
plt.xlabel(r'$\log p(x\mid f)/p(x\mid b)$')
plt.ylabel(r'$1-p_\mathrm{astro}$')
plt.legend(loc='best')
plt.yscale('log')
plt.savefig(path+'/p_astro.png')

p_astro = [1 + expm1(np.sort(log_pastros_imf)[::-1]), 1 + expm1(np.sort(log_pastros_lnm)[::-1])]
   
return p_astro 
