from pycbc.types import *
from pycbc.noise.gaussian import *
from pycbc.filter import *
from pycbc.waveform import *
from pycbc.vetoes import *
import pycbc.psd
from numpy.random import normal

sr = 4096.0
dt = 1.0/sr
bl = 256
df = 1.0/bl
N = int(sr * bl)
n = int(N/2 + 1)

psd = pycbc.psd.from_string("aLIGOZeroDetHighPower", n, df, 14)
strain = noise_from_psd(N, dt, psd, seed=0)

htildep, htildec = get_fd_waveform(approximant="TaylorF2", mass1=10, mass2=10, f_lower=15, delta_f=df)
htildep.resize(n)

snr, corr, norm = matched_filter_core(htildep, strain, psd, low_frequency_cutoff=15)

# TEST 1
bank_tilde1,_ = get_fd_waveform(approximant="TaylorF2", mass1=1.4, mass2=1.4, f_lower=15, delta_f=df)
bank_tilde2,_ = get_fd_waveform(approximant="TaylorF2", mass1=2.0, mass2=2.0, f_lower=15, delta_f=df)
bank_tilde3,_ = get_fd_waveform(approximant="TaylorF2", mass1=4.0, mass2=4.0, f_lower=15, delta_f=df)
bank_tilde4,_ = get_fd_waveform(approximant="TaylorF2", mass1=9.0, mass2=9.0, f_lower=15, delta_f=df)
bank_tilde5,_ = get_fd_waveform(approximant="TaylorF2", mass1=6.0, mass2=1.4, f_lower=15, delta_f=df)

bank_tilde1.resize(n)
bank_tilde2.resize(n)
bank_tilde3.resize(n)
bank_tilde4.resize(n)
bank_tilde5.resize(n)

bank_veto_bank = [bank_tilde1,bank_tilde2,bank_tilde3,bank_tilde4,bank_tilde5]
bank_veto_curr_overlaps = []
bank_snrs = []
bank_norms = []
for bank_template in bank_veto_bank:
    # For every bank veto template compute overlap between template
    # and the data
    curr_bank_snr,_,curr_bank_norm = matched_filter_core(bank_template,\
            strain,psd,low_frequency_cutoff=15)
    # SNR time series stored here
    bank_snrs.append(curr_bank_snr)
    # Template normalization factor stored here
    bank_norms.append(curr_bank_norm)
    bank_veto_curr_overlaps.append(overlap_cplx(htildep,\
             bank_template,psd=psd,\
             low_frequency_cutoff=15))

bank_veto = bank_chisq_from_filters(snr,norm,bank_snrs,bank_norms,bank_veto_curr_overlaps)
#numpy.savetxt('BV_TEST1.txt',bank_veto)

# TEST 2
bank_tilde1,_ = get_fd_waveform(approximant="TaylorF2", mass1=9.9, mass2=9.9, f_lower=15, delta_f=df)
bank_tilde2,_ = get_fd_waveform(approximant="TaylorF2", mass1=9.9, mass2=10., f_lower=15, delta_f=df)
bank_tilde3,_ = get_fd_waveform(approximant="TaylorF2", mass1=9.9, mass2=10.1, f_lower=15, delta_f=df)
bank_tilde4,_ = get_fd_waveform(approximant="TaylorF2", mass1=10., mass2=10.1, f_lower=15, delta_f=df)
bank_tilde5,_ = get_fd_waveform(approximant="TaylorF2", mass1=10.1, mass2=10.1, f_lower=15, delta_f=df)

bank_tilde1.resize(n)
bank_tilde2.resize(n)
bank_tilde3.resize(n)
bank_tilde4.resize(n)
bank_tilde5.resize(n)

bank_veto_bank = [bank_tilde1,bank_tilde2,bank_tilde3,bank_tilde4,bank_tilde5]
bank_veto_curr_overlaps = []
bank_snrs = []
bank_norms = []
for bank_template in bank_veto_bank:
    # For every bank veto template compute overlap between template
    # and the data
    curr_bank_snr,_,curr_bank_norm = matched_filter_core(bank_template,\
            strain,psd,low_frequency_cutoff=15)
    # SNR time series stored here
    bank_snrs.append(curr_bank_snr)
    # Template normalization factor stored here
    bank_norms.append(curr_bank_norm)
    bank_veto_curr_overlaps.append(overlap_cplx(htildep,\
             bank_template,psd=psd,\
             low_frequency_cutoff=15))
    test1,_,test2 = matched_filter_core(bank_template,\
            htildep,psd=psd,low_frequency_cutoff=15)
    sigmasq1 = sigmasq(htildep, psd, 15, None)
    sigmasq2 = sigmasq(bank_template, psd, 15, None)

bank_veto = bank_chisq_from_filters(snr,norm,bank_snrs,bank_norms,bank_veto_curr_overlaps)
#numpy.savetxt('BV_TEST2.txt',bank_veto)

