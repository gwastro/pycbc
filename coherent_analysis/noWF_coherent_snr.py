from pycbc.waveform import get_td_waveform
from pycbc.filter import matched_filter
from pycbc.filter import sigma
from pycbc.psd import aLIGOZeroDetHighPower
import pycbc.detector
import pycbc.types
import pycbc.fft
import math
import fstat 
import numpy
from strain_generator import strain_generator
from zeros_generator import zeros_generator
import sys


if len(sys.argv[1:]) == 0:
    strainType = "waveform"
elif sys.argv[1:][0] in ["waveform", "noise", "zeros"]:
    strainType = sys.argv[1:][0]
else:
    print "Please use strain type: waveform, noise, zeros"
    sys.exit()

print "Strain type used: " + str(strainType)


#The list of IFOs used in our analysis
ifoList = ['H1','L1']

#####################
# Define the filter and waveform parameters 
#####################
f_low = 30
sample_rate = 16384
pol = 0.0 #Choose a ploarisation angle
#Make up sky location and time
dec = 0.5
ra = 0.5
t_gps = 1000

"""
dist = 300
inc = math.pi/6
pol = 0.0 #math.pi/4
phase = 0.5
"""

###########
# Generate the filter
###########
# Generate a template (doesn't matter about distance, inclination, etc)
sp, sc = get_td_waveform(approximant="TaylorT4",
                mass1=1.4,
                mass2=1.4,
                f_lower=f_low,
                delta_t=1.0/sample_rate,
                distance=1)

####################
# Make our trigger
####################
# Strain is our trigger, later we will replace this with actual IFO data
if strainType == "waveform":
    strain = strain_generator()
elif strainType == "zeros":
    strain = zeros_generator()
#Add noise test

sp.prepend_zeros(len(strain['H1']) - len(sp))


#---------Now we have made the waveforms, we can do the data analysis stuff-----


#I'm not sure why, but we resize the template
newlen = int(2**(math.ceil(math.log(len(strain['H1']), 2))))
sp.resize(newlen)
sc.resize(newlen)

#####################
# Get the IFO response for this trigger
#####################

# Generate the aLIGO ZDHP PSD for each IFO
delta_f = 1.0 / sp.duration
flen = newlen/2 + 1
psd_dictionary = {}
for key in ifoList:
    psd_dictionary[key] = aLIGOZeroDetHighPower(flen, delta_f, f_low)

# filter the data 
#sig= sigma(sp, psd_dictionary['H1'], low_frequency_cutoff=30., high_frequency_cutoff=1000.)

#Get the antenna functions
ifo_dictionary = {}
for key in ifoList:
    ifo_dictionary[key] = pycbc.detector.Detector(key)


#Now get the antenna responses for our template for the sky location and time we just made up
Fp_dictionary = {}
Fc_dictionary = {}
for key in ifoList:
    Fp_dictionary[key], Fc_dictionary[key] = ifo_dictionary[key].antenna_pattern(ra, dec, pol, t_gps)

#############################################################
# Find the expected SNR in each detector to test our code (expected = coinc, our code = coherent, should be same for 2 detectors)
expected_snr_dictionary = {}
for key in ifoList:
    expected_snr_dictionary[key] = pycbc.filter.sigma(strain[key], psd=psd_dictionary[key], low_frequency_cutoff=f_low)
    print("Expected SNR for %s is %.4f" % (key, expected_snr_dictionary[key]) )
#############################################################

#Find sigma
sigma_dictionary = {}
for key in ifoList:
    sigma_dictionary[key] = pycbc.filter.sigma(sp, psd=psd_dictionary[key], low_frequency_cutoff=f_low, high_frequency_cutoff=None)

#########################################
#Used for tests
f_sig = numpy.array([sigma_dictionary['H1'] * numpy.array([Fp_dictionary['H1'],Fc_dictionary['H1']]),sigma_dictionary['L1'] * numpy.array([Fp_dictionary['L1'],Fc_dictionary['L1']])])
###########################################


#Change to DP frame
chi_dp = fstat.dominant_polarization(f_sig)[2]

"""
#Different way of calculating chi
tan_chi_top = 0
tan_chi_bottom = 0
for key in ifoList:
    tan_chi_top += 2*(sigma_dictionary[key] * Fp_dictionary[key])*(sigma_dictionary[key] * Fc_dictionary[key])
    tan_chi_bottom += ((sigma_dictionary[key] * Fp_dictionary[key])**2 - (sigma_dictionary[key] * Fc_dictionary[key])**2)
chi = numpy.arctan(tan_chi_top / tan_chi_bottom)/4 
chi_dp = chi
"""

Fp_dp_dictionary = {}
Fc_dp_dictionary = {}
for ifo in ifoList:
    Fp_dp_dictionary[ifo] = Fp_dictionary[ifo] * numpy.cos(2*chi_dp) + Fc_dictionary[ifo] * numpy.sin(2*chi_dp)
    Fc_dp_dictionary[ifo] = -Fp_dictionary[ifo] * numpy.sin(2*chi_dp) + Fc_dictionary[ifo] * numpy.cos(2*chi_dp)

#######################################
# Now calculate the variables needed for the rho_coh calulcation 
#######################################
#calculate f (The orthognal unit vectors in detector space which are used to calculate rho_coh)
fp_denominator0 = 0
fc_denominator0 = 0
for key in ifoList:
    fp_denominator0 += (sigma_dictionary[key] * Fp_dp_dictionary[key])**2
    fc_denominator0 += (sigma_dictionary[key] * Fc_dp_dictionary[key])**2
fp_denominator = math.sqrt(fp_denominator0)
fc_denominator = math.sqrt(fc_denominator0)

fp_dictionary = {}
fc_dictionary = {}
for key in ifoList:
    fp_dictionary[key] = (sigma_dictionary[key] * Fp_dp_dictionary[key]) / fp_denominator 
    fc_dictionary[key] = (sigma_dictionary[key] * Fc_dp_dictionary[key]) / fc_denominator 

#Check fc and fp

#print "fp dict is"
#print(fp_dictionary)
#print "fc dict is"
#print(fc_dictionary)

#Complex snrs for the two detectors (Used in formula for rho_coh)
cmplx_snr_dictionary = {}
for key in ifoList:
    cmplx_snr_dictionary[key] = matched_filter(sp, strain[key], psd=psd_dictionary[key], low_frequency_cutoff=f_low)

#######################################
#Calculate the coherent SNR
#######################################
rho_coh_sq = 0
for key1 in ifoList:
    for key2 in ifoList:
        rho_coh_sq += cmplx_snr_dictionary[key1].conj() * (fp_dictionary[key1]*fp_dictionary[key2] + fc_dictionary[key1]*fc_dictionary[key2]) * cmplx_snr_dictionary[key2]

rho_coh = numpy.sqrt(rho_coh_sq)

############################
# CHECKS!!
############################

imax = numpy.argmax(abs(rho_coh))
print("Maximum coherent SNR=%.4f" % max(abs(rho_coh)))
print("Argument of Maximum coherent SNR=%d" % imax)
for ifo in ifoList:
    print("SNR in ifo %s is %.4f" % (ifo, abs(cmplx_snr_dictionary[ifo][imax])) )

snr_array = numpy.array([cmplx_snr_dictionary['H1'][imax],cmplx_snr_dictionary['L1'][imax]])

recovered_a = fstat.snr_f_to_a(snr_array.conj(), f_sig)
d, cosi, psi, phi = fstat.a_to_params(recovered_a)

print("Recovered F-stat A values")
print fstat.snr_f_to_a(snr_array, f_sig)

print("Recovered parameters: D = %.2f Mpc; i = %.2f; psi = %.2f; phi = %.2f" 
    % (d, numpy.arccos(cosi), psi, phi))

