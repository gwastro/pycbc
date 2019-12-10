#! /usr/bin/env python
# Copyright (C) 2012  Alex Nitz
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
from time import sleep

import sys
from numpy import loadtxt,complex64,float32
from optparse import OptionParser
from glue.ligolw import utils as ligolw_utils
from glue.ligolw import table, lsctables
from math import pow

from scipy.interpolate import interp1d

from pycbc.utils import mass1_mass2_to_mchirp_eta
from pycbc.waveform import get_td_waveform, get_fd_waveform, td_approximants, fd_approximants
from pycbc import DYN_RANGE_FAC
from pycbc.types import FrequencySeries, TimeSeries, zeros, real_same_precision_as, complex_same_precision_as
from pycbc.filter import match, sigmasq, resample_to_delta_t
from pycbc.scheme import DefaultScheme, CUDAScheme, OpenCLScheme
from pycbc.fft import fft
from math import cos, sin
import pycbc.psd

def update_progress(progress):
    print('\r\r[{0}] {1:.2%}'.format('#'*(int(progress*100)/2)+' '*(50-int(progress*100)/2), progress), end=' ')
    if progress == 100:
        print("Done")
    sys.stdout.flush()

## Remove the need for these functions ########################################

def generate_fplus_fcross(latitude,longitude,polarization):
    f_plus = - (1.0/2.0) * (1.0 + cos(latitude)*cos(latitude)) * cos (2.0 * longitude) * cos (2.0 * polarization) - cos(latitude) * sin(2.0*longitude) * sin (2.0 * polarization)
    f_cross=  (1.0/2.0) * (1.0 + cos(latitude)*cos(latitude)) * cos (2.0 * longitude) * sin (2.0* polarization) - cos (latitude) * sin(2.0*longitude) * cos (2.0 * polarization)
    return f_plus, f_cross

def generate_detector_strain(template_params, h_plus, h_cross):
    latitude = 0
    longitude = 0
    polarization = 0

    if hasattr(template_params, 'latitude'):
        latitude = template_params.latitude
    if hasattr(template_params, 'longitude'):
        longitude = template_params.longitude
    if hasattr(template_params, 'polarization'):
        polarization = template_params.polarization

    f_plus, f_cross = generate_fplus_fcross(latitude, longitude, polarization)

    return (f_plus*h_plus+f_cross*h_cross)

def make_padded_frequency_series(vec,filter_N=None):
    """Pad a TimeSeries with a length of zeros greater than its length, such
    that the total length is the closest power of 2. This prevents the effects
    of wraparound.
    """
    if filter_N is None:
        power = ceil(log(len(vec),2))+1
        N = 2 ** power
    else:
        N = filter_N
    n = N/2+1


    if isinstance(vec,FrequencySeries):
        vectilde = FrequencySeries(zeros(n, dtype=complex_same_precision_as(vec)),
                                   delta_f=1.0,copy=False)
	if len(vectilde) < len(vec):
	    cplen = len(vectilde)
        else:
            cplen = len(vec)
        vectilde[0:cplen] = vec[0:cplen]
        delta_f = vec.delta_f


    if isinstance(vec,TimeSeries):
        vec_pad = TimeSeries(zeros(N),delta_t=vec.delta_t,
                         dtype=real_same_precision_as(vec))
        vec_pad[0:len(vec)] = vec
        delta_f = 1.0/(vec.delta_t*N)
        vectilde = FrequencySeries(zeros(n),delta_f=1.0,
                               dtype=complex_same_precision_as(vec))
        fft(vec_pad,vectilde)

    vectilde = FrequencySeries(vectilde * DYN_RANGE_FAC,delta_f=delta_f,dtype=complex64)
    return vectilde

def get_waveform(approximant, phase_order, amplitude_order, template_params, start_frequency, sample_rate, length):

    if approximant in td_approximants():
        hplus,hcross = get_td_waveform(template_params, approximant=approximant,
                                   phase_order=phase_order, delta_t=1.0 / sample_rate,
                                   f_lower=start_frequency, amplitude_order=amplitude_order)
        hvec = generate_detector_strain(template_params, hplus, hcross)

    elif approximant in fd_approximants():
        delta_f = sample_rate / length
        hvec = get_fd_waveform(template_params, approximant=approximant,
                               phase_order=phase_order, delta_f=delta_f,
                               f_lower=start_frequency, amplitude_order=amplitude_order)


    htilde = make_padded_frequency_series(hvec,filter_N)

    return htilde

###############################################################################




#File output Settings
parser = OptionParser()
parser.add_option("--match-file", dest="out_file", help="file to output match results", metavar="FILE")

#PSD Settings
parser.add_option("--asd-file", dest="asd_file", help="two-column ASCII file containing ASD data", metavar="FILE")
parser.add_option("--psd", dest="psd", help="Analytic PSD model from LALSimulation", choices=pycbc.psd.get_lalsim_psd_list())

aprs = list(set(td_approximants() + fd_approximants()))
#Template Settings
parser.add_option("--template-file", dest="bank_file", help="SimInspiral or SnglInspiral XML file containing the template parameters.", metavar="FILE")
parser.add_option("--template-approximant",help="Template Approximant Name: " + str(aprs), choices = aprs)
parser.add_option("--template-phase-order",help="PN order to use for the phase",default=-1,type=int)
parser.add_option("--template-amplitude-order",help="PN order to use for the amplitude",default=-1,type=int)
parser.add_option("--template-start-frequency",help="Starting frequency for injections",type=float)
parser.add_option("--template-sample-rate",help="Starting frequency for injections",type=float)

#Signal Settings
parser.add_option("--signal-file", dest="sim_file", help="SimInspiral or SnglInspiral XML file containing the signal parameters.", metavar="FILE")
parser.add_option("--signal-approximant",help="Signal Approximant Name: " + str(aprs), choices = aprs)
parser.add_option("--signal-phase-order",help="PN order to use for the phase",default=-1,type=int)
parser.add_option("--signal-amplitude-order",help="PN order to use for the amplitude",default=-1,type=int)
parser.add_option("--signal-start-frequency",help="Starting frequency for templates",type=float)
parser.add_option("--signal-sample-rate",help="Starting frequency for templates",type=float)

#Filtering Settings
parser.add_option('--filter-low-frequency-cutoff', metavar='FREQ', help='low frequency cutoff of matched filter', type=float)
parser.add_option("--filter-sample-rate",help="Filter Sample Rate [Hz]",type=float)
parser.add_option("--filter-signal-length",help="Length of signal for filtering, shoud be longer than all waveforms and include some padding",type=int)

#Hardware support
parser.add_option("--use-cuda",action="store_true")

#Restricted maximization
parser.add_option("--mchirp-window",type=float)
(options, args) = parser.parse_args()

template_sample_rate = options.filter_sample_rate
signal_sample_rate = options.filter_sample_rate

if options.template_sample_rate:
    template_sample_rate = options.template_sample_rate

if options.signal_sample_rate:
    template_sample_rate = options.signal_sample_rate

if options.psd and options.asd_file:
    parser.error("PSD and asd-file options are mututally exclusive")

if options.use_cuda:
    ctx = CUDAScheme()
else:
    ctx = DefaultScheme()

print("STARTING THE BANKSIM")

# Load in the template bank file
indoc = ligolw_utils.load_filename(options.bank_file, False)
try :
    template_table = table.get_table(indoc, lsctables.SnglInspiralTable.tableName)
except ValueError:
    template_table = table.get_table(indoc, lsctables.SimInspiralTable.tableName)

# open the output file where the max overlaps over the bank are stored
fout = open(options.out_file, "w")
fout2 = open(options.out_file+".found", "w")

print("Writing matches to " + options.out_file)
print("Writing recovered template in " + options.out_file+".found")

# Load in the simulation list
indoc = ligolw_utils.load_filename(options.sim_file, False)
try:
    signal_table = table.get_table(indoc, lsctables.SimInspiralTable.tableName)
except ValueError:
    signal_table = table.get_table(indoc, lsctables.SnglInspiralTable.tableName)

def outside_mchirp_window(template,signal,w):
    template_mchirp,et = mass1_mass2_to_mchirp_eta(template.mass1,template.mass2)
    signal_mchirp ,et  = mass1_mass2_to_mchirp_eta(signal.mass1,signal.mass2)
    if abs(signal_mchirp - template_mchirp) > (w*signal_mchirp) :
        return True
    else :
        False

filter_N = int(options.filter_signal_length * options.filter_sample_rate)
filter_n = filter_N / 2 + 1
filter_delta_f = 1.0 / float(options.filter_signal_length)

print("Number of Signal Waveforms: ",len(signal_table))
print("Number of Templates       : ",len(template_table))


print("Reading and Interpolating PSD")
if options.asd_file:
    psd = pycbc.psd.read.from_txt(options.asd_file, filter_n, filter_delta_f,
                           options.filter_low_frequency_cutoff)
elif options.psd:
    psd = pycbc.psd.analytic.from_string(options.psd, filter_n, filter_delta_f,
                           options.filter_low_frequency_cutoff)

psd *= DYN_RANGE_FAC **2
psd = FrequencySeries(psd,delta_f=psd.delta_f,dtype=float32)

with ctx:
    print("Pregenerating Signals")
    signals = []
    index = 0
    for signal_params in signal_table:
        index += 1
        update_progress(index/len(signal_table))
        stilde = get_waveform(options.signal_approximant,
                      options.signal_phase_order,
                      options.signal_amplitude_order,
                      signal_params,
                      options.signal_start_frequency,
                      options.filter_sample_rate,
                      filter_N)
        s_norm = sigmasq(stilde, psd=psd,
                          low_frequency_cutoff=options.filter_low_frequency_cutoff)
        stilde /= psd
        signals.append( (stilde, s_norm, [], signal_params) )

    print("Calculating Overlaps")
    index = 0
    # Calculate the overlaps
    for template_params in template_table:
        index += 1
        update_progress(float(index)/len(template_table))
        h_norm = htilde = None
        for stilde, s_norm, matches, signal_params in signals:
            # Check if we need to look at this template
            if options.mchirp_window and outside_mchirp_window(template_params,
                                        signal_params, options.mchirp_window):
                matches.append(-1)
                continue

            # Generate htilde if we haven't already done so
            if htilde is None:
                htilde = get_waveform(options.template_approximant,
                                      options.template_phase_order,
                                      options.template_amplitude_order,
                                      template_params,
                                      options.template_start_frequency,
                                      options.filter_sample_rate,
                                      filter_N)
                h_norm = sigmasq(htilde, psd=psd,
                       low_frequency_cutoff=options.filter_low_frequency_cutoff)

            o,i = match(htilde, stilde, h_norm=h_norm, s_norm=s_norm,
                     low_frequency_cutoff=options.filter_low_frequency_cutoff)
            matches.append(o)


#Find the maximum overlap in the bank and output to a file
for stilde, s_norm, matches, sim_template in signals:
    match_str= "%5.5f \n" % (max(matches))
    match_str2="  "+options.bank_file+" "+str(matches.index(max(matches)))+"\n"
    fout.write(match_str)
    fout2.write(match_str2)

