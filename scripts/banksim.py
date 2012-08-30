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
from pycbc.filter import get_cutoff_indices, match, sigmasq
from pycbc.scheme import DefaultScheme, CUDAScheme, OpenCLScheme
from pycbc.fft import fft
import pycbc.psd 

def update_progress(progress):
    print '\r\r[{0}] {1}%'.format('#'*(progress/2)+' '*(50-progress/2), progress),
    if progress == 100:
        print "Done"
    sys.stdout.flush()

## Remove the need for these functions ########################################
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

def get_waveform(approximant, order, template_params, start_frequency, sample_rate, length):

    if approximant in fd_approximants():
        delta_f = sample_rate / length
        hvec = get_fd_waveform(template_params, approximant=approximant,
                               phase_order=order, delta_f=delta_f,
                               f_lower=start_frequency, amplitude_order=order)     

    if approximant in td_approximants():
        hplus,hcross = get_td_waveform(template_params, approximant=approximant,
                                   phase_order=order, delta_t=1.0 / sample_rate,
                                   f_lower=start_frequency, amplitude_order=order) 
        hvec = hplus

    htilde = make_padded_frequency_series(hvec,filter_N)

    return htilde

###############################################################################


print "STARTING THE BANKSIM"

#File I/O Settings
parser = OptionParser()
parser.add_option("--signal-file", dest="sim_file", help="simulated waveform parameters", metavar="FILE")
parser.add_option("--template-file", dest="bank_file", help="template bank parameters", metavar="FILE")
parser.add_option("--asd-file", dest="asd_file",help="ASD data", metavar="FILE")
parser.add_option("--match-file", dest="out_file",help="file to output match results", metavar="FILE")

#Waveform generation Settings
parser.add_option("--template-approximant",help="Template Approximant Name")
parser.add_option("--template-order",help="PN order to use for the aproximant",default=-1,type=int) 
parser.add_option("--template-start-frequency",help="Starting frequency for injections",type=float) 

parser.add_option("--signal-approximant",help="Singal Approximant Name"  "[SPA]")
parser.add_option("--signal-order",help="PN order to use for the aproximant",default=-1,type=int)  
parser.add_option("--signal-start-frequency",help="Starting frequency for templates",type=float)  

#Filter Settings
parser.add_option('--filter-low-frequency-cutoff', metavar='FREQ', help='low frequency cutoff of matched filter', type=float)
parser.add_option("--filter-sample-rate",help="Filter Sample Rate [Hz]",type=float)
parser.add_option("--filter-signal-length",help="Length of signal for filtering, shoud be longer than all waveforms and include some padding",type=int)

parser.add_option("--cuda",action="store_true")
parser.add_option("--mchirp-window",type=float)               
(options, args) = parser.parse_args()   

if options.cuda:
    ctx = CUDAScheme()
else:
    ctx = DefaultScheme()

# Load in the template bank file
indoc = ligolw_utils.load_filename(options.bank_file, False)
try :
    template_table = table.get_table(indoc, lsctables.SnglInspiralTable.tableName) 
except ValueError:
    template_table = table.get_table(indoc, lsctables.SimInspiralTable.tableName)

# open the output file where the max overlaps over the bank are stored 
fout = open(options.out_file, "w")
fout2 = open(options.out_file+".found", "w")

print "Writing matches to " + options.out_file
print "Writing recovered template in " + options.out_file+".found"

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

print("Number of Signal Waveforms: ",len(signal_table))
print("Number of Templates       : ",len(template_table))

print("Pregenerating Signals")
signals = []
# Pregenerate the waveforms to simulate 

index = 0 
for signal_params in signal_table:
    index += 1
    update_progress(index*100/len(signal_table))
    stilde = get_waveform(options.signal_approximant, 
                  options.signal_order, 
                  signal_params, 
                  options.signal_start_frequency, 
                  options.filter_sample_rate, 
                  filter_N)
    signals.append((stilde,[],signal_params))


print("Reading and Interpolating PSD")
# Load the asd file
psd = pycbc.psd.from_txt(options.asd_file, len(stilde),  stilde.delta_f, options.filter_low_frequency_cutoff )
psd *= DYN_RANGE_FAC **2
psd = FrequencySeries(psd,delta_f=psd.delta_f,dtype=float32) 

print("Calculating Overlaps")
with ctx:
    index = 0 
    # Calculate the overlaps
    for template_params in template_table:
        index += 1
        update_progress(index*100/len(template_table))
        htilde = None
        for stilde,matches,signal_params in signals:
            # Check if we need to look at this template
            if options.mchirp_window and outside_mchirp_window(template_params,signal_params,options.mchirp_window):
                matches.append(-1)
                continue

            # Generate htilde if we haven't already done so 
            if htilde is None:
                htilde = get_waveform(options.template_approximant, 
                                      options.template_order, 
                                      template_params, 
                                      options.template_start_frequency, 
                                      options.filter_sample_rate, 
                                      filter_N)

            o,i = match(htilde, stilde, psd=psd, low_frequency_cutoff=options.filter_low_frequency_cutoff)         
            matches.append(o)
         

#Find the maximum overlap in the bank and output to a file
for stilde,matches,sim_template in signals:
    match_str= "%5.5f \n" % (max(matches))
    match_str2="  "+options.bank_file+" "+str(matches.index(max(matches)))+"\n"
    fout.write(match_str)
    fout2.write(match_str2)

