from pycbc.waveform import get_td_waveform
from pycbc.filter import match, overlap
from pycbc.psd import aLIGOZeroDetHighPower

# Buffer size in seconds. This is presumed to be
# longer than the longest waveform. 
time_buffer = 4

f_low = 30
sample_rate = 4096

# Length of corresponding time series and frequency series
tlen = sample_rate * time_buffer
flen = tlen / 2 + 1

delta_t = 1.0 / sample_rate
delta_f = 1.0 / time_buffer

print "Generating waveform 1"
hp, hc = get_td_waveform(approximant="EOBNRv2",
                         mass1=10,
                         mass2=10,
                         f_lower=f_low,
                         delta_t=1.0/4096)
print "waveform is %s seconds long" % hp.duration

print "Generating waveform 2"
sp, sc = get_td_waveform(approximant="TaylorT4",
                         mass1=10,
                         mass2=10,
                         f_lower=f_low, 
                         delta_t=1.0/4096)
                         
print "waveform is %s seconds long" % sp.duration

# Ensure that the waveforms are resized to the same length
sp.resize(tlen)
hp.resize(tlen)

print "Calculating analytic PSD"
psd = aLIGOZeroDetHighPower(flen, delta_f, f_low) 

print "Calculating match and overlap"
# Note: This takes a while the first time as an FFT plan is generated
# subsequent calls within the same program will be faster
m, i = match(hp, sp, psd=psd, low_frequency_cutoff=f_low)
o = overlap(hp, sp, psd=psd, low_frequency_cutoff=f_low)
print "Overlap %s" % o
print "Maximized Overlap %s" % m

