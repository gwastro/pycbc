import pycbc.noise
import pycbc.psd
import pycbc.filter
import pycbc.waveform
import pylab

# Generate some noise with an advanced ligo psd
flow = 30.0
delta_f = 1.0 / 16
flen = int(2048 / delta_f) + 1
psd = pycbc.psd.aLIGOZeroDetHighPower(flen, delta_f, flow)

# Generate 16 seconds of noise at 4096 Hz
delta_t = 1.0 / 4096
tsamples = int(16 / delta_t)
strain = pycbc.noise.noise_from_psd(tsamples, delta_t, psd, seed=127)
stilde = strain.to_frequencyseries()

# Use a waveform as a matched filter
hp, hc = pycbc.waveform.get_fd_waveform(approximant="SEOBNRv2_ROM_DoubleSpin",
                             mass1=25, mass2=25,
                             f_lower=flow, delta_f=stilde.delta_f)

hp.resize(len(stilde))
snr = pycbc.filter.matched_filter(hp, stilde, psd=psd,
                                      low_frequency_cutoff=flow)


pylab.plot(snr.sample_times, abs(snr))
pylab.ylabel('signal-to-noise ratio')
pylab.xlabel('time (s)')
pylab.show()
