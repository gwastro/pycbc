import pycbc.noise
import pycbc.psd
import pycbc.waveform
import pycbc.vetoes
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

# Calculate the power chisq (Bruce) chisq time series
hp, hc = pycbc.waveform.get_fd_waveform(approximant="SEOBNRv2_ROM_DoubleSpin",
                             mass1=25, mass2=25,
                             f_lower=flow, delta_f=stilde.delta_f)

hp.resize(len(stilde))
num_bins = 16
chisq = pycbc.vetoes.power_chisq(hp, stilde, num_bins, psd,
                                      low_frequency_cutoff=flow)

# convert to a reduced chisq
chisq /= (num_bins * 2) - 2

pylab.plot(chisq.sample_times, chisq)
pylab.legend()
pylab.ylabel('$\chi^2_r$')
pylab.xlabel('time (s)')
pylab.show()
