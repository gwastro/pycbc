import pycbc.noise
import pycbc.psd
import pycbc.filter
import pylab

# Generate some noise with an advanced ligo psd
flow = 5.0
delta_f = 1.0 / 16
flen = int(2048 / delta_f) + 1
psd = pycbc.psd.aLIGOZeroDetHighPower(flen, delta_f, flow)

# Generate 1 seconds of noise at 4096 Hz
delta_t = 1.0 / 4096
tsamples = int(1 / delta_t)
ts = pycbc.noise.noise_from_psd(tsamples, delta_t, psd, seed=127)
pylab.plot(ts.sample_times, ts, label='Original')

# Suppress the low frequencies below 30 Hz
ts = pycbc.filter.highpass(ts, 30.0)
pylab.plot(ts.sample_times, ts, label='Highpassed')

# Suppress the high frequencies
ts = pycbc.filter.lowpass_fir(ts, 1000.0, 8)
pylab.plot(ts.sample_times, ts, label='Highpassed + Lowpassed')

pylab.legend()
pylab.ylabel('Strain')
pylab.xlabel('Time (s)')
pylab.show()
