import matplotlib.pyplot as pp
import pycbc.noise
import pycbc.psd


# The color of the noise matches a PSD which you provide
flow = 30.0
delta_f = 1.0 / 64
flen = int(2048 / delta_f) + 1
psd = pycbc.psd.aLIGOZeroDetHighPower(flen, delta_f, flow)

# Here the noise is generated directly in the frequency domain so it matches
# the parameters of the PSD you give.
fs = pycbc.noise.frequency_noise_from_psd(psd, seed=127)

pp.loglog(fs.sample_frequencies, abs(fs**2.0), label='Magnitude')
pp.loglog(fs.sample_frequencies, abs(psd), label='Magnitude')
pp.xlim(flow, psd.sample_frequencies[-1])
pp.legend()
pp.ylabel('Strain^2 / Hz')
pp.xlabel('Frequency (Hz)')
pp.show()
