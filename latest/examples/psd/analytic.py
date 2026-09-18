import matplotlib.pyplot as pp
import pycbc.psd


# List the available analytic psds
print(pycbc.psd.get_lalsim_psd_list())

delta_f = 1.0 / 4
flen = int(1024 / delta_f)
low_frequency_cutoff = 30.0

# One can either call the psd generator by name
p1 = pycbc.psd.aLIGOZeroDetHighPower(flen, delta_f, low_frequency_cutoff)

# or by using the name as a string.
p2 = pycbc.psd.from_string('aLIGOZeroDetLowPower', flen, delta_f, low_frequency_cutoff)

pp.plot(p1.sample_frequencies, p1, label='HighPower')
pp.plot(p2.sample_frequencies, p2, label='LowPower')
pp.legend()
pp.show()
