import matplotlib.pyplot as pp
from pycbc import waveform


for phase_order in [2, 3, 4, 5, 6, 7]:
    hp, hc = waveform.get_td_waveform(approximant='SpinTaylorT4',
                                 mass1=10, mass2=10,
                                 phase_order=phase_order,
                                 delta_t=1.0/4096,
                                 f_lower=100)

    hp, hc = hp.trim_zeros(), hc.trim_zeros()
    amp = waveform.utils.amplitude_from_polarizations(hp, hc)
    f = waveform.utils.frequency_from_polarizations(hp, hc)

    pp.plot(f.sample_times, f, label="PN Order = %s" % phase_order)

pp.ylabel('Frequency (Hz)')
pp.xlabel('Time (s)')
pp.legend(loc='upper left')
pp.show()
