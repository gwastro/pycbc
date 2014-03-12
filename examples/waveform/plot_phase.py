import pylab
from pycbc import waveform

for apx in ['EOBNRv2', 'TaylorT4', 'IMRPhenomB']:
    hp, hc = waveform.get_td_waveform(approximant=apx,
                                 mass1=10,
                                 mass2=10,
                                 delta_t=1.0/4096,
                                 f_lower=40)
                                 
    hp, hc = hp.trim_zeros(), hc.trim_zeros()
    amp = waveform.utils.amplitude_from_polarizations(hp, hc)
    phase = waveform.utils.phase_from_polarizations(hp, hc)
    
    pylab.plot(phase, amp, label=apx)
    
pylab.ylabel('GW Strain Amplitude')
pylab.xlabel('GW Phase (radians)')
pylab.legend(loc='upper left')
pylab.show()
