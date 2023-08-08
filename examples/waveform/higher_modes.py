import matplotlib.pyplot as pp
from pycbc.waveform import get_td_waveform

# Let's plot what our new waveform looks like
pp.figure()

# You can select sets of modes or individual modes using the 'mode_array'
for mode_select in [None,  "22 33", "21", "22", "33", "32", "44",]: 
    hp, hc = get_td_waveform(approximant="IMRPhenomXPHM",
                         mass1=7,
                         mass2=40,
                         f_lower=20.0,
                         mode_array=mode_select,
                         inclination = 1.0,
                         delta_t=1.0/4096)
                         

    
    if mode_select is None:
        label = 'Full Waveform'
        a = hp.max()
    else:
        label = "l,m=" + mode_select
    
    (hp / a).plot(label=label)
   
pp.xlim(-2, 0.3)
pp.legend()
pp.xlabel('Time [s]')
pp.ylabel('Relative Strain')
pp.show()
