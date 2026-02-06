import matplotlib.pyplot as pp
from pycbc.waveform import get_td_waveform

# Let's plot what our new waveform looks like
pp.figure()

# You can select sets of modes or individual modes using the 'mode_array'
# The standard format is to provide a list of (l, m) modes, however
# a string format is also provided to aid use in population from config files.
# e.g. "22 33" is also acceptable to select these two modes.
# "None" will result in the waveform return its default which is usually
# to return all implemented modes.
for mode_select in [None,  
                    [(2, 2), (3, 3)], # Select two modes at once
                    [(2, 2)],
                    [(2, 1)], 
                    [(3, 2)],
                    [(4, 4)],
                   ]: 
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
        label = "l, m = " + '  '.join([f"{l}, {m}" for l, m in mode_select])
    
    (hp / a).plot(label=label)
   
pp.xlim(-1, 0.05)
pp.legend()
pp.xlabel('Time [s]')
pp.ylabel('Relative Strain')
pp.show()
