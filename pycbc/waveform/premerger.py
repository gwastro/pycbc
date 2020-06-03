""" Waveform approximants for the pre-merger detection of gravitational waves
"""

def premerger_taylorf2(**p):
    from pycbc.waveform import get_fd_waveform
    from pycbc.waveform.spa_tmplt import spa_length_in_time
    
    p.pop('approximant')
    hp, hc = get_fd_waveform(approximant="TaylorF2", **p)

    removed = spa_length_in_time(mass1=p['mass1'],
                                 mass2=p['mass2'],
                                 f_lower=p['f_final'],
                                 phase_order=-1)

    hp = hp.cyclic_time_shift(removed)
    hp.start_time += removed
    
    hc = hc.cyclic_time_shift(removed)
    hc.start_time += removed
    
    return hp, hc
