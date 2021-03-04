""" Waveform approximants for the pre-merger detection of gravitational waves
"""
import logging


def premerger_taylorf2(**p):
    """ Generate time-shifted TaylorF2"""
    from pycbc.waveform import get_fd_waveform
    from pycbc.waveform.spa_tmplt import spa_length_in_time
    from pycbc.waveform.utils import fd_taper

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

    logging.info("PreTaylorF2, m1=%.1f, m2=%.1f, fmax=%.1f, timeshift=%.1f",
                 p['mass1'], p['mass2'], p['f_final'], removed)
    kmin = int(p['f_lower'] / p['delta_f'])
    hp[0:kmin] = 0
    hc[0:kmin] = 0

    if 'final_taper' in p:
        taper_size = p['final_taper']
        hp = fd_taper(hp, p['f_final'] - taper_size, p['f_final'], side='right')
        hc = fd_taper(hc, p['f_final'] - taper_size, p['f_final'], side='right')

    hp.time_offset = removed
    hc.time_offset = removed

    return hp, hc
