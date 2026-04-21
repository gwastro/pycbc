'''
This file contains waveform data used in ensuring proper functionality of filter.matchedfilter.optimized_match as described in GitHub Issue #5144
'''

def load_match_testing_waveforms_from_hdf5():
    from pycbc.io.hdf import HFile

    with HFile('match_testing_waveforms.hdf5', 'r') as f:
        NL_fs = f['NL'][:]
        L_fs = f['L'][:]
        psd_fs = f['psd'][:]

    return L_fs, NL_fs, psd_fs