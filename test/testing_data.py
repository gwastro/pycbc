'''
This file contains waveform data used in ensuring proper functionality of filter.matchedfilter.optimized_match as described in GitHub Issue #5144
'''
import os

def load_match_testing_waveforms_from_hdf5():
    from pycbc.io.hdf import HFile

    # Get the directory where this file is located
    data_dir = os.path.dirname(__file__)
    hdf5_path = os.path.join(data_dir, 'match_testing_waveforms.hdf5')

    with HFile(hdf5_path, 'r') as f:
        NL_fs = f['NL'][:]
        L_fs = f['L'][:]
        psd_fs = f['psd'][:]

    return L_fs, NL_fs, psd_fs