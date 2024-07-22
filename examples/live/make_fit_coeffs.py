"""
Makes files which can be used as the fit_coeffs statistic.
These are not of any scientific use, but the code will accept them
and run properly
"""

import numpy as np
from pycbc.io.hdf import HFile

# Get number of templates from bank file
with HFile('template_bank.hdf', 'r') as bankf:
    n_templates = bankf['mass1'].size

for ifo in ['H1','L1','V1']:
    with HFile(f'{ifo}-fit_coeffs.hdf','w') as fits_f:
        fits_f.attrs['analysis_time'] = 430000
        fits_f.attrs['ifo'] = ifo
        fits_f.attrs['stat'] = f'{ifo}-fit_coeffs'
        fits_f.attrs['stat_threshold'] = 5

        fits_f['count_above_thresh'] = np.ones(n_templates) * 100
        fits_f['count_in_template'] = np.ones(n_templates) * 20000
        fits_f['fit_coeff'] = np.ones(n_templates) * 5.5
        fits_f['median_sigma'] = np.ones(n_templates) * 5800
        fits_f['template_id'] = np.arange(n_templates)

    
