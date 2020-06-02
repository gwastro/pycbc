#!/usr/bin/env python

import sys
import glob
import logging as log
import numpy as np
import h5py
import pycbc
from glue.ligolw import utils as ligolw_utils
from glue.ligolw import ligolw, table, lsctables


def close(a, b):
    return abs(a - b) < 1e-7


log.basicConfig(level=log.INFO, format='%(asctime)s %(message)s')

tested_detectors = {'H1', 'L1', 'V1'}
sim_gps_start = 1272790000
sim_gps_end =  1272790500
sim_f_lower = 18

with h5py.File('template_bank.hdf', 'r') as bankf:
    temp_mass1 = bankf['mass1'][:]
    temp_mass2 = bankf['mass2'][:]
    temp_s1z = bankf['spin1z'][:]
    temp_s2z = bankf['spin2z'][:]

detectors_with_trigs = set()
fail = False

log.info('Starting Test')

trig_paths = sorted(glob.glob('output/????_??_??/*.hdf'))
for trigfp in trig_paths:
    with h5py.File(trigfp, 'r') as trigf:
        for detector in tested_detectors:
            if detector not in trigf:
                continue
            group = trigf[detector]

            if 'psd' not in group:
                continue

            # check that PSD is sane
            psd = group['psd'][:] / pycbc.DYN_RANGE_FAC ** 2
            psd_df = group['psd'].attrs['delta_f']
            psd_f = np.arange(len(psd)) * psd_df
            psd_epoch = group['psd'].attrs['epoch']
            in_band_asd = psd[psd_f > sim_f_lower] ** 0.5
            if len(in_band_asd) == 0 or (in_band_asd < 1e-24).any() \
                    or (in_band_asd > 1e-20).any() \
                    or not np.isfinite(in_band_asd).all() \
                    or psd_epoch < sim_gps_start or psd_epoch > sim_gps_end:
                log.info('Invalid PSD in %s %s', trigfp, detector)
                fail = True

            if 'snr' not in group or len(group['snr']) == 0:
                continue

            detectors_with_trigs.add(detector)

            # check that SNR is non-negative and finite
            snr = group['snr'][:]
            if (snr < 0).any() or not np.isfinite(snr).all():
                log.error('Invalid SNR in %s %s', trigfp, detector)
                fail = True

            # check that Allen chi^2 is non-negative and finite
            chisq = group['chisq'][:]
            chisq_dof = group['chisq_dof'][:]
            if (chisq < 0).any() or not np.isfinite(chisq).all() \
                    or (chisq_dof < 0).any() or not np.isfinite(chisq_dof).all():
                log.error('Invalid Allen chi^2 in %s %s', trigfp, detector)
                fail = True

            # check that merger time is within the simulated time range
            end_time = group['end_time'][:]
            if (end_time < sim_gps_start).any() or (end_time > sim_gps_end).any():
                log.error('Invalid merger time in %s %s', trigfp, detector)
                fail = True

            # check that template parameters are consistent with the bank
            trig_mass1 = group['mass1'][:]
            trig_mass2 = group['mass2'][:]
            trig_s1z = group['spin1z'][:]
            trig_s2z = group['spin2z'][:]
            trig_temp_id = group['template_id'][:]
            test_mass1 = close(temp_mass1[trig_temp_id], trig_mass1)
            test_mass2 = close(temp_mass2[trig_temp_id], trig_mass2)
            test_s1z = close(temp_s1z[trig_temp_id], trig_s1z)
            test_s2z = close(temp_s2z[trig_temp_id], trig_s2z)
            test_all = np.logical_and.reduce((test_mass1, test_mass2,
                                              test_s1z, test_s2z))
            if not test_all.all():
                log.error('Invalid template parameters in %s %s',
                              trigfp, detector)
                fail = True

# check that triggers were produced in all detectors
if detectors_with_trigs != tested_detectors:
    missing = sorted(tested_detectors - detectors_with_trigs)
    log.error('No triggers found in %s', ', '.join(missing))
    fail = True
    
# dummy class needed for loading LIGOLW files
class LIGOLWContentHandler(ligolw.LIGOLWContentHandler):
    pass

lsctables.use_in(LIGOLWContentHandler)
    
# check properties of coincident triggers
coinc_trig_paths = sorted(glob.glob('output/coinc*.xml.gz'))
l=len(coinc_trig_paths)
if l==0:
    log.error('No coincident triggers detected')
    fail = True
elif l>=10: 
    log.error('Too many coincident triggers detected')
    fail = True
else: log.info(str(l)+' coincident trigger(s) detected')


for ctrigfp in coinc_trig_paths:
    xmldoc = ligolw_utils.load_filename(
            ctrigfp, False, contenthandler=LIGOLWContentHandler)
    sngl_inspiral_table = lsctables.SnglInspiralTable.get_table(xmldoc)
    log.info('acquired table')
    new_snr = sngl_inspiral_table.get_new_snr()  
    log.info('New SNR'+str(new_snr))    
    log.info('finished test')
  

if fail:
    log.error('Test Failed')
else:
    log.info('Test Passed')
 
sys.exit(1 if fail else 0)
