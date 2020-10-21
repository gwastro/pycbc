#!/usr/bin/env python

import sys
import glob
import logging as log
import numpy as np
import h5py
import pycbc
from pycbc.io import record
from glue.ligolw import utils as ligolw_utils
from glue.ligolw import ligolw, lsctables

# dummy class needed for loading LIGOLW files
class LIGOLWContentHandler(ligolw.LIGOLWContentHandler):
    pass
      
def close(a, b, c):
    return abs(a - b) <= c
   
def check_single_results(tested_detectors):
    single_fail = False
    with h5py.File('template_bank.hdf', 'r') as bankf:
        temp_mass1 = bankf['mass1'][:]
        temp_mass2 = bankf['mass2'][:]
        temp_s1z = bankf['spin1z'][:]
        temp_s2z = bankf['spin2z'][:]

    detectors_with_trigs = set()
      
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
                    single_fail = True

                if 'snr' not in group or len(group['snr']) == 0:
                    continue

                detectors_with_trigs.add(detector)

                # check that SNR is non-negative and finite
                snr = group['snr'][:]
                if (snr < 0).any() or not np.isfinite(snr).all():
                    log.error('Invalid SNR in %s %s', trigfp, detector)
                    single_fail = True

                # check that Allen chi^2 is non-negative and finite
                chisq = group['chisq'][:]
                chisq_dof = group['chisq_dof'][:]
                if (chisq < 0).any() or not np.isfinite(chisq).all() \
                        or (chisq_dof < 0).any() or not np.isfinite(chisq_dof).all():
                    log.error('Invalid Allen chi^2 in %s %s', trigfp, detector)
                    single_fail = True

                # check that merger time is within the simulated time range
                end_time = group['end_time'][:]
                if (end_time < sim_gps_start).any() or (end_time > sim_gps_end).any():
                    log.error('Invalid merger time in %s %s', trigfp, detector)
                    single_fail = True

                # check that template parameters are consistent with the bank
                trig_mass1 = group['mass1'][:]
                trig_mass2 = group['mass2'][:]
                trig_s1z = group['spin1z'][:]
                trig_s2z = group['spin2z'][:]
                trig_temp_id = group['template_id'][:]
                test_mass1 = close(temp_mass1[trig_temp_id], trig_mass1, 1e-7)
                test_mass2 = close(temp_mass2[trig_temp_id], trig_mass2, 1e-7)
                test_s1z = close(temp_s1z[trig_temp_id], trig_s1z, 1e-7)
                test_s2z = close(temp_s2z[trig_temp_id], trig_s2z, 1e-7)
                test_all = np.logical_and.reduce((test_mass1, test_mass2,
                                              test_s1z, test_s2z))
                if not test_all.all():
                    log.error('Invalid template parameters in %s %s',
                              trigfp, detector)
                    single_fail = True

    # check that triggers were produced in all detectors
    if detectors_with_trigs != tested_detectors:
        missing = sorted(tested_detectors - detectors_with_trigs)
        log.error('No triggers found in %s', ', '.join(missing))
        single_fail = True
      
    if single_fail:
        log.error('Single Trigger Test Failed')
    return single_fail


def check_coinc_results():
    coinc_fail = False    
    # gather coincident triggers
    coinc_trig_paths = sorted(glob.glob('output/coinc*.xml.gz'))
    n_coincs=len(coinc_trig_paths)
    if n_coincs==0:
        log.error('No coincident triggers detected')
        coinc_fail = True
    elif n_coincs>=10: 
        log.error('Too many coincident triggers detected')
        coinc_fail = True
    else: 
        log.info('%d coincident trigger(s) detected', n_coincs)
      
    injs = sorted(glob.glob('test_inj*.hdf'))
    n_injs = len(injs)
    inj_mass1 = np.empty(n_injs)
    inj_mass2 = np.empty(n_injs)
    inj_spin1z = np.empty(n_injs)
    inj_spin2z = np.empty(n_injs)
    inj_time = np.empty(n_injs)
    
    for idx, inj_path in enumerate(injs):
        with h5py.File(inj_path, 'r') as inj:
            inj_mass1[idx] = inj['mass1'][0]
            inj_mass2[idx] = inj['mass2'][0]
            inj_spin1z[idx] = inj['spin1z'][0]
            inj_spin2z[idx] = inj['spin2z'][0]
            inj_time[idx] = inj['tc'][0]
         
    
    if n_injs > n_coincs :
        log.error('More injections than coincident triggers')
        coinc_fail = True
              
      
    #create field array to store properties of triggers
    trig_props = record.FieldArray(n_coincs, dtype=[('mass1', float), \
                ('mass2', float), ('spin1z', float), ('spin2z', float), \
                ('tc', float), ('net_snr', float)])
      
    #store properties of coincident triggers
    for x, ctrigfp in enumerate(coinc_trig_paths):
        log.info('Checking trigger %s', ctrigfp)
        xmldoc = ligolw_utils.load_filename(
                  ctrigfp, False, contenthandler=LIGOLWContentHandler)
        sngl_inspiral_table = lsctables.SnglInspiralTable.get_table(xmldoc)
        
        trig_props['tc'][x]=sngl_inspiral_table.get_column('end_time')[0]                   
        trig_props['mass1'][x]=sngl_inspiral_table.get_column('mass1')[0]
        trig_props['mass2'][x]=sngl_inspiral_table.get_column('mass2')[0]    
        trig_props['spin1z'][x]=sngl_inspiral_table.get_column('spin1z')[0] 
        trig_props['spin2z'][x]=sngl_inspiral_table.get_column('spin2z')[0] 
        
        
        snr_list=sngl_inspiral_table.get_column('snr')
        network_snr = sum([snr ** 2 for snr in snr_list]) ** 0.5      
        trig_props['net_snr'][x]=network_snr
            
        
        log.info('IFO SNRs: '+str(snr_list))
        log.info('Network SNR: %f', network_snr)
        log.info('IFO End Time: %f', trig_props['tc'][x]) 
        log.info('Mass 1: %f', trig_props['mass1'][x])
        log.info('Mass 2: %f', trig_props['mass2'][x])
        log.info('Spin1z: %f', trig_props['spin1z'][x])
        log.info('Spin2z: %f', trig_props['spin2z'][x])
     

    # check if injections match trigger params
    for i in range(n_injs):
        has_match = False
        for j in range(n_coincs):
            if (close(inj_time[i],trig_props['tc'][j],1.0) and close(inj_mass1[i], trig_props['mass1'][j], 5e-7)
            and close(inj_mass2[i], trig_props['mass2'][j], 5e-7) and close(inj_spin1z[i], trig_props['spin1z'][j], 5e-7)
            and close(inj_spin2z[i], trig_props['spin2z'][j], 5e-7) and close(15.0, trig_props['net_snr'][j], 1.0)):
                has_match = True
                break
        
        if not has_match:
            coinc_fail = True
            log.error('Injection %i has no match',i)
            
    if coinc_fail:
        log.error('Coincident Trigger Test Failed')
    return coinc_fail

      

log.basicConfig(level=log.INFO, format='%(asctime)s %(message)s')

# gps times need to match those found in run.sh
sim_gps_start = 1272790000
sim_gps_end =  1272790500

# sim f lower needs to match the value in generate_injection.sh
sim_f_lower = 18


fail = False
lsctables.use_in(LIGOLWContentHandler)

       
tested_detectors = {'H1', 'L1', 'V1'}
   
single_fail = check_single_results(tested_detectors)
coinc_fail = check_coinc_results()
fail = single_fail or coinc_fail
   
if fail:
      log.error('Test Failed')
else:
      log.info('Test Passed')
 

sys.exit(1 if fail else 0)