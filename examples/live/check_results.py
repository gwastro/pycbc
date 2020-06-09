#!/usr/bin/env python

import sys
import glob
import logging as log
import numpy as np
import h5py
import pycbc
from glue.ligolw import utils as ligolw_utils
from glue.ligolw import ligolw, table, lsctables

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

def check_coinc_results(inj_table):
      coinc_fail = False
      print_template_parameters = False
      # gather coincident triggers
      coinc_trig_paths = sorted(glob.glob('output/coinc*.xml.gz'))
      l=len(coinc_trig_paths)
      if l==0:
            log.error('No coincident triggers detected')
            coinc_fail = True
      elif l>=10: 
            log.error('Too many coincident triggers detected')
            coinc_fail = True
      else: 
            log.info('%d coincident trigger(s) detected', l)
      
      #get properties of injection
      inj_mass1=inj_table.get_column('mass1')[0]
      inj_mass2=inj_table.get_column('mass2')[0]    
      inj_spin1z=inj_table.get_column('spin1z')[0] 
      inj_spin2z=inj_table.get_column('spin2z')[0] 
      inj_time=inj_table.get_column('end_time')[0]
      inj_snr=15
      
      #print properties of coincident triggers
      for ctrigfp in coinc_trig_paths:            
            xmldoc = ligolw_utils.load_filename(
                  ctrigfp, False, contenthandler=LIGOLWContentHandler)
            sngl_inspiral_table = lsctables.SnglInspiralTable.get_table(xmldoc)    
            end_time=sngl_inspiral_table.get_column('end_time')  
            snr_list=sngl_inspiral_table.get_column('snr')
            network_snr_squared=0
            for snr in snr_list:
                  network_snr_squared+=snr**2
            network_snr=np.sqrt(network_snr_squared)            
            mass1=sngl_inspiral_table.get_column('mass1')[0]
            mass2=sngl_inspiral_table.get_column('mass2')[0]    
            spin1z=sngl_inspiral_table.get_column('spin1z')[0] 
            spin2z=sngl_inspiral_table.get_column('spin2z')[0] 
            log.info('IFO SNRs: '+str(snr_list))
            log.info('Network SNR: '+str(network_snr))
            
            if print_template_parameters:
                  log.info('IFO End Times: '+str(end_time)) 
                  log.info('Mass 1: %f', mass1)
                  log.info('Mass 2: %f', mass2)
                  log.info('Spin1z: %f', spin1z)
                  log.info('Spin2z: %f', spin2z)
            
            #check if trigger properties match injection's properties
            for t in end_time:
                  if not close(t,inj_time,1.0):    
                        coinc_fail = True
                        log.error('Trigger time does not match injection time')  
            if not close(mass1, inj_mass1, 5e-7):
                  coinc_fail = True
                  log.error('Trigger mass 1 does not match injection mass 1')
            if not close(mass2, inj_mass2, 5e-7):
                  coinc_fail = True                  
                  log.error('Trigger mass 2 does not match injection mass 2')
            if not close(spin1z, inj_spin1z, 5e-7):
                  coinc_fail = True
                  log.error('Trigger spin1z does not match injection spin1z')
            if not close(spin2z, inj_spin2z, 5e-7):
                  coinc_fail = True
                  log.error('Trigger spin2z does not match injection spin1z')   
            if not close(network_snr, inj_snr, 0.1):
                  coinc_fail = True
                  log.error('Network SNR test failed')
            
      if coinc_fail:
            log.error('Coincident Trigger Test Failed')
      return coinc_fail
      

log.basicConfig(level=log.INFO, format='%(asctime)s %(message)s')
sim_gps_start = 1272790000
sim_gps_end =  1272790500
sim_f_lower = 18
fail = False
lsctables.use_in(LIGOLWContentHandler)
      
inj_xmls = sorted(glob.glob('./hwinjcbc*.xml.gz'))
if len(inj_xmls)==0:
      log.error('No inj_xml found when checking results')
      fail=True
else:
      path_to_inj_xml = inj_xmls[0]
      log.info('injection xml found')
      inj_doc = ligolw_utils.load_filename(
            path_to_inj_xml, False, contenthandler=LIGOLWContentHandler)
      inj_table = lsctables.SnglInspiralTable.get_table(inj_doc)
   
      tested_detectors = {'H1', 'L1', 'V1'}
   
      single_fail = check_single_results(tested_detectors)
      coinc_fail = check_coinc_results(inj_table)
      fail = single_fail or coinc_fail
   
if fail:
      log.error('Test Failed')
else:
      log.info('Test Passed')
 
sys.exit(1 if fail else 0)
