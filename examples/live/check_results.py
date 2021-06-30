#!/usr/bin/env python

import sys
import argparse
import glob
import logging as log
import numpy as np
import h5py
import pycbc
from pycbc.io import FieldArray
from glue.ligolw import utils as ligolw_utils
from glue.ligolw import ligolw, lsctables


# dummy class needed for loading LIGOLW files
class LIGOLWContentHandler(ligolw.LIGOLWContentHandler):
    pass

def close(a, b, c):
    return abs(a - b) <= c

def check_single_results(args):
    single_fail = False
    with h5py.File(args.bank, 'r') as bankf:
        temp_mass1 = bankf['mass1'][:]
        temp_mass2 = bankf['mass2'][:]
        temp_s1z = bankf['spin1z'][:]
        temp_s2z = bankf['spin2z'][:]

    detectors = set(args.detectors)
    detectors_with_trigs = set()

    trig_paths = sorted(glob.glob('output/????_??_??/*.hdf'))
    for trigfp in trig_paths:
        with h5py.File(trigfp, 'r') as trigf:
            for detector in detectors:
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
                in_band_asd = psd[psd_f > args.f_min] ** 0.5
                if len(in_band_asd) == 0 or (in_band_asd < 1e-24).any() \
                        or (in_band_asd > 1e-20).any() \
                        or not np.isfinite(in_band_asd).all() \
                        or psd_epoch < args.gps_start or psd_epoch > args.gps_end:
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
                if (end_time < args.gps_start).any() or (end_time > args.gps_end).any():
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
    missing = sorted(detectors - detectors_with_trigs)
    if missing:
        log.error('No triggers found in %s', ', '.join(missing))
        single_fail = True

    if single_fail:
        log.error('Single Trigger Test Failed')
    return single_fail


def check_coinc_results(args):
    coinc_fail = False
    # gather coincident triggers
    coinc_trig_paths = sorted(glob.glob('output/coinc*.xml.gz'))
    n_coincs = len(coinc_trig_paths)
    if n_coincs == 0:
        log.error('No coincident triggers detected')
        coinc_fail = True
    elif n_coincs >= 10:
        log.error('Too many coincident triggers detected')
        coinc_fail = True
    else:
        log.info('%d coincident trigger(s) detected', n_coincs)

    with h5py.File(args.injections, 'r') as injfile:
        inj_mass1 = injfile['mass1'][:]
        inj_mass2 = injfile['mass2'][:]
        inj_spin1z = injfile['spin1z'][:]
        inj_spin2z = injfile['spin2z'][:]
        inj_time = injfile['tc'][:]

    if len(inj_mass1) > n_coincs:
        log.error('More injections than coincident triggers')
        coinc_fail = True

    # create field array to store properties of triggers
    dtype = [('mass1', float), ('mass2', float),
             ('spin1z', float), ('spin2z', float),
             ('tc', float), ('net_snr', float)]
    trig_props = FieldArray(n_coincs, dtype=dtype)

    # store properties of coincident triggers
    for x, ctrigfp in enumerate(coinc_trig_paths):
        log.info('Checking trigger %s', ctrigfp)
        xmldoc = ligolw_utils.load_filename(
            ctrigfp, False, contenthandler=LIGOLWContentHandler)
        sngl_inspiral_table = lsctables.SnglInspiralTable.get_table(xmldoc)

        trig_props['tc'][x] = sngl_inspiral_table.get_column('end_time')[0]
        trig_props['mass1'][x] = sngl_inspiral_table.get_column('mass1')[0]
        trig_props['mass2'][x] = sngl_inspiral_table.get_column('mass2')[0]
        trig_props['spin1z'][x] = sngl_inspiral_table.get_column('spin1z')[0]
        trig_props['spin2z'][x] = sngl_inspiral_table.get_column('spin2z')[0]

        snr_list = sngl_inspiral_table.get_column('snr')
        trig_props['net_snr'][x] = sum(snr_list ** 2) ** 0.5

        log.info('IFO SNRs: %s', snr_list)
        log.info('Network SNR: %f', trig_props['net_snr'][x])
        log.info('IFO End Time: %f', trig_props['tc'][x])
        log.info('Mass 1: %f', trig_props['mass1'][x])
        log.info('Mass 2: %f', trig_props['mass2'][x])
        log.info('Spin1z: %f', trig_props['spin1z'][x])
        log.info('Spin2z: %f', trig_props['spin2z'][x])

    # check if injections match trigger params
    for i in range(len(inj_mass1)):
        has_match = False
        for j in range(n_coincs):
            # FIXME should calculate the optimal SNRs of the injections
            # and use those for checking net_snr
            if (close(inj_time[i], trig_props['tc'][j], 1.0)
                    and close(inj_mass1[i], trig_props['mass1'][j], 5e-7)
                    and close(inj_mass2[i], trig_props['mass2'][j], 5e-7)
                    and close(inj_spin1z[i], trig_props['spin1z'][j], 5e-7)
                    and close(inj_spin2z[i], trig_props['spin2z'][j], 5e-7)
                    and close(15.0, trig_props['net_snr'][j], 2.0)):
                has_match = True
                break

        if not has_match:
            coinc_fail = True
            log.error('Injection %i has no match', i)

    if coinc_fail:
        log.error('Coincident Trigger Test Failed')
    return coinc_fail


parser = argparse.ArgumentParser()
parser.add_argument('--gps-start', type=float, required=True)
parser.add_argument('--gps-end', type=float, required=True)
parser.add_argument('--f-min', type=float, required=True)
parser.add_argument('--bank', type=str, required=True)
parser.add_argument('--injections', type=str, required=True)
parser.add_argument('--detectors', type=str, required=True, nargs='+')
args = parser.parse_args()

log.basicConfig(level=log.INFO, format='%(asctime)s %(message)s')

lsctables.use_in(LIGOLWContentHandler)

single_fail = check_single_results(args)
coinc_fail = check_coinc_results(args)
fail = single_fail or coinc_fail

if fail:
    log.error('Test Failed')
else:
    log.info('Test Passed')

sys.exit(1 if fail else 0)
