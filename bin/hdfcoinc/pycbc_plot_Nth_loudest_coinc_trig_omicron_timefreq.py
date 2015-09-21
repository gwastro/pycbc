"""
Generates two plots, one for each IFO, that show the time-frequency trace of 
Nth loudest coincident foreground trigger overlaid on a background of 
Omicron triggers.
"""

import logging
import h5py
import numpy as np
import argparse
import glob
from glue.ligolw import ligolw
from glue.ligolw import lsctables
from glue.ligolw import table
from glue.ligolw import utils
import matplotlib
matplotlib.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from pycbc.workflow.segment import fromsegmentxml
import pycbc.pnutils
import pycbc.events

logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

class DefaultContentHandler(ligolw.LIGOLWContentHandler):
    pass
lsctables.use_in(DefaultContentHandler)

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--coinc-file', type=str, required=True,
        help='HDF file containing coincident CBC triggers')
parser.add_argument('--L1-trigs', type=str, required=True,
		help='HDF file containing L1 CBC triggers')
parser.add_argument('--H1-trigs', type=str, required=True,
        help='HDF file containing H1 CBC triggers')
parser.add_argument('--tmpltbank-file', type=str, required=True,
		help='HDF file containing template information for CBC search')
parser.add_argument('--L1-output-file', type=str, required=True,
        help='Full path to output file for L1 trigger')
parser.add_argument('--H1-output-file', type=str, required=True,
        help='Full path to output file for H1 trigger')
parser.add_argument('--loudest-event-number', type=int, required=True, default=1,
        help='Script will plot the Nth loudest coincident foreground trigger')
parser.add_argument('--L1-omicron-dir', type=str, required=True,
        help='Directory containing L1 Omicron triggers. Ex: /home/detchar/triggers/ER7/L1/')
parser.add_argument('--H1-omicron-dir', type=str, required=True,
        help='Directory containing H1 Omicron triggers. Ex: /home/detchar/triggers/ER7/H1/')
parser.add_argument('--omicron-snr-thresh', type=int, required=False, default=5,
        help='SNR threshold for choosing which Omicron triggers to plot.')
args = parser.parse_args()

logging.info('Reading HDF files')

coinc_trig_file = h5py.File(args.coinc_file,'r')
L1_trig_file = h5py.File(args.L1_trigs,'r')
H1_trig_file = h5py.File(args.H1_trigs,'r')
template_file = h5py.File(args.tmpltbank_file,'r')

logging.info('Parsing HDF files')

coinc_newsnr = coinc_trig_file['background_exc']['stat'][:]
Nth_loudest_idx = np.argsort(coinc_newsnr)[-args.loudest_event_number]

if coinc_trig_file.attrs['detector_1'] == 'H1':
    H1_idx = coinc_trig_file['background_exc']['trigger_id1'][Nth_loudest_idx]
    L1_idx = coinc_trig_file['background_exc']['trigger_id2'][Nth_loudest_idx]
elif coinc_trig_file.attrs['detector_1'] == 'L1':
    L1_idx = coinc_trig_file['background_exc']['trigger_id1'][Nth_loudest_idx]
    H1_idx = coinc_trig_file['background_exc']['trigger_id2'][Nth_loudest_idx]

# get info about single detector triggers that comprise loudest background event
# and calculate newSNR
H1_snr = H1_trig_file['H1']['snr'][H1_idx]
H1_chisq = H1_trig_file['H1']['chisq'][H1_idx]
H1_chisq_dof = H1_trig_file['H1']['chisq_dof'][H1_idx]
H1_reduced_chisq = H1_chisq/(2*H1_chisq_dof - 2)
H1_newsnr = pycbc.events.newsnr(H1_snr,H1_reduced_chisq)
H1_cbc_end_time = H1_trig_file['H1']['end_time'][H1_idx]
H1_template_id = H1_trig_file['H1']['template_id'][H1_idx]

L1_snr = L1_trig_file['L1']['snr'][L1_idx]
L1_chisq = L1_trig_file['L1']['chisq'][L1_idx]
L1_chisq_dof = L1_trig_file['L1']['chisq_dof'][L1_idx]
L1_reduced_chisq = L1_chisq/(2*L1_chisq_dof - 2)
L1_newsnr = pycbc.events.newsnr(L1_snr,L1_reduced_chisq)
L1_cbc_end_time = L1_trig_file['L1']['end_time'][L1_idx]
L1_template_id = L1_trig_file['L1']['template_id'][L1_idx]

# Sanity check: exact match means that the templates should be the same
if L1_template_id != H1_template_id:
    logging.info('Templates did not match between detectors!')

mass1 = template_file['mass1'][H1_template_id]
mass2 = template_file['mass2'][H1_template_id]

H1_hoft_chan = 'H1:GDS-CALIB_STRAIN'
L1_hoft_chan = 'L1:GDS-CALIB_STRAIN'
H1_omicron_start_time = H1_cbc_end_time - 100
H1_omicron_end_time = H1_cbc_end_time + 100
L1_omicron_start_time = L1_cbc_end_time - 100
L1_omicron_end_time = L1_cbc_end_time + 100

logging.info('Fetching omicron triggers')

# Generate list of directories to search over
H1_gps_era_start = str(H1_omicron_start_time)[:5]
H1_gps_era_end = str(H1_omicron_end_time)[:5]

H1_eras = map(str,range(int(H1_gps_era_start),int(H1_gps_era_end)))
if not H1_eras:
    H1_eras = [H1_gps_era_start]

L1_gps_era_start = str(L1_omicron_start_time)[:5]
L1_gps_era_end = str(L1_omicron_end_time)[:5]

L1_eras = map(str,range(int(L1_gps_era_start),int(L1_gps_era_end)))
if not L1_eras:
    L1_eras = [L1_gps_era_start]

# Grab all relevant Omicron trigger files
H1_omicron_times = []
H1_omicron_snr = []
H1_omicron_freq = []

for era in H1_eras:
    # Generate list of all Omicron SnglBurst xml trigger files
    H1_file_list = glob.glob(args.H1_omicron_dir + 
            '/GDS-CALIB_STRAIN_Omicron/%s/H1-GDS_CALIB_STRAIN_Omicron-*.xml.gz' %(era))
    
    # Parse trigger files into SNR, time, and frequency for Omicron triggers
    for file in H1_file_list:
        omicron_xml = utils.load_filename(file, contenthandler=DefaultContentHandler)
        H1snglburst_table = table.get_table(omicron_xml, lsctables.SnglBurstTable.tableName)

        for row in H1snglburst_table:
            if (row.snr > args.omicron_snr_thresh and 
                    H1_omicron_start_time < row.peak_time < H1_omicron_end_time):
                H1_omicron_times.append(row.peak_time + row.peak_time_ns * 10**(-9))
                H1_omicron_snr.append(row.snr)
                H1_omicron_freq.append(row.peak_frequency)


L1_omicron_times = []
L1_omicron_snr = []
L1_omicron_freq = []

for era in L1_eras:
    # Generate list of all Omicron SnglBurst xml trigger files
    L1_file_list = glob.glob(args.L1_omicron_dir + 
            '/GDS-CALIB_STRAIN_Omicron/%s/L1-GDS_CALIB_STRAIN_Omicron-*.xml.gz' %(era))
    
    # Parse trigger files into SNR, time, and frequency for Omicron triggers
    for file in L1_file_list:
        omicron_xml = utils.load_filename(file, contenthandler=DefaultContentHandler)
        L1snglburst_table = table.get_table(omicron_xml, lsctables.SnglBurstTable.tableName)

        for row in L1snglburst_table:
            if (row.snr > args.omicron_snr_thresh and
                    L1_omicron_start_time < row.peak_time < L1_omicron_end_time):
                L1_omicron_times.append(row.peak_time + row.peak_time_ns * 10**(-9))
                L1_omicron_snr.append(row.snr)
                L1_omicron_freq.append(row.peak_frequency)

f_low = 30
f_high = pycbc.pnutils.f_SchwarzISCO(mass1+mass2)
H1_inspiral_t, H1_inspiral_f = pycbc.pnutils.get_inspiral_tf(H1_cbc_end_time, mass1, mass2, f_low, f_high)
L1_inspiral_t, L1_inspiral_f = pycbc.pnutils.get_inspiral_tf(L1_cbc_end_time, mass1, mass2, f_low, f_high)

logging.info('Plotting')

plt.figure(0)
cm = plt.cm.get_cmap('jet')
plt.scatter(H1_omicron_times,H1_omicron_freq,c=H1_omicron_snr,s=30,cmap=cm,linewidth=0)
plt.grid(b=True, which='both')
cbar = plt.colorbar()
cbar.set_label('H1 Omicron trigger SNR')
plt.yscale('log')
plt.ylabel('Frequency (Hz)')
plt.xlabel('Time (s)')
plt.xlim(H1_omicron_start_time,H1_omicron_end_time)
plt.suptitle('H1 CBC trigger SNR = ' + format(H1_snr,'.2f') + ", newSNR = " + format(H1_newsnr,'.2f'),fontsize=12)
plt.title(format(mass1,'.2f') + " - " + format(mass2,'.2f') + " solar masses at GPS time " + format(H1_cbc_end_time,'.2f'),fontsize=12)
plt.hold(True)
plt.plot(H1_inspiral_t,H1_inspiral_f)
plt.savefig(args.H1_output_file)

plt.figure(1)
cm = plt.cm.get_cmap('jet')
plt.scatter(L1_omicron_times,L1_omicron_freq,c=L1_omicron_snr,s=30,cmap=cm,linewidth=0)
plt.grid(b=True, which='both')
cbar = plt.colorbar()
cbar.set_label('L1 Omicron trigger SNR')
plt.yscale('log')
plt.ylabel('Frequency (Hz)')
plt.xlabel('Time (s)')
plt.xlim(L1_omicron_start_time,L1_omicron_end_time)
plt.suptitle('L1 CBC trigger SNR = ' + format(L1_snr,'.2f') + ", newSNR = " + format(L1_newsnr,'.2f'),fontsize=12)
plt.title(format(mass1,'.2f') + " - " + format(mass2,'.2f') + " solar masses at GPS time " + format(L1_cbc_end_time,'.2f'),fontsize=12)
plt.hold(True)
plt.plot(L1_inspiral_t,L1_inspiral_f)
plt.savefig(args.L1_output_file)

logging.info('Done! Exiting script.')
