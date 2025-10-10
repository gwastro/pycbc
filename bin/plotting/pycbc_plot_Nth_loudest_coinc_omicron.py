"""
Generates a plot that shows the time-frequency trace of
Nth loudest coincident trigger overlaid on a background of
Omicron triggers.
"""

import logging
import numpy as np
import argparse
import glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from igwn_ligolw import lsctables, utils

import pycbc.events
from pycbc.waveform import (
    get_td_waveform, frequency_from_polarizations,
    amplitude_from_polarizations
)
from pycbc.io.ligolw import LIGOLWContentHandler
from pycbc.io.hdf import HFile


parser = argparse.ArgumentParser(description=__doc__)
pycbc.add_common_pycbc_options(parser)
parser.add_argument('--coinc-file', type=str, required=True,
                    help='HDF file containing coincident CBC triggers')
parser.add_argument('--single-ifo-trigs', type=str, required=True,
                    help='HDF file containing single IFO CBC triggers')
parser.add_argument('--ifo', type=str, required=True,
                    help='IFO, L1 or H1')
parser.add_argument('--tmpltbank-file', type=str, required=True,
                    help='HDF file containing template information for CBC search')
parser.add_argument('--output-file', type=str, required=True,
                    help='Full path to output file')
parser.add_argument('--loudest-event-number', type=int, required=True, default=1,
                    help='Script will plot the Nth loudest coincident trigger')
parser.add_argument('--omicron-dir', type=str, required=True,
                    help='Directory containing Omicron triggers. Ex: /home/detchar/triggers/ER7/')
parser.add_argument('--omicron-snr-thresh', type=int, required=False, default=5,
                    help='SNR threshold for choosing which Omicron triggers to plot.')
parser.add_argument('--plot-window', type=float, required=False, default=32,
                    help='Time window to plot around CBC trigger')
parser.add_argument('--omicron-channel',type=str, required=False, default='GDS-CALIB_STRAIN',
                    help='Channel to plot Omicron triggers for, do not include IFO')
parser.add_argument('--analysis-level', type=str, required=False, default='foreground',
                    choices = ['foreground','background','background_exc'],
                    help='Designates which level of the analysis output to search')
args = parser.parse_args()

pycbc.init_logging(args.verbose)

logging.info('Reading HDF files')

coinc_trig_file = HFile(args.coinc_file,'r')
single_trig_file = HFile(args.single_ifo_trigs,'r')
template_file = HFile(args.tmpltbank_file,'r')

logging.info('Parsing HDF files')

coinc_newsnr = coinc_trig_file[args.analysis_level]['stat'][:]
Nth_loudest_idx = np.argsort(coinc_newsnr)[-args.loudest_event_number]

if coinc_trig_file.attrs['detector_1'] == args.ifo:
    idx = coinc_trig_file[args.analysis_level]['trigger_id1'][Nth_loudest_idx]
else:
    idx = coinc_trig_file[args.analysis_level]['trigger_id2'][Nth_loudest_idx]

# get info about single detector triggers that comprise loudest background event
# and calculate newSNR
snr = single_trig_file[args.ifo]['snr'][idx]
chisq = single_trig_file[args.ifo]['chisq'][idx]
chisq_dof = single_trig_file[args.ifo]['chisq_dof'][idx]
reduced_chisq = chisq/(2*chisq_dof - 2)
newsnr = pycbc.events.ranking.newsnr(snr,reduced_chisq)
cbc_end_time = single_trig_file[args.ifo]['end_time'][idx]
template_id = single_trig_file[args.ifo]['template_id'][idx]

m1 = template_file['mass1'][template_id]
m2 = template_file['mass2'][template_id]
s1z = template_file['spin1z'][template_id]
s2z = template_file['spin2z'][template_id]

omicron_start_time = cbc_end_time - args.plot_window
omicron_end_time = cbc_end_time + args.plot_window

logging.info('Fetching omicron triggers')

# Generate list of directories to search over
gps_era_start = str(omicron_start_time)[:5]
gps_era_end = str(omicron_end_time)[:5]

eras = map(str,range(int(gps_era_start),int(gps_era_end)))
if not eras:
    eras = [gps_era_start]

# Grab all relevant Omicron trigger files
omicron_times = []
omicron_snr = []
omicron_freq = []

for era in eras:
    # Generate list of all Omicron SnglBurst xml trigger files
    file_list = glob.glob(args.omicron_dir +
            '/%s/%s_Omicron/%s/%s-%s_Omicron-*.xml.gz'
            %(args.ifo,args.omicron_channel,era,args.ifo,args.omicron_channel.replace('-','_')))

    # Parse trigger files into SNR, time, and frequency for Omicron triggers
    for file_name in file_list:
        omicron_xml = utils.load_filename(
                file_name, contenthandler=LIGOLWContentHandler)
        snglburst_table = lsctables.SnglBurstTable.get_table(omicron_xml)

        for row in snglburst_table:
            if (row.snr > args.omicron_snr_thresh and
                    omicron_start_time < row.peak_time < omicron_end_time):
                omicron_times.append(row.peak_time + row.peak_time_ns * 10**(-9))
                omicron_snr.append(row.snr)
                omicron_freq.append(row.peak_frequency)


# Generate inspiral waveform and calculate f(t) to plot on top of Omicron triggers
hp, hc = get_td_waveform(approximant='SEOBNRv2', mass1=m1, mass2=m2,
                 spin1x=0, spin1y=0, spin1z=s1z,
                 spin2x=0, spin2y=0, spin2z=s2z,
                 delta_t=(1./32768.), f_lower=30)


f = frequency_from_polarizations(hp, hc)
amp = amplitude_from_polarizations(hp, hc)
stop_idx = amp.abs_max_loc()[1]

f = f[:stop_idx]

freq = np.array(f.data)
times = np.array(f.sample_times) + cbc_end_time

logging.info('Plotting')

plt.figure(0)
cm = plt.cm.get_cmap('Reds')
plt.scatter(omicron_times,omicron_freq,c=omicron_snr,s=30,cmap=cm,linewidth=0)
plt.grid(b=True, which='both')
cbar = plt.colorbar()
cbar.set_label('%s Omicron trigger SNR' % (args.ifo))
plt.yscale('log')
plt.ylabel('Frequency (Hz)')
plt.xlabel('Time (s)')
plt.xlim(omicron_start_time,omicron_end_time)
plt.suptitle('%s CBC trigger SNR = ' % (args.ifo) + format(snr,'.2f') +
            ", newSNR = " + format(newsnr,'.2f'),fontsize=12)
plt.title(format(m1,'.2f') + " - " + format(m2,'.2f') +
            " solar masses at GPS time " + format(cbc_end_time,'.2f'),fontsize=12)
plt.hold(True)
plt.plot(times,freq)
plt.savefig(args.output_file)

logging.info('Done! Exiting script.')
