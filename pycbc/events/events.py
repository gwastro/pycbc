# Copyright (C) 2012  Alex Nitz
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
"""This modules defines functions for clustering and thresholding timeseries to 
produces event triggers
"""
import glue.ligolw.utils.process
import lal
import numpy

from pycbc.scheme import schemed
import numpy

complex64_subset = numpy.dtype([('val', numpy.complex64), ('loc', numpy.int64),])
complex128_subset = numpy.dtype([('val', numpy.complex128), ('loc', numpy.int64),])
float32_subset = numpy.dtype([('val', numpy.float32), ('loc', numpy.int64),])
float64_subset = numpy.dtype([('val', numpy.float64), ('loc', numpy.int64),])

_dmap = {numpy.dtype('float32'): float32_subset,
         numpy.dtype('float64'): float64_subset,
         numpy.dtype('complex64'): complex64_subset,
         numpy.dtype('complex128'): complex128_subset}

def subset_dtype(dtype):
    return _dmap[dtype]

@schemed("pycbc.events.threshold_")
def threshold(series, value, offset=0):
    """Return list of values and indices values over threshold in series. 
    """ 
    
@schemed("pycbc.events.threshold_")
def threshold_and_centered_window_cluster(series, threshold, window):
    """Return list of values and indices values over threshold in series. 
    """ 

def findchirp_cluster_over_window(events, window_length):
    clustered_events = numpy.zeros(len(events), dtype=events.dtype)
    j = 0 
    for i in range(len(events)):
        if i==0:
            clustered_events[0] = events[0]
        if events[i]['loc'] - clustered_events[j]['loc'] > window_length:
            j += 1
            clustered_events[j] = events[i]
        else:
            if abs(events[i]['val']) > abs(clustered_events[j]['val']):
                clustered_events[j] = events[i]
            else:
                continue
    return clustered_events[0:j+1]
    
def write_events(outname, all_events, bank, filter_params, opt):
    """ Write the found events to a sngl inspiral table 
    """
    outdoc = glue.ligolw.ligolw.Document()
    outdoc.appendChild(glue.ligolw.ligolw.LIGO_LW())

    proc_id = glue.ligolw.utils.process.register_to_xmldoc(outdoc, 
                    "pycbc_inspiral", opt.__dict__, comment="", ifos=[""],
                    version=glue.git_version.id, cvs_repository=glue.git_version.branch,
                    cvs_entry_time=glue.git_version.date).process_id
    
    sngl_table = glue.ligolw.lsctables.New(glue.ligolw.lsctables.SnglInspiralTable)
    outdoc.childNodes[0].appendChild(sngl_table)
    
    start_time = lal.LIGOTimeGPS(opt.gps_start_time)
    
    ind = 0
    for events, tmplt in zip(all_events, bank.table):
        snr_norm = filter_params['snr_norm'][ind]
        sigmasq = filter_params['sigmasq'][ind]
        
        for event in events:
            row = tmplt
            snr = event['val']
            idx = event['loc']
            end_time = start_time + float(idx) / opt.sample_rate
                      
            row.mass1 = tmplt.mass1
            row.mass2 = tmplt.mass2
            row.snr = abs(snr) * snr_norm
            row.end_time = int(end_time.gpsSeconds)
            row.end_time_ns = int(end_time.gpsNanoSeconds)
            row.process_id = proc_id
            row.coa_phase = numpy.angle(snr)
            row.sigmasq = sigmasq
            row.spin1x = 0
            row.spin1y = 0 
            row.spin1z = 0
            row.spin2x = 0 
            row.spin2y = 0
            row.spin2z = 0    
                  
     
            sngl_table.append(row)
        ind += 1
            
    glue.ligolw.utils.write_filename(outdoc, outname)             

__all__ = ['threshold_and_centered_window_cluster', 
           'findchirp_cluster_over_window', 'threshold', 
           'write_events', 'float32_subset', 'float64_subset',
           'complex64_subset', 'complex128_subset', 'subset_dtype']



