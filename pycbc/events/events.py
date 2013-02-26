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
import copy

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
def threshold(series, value):
    """Return list of values and indices values over threshold in series. 
    """ 
    
@schemed("pycbc.events.threshold_")
def threshold_and_centered_window_cluster(series, threshold, window):
    """Return list of values and indices values over threshold in series. 
    """ 

def findchirp_cluster_over_window(events, window_length):
    clustered_events = numpy.zeros(len(events), dtype=events.dtype)
    indices = numpy.zeros(len(events), dtype=int)
    j = 0 
    for i in range(len(events)):
        if i==0:
            clustered_events[0] = events[0]
        if events[i]['loc'] - clustered_events[j]['loc'] > window_length:
            j += 1
            clustered_events[j] = events[i]
            indices[j] = i
        else:
            if abs(events[i]['val']) > abs(clustered_events[j]['val']):
                clustered_events[j] = events[i]
                indices[j] = i
            else:
                continue
    return clustered_events[0:j+1], indices[0:j+1]
    
class EventManager(object):
    def __init__(self, columns, opt, **kwds):
        self.opt = opt
        self.global_params = kwds
        self.events = {}
        self.template_params = []
        for column in columns:
            self.events[column] = []
        self.output_file = "out.xml"
        
    def add_template_events(self, columns, vectors):
        """ Add a vector indexed values
        """
        for column, vector in zip(columns, vectors):
            data = self.events[column][-1]
            if self.events[column][-1] is not None:
                self.events[column][-1] = numpy.append(data, vector)
            else:
                self.events[column][-1] = vector
                
     
    def cluster_template_events(self, name, window_size):
        """ Cluster the internal events over the named column
        """
        cevents, indices = findchirp_cluster_over_window(self.events[name][-1], window_size)
        
        for column in self.events:
            self.events[column][-1] = numpy.take(self.events[column][-1], indices)
        
    def new_template(self, **kwds):
        self.template_params.append(kwds)
        for column in self.events:
            self.events[column].append(None)
    
    def add_template_params(self, **kwds):
            self.template_params[-1].update(kwds)           
        
    def write_events(self):
        """ Write the found events to a sngl inspiral table 
        """
        outdoc = glue.ligolw.ligolw.Document()
        outdoc.appendChild(glue.ligolw.ligolw.LIGO_LW())

        proc_id = glue.ligolw.utils.process.register_to_xmldoc(outdoc, 
                        "pycbc_inspiral", self.opt.__dict__, comment="", ifos=[""],
                        version=glue.git_version.id, cvs_repository=glue.git_version.branch,
                        cvs_entry_time=glue.git_version.date).process_id
        
        sngl_table = glue.ligolw.lsctables.New(glue.ligolw.lsctables.SnglInspiralTable)
        outdoc.childNodes[0].appendChild(sngl_table)
        
        start_time = lal.LIGOTimeGPS(self.opt.gps_start_time)
        
        for tind in range(len(self.template_params)):
            tmplt = self.template_params[tind]['tmplt']
            snr_norm = self.template_params[tind]['snr_norm']
            sigmasq = self.template_params[tind]['sigmasq']
            
            for eind in range(len(self.events['snr'][tind])):
                print type(tmplt), "HUH?", len(self.template_params), len(self.events['snr'][tind])
                row = copy.deepcopy(tmplt)
                    
                snr = self.events['snr'][tind][eind]['val']
                idx = self.events['snr'][tind][eind]['loc']
                end_time = start_time + float(idx) / self.opt.sample_rate
                        
                row.chisq = self.events['chisq'][tind][eind]
                row.snr = abs(snr) * snr_norm
                row.end_time = int(end_time.gpsSeconds)
                row.end_time_ns = int(end_time.gpsNanoSeconds)
                row.process_id = proc_id
                row.coa_phase = numpy.angle(snr)
                row.sigmasq = sigmasq
                
                sngl_table.append(row)
                
        glue.ligolw.utils.write_filename(outdoc, self.output_file)     

__all__ = ['threshold_and_centered_window_cluster', 
           'findchirp_cluster_over_window', 'threshold', 
           'EventManager', 'float32_subset', 'float64_subset',
           'complex64_subset', 'complex128_subset', 'subset_dtype']



