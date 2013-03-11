# Copyright (C) 2012  Alex Nitz
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# self.option) any later version.
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
        if self.events[name][-1] is not None:
            cevents, indices = findchirp_cluster_over_window(self.events[name][-1], window_size)
        else:
            return
        
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
        
        # Create sngl_inspiral table ###########################################
        sngl_table = glue.ligolw.lsctables.New(glue.ligolw.lsctables.SnglInspiralTable)
        outdoc.childNodes[0].appendChild(sngl_table)
        
        start_time = lal.LIGOTimeGPS(self.opt.gps_start_time)
        ifo = self.opt.channel_name[0:2]
        
        if self.opt.trig_start_time:
            tstart_time = self.opt.trig_start_time
        else:
            tstart_time = self.opt.gps_start_time + 64
            
        if self.opt.trig_end_time:
            tend_time = self.opt.trig_end_time
        else:
            tend_time = self.opt.gps_end_time - 64
             
        
        for tind in range(len(self.template_params)):
            tmplt = self.template_params[tind]['tmplt']
            snr_norm = self.template_params[tind]['snr_norm']
            sigmasq = self.template_params[tind]['sigmasq']
            
            for eind in range(len(self.events['snr'][tind])):
                row = copy.deepcopy(tmplt)
                    
                snr = self.events['snr'][tind][eind]['val']
                idx = self.events['snr'][tind][eind]['loc']
                end_time = start_time + float(idx) / self.opt.sample_rate
                        
                row.channel = self.opt.channel_name[3:]
                row.ifo = ifo
                
                if self.opt.chisq_bins != 0:
                    row.chisq_dof = self.opt.chisq_bins
                    row.chisq = self.events['chisq'][tind][eind]
                
                row.snr = abs(snr) * snr_norm
                row.end_time = int(end_time.gpsSeconds)
                row.end_time_ns = int(end_time.gpsNanoSeconds)
                row.process_id = proc_id
                row.coa_phase = numpy.angle(snr)
                row.sigmasq = sigmasq
                
                sngl_table.append(row)
                
        # Create Search Summary Table ########################################
        search_summary_table = glue.ligolw.lsctables.New(glue.ligolw.lsctables.SearchSummaryTable)
        outdoc.childNodes[0].appendChild(search_summary_table)
        
        row = glue.ligolw.lsctables.SearchSummary()
        row.nevents = len(sngl_table)
        row.process_id = proc_id
        row.shared_object = ""
        row.lalwrapper_cvs_tag = ""
        row.lal_cvs_tag = ""
        row.comment = ""
        row.ifos = ifo
        row.in_start_time = self.opt.gps_start_time - self.opt.pad_data
        row.in_start_time_ns = 0
        row.in_end_time = self.opt.gps_end_time + self.opt.pad_data
        row.in_end_time_ns = 0
        row.out_start_time = tstart_time
        row.out_start_time_ns = 0
        row.out_end_time = tend_time
        row.out_end_time_ns = 0
        row.nnodes = 1
        
        search_summary_table.append(row)
        
        # Create Filter Table ########################################
        filter_table = glue.ligolw.lsctables.New(glue.ligolw.lsctables.FilterTable)
        outdoc.childNodes[0].appendChild(filter_table)
        
        row = glue.ligolw.lsctables.Filter()
        row.process_id = proc_id
        row.program = "PyCBC_INSPIRAL"
        row.start_time = self.opt.gps_start_time
        row.filter_name = self.opt.approximant
        row.param_set = 0
        row.comment = ""
        row.filter_id = str(glue.ligolw.lsctables.FilterID(0))
        
        filter_table.append(row)
        
        # SumVars Table ########################################
        search_summvars_table = glue.ligolw.lsctables.New(glue.ligolw.lsctables.SearchSummVarsTable)
        outdoc.childNodes[0].appendChild(search_summvars_table)
        
        row = glue.ligolw.lsctables.SearchSummVars()
        row.process_id = proc_id
        row.name = "raw data sample rate"
        row.string = ""
        row.value = 1.0 /16384
        row.search_summvar_id = str(glue.ligolw.lsctables.SearchSummVarsID(0))       
        search_summvars_table.append(row)
        
        row = glue.ligolw.lsctables.SearchSummVars()
        row.process_id = proc_id
        row.name = "filter data sample rate"
        row.string = ""
        row.value = 1.0 / self.opt.sample_rate
        row.search_summvar_id = str(glue.ligolw.lsctables.SearchSummVarsID(1))       
        search_summvars_table.append(row)
        
        # SumValue Table ########################################
        summ_val_columns = ['program', 'process_id', 'start_time', 'start_time_ns',
                           'end_time', 'end_time_ns', 'ifo', 'name', 'value', 'comment',
                           'summ_value_id']
        summ_value_table = glue.ligolw.lsctables.New(glue.ligolw.lsctables.SummValueTable, columns = summ_val_columns)
        outdoc.childNodes[0].appendChild(summ_value_table)
        
        row = glue.ligolw.lsctables.SummValue()
        row.process_id = proc_id
        row.start_time = tstart_time
        row.start_time_ns = 0
        row.end_time = tend_time
        row.end_time_ns = 0
        row.ifo = ifo
        row.frameset_group = ""
        row.program = "PyCBC-INSPIRAL"
        row.error = 0
        row.intvalue = 0
        
        
        row1 = copy.deepcopy(row)
        row2 = copy.deepcopy(row)
        row3 = copy.deepcopy(row)
        row1.name = "inspiral_effective_distance"
        row1.value = 400
        row1.comment = "1.4_1.4_8"
        row1.summ_value_id = str(glue.ligolw.lsctables.SummValueID(0))       
        summ_value_table.append(row1)
        
        row2.name = "calibration alpha"
        row2.value = 0
        row2.comment = "analysis"
        row2.summ_value_id = str(glue.ligolw.lsctables.SummValueID(1))       
        summ_value_table.append(row2)
        
        row3.name = "calibration alphabeta"
        row3.value = 0
        row3.comment = "analysis"
        row3.summ_value_id = str(glue.ligolw.lsctables.SummValueID(2))       
        summ_value_table.append(row3)
        
        
        
        # Write out file #####################################################
        duration = str(int(self.opt.gps_end_time - self.opt.gps_start_time))
        out_name = ifo + "-" + "INSPIRAL_" + self.opt.ifo_tag + "_" + self.opt.user_tag + "-" + str(self.opt.gps_start_time) + "-" + duration + ".xml.gz"
                
        glue.ligolw.utils.write_filename(outdoc, out_name, gz=True)     

__all__ = ['threshold_and_centered_window_cluster', 
           'findchirp_cluster_over_window', 'threshold', 
           'EventManager', 'float32_subset', 'float64_subset',
           'complex64_subset', 'complex128_subset', 'subset_dtype']



