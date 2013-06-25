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
import pycbc

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
def threshold_and_cluster(series, threshold, window):
    """Return list of values and indices values over threshold in series. 
    """ 

def findchirp_cluster_over_window(times, values, window_length):
    indices = numpy.zeros(len(times), dtype=int)
    j = 0 
    for i in range(len(times)):
        if i==0:
            indices[0] = 0
        if times[i] - times[indices[j]] > window_length:
            j += 1
            indices[j] = i
        else:
            if abs(values[i]) > abs(values[indices[j]]):
                indices[j] = i
            else:
                continue
    return indices[0:j+1]
    
class EventManager(object):
    def __init__(self, opt, column, column_types, **kwds):
        self.opt = opt
        self.global_params = kwds

        self.event_dtype = [ ('template_id', int) ]
        for column, coltype in zip (column, column_types):
            self.event_dtype.append( (column, coltype) )
        
        self.events = numpy.events = numpy.array([], dtype=self.event_dtype)
        self.template_params = []
        self.template_index = -1
        self.template_events = numpy.array([], dtype=self.event_dtype)
                  
    def chisq_threshold(self, value, num_bins, delta=0):
        remove = []
        for i in range(len(self.events)):
            event = self.events[i]
            tind = event['template_id']
            norm = self.template_params[tind]['snr_norm'] ** 2.0
            xi = event['chisq'] / (num_bins + delta * event['snr'].conj() * event['snr'] * norm)
            if xi > value:
                remove.append(i)
        self.events = numpy.delete(self.events, remove)          
   
    def maximize_over_bank(self, tcolumn, column, window):
        pass
#        self.events.sort(order=[tcolumn])
##        nbins = numpy.ceil(self.events[-1][tcolumn]/float(window))
#        edges = numpy.arange(0, nbins, 1) * float(window)
#        indices = numpy.searchsorted(self.events[tcolumn], edges)
#        indices = numpy.append(indices, len(self.events))
#        
#        maxes = []
#        locs = []
#        row = []
#        remove = []
##        #print nbins, edges, indices, self.events[-1][tcolumn], window
#        for i in range(len(indices)-1):
#            l = indices[i]
#            r = indices[i+1]
#            if l == r:
#                continue
#            maxid = abs(self.events[l:r][column]).argmax()
#            maxes.append(self.events[l:r][column][maxid])
#            locs.append(self.events[l:r][tcolumn][maxid])
##            rows.append(self.events[l:r][maxid])
#
#        for i in range(len(maxes)):
#            if i == 0:
#                lm = 0;
               
        
        
    def add_template_events(self, columns, vectors):
        """ Add a vector indexed """
        new_events = numpy.zeros(len(vectors[0]), dtype=self.event_dtype)
        new_events['template_id'] = self.template_index
        for c, v in zip(columns, vectors):
            if v is not None:
                new_events[c] = v  
        self.template_events = numpy.append(self.template_events, new_events)                     
     
    def cluster_template_events(self, tcolumn, column, window_size):
        """ Cluster the internal events over the named column
        """
        cvec = self.template_events[column]
        tvec = self.template_events[tcolumn]
        indices = findchirp_cluster_over_window(tvec, cvec, window_size)
        self.template_events = numpy.take(self.template_events, indices)       
        
    def new_template(self, **kwds):
        self.template_params.append(kwds)
        self.template_index += 1
    
    def add_template_params(self, **kwds):
        self.template_params[-1].update(kwds)  
            
    def finalize_template_events(self):
        self.events = numpy.append(self.events, self.template_events)
        self.template_events = numpy.array([], dtype=self.event_dtype)        
        
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
            
        event_num = 0      
        for event in self.events:
            tind = event['template_id']
        
            tmplt = self.template_params[tind]['tmplt']
            snr_norm = self.template_params[tind]['snr_norm']
            sigmasq = self.template_params[tind]['sigmasq']
            
            row = copy.deepcopy(tmplt)
                
            snr = event['snr']
            idx = event['time_index']
            end_time = start_time + float(idx) / self.opt.sample_rate
                    
            row.channel = self.opt.channel_name[3:]
            row.ifo = ifo
            
            if self.opt.chisq_bins != 0:
                # FIXME: This is *not* the dof!!!
                row.chisq_dof = self.opt.chisq_bins
                row.chisq = event['chisq']

            if self.opt.bank_veto_bank_file:
                row.bank_chisq_dof = self.global_params['num_bank_templates']
                row.bank_chisq = event['bank_chisq']
            
            row.eff_distance = sigmasq ** (0.5) / abs(snr * snr_norm) * pycbc.DYN_RANGE_FAC
            row.snr = abs(snr) * snr_norm
            row.end_time = int(end_time.gpsSeconds)
            row.end_time_ns = int(end_time.gpsNanoSeconds)
            row.process_id = proc_id
            row.coa_phase = numpy.angle(snr)
            row.sigmasq = sigmasq
            
            row.event_id = glue.ligolw.lsctables.SnglInspiralID(event_num)
            event_num += 1
            
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


__all__ = ['threshold_and_cluster', 
           'findchirp_cluster_over_window', 'threshold', 
           'EventManager', 'float32_subset', 'float64_subset',
           'complex64_subset', 'complex128_subset', 'subset_dtype']



