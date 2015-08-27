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
import lal, numpy, copy, os.path

from pycbc.types import Array
from pycbc.types import convert_to_process_params_dict
from pycbc.scheme import schemed
from pycbc.detector import Detector

from . import coinc

@schemed("pycbc.events.threshold_")
def threshold(series, value):
    """Return list of values and indices values over threshold in series.
    """
@schemed("pycbc.events.threshold_")
def threshold_only(series, value):
    """Return list of values and indices whose values in series are
       larger (in absolute value) than value
    """

@schemed("pycbc.events.threshold_")
def threshold_and_cluster(series, threshold, window):
    """Return list of values and indices values over threshold in series.
    """

@schemed("pycbc.events.threshold_")
def _threshold_cluster_factory(series):
    pass

class ThresholdCluster(object):
    """Create a threshold and cluster engine

    Parameters
    ----------
    series : complex64
      Input pycbc.types.Array (or subclass); it will be searched for
      points above threshold that are then clustered
    """
    def __new__(cls, *args, **kwargs):
        real_cls = _threshold_cluster_factory(*args, **kwargs)
        return real_cls(*args, **kwargs)

# The class below should serve as the parent for all schemed classes.
# The intention is that this class serves simply as the location for
# all documentation of the class and its methods, though that is not
# yet implemented.  Perhaps something along the lines of:
#
#    http://stackoverflow.com/questions/2025562/inherit-docstrings-in-python-class-inheritance
#
# will work? Is there a better way?
class _BaseThresholdCluster(object):
    def threshold_and_cluster(self, threshold, window):
        """
        Threshold and cluster the memory specified at instantiation with the
        threshold specified at creation and the window size specified at creation.

        Parameters:
        -----------
        threshold : float32
          The minimum absolute value of the series given at object initialization
          to return when thresholding and clustering.
        window : uint32
          The size (in number of samples) of the window over which to cluster

        Returns:
        --------
        event_vals : complex64
          Numpy array, complex values of the clustered events
        event_locs : uint32
          Numpy array, indices into series of location of events
        """
        pass


def findchirp_cluster_over_window(times, values, window_length):
    """ Reduce the events by clustering over a window using
    the FindChirp clustering algorithm

    Parameters
    -----------
    indices: Array
        The list of indices of the SNR values
    snr: Array
        The list of SNR value
    window_size: int
        The size of the window in integer samples. Must be positive.

    Returns
    -------
    indices: Array
        The reduced list of indices of the SNR values
    """
    assert window_length > 0, 'Clustering window length is not positive'

    from scipy.weave import inline
    indices = numpy.zeros(len(times), dtype=int)
    tlen = len(times)
    k = numpy.zeros(1, dtype=int)
    absvalues = abs(values)
    times = times.astype(int)
    code = """
        int j = 0;
        for (int i=0; i < tlen; i++){
            if ((times[i] - times[indices[j]]) > window_length){
                j += 1;
                indices[j] = i;
            }
            else if (absvalues[i] > absvalues[indices[j]]){
                indices[j] = i;
            }
        }
        k[0] = j;
    """
    inline(code, ['times', 'absvalues', 'window_length', 'indices', 'tlen', 'k'],
                 extra_compile_args=['-march=native -O3 -w'])
    return indices[0:k[0]+1]

def cluster_reduce(idx, snr, window_size):
    """ Reduce the events by clustering over a window

    Parameters
    -----------
    indices: Array
        The list of indices of the SNR values
    snr: Array
        The list of SNR value
    window_size: int
        The size of the window in integer samples.

    Returns
    -------
    indices: Array
        The list of indices of the SNR values
    snr: Array
        The list of SNR values
    """
    ind = findchirp_cluster_over_window(idx, snr, window_size)
    return idx.take(ind), snr.take(ind)

def newsnr(snr, reduced_x2, q=6., n=2.):
    """Calculate the re-weighted SNR statistic ('newSNR') from given SNR and
    reduced chi-squared values. See http://arxiv.org/abs/1208.3491 for
    definition. Previous implementation in glue/ligolw/lsctables.py
    """
    newsnr = numpy.array(snr, ndmin=1, dtype=numpy.float64)
    reduced_x2 = numpy.array(reduced_x2, ndmin=1, dtype=numpy.float64)

    # newsnr is only different from snr if reduced chisq > 1
    ind = numpy.where(reduced_x2 > 1.)[0]
    newsnr[ind] *= ( 0.5 * (1. + reduced_x2[ind] ** (q/n)) ) ** (-1./q)

    if len(newsnr) > 1:
        return newsnr
    else:
        return newsnr[0]

def effsnr(snr, reduced_x2, fac=250.):
    """Calculate the effective SNR statistic. See (S5y1 paper) for definition.
    Previous implementation in glue/ligolw/lsctables.py
    """
    snr = numpy.array(snr, ndmin=1, dtype=numpy.float64)
    rchisq = numpy.array(reduced_x2, ndmin=1, dtype=numpy.float64)
    effsnr = snr / (1 + snr ** 2 / fac) ** 0.25 / rchisq ** 0.25

    if len(effsnr) > 1:
        return effsnr
    else:
        return effsnr[0]

class EventManager(object):
    def __init__(self, opt, column, column_types, **kwds):
        self.opt = opt
        self.global_params = kwds

        self.event_dtype = [ ('template_id', int) ]
        for column, coltype in zip (column, column_types):
            self.event_dtype.append( (column, coltype) )

        self.events = numpy.array([], dtype=self.event_dtype)
        self.template_params = []
        self.template_index = -1
        self.template_events = numpy.array([], dtype=self.event_dtype)

    @classmethod
    def from_multi_ifo_interface(cls, opt, ifo, column, column_types, **kwds):
        """
        To use this for a single ifo from the multi ifo interface requires
        some small fixing of the opt structure. This does that. As we edit the
        opt structure the process_params table will not be correct.
        """
        opt = copy.deepcopy(opt)
        opt_dict = vars(opt)
        for arg, value in opt_dict.items():
            if isinstance(value, dict):
                setattr(opt, arg, getattr(opt, arg)[ifo])
        return cls(opt, column, column_types, **kwds)

    def chisq_threshold(self, value, num_bins, delta=0):
        remove = []
        for i, event in enumerate(self.events):
            xi = event['chisq'] / (event['chisq_dof'] / 2 + 1 + delta * event['snr'].conj() * event['snr'])
            if xi > value:
                remove.append(i)
        self.events = numpy.delete(self.events, remove)

    def newsnr_threshold(self, threshold):
        """ Remove events with newsnr smaller than given threshold
        """
        if not self.opt.chisq_bins:
            raise RuntimeError('Chi-square test must be enabled in order to use newsnr threshold')

        remove = [i for i, e in enumerate(self.events) if \
            newsnr(abs(e['snr']), e['chisq'] / e['chisq_dof']) < threshold]
        self.events = numpy.delete(self.events, remove)
        
    def keep_near_injection(self, window, injections):
        from pycbc.events.veto import indices_within_times
        if len(self.events) == 0:
            return
        
        inj_time = numpy.array(injections.end_times())
        gpstime = self.events['time_index'].astype(numpy.float64)
        gpstime = gpstime / self.opt.sample_rate + self.opt.gps_start_time
        i = indices_within_times(gpstime, inj_time - window, inj_time + window)
        self.events = self.events[i]

    def keep_loudest_in_interval(self, window, num_keep):
        if len(self.events) == 0:
            return
        
        e = self.events
        stat = newsnr(abs(e['snr']), e['chisq'] / e['chisq_dof'])
        time = e['time_index']
        
        wtime = (time / window).astype(numpy.int32)
        bins = numpy.unique(wtime)
        
        keep = []
        for b in bins:
            bloc = numpy.where((wtime == b))[0]
            bloudest = stat[bloc].argsort()[-num_keep:]
            keep.append(bloc[bloudest])
        keep = numpy.concatenate(keep)
        self.events = e[keep]

    def maximize_over_bank(self, tcolumn, column, window):
        if len(self.events) == 0:
            return

        self.events = numpy.sort(self.events, order=tcolumn)
        cvec = self.events[column]
        tvec = self.events[tcolumn]

        indices = []
#        mint = tvec.min()
#        maxt = tvec.max()
#        edges = numpy.arange(mint, maxt, window)

#        # Get the location of each time bin
#        bins = numpy.searchsorted(tvec, edges)
#        bins = numpy.append(bins, len(tvec))
#        for i in range(len(bins)-1):
#            kmin = bins[i]
#            kmax = bins[i+1]
#            if kmin == kmax:
#                continue
#            event_idx = numpy.argmax(cvec[kmin:kmax]) + kmin
#            indices.append(event_idx)

        # This algorithm is confusing, but it is what lalapps_inspiral does
        # REMOVE ME!!!!!!!!!!!
        gps = tvec.astype(numpy.float64) / self.opt.sample_rate + self.opt.gps_start_time
        gps_sec  = numpy.floor(gps)
        gps_nsec = (gps - gps_sec) * 1e9

        wnsec = int(window * 1e9 / self.opt.sample_rate)
        win = gps_nsec.astype(int) / wnsec

        indices.append(0)
        for i in xrange(len(tvec)):
            if gps_sec[i] == gps_sec[indices[-1]] and  win[i] == win[indices[-1]]:
                    if abs(cvec[i]) > abs(cvec[indices[-1]]):
                        indices[-1] = i
            else:
                indices.append(i)

        self.events = numpy.take(self.events, indices)

    def add_template_events(self, columns, vectors):
        """ Add a vector indexed """
        # initialize with zeros - since vectors can be None, look for the
        # first one that isn't
        new_events = None
        for v in vectors:
            if v is not None:
                new_events = numpy.zeros(len(v), dtype=self.event_dtype)
                break
        # they shouldn't all be None
        assert new_events is not None
        new_events['template_id'] = self.template_index
        for c, v in zip(columns, vectors):
            if v is not None:
                if isinstance(v, Array):
                    new_events[c] = v.numpy()
                else:
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

    def make_output_dir(self, outname):
        path = os.path.dirname(outname)
        if path != '':
            if not os.path.exists(path) and path is not None:
                os.makedirs(path)

    def write_events(self, outname):
        """ Write the found events to a sngl inspiral table
        """
        self.make_output_dir(outname)

        if '.xml' in outname:
            self.write_to_xml(outname)
        elif '.hdf' in outname:
            self.write_to_hdf(outname)
        else:
            raise ValueError('Cannot write to this format')

    def write_to_hdf(self, outname):
        class fw(object):
            def __init__(self, name, prefix):
                import h5py
                self.f = h5py.File(name, 'w')
                self.prefix = prefix

            def __setitem__(self, name, data):
                col = self.prefix + '/' + name
                self.f.create_dataset(col, data=data,
                                      compression='gzip',
                                      compression_opts=9,
                                      shuffle=True)

        self.events.sort(order='template_id')

        # Template id hack
        m1 = numpy.array([p['tmplt'].mass1 for p in self.template_params], dtype=numpy.float32)
        m2 = numpy.array([p['tmplt'].mass2 for p in self.template_params], dtype=numpy.float32)
        s1 = numpy.array([p['tmplt'].spin1z for p in self.template_params], dtype=numpy.float32)
        s2 = numpy.array([p['tmplt'].spin2z for p in self.template_params], dtype=numpy.float32)
        th = numpy.zeros(len(m1), dtype=int)
        for j, v in enumerate(zip(m1, m2, s1, s2)):
            th[j] = hash(v)

        tid = self.events['template_id']
        f = fw(outname, self.opt.channel_name[0:2])

        if len(self.events):
            f['snr'] = abs(self.events['snr'])
            f['coa_phase'] = numpy.angle(self.events['snr'])
            f['chisq'] = self.events['chisq']
            f['bank_chisq'] = self.events['bank_chisq']
            f['bank_chisq_dof'] = self.events['bank_chisq_dof']
            f['cont_chisq'] = self.events['cont_chisq']
            f['end_time'] = self.events['time_index'] / float(self.opt.sample_rate) + self.opt.gps_start_time

            template_sigmasq = numpy.array([t['sigmasq'] for t in self.template_params], dtype=numpy.float32)
            f['sigmasq'] = template_sigmasq[tid]

            template_durations = [p['tmplt'].template_duration for p in self.template_params]
            f['template_duration'] = numpy.array(template_durations, dtype=numpy.float32)[tid]

            # FIXME: Can we get this value from the autochisq instance?
            cont_dof = self.opt.autochi_number_points
            if self.opt.autochi_onesided is None:
                cont_dof = cont_dof * 2
            if self.opt.autochi_two_phase:
                cont_dof = cont_dof * 2
            if self.opt.autochi_max_valued_dof:
                cont_dof = self.opt.autochi_max_valued_dof
            f['cont_chisq_dof'] = numpy.repeat(cont_dof, len(self.events))

            if 'chisq_dof' in self.events.dtype.names:
                f['chisq_dof'] = self.events['chisq_dof'] / 2 + 1
            else:
                f['chisq_dof'] = numpy.zeros(len(self.events))

            f['template_hash'] = th[tid]

        if self.opt.trig_start_time:
            f['search/start_time'] = numpy.array([self.opt.trig_start_time])
        else:
            f['search/start_time'] = numpy.array([self.opt.gps_start_time + self.opt.segment_start_pad])

        if self.opt.trig_end_time:
            f['search/end_time'] = numpy.array([self.opt.trig_end_time])
        else:
            f['search/end_time'] = numpy.array([self.opt.gps_end_time - self.opt.segment_end_pad])

    def write_to_xml(self, outname):
        """ Write the found events to a sngl inspiral table
        """
        outdoc = glue.ligolw.ligolw.Document()
        outdoc.appendChild(glue.ligolw.ligolw.LIGO_LW())

        ifo = self.opt.channel_name[0:2]

        proc_id = glue.ligolw.utils.process.register_to_xmldoc(outdoc,
                        "inspiral", self.opt.__dict__, comment="", ifos=[ifo],
                        version=glue.git_version.id,
                        cvs_repository=glue.git_version.branch,
                        cvs_entry_time=glue.git_version.date).process_id

        # Create sngl_inspiral table ###########################################
        sngl_table = glue.ligolw.lsctables.New(\
                                       glue.ligolw.lsctables.SnglInspiralTable)
        self._add_sngls_to_output(sngl_table, proc_id)
        outdoc.childNodes[0].appendChild(sngl_table)

        # Create Search Summary Table ########################################
        search_summary_table = self._create_search_summary_table(proc_id,
                                                               len(sngl_table))
        outdoc.childNodes[0].appendChild(search_summary_table)

        # Create Filter Table ########################################
        filter_table = self._create_filter_table(proc_id)
        outdoc.childNodes[0].appendChild(filter_table)

        # SumVars Table ########################################
        search_summvars_table = self._create_search_summvars_table(proc_id)
        outdoc.childNodes[0].appendChild(search_summvars_table)

        # SumValue Table ########################################
        summ_value_table = self._create_summ_val_table(proc_id)
        outdoc.childNodes[0].appendChild(summ_value_table)

        # Write out file #####################################################
        glue.ligolw.utils.write_filename(outdoc, outname,
                                         gz=outname.endswith('gz'))

    def _add_sngls_to_output(self, sngl_table, proc_id, ifo=None, channel=None,
                             start_time=None, sample_rate=None,
                             multi_ifo=False):
        """
        Add events to sngl inspiral table.
        """
        if multi_ifo and ifo is not None:
            err_msg = "If using multiple ifos you cannot supply the ifo kwarg "
            err_msg += "to _add_sngls_to_output"
            raise ValueError(err_msg)

        if start_time is None:
            start_time = lal.LIGOTimeGPS(self.opt.gps_start_time)
        if sample_rate is None:
            sample_rate = self.opt.sample_rate
        if ifo is None and not multi_ifo:
            ifo = self.opt.channel_name[0:2]
        if channel is None:
            if multi_ifo:
                channel = {}
                for ifo in self.ifos:
                    channel[ifo] = self.opt.channel_name[ifo].split(':')[1]
            else:
                channel = self.opt.channel_name.split(':')[1]

        for event_num, event in enumerate(self.events):
            tind = event['template_id']

            tmplt = self.template_params[tind]['tmplt']

            row = copy.deepcopy(tmplt)

            snr = event['snr']
            idx = event['time_index']
            end_time = start_time + float(idx) / sample_rate

            if multi_ifo:
                sigmasq = self.template_params[tind]['sigmasq'][ifo]
                ifo = self.ifo_reverse[event['ifo']]
                row.channel = channel[ifo]
            else:
                sigmasq = self.template_params[tind]['sigmasq']
                row.channel = channel
            row.ifo = ifo

            row.chisq = event['chisq']
            # FIXME: This is *not* the dof!!!
            # but is needed for later programs not to fail
            if 'chisq_dof' in event.dtype.names:
                row.chisq_dof = event['chisq_dof'] / 2 + 1
            else:
                row.chisq_dof = 0

            if hasattr(self.opt, 'bank_veto_bank_file')\
                    and self.opt.bank_veto_bank_file:
                row.bank_chisq = event['bank_chisq']
                row.bank_chisq_dof = event['bank_chisq_dof']
            else:
                row.bank_chisq_dof = 0
                row.bank_chisq = 0

            if hasattr(self.opt, 'autochi_number_points')\
                    and self.opt.autochi_number_points > 0:
                row.cont_chisq = event['cont_chisq']
                # FIXME: Can this come from the autochisq instance?
                cont_dof = self.opt.autochi_number_points
                if self.opt.autochi_onesided is None:
                    cont_dof = cont_dof * 2
                if self.opt.autochi_two_phase:
                    cont_dof = cont_dof * 2
                if self.opt.autochi_max_valued_dof:
                    cont_dof = self.opt.autochi_max_valued_dof
                row.cont_chisq_dof = cont_dof

            row.eff_distance = sigmasq ** (0.5) / abs(snr)
            row.snr = abs(snr)
            row.end_time = int(end_time.gpsSeconds)
            row.end_time_ns = int(end_time.gpsNanoSeconds)
            row.process_id = proc_id
            row.coa_phase = numpy.angle(snr)
            row.sigmasq = sigmasq

            row.event_id = glue.ligolw.lsctables.SnglInspiralID(event_num)

            sngl_table.append(row)

    def _create_search_summary_table(self, proc_id, nevents,
                                     ifo=None, start_time=None, end_time=None,
                                     trig_start_time=None, trig_end_time=None):
        if ifo is None:
            ifo = self.opt.channel_name[0:2]
        if start_time is None:
            start_time = self.opt.gps_start_time - self.opt.pad_data
        if end_time is None:
            end_time = self.opt.gps_end_time + self.opt.pad_data
        if trig_start_time is None:
            if self.opt.trig_start_time:
                trig_start_time = self.opt.trig_start_time
            else:
                trig_start_time = self.opt.gps_start_time +\
                              self.opt.segment_start_pad
        if trig_end_time is None:
            if self.opt.trig_end_time:
                trig_end_time = self.opt.trig_end_time
            else:
                trig_end_time = self.opt.gps_end_time - self.opt.segment_end_pad

        search_summary_table = glue.ligolw.lsctables.New(\
                                      glue.ligolw.lsctables.SearchSummaryTable)
        row = glue.ligolw.lsctables.SearchSummary()
        row.nevents = nevents
        row.process_id = proc_id
        row.shared_object = ""
        row.lalwrapper_cvs_tag = ""
        row.lal_cvs_tag = ""
        row.comment = ""
        row.ifos = ifo
        row.in_start_time = start_time
        row.in_start_time_ns = 0
        row.in_end_time = end_time
        row.in_end_time_ns = 0
        row.out_start_time = trig_start_time
        row.out_start_time_ns = 0
        row.out_end_time = trig_end_time
        row.out_end_time_ns = 0
        row.nnodes = 1

        search_summary_table.append(row)
        return search_summary_table

    def _create_filter_table(self, proc_id, start_time=None, approximant=None):
        if start_time is None:
            start_time = self.opt.gps_start_time
        if approximant is None:
            approximant = self.opt.approximant

        filter_table = glue.ligolw.lsctables.New(\
                                             glue.ligolw.lsctables.FilterTable)

        row = glue.ligolw.lsctables.Filter()
        row.process_id = proc_id
        row.program = "PyCBC_INSPIRAL"
        row.start_time = start_time
        row.filter_name = approximant
        row.param_set = 0
        row.comment = ""
        row.filter_id = str(glue.ligolw.lsctables.FilterID(0))

        filter_table.append(row)
        return filter_table

    def _create_search_summvars_table(self, proc_id, sample_rate=None):
        if sample_rate is None:
            sample_rate = self.opt.sample_rate

        search_summvars_table = glue.ligolw.lsctables.New(\
                                     glue.ligolw.lsctables.SearchSummVarsTable)

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
        row.value = 1.0 / sample_rate
        row.search_summvar_id = str(glue.ligolw.lsctables.SearchSummVarsID(1))
        search_summvars_table.append(row)
        return search_summvars_table

    def _create_summ_val_table(self, proc_id, ifo=None, trig_start_time=None,
                                trig_end_time=None, low_frequency_cutoff=None):

        if ifo is None:
            ifo = self.opt.channel_name[0:2]

        if trig_start_time is None:
            if self.opt.trig_start_time:
                trig_start_time = self.opt.trig_start_time
            else:
                trig_start_time = self.opt.gps_start_time +\
                              self.opt.segment_start_pad
        if trig_end_time is None:
            if self.opt.trig_end_time:
                trig_end_time = self.opt.trig_end_time
            else:
                trig_end_time = self.opt.gps_end_time - self.opt.segment_end_pad
        if low_frequency_cutoff is None:
            low_frequency_cutoff = self.opt.low_frequency_cutoff

        summ_val_columns = ['program', 'process_id', 'start_time',
                            'start_time_ns', 'end_time', 'end_time_ns', 'ifo',
                            'name', 'value', 'comment', 'summ_value_id']
        summ_value_table = glue.ligolw.lsctables.New(\
                glue.ligolw.lsctables.SummValueTable, columns=summ_val_columns)

        row = glue.ligolw.lsctables.SummValue()
        row.process_id = proc_id
        row.start_time = trig_start_time
        row.start_time_ns = 0
        row.end_time = trig_end_time
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

        psd = self.global_params['psd']
        from pycbc.waveform.spa_tmplt import spa_distance
        from pycbc import DYN_RANGE_FAC
        # FIXME: Lalapps did this "right" for non-spa waveforms.
        #        Should also be right here (maybe covering a range of masses)
        row1.value = spa_distance(psd, 1.4, 1.4, self.opt.low_frequency_cutoff,
                                  snr=8) * DYN_RANGE_FAC
        row1.comment = "1.4_1.4_8"
        row1.summ_value_id = str(glue.ligolw.lsctables.SummValueID(0))
        summ_value_table.append(row1)

        # FIXME: We haven't run on uncalibrated data since S4(?)
        #        Do we really still need this?
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

        return summ_value_table

class EventManagerMultiDet(EventManager):
    def __init__(self, opt, ifos, column, column_types, psd=None, **kwargs):
        self.opt = opt
        self.ifos = ifos
        self.global_params = kwargs
        if psd is not None:
            self.global_params['psd'] = psd[ifos[0]]

        # The events array does not like holding the ifo as string,
        # so create a mapping dict and hold as an int
        self.ifo_dict = {}
        self.ifo_reverse = {}
        for i, ifo in enumerate(ifos):
            self.ifo_dict[ifo] = i
            self.ifo_reverse[i] = ifo

        self.event_dtype = [ ('template_id', int), ('event_id', int) ]
        for column, coltype in zip (column, column_types):
            self.event_dtype.append( (column, coltype) )

        self.events = numpy.array([], dtype=self.event_dtype)
        self.event_id_map = {}
        self.event_index = 0
        self.template_params = []
        self.template_index = -1
        self.template_event_dict = {}
        self.coinc_list = []
        for ifo in ifos:
            self.template_event_dict[ifo] = numpy.array([],
                                                        dtype=self.event_dtype)

    def add_template_events_to_ifo(self, ifo, columns, vectors):
        """ Add a vector indexed """
        # Just call through to the standard function
        self.template_events = self.template_event_dict[ifo]
        self.add_template_events(columns, vectors)
        self.template_event_dict[ifo] = self.template_events
        self.template_events = None

    def cluster_template_events_single_ifo(self, tcolumn, column, window_size,
                                          ifo):
        """ Cluster the internal events over the named column
        """
        # Just call through to the standard function
        self.template_events = self.template_event_dict[ifo]
        self.cluster_template_events(tcolumn, column, window_size)
        self.template_event_dict[ifo] = self.template_events
        self.template_events = None

    def finalize_template_events(self, perform_coincidence=True,
                                 coinc_window=0.0):
        # Set ids
        for ifo in self.ifos:
            num_events = len(self.template_event_dict[ifo])
            new_event_ids = numpy.arange(self.event_index,
                                                   self.event_index+num_events)
            self.template_event_dict[ifo]['event_id'] = new_event_ids
            self.event_index = self.event_index+num_events

        if perform_coincidence:
            if not len(self.ifos) == 2:
                err_msg = "Coincidence currently only supported for 2 ifos."
                raise ValueError(err_msg)
            ifo1 = self.ifos[0]
            ifo2 = self.ifos[1]
            end_times1 = self.template_event_dict[ifo1]['time_index']  /\
              float(self.opt.sample_rate[ifo1]) + self.opt.gps_start_time[ifo1]
            end_times2 = self.template_event_dict[ifo2]['time_index']  /\
              float(self.opt.sample_rate[ifo2]) + self.opt.gps_start_time[ifo2]
            light_travel_time = Detector(ifo1).light_travel_time_to_detector(\
                                                                Detector(ifo2))
            coinc_window = coinc_window + light_travel_time
            # FIXME: Remove!!!
            coinc_window = 2.0
            if len(end_times1) and len(end_times2):
                idx_list1, idx_list2, _ = \
                                 coinc.time_coincidence(end_times1, end_times2,
                                                        coinc_window)
                if len(idx_list1):
                    for idx1, idx2 in zip(idx_list1, idx_list2):
                        event1 = self.template_event_dict[ifo1][idx1]
                        event2 = self.template_event_dict[ifo2][idx2]
                        self.coinc_list.append((event1, event2))
        for ifo in self.ifos:
            self.events = numpy.append(self.events,
                                                 self.template_event_dict[ifo])
            self.template_event_dict[ifo] = numpy.array([],
                                                        dtype=self.event_dtype)

    def write_events(self, outname):
        """ Write the found events to a sngl inspiral table
        """
        self.make_output_dir(outname)
        outdoc = glue.ligolw.ligolw.Document()
        outdoc.appendChild(glue.ligolw.ligolw.LIGO_LW())

        ifostring = ''.join(self.ifos)
        ifo_ex = self.ifos[0]
        start_time = self.opt.gps_start_time[ifo_ex]
        start_time_gps = lal.LIGOTimeGPS(start_time)
        start_time_padded = self.opt.gps_start_time[ifo_ex] \
                             - self.opt.pad_data[ifo_ex]
        end_time_padded = self.opt.gps_end_time[ifo_ex] \
                           + self.opt.pad_data[ifo_ex]
        if self.opt.trig_start_time[ifo_ex]:
            trig_start_time = self.opt.trig_start_time[ifo_ex]
        else:
            trig_start_time = self.opt.gps_start_time[ifo_ex] \
                               + self.opt.segment_start_pad[ifo_ex]
        if self.opt.trig_end_time[ifo_ex]:
            trig_end_time = self.opt.trig_end_time[ifo_ex]
        else:
            trig_end_time = self.opt.gps_end_time[ifo_ex] \
                             + self.opt.segment_end_pad[ifo_ex]
        sample_rate = self.opt.sample_rate[ifo_ex]
        approximant = self.opt.approximant
        low_frequency_cutoff = self.opt.low_frequency_cutoff

        proc_id = glue.ligolw.utils.process.register_to_xmldoc(outdoc,
                        "inspiral", convert_to_process_params_dict(self.opt),
                        comment="", ifos=[ifostring],
                        version=glue.git_version.id,
                        cvs_repository=glue.git_version.branch,
                        cvs_entry_time=glue.git_version.date).process_id

        # Create sngl_inspiral table ###########################################
        sngl_table = glue.ligolw.lsctables.New(\
                                       glue.ligolw.lsctables.SnglInspiralTable)
        self._add_sngls_to_output(sngl_table, proc_id,
                                  start_time=start_time_gps,
                                  sample_rate=sample_rate, multi_ifo=True)
        outdoc.childNodes[0].appendChild(sngl_table)

        # Create the coincidence tables ######################################
        coinc_def_table = glue.ligolw.lsctables.New(\
                                       glue.ligolw.lsctables.CoincDefTable)
        coinc_event_table = glue.ligolw.lsctables.New(\
                                         glue.ligolw.lsctables.CoincTable)
        coinc_event_map_table = glue.ligolw.lsctables.New(\
                                      glue.ligolw.lsctables.CoincMapTable)
        time_slide_table = glue.ligolw.lsctables.New(\
                                          glue.ligolw.lsctables.TimeSlideTable)
        coinc_inspiral_table = glue.ligolw.lsctables.New(\
                                      glue.ligolw.lsctables.CoincInspiralTable)
        self._add_coincs_to_output(coinc_def_table, coinc_event_table,
                                   coinc_event_map_table, time_slide_table,
                                   coinc_inspiral_table, sngl_table, proc_id)
        outdoc.childNodes[0].appendChild(coinc_def_table)
        outdoc.childNodes[0].appendChild(coinc_event_table)
        outdoc.childNodes[0].appendChild(coinc_event_map_table)
        outdoc.childNodes[0].appendChild(time_slide_table)
        outdoc.childNodes[0].appendChild(coinc_inspiral_table)

        # Create Search Summary Table ########################################
        search_summary_table = self._create_search_summary_table(proc_id,
                           len(sngl_table),
                           ifo=ifostring, start_time=start_time_padded,
                           end_time=end_time_padded,
                           trig_start_time=trig_start_time,
                           trig_end_time=trig_end_time)
        outdoc.childNodes[0].appendChild(search_summary_table)

        # Create Filter Table ########################################
        filter_table = self._create_filter_table(proc_id, start_time=start_time,
                                                 approximant=approximant)
        outdoc.childNodes[0].appendChild(filter_table)

        # SumVars Table ########################################
        search_summvars_table = self._create_search_summvars_table(proc_id,
                                                       sample_rate=sample_rate)
        outdoc.childNodes[0].appendChild(search_summvars_table)

        # SumValue Table ########################################
        summ_value_table = self._create_summ_val_table(proc_id, ifo=ifostring,
                                     trig_start_time=trig_start_time,
                                     trig_end_time=trig_end_time,
                                     low_frequency_cutoff=low_frequency_cutoff)
        outdoc.childNodes[0].appendChild(summ_value_table)

        # Write out file #####################################################
        glue.ligolw.utils.write_filename(outdoc, outname,
                                         gz=outname.endswith('gz'))

    def _add_coincs_to_output(self, coinc_def_table, coinc_event_table,
                              coinc_event_map_table, time_slide_table,
                              coinc_inspiral_table, sngl_table, proc_id):
        # FIXME: This shouldn't live here
        # FIXME: More choices would be good
        magic_number = 6.0
        def get_weighted_snr(self, fac):
            rchisq = self.chisq/(2*self.chisq_dof - 2)
            nhigh = 2.
            if rchisq > 1.:
                return self.snr/((1+rchisq**(fac/nhigh))/2)**(1./fac)
            else:
                return self.snr

        # Define global IDs up front:
        coinc_def_id = glue.ligolw.lsctables.CoincDefID(0)
        # FIXME: Add support for multiple slides
        time_slide_id = glue.ligolw.lsctables.TimeSlideID(0)
        for ifo in self.ifos:
            time_slide_row = glue.ligolw.lsctables.TimeSlide()
            time_slide_row.instrument = ifo
            time_slide_row.time_slide_id = time_slide_id
            time_slide_row.offset = 0
            time_slide_row.process_id = proc_id
            time_slide_table.append(time_slide_row)
        time_slide_dict = time_slide_table.as_dict()

        ifostring = ''.join(self.ifos)

        count = 0
        for coinc in self.coinc_list:
            # Check that all sngls are present
            coinc_removed_flag=0
            for sngl in coinc:
                if not self.event_id_map.has_key(sngl['event_id']):
                    # If not event_id then one of the sngls is not in the sngl
                    # table because it was removed at some point after testing
                    # coincidence. Therefore this is not still a coincident
                    # event.
                    coinc_removed_flag=1
                    break
            if coinc_removed_flag:
                continue
            coinc_id = glue.ligolw.lsctables.CoincID(count)
            count = count+1
            # Create the coinc map entry
            sngl_xmls = []
            for sngl in coinc:
                coinc_map_row = glue.ligolw.lsctables.CoincMap()
                # I really need this .... every time?!
                coinc_map_row.table_name = 'sngl_inspiral'
                coinc_map_row.coinc_event_id = coinc_id
                sngl_id_num = self.event_id_map[sngl['event_id']]
                sngl_id = glue.ligolw.lsctables.SnglInspiralID(sngl_id_num)
                coinc_map_row.event_id = sngl_id
                coinc_event_map_table.append(coinc_map_row)
                # NOTE: This now assumes event_ids are ordered in sngl_inspiral
                #       table.
                sngl_xmls.append(sngl_table[sngl_id_num])

            # Now construct the coinc_inspiral, which is actually *two* tables
            coinc_event_row = glue.ligolw.lsctables.Coinc()
            coinc_inspiral_row = glue.ligolw.lsctables.CoincInspiral()
            # Fill the joining/meta columns
            coinc_event_row.coinc_def_id = coinc_def_id
            coinc_event_row.nevents = len(coinc)
            coinc_event_row.instruments = ifostring
            coinc_inspiral_row.set_ifos(self.ifos)
            coinc_event_row.time_slide_id = time_slide_id
            coinc_event_row.process_id = proc_id
            coinc_event_row.coinc_event_id = coinc_id
            coinc_inspiral_row.coinc_event_id = coinc_id

            # Meaningful rows
            coinc_inspiral_row.mchirp = sum(sngl.mchirp for sngl in sngl_xmls)\
                                                               / len(sngl_xmls)
            coinc_inspiral_row.minimum_duration = \
                              min(sngl.template_duration for sngl in sngl_xmls)
            coinc_inspiral_row.mass = sum(sngl.mass1 + sngl.mass2 \
                                        for sngl in sngl_xmls) / len(sngl_xmls)
            # End time is chosen as the unslid time of the first ifo in the
            # coincidence, where "first" ifo is chosen alphabetically.
            first_xml = min(sngl_xmls, key = lambda sngl: sngl.ifo)
            end_time = first_xml.get_end() + \
                                    timeslid_dict[time_slide_id][first_xml.ifo]
            coinc_inspiral_row.set_end(end_time)
            coinc_inspiral_row.snr = numpy.sqrt( sum( \
                                  get_weighted_snr(sngl, fac=magic_number)**2 \
                                                        for sngl in sngl_xmls))

            # Rows that are populated later
            coinc_event_row.likelihood = 0.
            coinc_inspiral_row.false_alarm_rate = 0.
            coinc_inspiral_row.combined_far = 0.

            # Add new row
            coinc_event_table.append(coinc_event_row)
            coinc_inspiral_table.append(coinc_inspiral_row)

        # Create coinc_definer table
        coinc_def_row = glue.ligolw.lsctables.CoincDef()
        coinc_def_row.search = "inspiral"
        coinc_def_row.description = "sngl_inspiral-sngl_inspiral coincidences"
        coinc_def_row.coinc_def_id = coinc_def_id
        coinc_def_row.search_coinc_type = 0
        coinc_def_table.append(coinc_def_row)

__all__ = ['threshold_and_cluster', 'newsnr', 'effsnr',
           'findchirp_cluster_over_window',
           'threshold', 'cluster_reduce', 'ThresholdCluster',
           'EventManager', 'EventManagerMultiDet']
