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

from pycbc import WEAVE_FLAGS
from pycbc.types import Array
from pycbc.types import convert_to_process_params_dict
from pycbc.scheme import schemed
from pycbc.detector import Detector

from . import coinc

@schemed("pycbc.events.threshold_")
def threshold(series, value):
    """Return list of values and indices values over threshold in series.
    """
    return
    
@schemed("pycbc.events.threshold_")
def threshold_only(series, value):
    """Return list of values and indices whose values in series are
       larger (in absolute value) than value
    """
    return

#FIXME: This should be under schemed, but I don't understand that yet!
def threshold_real_numpy(series, value):
    arr = series.data
    locs = numpy.where(arr > value)[0]
    vals = arr[locs]
    return locs, vals

@schemed("pycbc.events.threshold_")
def threshold_and_cluster(series, threshold, window):
    """Return list of values and indices values over threshold in series.
    """
    return

@schemed("pycbc.events.threshold_")
def _threshold_cluster_factory(series):
    return

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
        int curr_ind = 0;
        for (int i=0; i < tlen; i++){
            if ((times[i] - times[curr_ind]) > window_length){
                j += 1;
                indices[j] = i;
                curr_ind = i;
            }
            else if (absvalues[i] > absvalues[curr_ind]){
                indices[j] = i;
                curr_ind = i;
            }
        }
        k[0] = j;
    """
    inline(code, ['times', 'absvalues', 'window_length', 'indices', 'tlen', 'k'],
                 extra_compile_args=[WEAVE_FLAGS])
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
        self.write_performance = False

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

    def save_performance(self, ncores, nfilters, ntemplates, run_time,
                         setup_time):
        """
        Calls variables from pycbc_inspiral to be used in a timing calculation
        """
        self.run_time = run_time
        self.setup_time = setup_time
        self.ncores = ncores
        self.nfilters = nfilters
        self.ntemplates = ntemplates
        self.write_performance = True

    def write_events(self, outname):
        """ Write the found events to a sngl inspiral table
        """
        self.make_output_dir(outname)

        if '.hdf' in outname:
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
        th = numpy.array([p['tmplt'].template_hash for p in self.template_params])
        tid = self.events['template_id']
        f = fw(outname, self.opt.channel_name[0:2])

        if len(self.events):
            f['snr'] = abs(self.events['snr'])
            try:
                # Precessing
                f['u_vals'] = self.events['u_vals']
                f['coa_phase'] = self.events['coa_phase']
                f['hplus_cross_corr'] = self.events['hplus_cross_corr']
            except:
                # Not precessing
                f['coa_phase'] = numpy.angle(self.events['snr'])
            f['chisq'] = self.events['chisq']
            f['bank_chisq'] = self.events['bank_chisq']
            f['bank_chisq_dof'] = self.events['bank_chisq_dof']
            f['cont_chisq'] = self.events['cont_chisq']
            f['end_time'] = self.events['time_index'] / float(self.opt.sample_rate) + self.opt.gps_start_time
            try:
                # Precessing
                template_sigmasq_plus = numpy.array([t['sigmasq_plus'] for t in self.template_params], dtype=numpy.float32)
                f['sigmasq_plus'] = template_sigmasq_plus[tid]
                template_sigmasq_cross = numpy.array([t['sigmasq_cross'] for t in self.template_params], dtype=numpy.float32)
                f['sigmasq_cross'] = template_sigmasq_cross[tid]
                # FIXME: I want to put something here, but I haven't yet
                #        figured out what it should be. I think we would also
                #        need information from the plus and cross correlation
                #        (both real and imaginary(?)) to get this.
                f['sigmasq'] = template_sigmasq_plus[tid]
            except:
                # Not precessing
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
            search_start_time = float(self.opt.trig_start_time)
        else:
            f['search/start_time'] = numpy.array([self.opt.gps_start_time + self.opt.segment_start_pad])
            search_start_time = float(self.opt.gps_start_time + self.opt.segment_start_pad)
        if self.opt.trig_end_time:
            f['search/end_time'] = numpy.array([self.opt.trig_end_time])
            search_end_time = float(self.opt.trig_end_time)
        else:
            f['search/end_time'] = numpy.array([self.opt.gps_end_time - self.opt.segment_end_pad])
            search_end_time = float(self.opt.gps_end_time - self.opt.segment_end_pad)

        if self.write_performance:
            self.analysis_time = search_end_time - search_start_time
            time_ratio = numpy.array([float(self.analysis_time) / float(self.run_time)])
            temps_per_core = float(self.ntemplates) / float(self.ncores)
            filters_per_core = float(self.nfilters) / float(self.ncores)
            f['search/templates_per_core'] = \
                numpy.array([float(temps_per_core) * float(time_ratio)])
            f['search/filter_rate_per_core'] = \
                numpy.array([filters_per_core / float(self.run_time)])
            f['search/setup_time_fraction'] = \
                numpy.array([float(self.setup_time) / float(self.run_time)])

        if 'gating_info' in self.global_params:
            gating_info = self.global_params['gating_info']
            for gate_type in ['file', 'auto']:
                if gate_type in gating_info:
                    f['gating/' + gate_type + '/time'] = \
                            numpy.array([float(g[0]) for g in gating_info[gate_type]])
                    f['gating/' + gate_type + '/width'] = \
                            numpy.array([g[1] for g in gating_info[gate_type]])
                    f['gating/' + gate_type + '/pad'] = \
                            numpy.array([g[2] for g in gating_info[gate_type]])

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
        raise NotImplementedError


__all__ = ['threshold_and_cluster', 'newsnr', 'effsnr',
           'findchirp_cluster_over_window',
           'threshold', 'cluster_reduce', 'ThresholdCluster',
           'threshold_real_numpy',
           'EventManager', 'EventManagerMultiDet']
