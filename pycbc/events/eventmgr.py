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
from __future__ import absolute_import
import numpy, copy, os.path

from pycbc import WEAVE_FLAGS
from pycbc.types import Array
from pycbc.scheme import schemed
from pycbc.detector import Detector

from . import coinc, ranking

@schemed("pycbc.events.threshold_")
def threshold(series, value):
    """Return list of values and indices values over threshold in series.
    """
    err_msg = "This function is a stub that should be overridden using the "
    err_msg += "scheme. You shouldn't be seeing this error!"
    raise ValueError(err_msg)

@schemed("pycbc.events.threshold_")
def threshold_only(series, value):
    """Return list of values and indices whose values in series are
       larger (in absolute value) than value
    """
    err_msg = "This function is a stub that should be overridden using the "
    err_msg += "scheme. You shouldn't be seeing this error!"
    raise ValueError(err_msg)

# FIXME: This should be under schemed, but I don't understand that yet!
def threshold_real_numpy(series, value):
    arr = series.data
    locs = numpy.where(arr > value)[0]
    vals = arr[locs]
    return locs, vals

@schemed("pycbc.events.threshold_")
def threshold_and_cluster(series, threshold, window):
    """Return list of values and indices values over threshold in series.
    """
    err_msg = "This function is a stub that should be overridden using the "
    err_msg += "scheme. You shouldn't be seeing this error!"
    raise ValueError(err_msg)

@schemed("pycbc.events.threshold_")
def _threshold_cluster_factory(series):
    err_msg = "This class is a stub that should be overridden using the "
    err_msg += "scheme. You shouldn't be seeing this error!"
    raise ValueError(err_msg)


class ThresholdCluster(object):
    """Create a threshold and cluster engine

    Parameters
    -----------
    series : complex64
      Input pycbc.types.Array (or subclass); it will be searched for
      points above threshold that are then clustered
    """
    def __new__(cls, *args, **kwargs):
        real_cls = _threshold_cluster_factory(*args, **kwargs)
        return real_cls(*args, **kwargs) # pylint:disable=not-callable


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
        threshold and window size specified at creation.

        Parameters
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

    from weave import inline
    indices = numpy.zeros(len(times), dtype=int)
    tlen = len(times) # pylint:disable=unused-variable
    k = numpy.zeros(1, dtype=int)
    absvalues = abs(values) # pylint:disable=unused-variable
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


class EventManager(object):
    def __init__(self, opt, column, column_types, **kwds):
        self.opt = opt
        self.global_params = kwds

        self.event_dtype = [('template_id', int)]
        for col, coltype in zip(column, column_types):
            self.event_dtype.append((col, coltype))

        self.events = numpy.array([], dtype=self.event_dtype)
        self.accumulate = [self.events]
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
            xi = event['chisq'] / (event['chisq_dof'] / 2 + 1 +
                                   delta * event['snr'].conj() * event['snr'])
            if xi > value:
                remove.append(i)
        self.events = numpy.delete(self.events, remove)

    def newsnr_threshold(self, threshold):
        """ Remove events with newsnr smaller than given threshold
        """
        if not self.opt.chisq_bins:
            raise RuntimeError('Chi-square test must be enabled in order to '
                               'use newsnr threshold')

        remove = [i for i, e in enumerate(self.events) if
                  ranking.newsnr(abs(e['snr']), e['chisq'] / e['chisq_dof'])
                  < threshold]
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

    def keep_loudest_in_interval(self, window, num_keep, statname="newsnr",
                                 log_chirp_width=None):
        if len(self.events) == 0:
            return

        from pycbc.events import stat
        e_copy = self.events.copy()

        # Here self.events['snr'] is the complex SNR
        e_copy['snr'] = abs(e_copy['snr'])
        # Messy step because pycbc inspiral's internal 'chisq_dof' is 2p-2
        # but stat.py / ranking.py functions use 'chisq_dof' = p
        e_copy['chisq_dof'] = e_copy['chisq_dof'] / 2 + 1
        # Initialize statclass with an empty file list
        stat_instance = stat.sngl_statistic_dict[statname]([])
        statv = stat_instance.single(e_copy)

        # Convert trigger time to integer bin number
        # NB time_index and window are in units of samples
        wtime = (e_copy['time_index'] / window).astype(numpy.int32)
        bins = numpy.unique(wtime)

        if log_chirp_width:
            from pycbc.conversions import mchirp_from_mass1_mass2
            m1 = numpy.array([p['tmplt'].mass1 for p in self.template_params])
            m2 = numpy.array([p['tmplt'].mass2 for p in self.template_params])
            mc = mchirp_from_mass1_mass2(m1, m2)[e_copy['template_id']]

            # convert chirp mass to integer bin number
            imc = (numpy.log(mc) / log_chirp_width).astype(numpy.int32)
            cbins = numpy.unique(imc)

        keep = []
        for b in bins:
            if log_chirp_width:
                for b2 in cbins:
                    bloc = numpy.where((wtime == b) & (imc == b2))[0]
                    bloudest = statv[bloc].argsort()[-num_keep:]
                    keep.append(bloc[bloudest])
            else:
                bloc = numpy.where((wtime == b))[0]
                bloudest = statv[bloc].argsort()[-num_keep:]
                keep.append(bloc[bloudest])

        keep = numpy.concatenate(keep)
        self.events = self.events[keep]

    def add_template_events(self, columns, vectors):
        """ Add a vector indexed """
        # initialize with zeros - since vectors can be None, look for the
        # longest one that isn't
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
        if window_size == 0:
            indices = numpy.arange(len(tvec))
        else:
            indices = findchirp_cluster_over_window(tvec, cvec, window_size)
        self.template_events = numpy.take(self.template_events, indices)

    def new_template(self, **kwds):
        self.template_params.append(kwds)
        self.template_index += 1

    def add_template_params(self, **kwds):
        self.template_params[-1].update(kwds)

    def finalize_template_events(self):
        self.accumulate.append(self.template_events)
        self.template_events = numpy.array([], dtype=self.event_dtype)

    def finalize_events(self):
        self.events = numpy.concatenate(self.accumulate)

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
        th = numpy.array([p['tmplt'].template_hash for p in
                          self.template_params])
        tid = self.events['template_id']
        f = fw(outname, self.opt.channel_name[0:2])

        if len(self.events):
            f['snr'] = abs(self.events['snr'])
            try:
                # Precessing
                f['u_vals'] = self.events['u_vals']
                f['coa_phase'] = self.events['coa_phase']
                f['hplus_cross_corr'] = self.events['hplus_cross_corr']
            except Exception:
                # Not precessing
                f['coa_phase'] = numpy.angle(self.events['snr'])
            f['chisq'] = self.events['chisq']
            f['bank_chisq'] = self.events['bank_chisq']
            f['bank_chisq_dof'] = self.events['bank_chisq_dof']
            f['cont_chisq'] = self.events['cont_chisq']
            f['end_time'] = self.events['time_index'] / \
                              float(self.opt.sample_rate) \
                            + self.opt.gps_start_time
            try:
                # Precessing
                template_sigmasq_plus = numpy.array(
                             [t['sigmasq_plus'] for t in self.template_params],
                                                    dtype=numpy.float32)
                f['sigmasq_plus'] = template_sigmasq_plus[tid]
                template_sigmasq_cross = numpy.array(
                            [t['sigmasq_cross'] for t in self.template_params],
                                                     dtype=numpy.float32)
                f['sigmasq_cross'] = template_sigmasq_cross[tid]
                # FIXME: I want to put something here, but I haven't yet
                #        figured out what it should be. I think we would also
                #        need information from the plus and cross correlation
                #        (both real and imaginary(?)) to get this.
                f['sigmasq'] = template_sigmasq_plus[tid]
            except Exception:
                # Not precessing
                template_sigmasq = numpy.array(
                                  [t['sigmasq'] for t in self.template_params],
                                               dtype=numpy.float32)
                f['sigmasq'] = template_sigmasq[tid]

            template_durations = [p['tmplt'].template_duration for p in
                                  self.template_params]
            f['template_duration'] = numpy.array(template_durations,
                                                 dtype=numpy.float32)[tid]

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

            if 'sg_chisq' in self.events.dtype.names:
                f['sg_chisq'] = self.events['sg_chisq']

            if self.opt.psdvar_short_segment is not None:
                f['psd_var_val'] = self.events['psd_var_val']

        if self.opt.trig_start_time:
            f['search/start_time'] = numpy.array([self.opt.trig_start_time])
            search_start_time = float(self.opt.trig_start_time)
        else:
            f['search/start_time'] = numpy.array([self.opt.gps_start_time +
                                                  self.opt.segment_start_pad])
            search_start_time = float(self.opt.gps_start_time +
                                      self.opt.segment_start_pad)
        if self.opt.trig_end_time:
            f['search/end_time'] = numpy.array([self.opt.trig_end_time])
            search_end_time = float(self.opt.trig_end_time)
        else:
            f['search/end_time'] = numpy.array([self.opt.gps_end_time -
                                                self.opt.segment_end_pad])
            search_end_time = float(self.opt.gps_end_time -
                                    self.opt.segment_end_pad)

        if self.write_performance:
            self.analysis_time = search_end_time - search_start_time
            time_ratio = numpy.array(
                [float(self.analysis_time) / float(self.run_time)])
            temps_per_core = float(self.ntemplates) / float(self.ncores)
            filters_per_core = float(self.nfilters) / float(self.ncores)
            f['search/templates_per_core'] = \
                numpy.array([float(temps_per_core) * float(time_ratio)])
            f['search/filter_rate_per_core'] = \
                numpy.array([filters_per_core / float(self.run_time)])
            f['search/setup_time_fraction'] = \
                numpy.array([float(self.setup_time) / float(self.run_time)])
            f['search/run_time'] = numpy.array([float(self.run_time)])

        if 'q_trans' in self.global_params:
            qtrans = self.global_params['q_trans']
            for key in qtrans:
                if key == 'qtiles':
                    for seg in qtrans[key]:
                        for q in qtrans[key][seg]:
                            f['qtransform/%s/%s/%s' % (key, seg, q)] = \
                                                            qtrans[key][seg][q]
                elif key == 'qplanes':
                    for seg in qtrans[key]:
                        f['qtransform/%s/%s' % (key, seg)] = qtrans[key][seg]

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


class EventManagerMultiDetBase(EventManager):
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

        self.event_dtype = [('template_id', int), ('event_id', int)]
        for col, coltype in zip(column, column_types):
            self.event_dtype.append((col, coltype))

        self.events = numpy.array([], dtype=self.event_dtype)
        self.event_id_map = {}
        self.template_params = []
        self.template_index = -1
        self.template_event_dict = {}
        self.coinc_list = []
        self.write_performance = False
        for ifo in ifos:
            self.template_event_dict[ifo] = \
                numpy.array([], dtype=self.event_dtype)

    def add_template_events_to_ifo(self, ifo, columns, vectors):
        """ Add a vector indexed """
        # Just call through to the standard function
        self.template_events = self.template_event_dict[ifo]
        self.add_template_events(columns, vectors)
        self.template_event_dict[ifo] = self.template_events
        self.template_events = None


class EventManagerCoherent(EventManagerMultiDetBase):
    def __init__(self, opt, ifos, column, column_types, network_column,
                 network_column_types, psd=None, **kwargs):
        super(EventManagerCoherent, self).__init__(
            opt, ifos, column, column_types, psd=None, **kwargs)
        self.network_event_dtype = \
            [(ifo + '_event_id', int) for ifo in self.ifos]
        self.network_event_dtype.append(('template_id', int))
        self.network_event_dtype.append(('event_id', int))
        for col, coltype in zip(network_column, network_column_types):
            self.network_event_dtype.append((col, coltype))
        self.network_events = numpy.array([], dtype=self.network_event_dtype)
        self.event_index = {}
        for ifo in self.ifos:
            self.event_index[ifo] = 0
        self.event_index['network'] = 0
        self.template_event_dict['network'] = \
                                numpy.array([], dtype=self.network_event_dtype)

    def add_template_network_events(self, columns, vectors):
        """ Add a vector indexed """
        # initialize with zeros - since vectors can be None, look for the
        # longest one that isn't
        new_events = None
        new_events = numpy.zeros(
            max([len(v) for v in vectors if v is not None]),
            dtype=self.network_event_dtype
        )
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

    def add_template_events_to_network(self, columns, vectors):
        """ Add a vector indexed """
        # Just call through to the standard function
        self.template_events = self.template_event_dict['network']
        self.add_template_network_events(columns, vectors)
        self.template_event_dict['network'] = self.template_events
        self.template_events = None

    def write_to_hdf(self, outname):
        class fw(object):
            def __init__(self, name):
                import h5py
                self.f = h5py.File(name, 'w')

            def __setitem__(self, name, data):
                col = self.prefix + '/' + name
                self.f.create_dataset(col, data=data,
                                      compression='gzip',
                                      compression_opts=9,
                                      shuffle=True)

        self.events.sort(order='template_id')
        th = numpy.array([p['tmplt'].template_hash for p in
                          self.template_params])
        tid = self.events['template_id']
        f = fw(outname)
        # Output network stuff
        f.prefix = 'network'
        network_events = numpy.array([e for e in self.network_events],
                                     dtype=self.network_event_dtype)
        f['event_id'] = network_events['event_id']
        f['network_snr'] = network_events['network_snr']
        f['null_snr'] = network_events['null_snr']
        f['end_time_gc'] = network_events['time_index'] / \
                float(self.opt.sample_rate[self.ifos[0].lower()]) + \
                        self.opt.gps_start_time[self.ifos[0].lower()]
        f['nifo'] = network_events['nifo']
        f['latitude'] = network_events['latitude']
        f['longitude'] = network_events['longitude']
        f['template_id'] = network_events['template_id']
        for ifo in self.ifos:
            # First add the ifo event ids to the network branch
            f[ifo + '_event_id'] = network_events[ifo + '_event_id']
        # Individual ifo stuff
        for i, ifo in enumerate(self.ifos):
            f.prefix = ifo
            ifo_events = numpy.array([e for e in self.events
                    if e['ifo'] == self.ifo_dict[ifo]], dtype=self.event_dtype)
            if len(ifo_events):
                ifo_str = ifo.lower()[0] if ifo != 'H1' else ifo.lower()
                f['snr_%s' % ifo_str] = abs(ifo_events['snr'])
                f['event_id'] = ifo_events['event_id']
                try:
                    # Precessing
                    f['u_vals'] = ifo_events['u_vals']
                    f['coa_phase'] = ifo_events['coa_phase']
                    f['hplus_cross_corr'] = ifo_events['hplus_cross_corr']
                except Exception:
                    f['coa_phase'] = numpy.angle(ifo_events['snr'])
                f['chisq'] = ifo_events['chisq']
                f['bank_chisq'] = ifo_events['bank_chisq']
                f['bank_chisq_dof'] = ifo_events['bank_chisq_dof']
                f['cont_chisq'] = ifo_events['cont_chisq']
                f['end_time'] = ifo_events['time_index'] / \
                        float(self.opt.sample_rate[ifo_str]) + \
                        self.opt.gps_start_time[ifo_str]
                f['time_index'] = ifo_events['time_index']
                try:
                    # Precessing
                    template_sigmasq_plus = numpy.array(
                        [t['sigmasq_plus'] for t in self.template_params],
                        dtype=numpy.float32
                    )
                    f['sigmasq_plus'] = template_sigmasq_plus[tid]
                    template_sigmasq_cross = numpy.array(
                        [t['sigmasq_cross'] for t in self.template_params],
                        dtype=numpy.float32
                    )
                    f['sigmasq_cross'] = template_sigmasq_cross[tid]
                    # FIXME: I want to put something here, but I haven't yet
                    #      figured out what it should be. I think we would also
                    #      need information from the plus and cross correlation
                    #      (both real and imaginary(?)) to get this.
                    f['sigmasq'] = template_sigmasq_plus[tid]
                except Exception:
                    # Not precessing
                    template_sigmasq = numpy.array(
                             [t['sigmasq'][ifo] for t in self.template_params],
                                                   dtype=numpy.float32)
                    f['sigmasq'] = template_sigmasq[tid]

                template_durations = [p['tmplt'].template_duration for p in
                                      self.template_params]
                f['template_duration'] = numpy.array(template_durations,
                                                     dtype=numpy.float32)[tid]

                # FIXME: Can we get this value from the autochisq instance?
                # cont_dof = self.opt.autochi_number_points
                # if self.opt.autochi_onesided is None:
                #     cont_dof = cont_dof * 2
                # if self.opt.autochi_two_phase:
                #     cont_dof = cont_dof * 2
                # if self.opt.autochi_max_valued_dof:
                #     cont_dof = self.opt.autochi_max_valued_dof
                # f['cont_chisq_dof'] = numpy.repeat(cont_dof, len(ifo_events))

                if 'chisq_dof' in ifo_events.dtype.names:
                    f['chisq_dof'] = ifo_events['chisq_dof'] / 2 + 1
                else:
                    f['chisq_dof'] = numpy.zeros(len(ifo_events))

                f['template_hash'] = th[tid][self.events['ifo'] == i]

            if self.opt.trig_start_time:
                f['search/start_time'] = numpy.array([
                             self.opt.trig_start_time[ifo]], dtype=numpy.int32)
                search_start_time = float(self.opt.trig_start_time[ifo])
            else:
                f['search/start_time'] = numpy.array([
                    self.opt.gps_start_time[ifo] +
                    self.opt.segment_start_pad[ifo]], dtype=numpy.int32)
                search_start_time = float(self.opt.gps_start_time[ifo] +
                                          self.opt.segment_start_pad[ifo])
            if self.opt.trig_end_time:
                f['search/end_time'] = numpy.array([
                        self.opt.trig_end_time[ifo]], dtype=numpy.int32)
                search_end_time = float(self.opt.trig_end_time[ifo])
            else:
                f['search/end_time'] = numpy.array(
                    [self.opt.gps_end_time[ifo] -
                     self.opt.segment_end_pad[ifo]], dtype=numpy.int32)
                search_end_time = float(self.opt.gps_end_time[ifo] -
                                        self.opt.segment_end_pad[ifo])

            if self.write_performance:
                self.analysis_time = search_end_time - search_start_time
                time_ratio = numpy.array([float(self.analysis_time) /
                                          float(self.run_time)])
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
                        f['gating/' + gate_type + '/time'] = numpy.array(
                                 [float(g[0]) for g in gating_info[gate_type]])
                        f['gating/' + gate_type + '/width'] = numpy.array(
                                        [g[1] for g in gating_info[gate_type]])
                        f['gating/' + gate_type + '/pad'] = numpy.array(
                                        [g[2] for g in gating_info[gate_type]])

    def finalize_template_events(self):
        # Check that none of the template events have the same time index as an
        # existing event in events. I.e. don't list the same ifo event multiple
        # times when looping over sky points and time slides.
        existing_times = {}
        new_times = {}
        existing_template_id = {}
        new_template_id = {}
        existing_events_mask = {}
        new_template_event_mask = {}
        existing_template_event_mask = {}
        for i, ifo in enumerate(self.ifos):
            ifo_events = numpy.where(self.events['ifo'] == i)
            existing_times[ifo] = self.events['time_index'][ifo_events]
            new_times[ifo] = self.template_event_dict[ifo]['time_index']
            existing_template_id[ifo] = self.events['template_id'][ifo_events]
            new_template_id[ifo] = self.template_event_dict[ifo]['template_id']
            # This is true for each existing event that has the same time index
            # and template id as a template trigger.
            existing_events_mask[ifo] = numpy.argwhere(
                numpy.logical_and(
                    numpy.isin(existing_times[ifo], new_times[ifo]),
                    numpy.isin(existing_template_id[ifo], new_template_id[ifo])
                                 )).reshape(-1,)
            # This is true for each template event that has either a new
            # trigger time or a new template id.
            new_template_event_mask[ifo] = numpy.argwhere(
                numpy.logical_or(
                   ~numpy.isin(new_times[ifo], existing_times[ifo]),
                   ~numpy.isin(new_template_id[ifo], existing_template_id[ifo])
                                )).reshape(-1,)
            # This is true for each template event that has the same time index
            # and template id as an exisitng event trigger.
            existing_template_event_mask[ifo] = numpy.argwhere(
                numpy.logical_and(
                    numpy.isin(new_times[ifo], existing_times[ifo]),
                    numpy.isin(new_template_id[ifo], existing_template_id[ifo])
                                 )).reshape(-1,)
            # Set ids (These show how each trigger in the single ifo trigger
            # list correspond to the network triggers)
            num_events = len(new_template_event_mask[ifo])
            new_event_ids = numpy.arange(self.event_index[ifo],
                                         self.event_index[ifo] + num_events)
            # Every template event that corresponds to a new trigger gets a new
            # id. Triggers that have been found before are not saved.
            self.template_event_dict[ifo]['event_id'][
                new_template_event_mask[ifo]] = new_event_ids
            self.template_event_dict['network'][ifo + '_event_id'][
                new_template_event_mask[ifo]] = new_event_ids
            # Template events that have been found before get the event id of
            # the first time they were found.
            self.template_event_dict['network'][ifo + '_event_id'][
                  existing_template_event_mask[ifo]] = \
                self.events[self.events['ifo'] == i][
                  existing_events_mask[ifo]]['event_id']
            self.event_index[ifo] = self.event_index[ifo] + num_events

        num_events = len(self.template_event_dict['network'])
        new_event_ids = numpy.arange(self.event_index['network'],
                                     self.event_index['network'] + num_events)
        self.template_event_dict['network']['event_id'] = new_event_ids
        # Move template events for each ifo to the events list
        for ifo in self.ifos:
            self.events = numpy.append(
                self.events,
                self.template_event_dict[ifo][new_template_event_mask[ifo]]
            )
            self.template_event_dict[ifo] = \
                                        numpy.array([], dtype=self.event_dtype)
        # Move the template events for the network to the network events list
        self.network_events = numpy.append(self.network_events,
                                   self.template_event_dict['network'])
        self.template_event_dict['network'] = \
                                numpy.array([], dtype=self.network_event_dtype)


class EventManagerMultiDet(EventManagerMultiDetBase):
    def __init__(self, opt, ifos, column, column_types, psd=None, **kwargs):
        super(EventManagerMultiDet, self).__init__(
            opt, ifos, column, column_types, psd=None, **kwargs)
        self.event_index = 0

    def cluster_template_events_single_ifo(
            self, tcolumn, column, window_size, ifo):
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
            end_times1 = self.template_event_dict[ifo1]['time_index'] /\
              float(self.opt.sample_rate[ifo1]) + self.opt.gps_start_time[ifo1]
            end_times2 = self.template_event_dict[ifo2]['time_index'] /\
              float(self.opt.sample_rate[ifo2]) + self.opt.gps_start_time[ifo2]
            light_travel_time = Detector(ifo1).light_travel_time_to_detector(
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

        if '.hdf' in outname:
            self.write_to_hdf(outname)
        else:
            raise ValueError('Cannot write to this format')

    def write_to_hdf(self, outname):
        class fw(object):
            def __init__(self, name):
                import h5py
                self.f = h5py.File(name, 'w')

            def __setitem__(self, name, data):
                col = self.prefix + '/' + name
                self.f.create_dataset(col, data=data,
                                      compression='gzip',
                                      compression_opts=9,
                                      shuffle=True)

        self.events.sort(order='template_id')
        th = numpy.array([p['tmplt'].template_hash for p in
                                                         self.template_params])
        tid = self.events['template_id']
        f = fw(outname)
        for ifo in self.ifos:
            f.prefix = ifo
            ifo_events = numpy.array([e for e in self.events if
                                      e['ifo'] == self.ifo_dict[ifo]],
                                     dtype=self.event_dtype)
            if len(ifo_events):
                ifo_str = ifo.lower()[0] if ifo != 'H1' else ifo.lower()
                f['snr_%s' % ifo_str] = abs(ifo_events['snr'])
                try:
                    # Precessing
                    f['u_vals'] = ifo_events['u_vals']
                    f['coa_phase'] = ifo_events['coa_phase']
                    f['hplus_cross_corr'] = ifo_events['hplus_cross_corr']
                except Exception:
                    f['coa_phase'] = numpy.angle(ifo_events['snr'])
                f['chisq'] = ifo_events['chisq']
                f['bank_chisq'] = ifo_events['bank_chisq']
                f['bank_chisq_dof'] = ifo_events['bank_chisq_dof']
                f['cont_chisq'] = ifo_events['cont_chisq']
                f['end_time'] = ifo_events['time_index'] / \
                        float(self.opt.sample_rate[ifo_str]) + \
                          self.opt.gps_start_time[ifo_str]
                try:
                    # Precessing
                    template_sigmasq_plus = numpy.array([t['sigmasq_plus'] for
                               t in self.template_params], dtype=numpy.float32)
                    f['sigmasq_plus'] = template_sigmasq_plus[tid]
                    template_sigmasq_cross = numpy.array([t['sigmasq_cross']
                           for t in self.template_params], dtype=numpy.float32)
                    f['sigmasq_cross'] = template_sigmasq_cross[tid]
                    # FIXME: I want to put something here, but I haven't yet
                    #      figured out what it should be. I think we would also
                    #      need information from the plus and cross correlation
                    #      (both real and imaginary(?)) to get this.
                    f['sigmasq'] = template_sigmasq_plus[tid]
                except Exception:
                    # Not precessing
                    template_sigmasq = numpy.array([t['sigmasq'][ifo] for t in
                                                    self.template_params],
                                                   dtype=numpy.float32)
                    f['sigmasq'] = template_sigmasq[tid]

                template_durations = [p['tmplt'].template_duration for p in
                                      self.template_params]
                f['template_duration'] = \
                      numpy.array(template_durations, dtype=numpy.float32)[tid]

                # FIXME: Can we get this value from the autochisq instance?
                cont_dof = self.opt.autochi_number_points
                if self.opt.autochi_onesided is None:
                    cont_dof = cont_dof * 2
                # if self.opt.autochi_two_phase:
                #     cont_dof = cont_dof * 2
                # if self.opt.autochi_max_valued_dof:
                #     cont_dof = self.opt.autochi_max_valued_dof
                f['cont_chisq_dof'] = numpy.repeat(cont_dof, len(ifo_events))

                if 'chisq_dof' in ifo_events.dtype.names:
                    f['chisq_dof'] = ifo_events['chisq_dof'] / 2 + 1
                else:
                    f['chisq_dof'] = numpy.zeros(len(ifo_events))

                f['template_hash'] = th[tid]

                if self.opt.psdvar_short_segment is not None:
                    f['psd_var_val'] = ifo_events['psd_var_val']

            if self.opt.trig_start_time:
                f['search/start_time'] = numpy.array(
                            [self.opt.trig_start_time[ifo]], dtype=numpy.int32)
                search_start_time = float(self.opt.trig_start_time[ifo])
            else:
                f['search/start_time'] = numpy.array(
                          [self.opt.gps_start_time[ifo] +
                           self.opt.segment_start_pad[ifo]], dtype=numpy.int32)
                search_start_time = float(self.opt.gps_start_time[ifo] +
                                          self.opt.segment_start_pad[ifo])
            if self.opt.trig_end_time:
                f['search/end_time'] = numpy.array(
                              [self.opt.trig_end_time[ifo]], dtype=numpy.int32)
                search_end_time = float(self.opt.trig_end_time[ifo])
            else:
                f['search/end_time'] = numpy.array(
                            [self.opt.gps_end_time[ifo] -
                             self.opt.segment_end_pad[ifo]], dtype=numpy.int32)
                search_end_time = float(self.opt.gps_end_time[ifo] -
                                        self.opt.segment_end_pad[ifo])

            if self.write_performance:
                self.analysis_time = search_end_time - search_start_time
                time_ratio = numpy.array(
                            [float(self.analysis_time) / float(self.run_time)])
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
                        f['gating/' + gate_type + '/time'] = numpy.array(
                                 [float(g[0]) for g in gating_info[gate_type]])
                        f['gating/' + gate_type + '/width'] = numpy.array(
                                        [g[1] for g in gating_info[gate_type]])
                        f['gating/' + gate_type + '/pad'] = numpy.array(
                                        [g[2] for g in gating_info[gate_type]])


__all__ = ['threshold_and_cluster', 'findchirp_cluster_over_window',
           'threshold', 'cluster_reduce', 'ThresholdCluster',
           'threshold_real_numpy', 'threshold_only',
           'EventManager', 'EventManagerMultiDet', 'EventManagerCoherent']
