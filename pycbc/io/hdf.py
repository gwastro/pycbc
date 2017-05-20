# convenience classes for accessing hdf5 trigger files
# the 'get_column()' method is implemented parallel to 
# the existing pylal.SnglInspiralUtils functions

import h5py
import numpy as np
import logging
import inspect

from lal import LIGOTimeGPS, YRJUL_SI

from pycbc_glue.ligolw import ligolw
from pycbc_glue.ligolw import table
from pycbc_glue.ligolw import lsctables
from pycbc_glue.ligolw import ilwd
from pycbc_glue.ligolw import utils as ligolw_utils
from pycbc_glue.ligolw.utils import process as ligolw_process

from pycbc import version as pycbc_version
from pycbc.tmpltbank import return_search_summary
from pycbc.tmpltbank import return_empty_sngl
from pycbc import events, pnutils

class HFile(h5py.File):   
    """ Low level extensions to the capabilities of reading an hdf5 File
    """

    def select(self, fcn, *args, **kwds):
        """ Return arrays from an hdf5 file that satisfy the given function

        Parameters
        ----------
        fcn : a function
            A function that accepts the same number of argument as keys given
            and returns a boolean array of the same length.

        args : strings
            A variable number of strings that are keys into the hdf5. These must
            refer to arrays of equal length.

        chunksize : {1e6, int}, optional
            Number of elements to read and process at a time.

        return_indices : bool, optional
            If True, also return the indices of elements passing the function.

        Returns
        -------
        values : np.ndarrays
            A variable number of arrays depending on the number of keys into
            the hdf5 file that are given. If return_indices is True, the first
            element is an array of indices of elements passing the function.

        >>> f = HFile(filename)
        >>> snr = f.select(lambda snr: snr > 6, 'H1/snr')
        """

        # get references to each array
        refs = {}
        data = {}
        for arg in args:
            refs[arg] = self[arg]
            data[arg] = []

        return_indices = kwds.get('return_indices', False)
        indices = np.array([], dtype=np.uint64)

        # To conserve memory read the array in chunks
        chunksize = kwds.get('chunksize', int(1e6))
        size = len(refs[arg])

        i = 0
        while i < size:
            r = i + chunksize if i + chunksize < size else size

            #Read each chunks worth of data and find where it passes the function
            partial = [refs[arg][i:r] for arg in args]
            keep = fcn(*partial)
            if return_indices:
                indices = np.concatenate([indices, np.flatnonzero(keep) + i])

            #store only the results that pass the function
            for arg, part in zip(args, partial):
                data[arg].append(part[keep])

            i += chunksize

        # Combine the partial results into full arrays
        if len(args) == 1:
            res = np.concatenate(data[args[0]])
            if return_indices:
                return indices, res
            else:
                return res
        else:
            res = tuple(np.concatenate(data[arg]) for arg in args)
            if return_indices:
                return (indices,) + res
            else:
                return res


class DictArray(object):
    """ Utility for organizing sets of arrays of equal length. 
    
    Manages a dictionary of arrays of equal length. This can also
    be instantiated with a set of hdf5 files and the key values. The full
    data is always in memory and all operations create new instances of the
    DictArray.
    """
    def __init__(self, data=None, files=None, groups=None):
        """ Create a DictArray
        
        Parameters
        ----------
        data: dict, optional
            Dictionary of equal length numpy arrays
        files: list of filenames, optional
            List of hdf5 file filenames. Incompatibile with the `data` option.
        groups: list of strings
            List of keys into each file. Required by the files options.
        """
        self.data = data
        
        if files:
            self.data = {}
            for g in groups:
                self.data[g] = []
            
            for f in files:
                d = HFile(f)
                for g in groups:
                    if g in d:
                        self.data[g].append(d[g][:])
                d.close()
                    
            for k in self.data:
                if not len(self.data[k]) == 0:
                    self.data[k] = np.concatenate(self.data[k])

        for k in self.data:
            setattr(self, k, self.data[k])

    def _return(self, data):
        return self.__class__(data=data)

    def __len__(self):
        return len(self.data[self.data.keys()[0]])

    def __add__(self, other):
        data = {}
        for k in self.data:
            data[k] = np.concatenate([self.data[k], other.data[k]])
        return self._return(data=data)

    def select(self, idx):
        """ Return a new DictArray containing only the indexed values
        """
        data = {}
        for k in self.data:
            data[k] = self.data[k][idx]
        return self._return(data=data)
   
    def remove(self, idx):
        """ Return a new DictArray that does not contain the indexed values
        """
        data = {}
        for k in self.data:
            data[k] = np.delete(self.data[k], idx)
        return self._return(data=data)


class StatmapData(DictArray):
    def __init__(self, data=None, seg=None, attrs=None,
                       files=None):
        groups = ['stat', 'time1', 'time2', 'trigger_id1', 'trigger_id2', 
               'template_id', 'decimation_factor', 'timeslide_id']
        super(StatmapData, self).__init__(data=data, files=files, groups=groups)
        
        if data:
            self.seg=seg
            self.attrs=attrs
        elif files:
            f = HFile(files[0], "r")
            self.seg = f['segments']
            self.attrs = f.attrs

    def _return(self, data):
        return self.__class__(data=data, attrs=self.attrs, seg=self.seg)

    def cluster(self, window):
        """ Cluster the dict array, assuming it has the relevant Coinc colums,
        time1, time2, stat, and timeslide_id
        """
        # If no events, do nothing
        if len(self.time1) == 0 or len(self.time2) == 0:
            return self
        from pycbc.events import cluster_coincs
        interval = self.attrs['timeslide_interval']
        cid = cluster_coincs(self.stat, self.time1, self.time2,
                                 self.timeslide_id, interval, window)
        return self.select(cid) 

    def save(self, outname):
        f = HFile(outname, "w")
        for k in self.attrs:
            f.attrs[k] = self.attrs[k]
            
        for k in self.data:
            f.create_dataset(k, data=self.data[k], 
                      compression='gzip',
                      compression_opts=9,
                      shuffle=True)
            

        for key in self.seg.keys():
            f['segments/%s/start' % key] = self.seg[key]['start'][:]
            f['segments/%s/end' % key] = self.seg[key]['end'][:]
        f.close()

class FileData(object):

    def __init__(self, fname, group=None, columnlist=None, filter_func=None):
        """
        Parameters
        ----------
        group : string
            Name of group to be read from the file
        columnlist : list of strings
            Names of columns to be read; if None, use all existing columns 
        filter_func : string 
            String should evaluate to a Boolean expression using attributes
            of the class instance derived from columns: ex. 'self.snr < 6.5'
        """
        if not fname: raise RuntimeError("Didn't get a file!")
        self.fname = fname
        self.h5file = HFile(fname, "r")
        if group is None:
            if len(self.h5file.keys()) == 1:
                group = self.h5file.keys()[0]
            else:
                raise RuntimeError("Didn't get a group!")
        self.group_key = group
        self.group = self.h5file[group]
        self.columns = columnlist if columnlist is not None \
                       else self.group.keys()
        self.filter_func = filter_func
        self._mask = None

    def close(self):
        self.h5file.close()

    @property
    def mask(self):
        """
        Create a mask implementing the requested filter on the datasets

        Returns
        -------
        array of Boolean
            True for dataset indices to be returned by the get_column method
        """
        if self.filter_func is None:
            raise RuntimeError("Can't get a mask without a filter function!")
        else:
            # only evaluate if no previous calculation was done
            if self._mask is None:
                # get required columns into the namespace as numpy arrays
                for column in self.columns:
                    if column in self.filter_func:
                        setattr(self, column, self.group[column][:])
                self._mask = eval(self.filter_func)
            return self._mask

    def get_column(self, col):
        """
        Parameters
        ----------
        col : string
            Name of the dataset to be returned

        Returns
        -------
        numpy array
            Values from the dataset, filtered if requested
        """
        # catch corner case with an empty file (group with no datasets)
        if not len(self.group.keys()):
            return np.array([])
        vals = self.group[col]
        if self.filter_func:
            return vals[self.mask]
        else:
            return vals[:]


class DataFromFiles(object):

    def __init__(self, filelist, group=None, columnlist=None, filter_func=None):
        self.files = filelist
        self.group = group
        self.columns = columnlist
        self.filter_func = filter_func

    def get_column(self, col):
        """
        Loop over files getting the requested dataset values from each

        Parameters
        ----------
        col : string
            Name of the dataset to be returned

        Returns
        -------
        numpy array
            Values from the dataset, filtered if requested and
            concatenated in order of file list
        """
        logging.info('getting %s' % col)
        vals = []
        for f in self.files:
            d = FileData(f, group=self.group, columnlist=self.columns,
                         filter_func=self.filter_func)
            vals.append(d.get_column(col))
            # Close each file since h5py has an upper limit on the number of
            # open file objects (approx. 1000)
            d.close()
        logging.info('- got %i values' % sum(len(v) for v in vals))
        return np.concatenate(vals)

class SingleDetTriggers(object):
    """
    Provides easy access to the parameters of single-detector CBC triggers.
    """
    # FIXME: Some of these are optional and should be kwargs.
    def __init__(self, trig_file, bank_file, veto_file, segment_name, filter_func, detector):
        logging.info('Loading triggers')
        self.trigs_f = HFile(trig_file, 'r')
        self.trigs = self.trigs_f[detector]
        if bank_file:
            logging.info('Loading bank')
            self.bank = HFile(bank_file, 'r')
        else:
            logging.info('No bank file given')
            # empty dict in place of non-existent hdf file
            self.bank = {}

        if veto_file:
            logging.info('Applying veto segments')
            # veto_mask is an array of indices into the trigger arrays
            # giving the surviving triggers 
            logging.info('%i triggers before vetoes',
                          len(self.trigs['end_time'][:]))
            self.veto_mask, _ = events.veto.indices_outside_segments(
                self.trigs['end_time'][:], [veto_file],
                ifo=detector, segment_name=segment_name)
            logging.info('%i triggers remain after vetoes',
                          len(self.veto_mask))
        else:
            self.veto_mask = np.arange(len(self.trigs['end_time']))

        if filter_func:
            # get required columns into the namespace with dummy attribute
            # names to avoid confusion with other class properties
            for c in self.trigs.keys():
                if c in filter_func:
                    setattr(self, '_'+c, self.trigs[c][:])
            for c in self.bank.keys():
                if c in filter_func:
                    # get template parameters corresponding to triggers
                    setattr(self, '_'+c,
                          np.array(self.bank[c])[self.trigs['template_id'][:]])
            self.filter_mask = eval(filter_func.replace('self.', 'self._'))
            # remove the dummy attributes
            for c in self.trigs.keys() + self.bank.keys():
                if c in filter_func: delattr(self, '_'+c)
            self.boolean_veto = np.in1d(np.arange(len(self.trigs['end_time'])),
                  self.veto_mask, assume_unique=True)
            self.mask = np.logical_and(self.boolean_veto, self.filter_mask)
            logging.info('%i triggers remain after cut on %s',
                          len(self.trigs['end_time'][self.mask]), filter_func)
        else:
            self.mask = self.veto_mask

    def checkbank(self, param):
        if self.bank == {}:
            return RuntimeError("Can't get %s values without a bank file" 
                                                                       % param)

    @classmethod
    def get_param_names(cls):
        """Returns a list of plottable CBC parameter variables"""
        return [m[0] for m in inspect.getmembers(cls) \
            if type(m[1]) == property]

    def mask_to_n_loudest_clustered_events(self, n_loudest=10,
                                           ranking_statistic="newsnr",
                                           cluster_window=10):
        """Edits the mask property of the class to point to the N loudest
        single detector events as ranked by ranking statistic. Events are
        clustered so that no more than 1 event within +/- cluster-window will
        be considered."""
        # If this becomes memory intensive we can optimize
        if ranking_statistic == "newsnr":
            stat = self.newsnr
            # newsnr doesn't return an array if len(stat) == 1
            if len(self.snr) == 1:
                stat = np.array([stat])
            self.stat_name = "Reweighted SNR"
        elif ranking_statistic == "newsnr_sgveto":
            stat = self.newsnr_sgveto
            # newsnr doesn't return an array if len(stat) == 1
            if len(self.snr) == 1:
                stat = np.array([stat])
            self.stat_name = "Reweighted SNR (+sgveto)"
        elif ranking_statistic == "snr":
            stat = self.snr
            self.stat_name = "SNR"
        else:
            err_msg = "Don't recognize statistic %s." % (ranking_statistic)
            raise ValueError(err_msg)
        times = self.end_time
        index = stat.argsort()[::-1]
        new_times = []
        new_index = []
        for curr_idx in index:
            curr_time = times[curr_idx] 
            for time in new_times:
                if abs(curr_time - time) < cluster_window:
                    break
            else:
                # Only get here if no other triggers within cluster window
                new_index.append(curr_idx)
                new_times.append(curr_time)
            if len(new_index) >= n_loudest:
                break
        index = np.array(new_index)
        self.stat = stat[index]
        if self.mask.dtype == 'bool':
            orig_indices = self.mask.nonzero()[0][index]
            self.mask = np.in1d(np.arange(len(self.mask)), orig_indices,
                                assume_unique=True)
        else:
            self.mask = self.mask[index]

    @property
    def template_id(self):
        return np.array(self.trigs['template_id'])[self.mask]

    @property
    def mass1(self):
        self.checkbank('mass1')
        return np.array(self.bank['mass1'])[self.template_id]

    @property
    def mass2(self):
        self.checkbank('mass2')
        return np.array(self.bank['mass2'])[self.template_id]

    @property
    def spin1z(self):
        self.checkbank('spin1z')
        return np.array(self.bank['spin1z'])[self.template_id]

    @property
    def spin2z(self):
        self.checkbank('spin2z')
        return np.array(self.bank['spin2z'])[self.template_id]

    @property
    def spin2x(self):
        self.checkbank('spin2x')
        return np.array(self.bank['spin2x'])[self.template_id]

    @property
    def spin2y(self):
        self.checkbank('spin2y')
        return np.array(self.bank['spin2y'])[self.template_id]

    @property
    def spin1x(self):
        self.checkbank('spin1x')
        return np.array(self.bank['spin1x'])[self.template_id]

    @property
    def spin1y(self):
        self.checkbank('spin1y')
        return np.array(self.bank['spin1y'])[self.template_id]

    @property
    def inclination(self):
        self.checkbank('inclination')
        return np.array(self.bank['inclination'])[self.template_id]

    @property
    def f_lower(self):
        self.checkbank('f_lower')
        return np.array(self.bank['f_lower'])[self.template_id]

    @property
    def mtotal(self):
        return self.mass1 + self.mass2

    @property
    def mchirp(self):
        mchirp, _ = pnutils.mass1_mass2_to_mchirp_eta(
            self.mass1, self.mass2)
        return mchirp

    @property
    def eta(self):
        _, eta = pnutils.mass1_mass2_to_mchirp_eta(
            self.mass1, self.mass2)
        return eta

    @property
    def effective_spin(self):
        # FIXME assumes aligned spins
        return pnutils.phenomb_chi(self.mass1, self.mass2,
                                   self.spin1z, self.spin2z)

    # IMPROVEME: would like to have a way to access all get_freq and/or
    # other pnutils.* names rather than hard-coding each one
    # - eg make this part of a fancy interface to the bank file ?
    @property
    def f_seobnrv2_peak(self):
        return pnutils.get_freq('fSEOBNRv2Peak', self.mass1, self.mass2,
                                self.spin1z, self.spin2z)

    @property
    def f_seobnrv4_peak(self):
        return pnutils.get_freq('fSEOBNRv4Peak', self.mass1, self.mass2,
                                self.spin1z, self.spin2z)

    @property
    def end_time(self):
        return np.array(self.trigs['end_time'])[self.mask]

    @property
    def template_duration(self):
        return np.array(self.trigs['template_duration'])[self.mask]

    @property
    def snr(self):
        return np.array(self.trigs['snr'])[self.mask]

    @property
    def sgchisq(self):
        return np.array(self.trigs['sg_chisq'])[self.mask]

    @property
    def u_vals(self):
        return np.array(self.trigs['u_vals'])[self.mask]

    @property
    def rchisq(self):
        return np.array(self.trigs['chisq'])[self.mask] \
            / (np.array(self.trigs['chisq_dof'])[self.mask] * 2 - 2)

    @property
    def newsnr(self):
        return events.newsnr(self.snr, self.rchisq)

    @property
    def newsnr_sgveto(self):
        return events.newsnr_sgveto(self.snr, self.rchisq, self.sgchisq)

    def get_column(self, cname):
        if hasattr(self, cname):
            return getattr(self, cname)
        else:
            return np.array(self.trigs[cname])[self.mask]


class ForegroundTriggers(object):
    # FIXME: A lot of this is hardcoded to expect two ifos
    def __init__(self, coinc_file, bank_file, sngl_files=None, n_loudest=None,
                     group='foreground'):
        self.coinc_file = FileData(coinc_file, group=group)
        self.sngl_files = {}
        if sngl_files is not None:
            for file in sngl_files:
                curr_dat = FileData(file)
                curr_ifo = curr_dat.group_key
                self.sngl_files[curr_ifo] = curr_dat
        self.bank_file = HFile(bank_file, "r")
        self.n_loudest = n_loudest

        self._sort_arr = None
        self._template_id = None
        self._trig_ids = None

    @property
    def sort_arr(self):
        if self._sort_arr is None:
            ifar = self.coinc_file.get_column('ifar')
            sorting = ifar.argsort()[::-1]
            if self.n_loudest:
                sorting = sorting[:self.n_loudest]
            self._sort_arr = sorting
        return self._sort_arr

    @property
    def template_id(self):
        if self._template_id is None:
            template_id = self.get_coincfile_array('template_id')
            self._template_id = template_id
        return self._template_id

    @property
    def trig_id(self):
        if self._trig_ids is not None:
            return self._trig_ids
        self._trig_ids = {}
        # FIXME: There is no clear mapping from trig_id to ifo. This is bad!!!
        #        for now a hack is in place.
        ifo1 = self.coinc_file.h5file.attrs['detector_1']
        ifo2 = self.coinc_file.h5file.attrs['detector_2']
        trigid1 = self.get_coincfile_array('trigger_id1')
        trigid2 = self.get_coincfile_array('trigger_id2')
        self._trig_ids[ifo1] = trigid1
        self._trig_ids[ifo2] = trigid2
        return self._trig_ids

    def get_coincfile_array(self, variable):
        return self.coinc_file.get_column(variable)[self.sort_arr]

    def get_bankfile_array(self, variable):
        try:
            return self.bank_file[variable][:][self.template_id]
        except IndexError:
            if len(self.template_id) == 0:
                return np.array([])
            raise

    def get_snglfile_array_dict(self, variable):
        return_dict = {}
        for ifo in self.sngl_files.keys():
            try:
                curr = self.sngl_files[ifo].get_column(variable)[\
                                                             self.trig_id[ifo]]
            except IndexError:
                if len(self.trig_id[ifo]) == 0:
                    curr = np.array([])
                else:
                    raise
            return_dict[ifo] = curr
        return return_dict

    def to_coinc_xml_object(self, file_name):
        # FIXME: This function will only work with two ifos!!

        outdoc = ligolw.Document()
        outdoc.appendChild(ligolw.LIGO_LW())

        ifos = [ifo for ifo in self.sngl_files.keys()]
        proc_id = ligolw_process.register_to_xmldoc(outdoc, 'pycbc',
                     {}, ifos=ifos, comment='', version=pycbc_version.git_hash,
                     cvs_repository='pycbc/'+pycbc_version.git_branch,
                     cvs_entry_time=pycbc_version.date).process_id

        search_summ_table = lsctables.New(lsctables.SearchSummaryTable)
        coinc_h5file = self.coinc_file.h5file
        start_time = coinc_h5file['segments']['coinc']['start'][:].min()
        end_time = coinc_h5file['segments']['coinc']['end'][:].max()
        num_trigs = len(self.sort_arr)
        search_summary = return_search_summary(start_time, end_time,
                                                 num_trigs, ifos)
        search_summ_table.append(search_summary)
        outdoc.childNodes[0].appendChild(search_summ_table)

        sngl_inspiral_table = lsctables.New(lsctables.SnglInspiralTable)
        coinc_def_table = lsctables.New(lsctables.CoincDefTable)
        coinc_event_table = lsctables.New(lsctables.CoincTable)
        coinc_inspiral_table = lsctables.New(lsctables.CoincInspiralTable)
        coinc_event_map_table = lsctables.New(lsctables.CoincMapTable)
        time_slide_table = lsctables.New(lsctables.TimeSlideTable)

        # Set up time_slide table
        time_slide_id = lsctables.TimeSlideID(0)
        for ifo in ifos:
            time_slide_row = lsctables.TimeSlide()
            time_slide_row.instrument = ifo
            time_slide_row.time_slide_id = time_slide_id
            time_slide_row.offset = 0
            time_slide_row.process_id = proc_id
            time_slide_table.append(time_slide_row)

        # Set up coinc_definer table
        coinc_def_id = lsctables.CoincDefID(0)
        coinc_def_row = lsctables.CoincDef()
        coinc_def_row.search = "inspiral"
        coinc_def_row.description = "sngl_inspiral<-->sngl_inspiral coincidences"
        coinc_def_row.coinc_def_id = coinc_def_id
        coinc_def_row.search_coinc_type = 0
        coinc_def_table.append(coinc_def_row)

        bank_col_names = ['mass1', 'mass2', 'spin1z', 'spin2z']
        bank_col_vals = {}
        for name in bank_col_names:
            bank_col_vals[name] = self.get_bankfile_array(name)

        coinc_event_names = ['ifar', 'time1', 'fap', 'stat']
        coinc_event_vals = {}
        for name in coinc_event_names:
            coinc_event_vals[name] = self.get_coincfile_array(name)

        sngl_col_names = ['snr', 'chisq', 'chisq_dof', 'bank_chisq',
                          'bank_chisq_dof', 'cont_chisq', 'cont_chisq_dof',
                          'end_time', 'template_duration', 'coa_phase',
                          'sigmasq']
        sngl_col_vals = {}
        for name in sngl_col_names:
            sngl_col_vals[name] = self.get_snglfile_array_dict(name)

        for idx in xrange(len(self.sort_arr)):
            # Set up IDs and mapping values
            coinc_id = lsctables.CoincID(idx)

            # Set up sngls
            # FIXME: As two-ifo is hardcoded loop over all ifos
            sngl_combined_mchirp = 0
            sngl_combined_mtot = 0
            for ifo in ifos:
                sngl_id = self.trig_id[ifo][idx]
                event_id = lsctables.SnglInspiralID(sngl_id)
                sngl = return_empty_sngl()
                sngl.event_id = event_id
                sngl.ifo = ifo
                for name in sngl_col_names:
                    val = sngl_col_vals[name][ifo][idx]
                    if name == 'end_time':
                        sngl.set_end(LIGOTimeGPS(val))
                    else:
                        setattr(sngl, name, val)
                for name in bank_col_names:
                    val = bank_col_vals[name][idx]
                    setattr(sngl, name, val)
                sngl.mtotal, sngl.eta = pnutils.mass1_mass2_to_mtotal_eta(
                        sngl.mass1, sngl.mass2)
                sngl.mchirp, _ = pnutils.mass1_mass2_to_mchirp_eta(
                        sngl.mass1, sngl.mass2)
                sngl.eff_distance = (sngl.sigmasq)**0.5 / sngl.snr
                sngl_combined_mchirp += sngl.mchirp
                sngl_combined_mtot += sngl.mtotal

                sngl_inspiral_table.append(sngl)

                # Set up coinc_map entry
                coinc_map_row = lsctables.CoincMap()
                coinc_map_row.table_name = 'sngl_inspiral'
                coinc_map_row.coinc_event_id = coinc_id
                coinc_map_row.event_id = event_id
                coinc_event_map_table.append(coinc_map_row)

            sngl_combined_mchirp = sngl_combined_mchirp / len(ifos)
            sngl_combined_mtot = sngl_combined_mtot / len(ifos)

            # Set up coinc inspiral and coinc event tables
            coinc_event_row = lsctables.Coinc()
            coinc_inspiral_row = lsctables.CoincInspiral()
            coinc_event_row.coinc_def_id = coinc_def_id
            coinc_event_row.nevents = len(ifos)
            coinc_event_row.instruments = ','.join(ifos)
            coinc_inspiral_row.set_ifos(ifos)
            coinc_event_row.time_slide_id = time_slide_id
            coinc_event_row.process_id = proc_id
            coinc_event_row.coinc_event_id = coinc_id
            coinc_inspiral_row.coinc_event_id = coinc_id
            coinc_inspiral_row.mchirp = sngl_combined_mchirp
            coinc_inspiral_row.mass = sngl_combined_mtot
            coinc_inspiral_row.set_end(\
                                   LIGOTimeGPS(coinc_event_vals['time1'][idx]))
            coinc_inspiral_row.snr = coinc_event_vals['stat'][idx]
            coinc_inspiral_row.false_alarm_rate = coinc_event_vals['fap'][idx]
            coinc_inspiral_row.combined_far = 1./coinc_event_vals['ifar'][idx]
            # Transform to Hz
            coinc_inspiral_row.combined_far = \
                                    coinc_inspiral_row.combined_far / YRJUL_SI
            coinc_event_row.likelihood = 0.
            coinc_inspiral_row.minimum_duration = 0.
            coinc_event_table.append(coinc_event_row)
            coinc_inspiral_table.append(coinc_inspiral_row)

        outdoc.childNodes[0].appendChild(coinc_def_table)
        outdoc.childNodes[0].appendChild(coinc_event_table)
        outdoc.childNodes[0].appendChild(coinc_event_map_table)
        outdoc.childNodes[0].appendChild(time_slide_table)
        outdoc.childNodes[0].appendChild(coinc_inspiral_table)
        outdoc.childNodes[0].appendChild(sngl_inspiral_table)

        ligolw_utils.write_filename(outdoc, file_name)

chisq_choices = ['traditional', 'cont', 'bank', 'max_cont_trad', 'sg',
                 'max_bank_cont', 'max_bank_trad', 'max_bank_cont_trad']

def get_chisq_from_file_choice(hdfile, chisq_choice):
    f = hdfile
    if chisq_choice in ['traditional','max_cont_trad', 'max_bank_trad',
                             'max_bank_cont_trad']:
        trad_chisq = f['chisq'][:]
        # We now need to handle the case where chisq is not actually calculated
        # 0 is used as a sentinel value
        trad_chisq_dof = f['chisq_dof'][:]
        trad_chisq /= (trad_chisq_dof * 2 - 2)
    if chisq_choice in ['cont', 'max_cont_trad', 'max_bank_cont',
                             'max_bank_cont_trad']:
        cont_chisq = f['cont_chisq'][:]
        cont_chisq_dof = f['cont_chisq_dof'][:]
        cont_chisq /= cont_chisq_dof
    if chisq_choice in ['bank', 'max_bank_cont', 'max_bank_trad',
                             'max_bank_cont_trad']:
        bank_chisq = f['bank_chisq'][:]
        bank_chisq_dof = f['bank_chisq_dof'][:]
        bank_chisq /= bank_chisq_dof
    if chisq_choice == 'sg':
        chisq = f['sg_chisq'][:]
    elif chisq_choice == 'traditional':
        chisq = trad_chisq
    elif chisq_choice == 'cont':
        chisq = cont_chisq
    elif chisq_choice == 'bank':
        chisq = bank_chisq
    elif chisq_choice == 'max_cont_trad':
        chisq = np.maximum(trad_chisq, cont_chisq)
    elif chisq_choice == 'max_bank_cont':
        chisq = np.maximum(bank_chisq, cont_chisq)
    elif chisq_choice == 'max_bank_trad':
        chisq = np.maximum(bank_chisq, trad_chisq)
    elif chisq_choice == 'max_bank_cont_trad':
        chisq = np.maximum(np.maximum(bank_chisq, cont_chisq), trad_chisq)
    else:
        err_msg="Do not recognized --chisq-choice %s" % chisq_choice
        raise ValueError(err_msg)
    return chisq

