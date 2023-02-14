# Convenience classes for accessing hdf5 trigger files
# The 'get_column()' method is implemented parallel to
# legacy pylal.SnglInspiralUtils functions

import h5py
import numpy as np
import logging
import inspect
import pickle

from itertools import chain
from io import BytesIO
from lal import LIGOTimeGPS, YRJUL_SI

from ligo.lw import ligolw
from ligo.lw import lsctables
from ligo.lw import utils as ligolw_utils
from ligo.lw.utils import process as ligolw_process

from pycbc import version as pycbc_version
from pycbc.io.ligolw import return_search_summary, return_empty_sngl
from pycbc import events, conversions, pnutils
from pycbc.events import ranking, veto
from pycbc.events import mean_if_greater_than_zero
from pycbc.pnutils import mass1_mass2_to_mchirp_eta


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

            # Read each chunk's worth of data and find where it passes the function
            partial = [refs[arg][i:r] for arg in args]
            keep = fcn(*partial)
            if return_indices:
                indices = np.concatenate([indices, np.flatnonzero(keep) + i])

            # Store only the results that pass the function
            for arg, part in zip(args, partial):
                data[arg].append(part[keep])

            i += chunksize

        # Combine the partial results into full arrays
        if len(args) == 1:
            res = np.concatenate(data[args[0]])
            if return_indices:
                return indices.astype(np.uint64), res
            else:
                return res
        else:
            res = tuple(np.concatenate(data[arg]) for arg in args)
            if return_indices:
                return (indices.astype(np.uint64),) + res
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
            List of keys into each file. Required by the files option.
        """
        # Check that input fits with how the DictArray is set up
        if data and files:
            raise RuntimeError('DictArray can only have data or files as '
                               'input, not both.')
        if files and not groups:
            raise RuntimeError('If files are given then need groups.')

        self.data = data
        self.groups = groups
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
        return len(self.data[tuple(self.data.keys())[0]])

    def __add__(self, other):
        data = {}
        for k in self.data:
            try:
                data[k] = np.concatenate([self.data[k], other.data[k]])
            except KeyError:
                logging.info('%s does not exist in other data' % k)
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
            data[k] = np.delete(self.data[k], np.array(idx, dtype=int))
        return self._return(data=data)

    def save(self, outname):
        f = HFile(outname, "w")
        for k in self.attrs:
            f.attrs[k] = self.attrs[k]

        for k in self.data:
            f.create_dataset(k, data=self.data[k],
                      compression='gzip',
                      compression_opts=9,
                      shuffle=True)
        f.close()


class StatmapData(DictArray):
    def __init__(self, data=None, seg=None, attrs=None, files=None,
                 groups=('stat', 'time1', 'time2', 'trigger_id1',
                         'trigger_id2', 'template_id', 'decimation_factor',
                         'timeslide_id')):
        super(StatmapData, self).__init__(data=data, files=files,
                                          groups=groups)

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
        super(StatmapData, self).save(outname)
        with HFile(outname, "w") as f:
            for key in self.seg.keys():
                f['segments/%s/start' % key] = self.seg[key]['start'][:]
                f['segments/%s/end' % key] = self.seg[key]['end'][:]


class MultiifoStatmapData(StatmapData):
    def __init__(self, data=None, seg=None, attrs=None,
                       files=None, ifos=None):
        groups = ['decimation_factor', 'stat', 'template_id', 'timeslide_id']
        for ifo in ifos:
            groups += ['%s/time' % ifo]
            groups += ['%s/trigger_id' % ifo]

        super(MultiifoStatmapData, self).__init__(data=data, files=files,
                                                  groups=groups, attrs=attrs,
                                                  seg=seg)

    def _return(self, data):
        ifolist = self.attrs['ifos'].split(' ')
        return self.__class__(data=data, attrs=self.attrs, seg=self.seg,
                              ifos=ifolist)

    def cluster(self, window):
        """ Cluster the dict array, assuming it has the relevant Coinc colums,
        time1, time2, stat, and timeslide_id
        """
        # If no events, do nothing
        pivot_ifo = self.attrs['pivot']
        fixed_ifo = self.attrs['fixed']
        if len(self.data['%s/time' % pivot_ifo]) == 0 or len(self.data['%s/time' % fixed_ifo]) == 0:
            return self
        from pycbc.events import cluster_coincs
        interval = self.attrs['timeslide_interval']
        cid = cluster_coincs(self.stat,
                             self.data['%s/time' % pivot_ifo],
                             self.data['%s/time' % fixed_ifo],
                             self.timeslide_id,
                             interval,
                             window)
        return self.select(cid)


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
                group, = self.h5file.keys()
            else:
                raise RuntimeError("Didn't get a group!")
        self.group_key = group
        self.group = self.h5file[group]
        self.columns = columnlist if columnlist is not None \
                       else list(self.group.keys())
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
    def __init__(self, trig_file, bank_file, veto_file,
                 segment_name, filter_func, detector, premask=None):
        logging.info('Loading triggers')
        self.trigs_f = HFile(trig_file, 'r')
        self.trigs = self.trigs_f[detector]
        self.ifo = detector  # convenience attributes
        self.detector = detector
        if bank_file:
            logging.info('Loading bank')
            self.bank = HFile(bank_file, 'r')
        else:
            logging.info('No bank file given')
            # empty dict in place of non-existent hdf file
            self.bank = {}

        if premask is not None:
            self.mask = premask
        else:
            self.mask = np.ones(len(self.trigs['end_time']), dtype=bool)

        if veto_file:
            logging.info('Applying veto segments')
            # veto_mask is an array of indices into the trigger arrays
            # giving the surviving triggers
            logging.info('%i triggers before vetoes', self.mask.sum())
            self.veto_mask, _ = events.veto.indices_outside_segments(
                self.end_time, [veto_file],
                ifo=detector, segment_name=segment_name)

            idx = np.flatnonzero(self.mask)[self.veto_mask]
            self.mask[:] = False
            self.mask[idx] = True
            logging.info('%i triggers remain after vetoes',
                          len(self.veto_mask))

        # FIXME this should use the hfile select interface to avoid
        # memory and processing limitations.
        if filter_func:
            # get required columns into the namespace with dummy attribute
            # names to avoid confusion with other class properties
            logging.info('Setting up filter function')
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
            for c in chain(self.trigs.keys(), self.bank.keys()):
                if c in filter_func: delattr(self, '_'+c)

            self.mask = self.mask & self.filter_mask
            logging.info('%i triggers remain after cut on %s',
                         sum(self.mask), filter_func)

    def __getitem__(self, key):
        # Is key in the TRIGGER_MERGE file?
        try:
            return self.get_column(key)
        except KeyError:
            pass

        # Is key in the bank file?
        try:
            self.checkbank(key)
            return self.bank[key][:][self.template_id]
        except (RuntimeError, KeyError) as exc:
            err_msg = "Cannot find {} in input files".format(key)
            raise ValueError(err_msg) from exc

    def checkbank(self, param):
        if self.bank == {}:
            return RuntimeError("Can't get %s values without a bank file"
                                                                       % param)

    def trig_dict(self):
        """Returns dict of the masked trigger valuse """
        mtrigs = {}
        for k in self.trigs:
            if len(self.trigs[k]) == len(self.trigs['end_time']):
                if self.mask is not None:
                    mtrigs[k] = self.trigs[k][self.mask]
                else:
                    mtrigs[k] = self.trigs[k][:]
        return mtrigs

    @classmethod
    def get_param_names(cls):
        """Returns a list of plottable CBC parameter variables"""
        return [m[0] for m in inspect.getmembers(cls) \
            if type(m[1]) == property]

    def apply_mask(self, logic_mask):
        """Apply a boolean array to the set of triggers"""
        if hasattr(self.mask, 'dtype') and (self.mask.dtype == 'bool'):
            orig_indices = self.mask.nonzero()[0][logic_mask]
            self.mask[:] = False
            self.mask[orig_indices] = True
        else:
            self.mask = list(np.array(self.mask)[logic_mask])

    def mask_to_n_loudest_clustered_events(self, rank_method,
                                           n_loudest=10,
                                           cluster_window=10):
        """Edits the mask property of the class to point to the N loudest
        single detector events as ranked by ranking statistic. Events are
        clustered so that no more than 1 event within +/- cluster-window will
        be considered."""

        # If this becomes memory intensive we can optimize
        stat = rank_method.rank_stat_single((self.ifo, self.trig_dict()))
        if len(stat) == 0:
            # No triggers, so just return here
            self.stat = np.array([])
            return

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
        index.sort()
        self.stat = stat[index]
        if hasattr(self.mask, 'dtype') and self.mask.dtype == 'bool':
            orig_indices = np.flatnonzero(self.mask)[index]
            self.mask = list(orig_indices)
        elif isinstance(self.mask, list):
            self.mask = list(np.array(self.mask)[index])

    @property
    def template_id(self):
        return self.get_column('template_id').astype(int)

    @property
    def mass1(self):
        self.checkbank('mass1')
        return self.bank['mass1'][:][self.template_id]

    @property
    def mass2(self):
        self.checkbank('mass2')
        return self.bank['mass2'][:][self.template_id]

    @property
    def spin1z(self):
        self.checkbank('spin1z')
        return self.bank['spin1z'][:][self.template_id]

    @property
    def spin2z(self):
        self.checkbank('spin2z')
        return self.bank['spin2z'][:][self.template_id]

    @property
    def spin2x(self):
        self.checkbank('spin2x')
        return self.bank['spin2x'][:][self.template_id]

    @property
    def spin2y(self):
        self.checkbank('spin2y')
        return self.bank['spin2y'][:][self.template_id]

    @property
    def spin1x(self):
        self.checkbank('spin1x')
        return self.bank['spin1x'][:][self.template_id]

    @property
    def spin1y(self):
        self.checkbank('spin1y')
        return self.bank['spin1y'][:][self.template_id]

    @property
    def inclination(self):
        self.checkbank('inclination')
        return self.bank['inclination'][:][self.template_id]

    @property
    def f_lower(self):
        self.checkbank('f_lower')
        return self.bank['f_lower'][:][self.template_id]

    @property
    def approximant(self):
        self.checkbank('approximant')
        return self.bank['approximant'][:][self.template_id]

    @property
    def mtotal(self):
        return self.mass1 + self.mass2

    @property
    def mchirp(self):
        return conversions.mchirp_from_mass1_mass2(self.mass1, self.mass2)

    @property
    def eta(self):
        return conversions.eta_from_mass1_mass2(self.mass1, self.mass2)

    @property
    def effective_spin(self):
        # FIXME assumes aligned spins
        return conversions.chi_eff(self.mass1, self.mass2,
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
        return self.get_column('end_time')

    @property
    def template_duration(self):
        return self.get_column('template_duration')

    @property
    def snr(self):
        return self.get_column('snr')

    @property
    def sgchisq(self):
        return self.get_column('sg_chisq')

    @property
    def u_vals(self):
        return self.get_column('u_vals')

    @property
    def rchisq(self):
        return self.get_column('chisq') \
            / (self.get_column('chisq_dof') * 2 - 2)

    @property
    def psd_var_val(self):
        return self.get_column('psd_var_val')

    @property
    def newsnr(self):
        return ranking.newsnr(self.snr, self.rchisq)

    @property
    def newsnr_sgveto(self):
        return ranking.newsnr_sgveto(self.snr, self.rchisq, self.sgchisq)

    @property
    def newsnr_sgveto_psdvar(self):
        return ranking.newsnr_sgveto_psdvar(self.snr, self.rchisq,
                                           self.sgchisq, self.psd_var_val)

    @property
    def newsnr_sgveto_psdvar_threshold(self):
        return ranking.newsnr_sgveto_psdvar_threshold(self.snr, self.rchisq,
                                           self.sgchisq, self.psd_var_val)

    def get_ranking(self, rank_name, **kwargs):
        return ranking.get_sngls_ranking_from_trigs(self, rank_name, **kwargs)

    def get_column(self, cname):
        # Fiducial value that seems to work, not extensively tuned.
        MFRAC = 0.3

        # If the mask accesses few enough elements then directly use it
        # This can be slower than reading in all the elements if most of them
        # will be read.
        if self.mask is not None and (isinstance(self.mask, list) or \
                (len(self.mask.nonzero()[0]) < (len(self.mask) * MFRAC))):
            return self.trigs[cname][self.mask]

        # We have a lot of elements to read so we resort to readin the entire
        # array before masking.
        elif self.mask is not None:
            return self.trigs[cname][:][self.mask]
        else:
            return self.trigs[cname][:]


class ForegroundTriggers(object):

    # Injection files are expected to only have 'exclusive' IFAR/FAP values,
    # should use has_inc=False for these.
    def __init__(self, coinc_file, bank_file, sngl_files=None, n_loudest=None,
                 group='foreground', has_inc=True):
        self.coinc_file = FileData(coinc_file, group=group)
        if 'ifos' in self.coinc_file.h5file.attrs:
            self.ifos = self.coinc_file.h5file.attrs['ifos'].split(' ')
        else:
            raise ValueError("File doesn't have an 'ifos' attribute!",
                             coinc_file)
        self.sngl_files = {}
        if sngl_files is not None:
            for sngl_file in sngl_files:
                curr_dat = FileData(sngl_file)
                curr_ifo = curr_dat.group_key
                self.sngl_files[curr_ifo] = curr_dat

        if not all([ifo in self.sngl_files.keys() for ifo in self.ifos]):
            print("sngl_files: {}".format(sngl_files))
            print("self.ifos: {}".format(self.ifos))
            raise RuntimeError("IFOs in statmap file not all represented "
                               "by single-detector trigger files.")
        if not sorted(self.sngl_files.keys()) == sorted(self.ifos):
            logging.warning("WARNING: Single-detector trigger files "
                            "given for IFOs not in the statmap file")

        self.bank_file = HFile(bank_file, "r")
        self.n_loudest = n_loudest

        self._inclusive = has_inc
        self._sort_arr = None
        self._template_id = None
        self._trig_ids = None
        self.get_active_segments()

    @property
    def sort_arr(self):
        if self._sort_arr is None:
            if self._inclusive:
                try:
                    ifar = self.coinc_file.get_column('ifar')
                except KeyError:
                    logging.warning("WARNING: Can't find inclusive IFAR!"
                                    "Using exclusive IFAR instead ...")
                    ifar = self.coinc_file.get_column('ifar_exc')
                    self._inclusive = False
            else:
                ifar = self.coinc_file.get_column('ifar_exc')
            sorting = ifar.argsort()[::-1]
            if self.n_loudest:
                sorting = sorting[:self.n_loudest]
            self._sort_arr = sorting
        return self._sort_arr

    @property
    def template_id(self):
        if self._template_id is None:
            template_id = self.get_coincfile_array('template_id')
            self._template_id = template_id.astype(int)
        return self._template_id

    @property
    def trig_id(self):
        if self._trig_ids is not None:
            return self._trig_ids

        self._trig_ids = {}
        for ifo in self.ifos:
            self._trig_ids[ifo] = self.get_coincfile_array(ifo + '/trigger_id')
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
        for ifo in self.ifos:
            try:
                tid = self.trig_id[ifo]
                lgc = tid == -1
                # Put in *some* value for the invalid points to avoid failure
                # Make sure this doesn't change the cached internal array!
                tid = np.copy(tid)
                tid[lgc] = 0
                # If small number of points don't read the full file
                if len(tid) < 1000:
                    curr = []
                    hdf_dataset = self.sngl_files[ifo].group[variable]
                    for idx in tid:
                        curr.append(hdf_dataset[idx])
                    curr = np.array(curr)
                else:
                    curr = self.sngl_files[ifo].get_column(variable)[tid]
            except IndexError:
                if len(self.trig_id[ifo]) == 0:
                    curr = np.array([])
                    lgc = curr == 0
                else:
                    raise
            return_dict[ifo] = (curr, np.logical_not(lgc))
        return return_dict

    def get_active_segments(self):
        self.active_segments = {}
        for ifo in self.ifos:
            starts = self.sngl_files[ifo].get_column('search/start_time')
            ends = self.sngl_files[ifo].get_column('search/end_time')
            self.active_segments[ifo] = veto.start_end_to_segments(starts,
                                                                   ends)

    def get_end_time(self):
        times_gen = (self.get_coincfile_array('{}/time'.format(ifo))
                     for ifo in self.ifos)
        ref_times = np.array([mean_if_greater_than_zero(t)[0]
                              for t in zip(*times_gen)])
        return ref_times

    def get_ifos(self):
        """
        Returns
        -------
        ifos_list
             List of lists of ifo names involved in each foreground event.
             Ifos will be listed in the same order as self.ifos
        """
        # Ian thinks this could be coded more simply and efficiently
        # Note also that effectively the same thing is done as part of the
        # to_coinc_hdf_object method
        ifo_or_minus = []
        for ifo in self.ifos:
            ifo_trigs = np.where(self.get_coincfile_array(ifo + '/time') < 0,
                                 '-', ifo)
            ifo_or_minus.append(ifo_trigs)
        ifos_list = [list(trig[trig != '-'])
                     for trig in iter(np.array(ifo_or_minus).T)]
        return ifos_list

    def to_coinc_xml_object(self, file_name):
        outdoc = ligolw.Document()
        outdoc.appendChild(ligolw.LIGO_LW())

        ifos = list(self.sngl_files.keys())
        proc_id = ligolw_process.register_to_xmldoc(outdoc, 'pycbc',
                     {}, instruments=ifos, comment='', version=pycbc_version.git_hash,
                     cvs_repository='pycbc/'+pycbc_version.git_branch,
                     cvs_entry_time=pycbc_version.date).process_id

        search_summ_table = lsctables.New(lsctables.SearchSummaryTable)
        coinc_h5file = self.coinc_file.h5file
        try:
            start_time = coinc_h5file['segments']['coinc']['start'][:].min()
            end_time = coinc_h5file['segments']['coinc']['end'][:].max()
        except KeyError:
            start_times = []
            end_times = []
            for ifo_comb in coinc_h5file['segments']:
                if ifo_comb == 'foreground_veto':
                    continue
                seg_group = coinc_h5file['segments'][ifo_comb]
                start_times.append(seg_group['start'][:].min())
                end_times.append(seg_group['end'][:].max())
            start_time = min(start_times)
            end_time = max(end_times)
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
        coinc_def_row.description = \
            "sngl_inspiral<-->sngl_inspiral coincidences"
        coinc_def_row.coinc_def_id = coinc_def_id
        coinc_def_row.search_coinc_type = 0
        coinc_def_table.append(coinc_def_row)

        bank_col_names = ['mass1', 'mass2', 'spin1z', 'spin2z']
        bank_col_vals = {}
        for name in bank_col_names:
            bank_col_vals[name] = self.get_bankfile_array(name)

        coinc_event_names = ['ifar', 'time', 'fap', 'stat']
        coinc_event_vals = {}
        for name in coinc_event_names:
            if name == 'time':
                coinc_event_vals[name] = self.get_end_time()
            else:
                coinc_event_vals[name] = self.get_coincfile_array(name)

        sngl_col_names = ['snr', 'chisq', 'chisq_dof', 'bank_chisq',
                          'bank_chisq_dof', 'cont_chisq', 'cont_chisq_dof',
                          'end_time', 'template_duration', 'coa_phase',
                          'sigmasq']
        sngl_col_vals = {}
        for name in sngl_col_names:
            sngl_col_vals[name] = self.get_snglfile_array_dict(name)

        sngl_event_count = 0
        for idx in range(len(self.sort_arr)):
            # Set up IDs and mapping values
            coinc_id = lsctables.CoincID(idx)

            # Set up sngls
            sngl_mchirps = []
            sngl_mtots = []
            net_snrsq = 0
            triggered_ifos = []
            for ifo in ifos:
                # If this ifo is not participating in this coincidence then
                # ignore it and move on.
                if not sngl_col_vals['snr'][ifo][1][idx]:
                    continue
                triggered_ifos += [ifo]
                event_id = lsctables.SnglInspiralID(sngl_event_count)
                sngl_event_count += 1
                sngl = return_empty_sngl()
                sngl.event_id = event_id
                sngl.ifo = ifo
                net_snrsq += sngl_col_vals['snr'][ifo][0][idx]**2
                for name in sngl_col_names:
                    val = sngl_col_vals[name][ifo][0][idx]
                    if name == 'end_time':
                        sngl.end = LIGOTimeGPS(val)
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
                # If exact match is not used, get masses from single triggers
                sngl_mchirps += [sngl.mchirp]
                sngl_mtots += [sngl.mtotal]

                sngl_inspiral_table.append(sngl)

                # Set up coinc_map entry
                coinc_map_row = lsctables.CoincMap()
                coinc_map_row.table_name = 'sngl_inspiral'
                coinc_map_row.coinc_event_id = coinc_id
                coinc_map_row.event_id = event_id
                coinc_event_map_table.append(coinc_map_row)

            # Take the mean if exact match is not used
            sngl_combined_mchirp = np.mean(sngl_mchirps)
            sngl_combined_mtot = np.mean(sngl_mtots)

            # Set up coinc inspiral and coinc event tables
            coinc_event_row = lsctables.Coinc()
            coinc_inspiral_row = lsctables.CoincInspiral()
            coinc_event_row.coinc_def_id = coinc_def_id
            coinc_event_row.nevents = len(triggered_ifos)
            # NB, `coinc_event_row.instruments = triggered_ifos does not give a
            # correct result with ligo.lw 1.7.1
            coinc_event_row.instruments = ','.join(sorted(triggered_ifos))
            coinc_inspiral_row.instruments = triggered_ifos
            coinc_event_row.time_slide_id = time_slide_id
            coinc_event_row.process_id = proc_id
            coinc_event_row.coinc_event_id = coinc_id
            coinc_inspiral_row.coinc_event_id = coinc_id
            coinc_inspiral_row.mchirp = sngl_combined_mchirp
            coinc_inspiral_row.mass = sngl_combined_mtot
            coinc_inspiral_row.end = LIGOTimeGPS(coinc_event_vals['time'][idx])
            coinc_inspiral_row.snr = net_snrsq**0.5
            coinc_inspiral_row.false_alarm_rate = coinc_event_vals['fap'][idx]
            coinc_inspiral_row.combined_far = 1./coinc_event_vals['ifar'][idx]
            # Transform to Hz
            coinc_inspiral_row.combined_far = \
                                    coinc_inspiral_row.combined_far / YRJUL_SI
            coinc_event_row.likelihood = coinc_event_vals['stat'][idx]
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

    def to_coinc_hdf_object(self, file_name):
        ofd = h5py.File(file_name,'w')

        # Some fields are special cases
        logging.info("Outputting search results")
        time = self.get_end_time()
        ofd.create_dataset('time', data=time, dtype=np.float32)

        if self._inclusive:
            ifar = self.get_coincfile_array('ifar')
            ofd.create_dataset('ifar', data=ifar, dtype=np.float32)

        ifar_exc = self.get_coincfile_array('ifar_exc')
        ofd.create_dataset('ifar_exclusive', data=ifar_exc,
                           dtype=np.float32)

        if self._inclusive:
            fap = self.get_coincfile_array('fap')
            ofd.create_dataset('p_value', data=fap,
                               dtype=np.float32)

        fap_exc = self.get_coincfile_array('fap_exc')
        ofd.create_dataset('p_value_exclusive', data=fap_exc,
                           dtype=np.float32)

        # Coinc fields
        for field in ['stat']:
            vals = self.get_coincfile_array(field)
            ofd.create_dataset(field, data=vals, dtype=np.float32)

        logging.info("Outputting template information")
        # Bank fields
        for field in ['mass1','mass2','spin1z','spin2z']:
            vals = self.get_bankfile_array(field)
            ofd.create_dataset(field, data=vals, dtype=np.float32)

        mass1 = self.get_bankfile_array('mass1')
        mass2 = self.get_bankfile_array('mass2')
        mchirp, _ = mass1_mass2_to_mchirp_eta(mass1, mass2)

        ofd.create_dataset('chirp_mass', data=mchirp, dtype=np.float32)

        logging.info("Outputting single-trigger information")
        logging.info("reduced chisquared")
        chisq_vals_valid = self.get_snglfile_array_dict('chisq')
        chisq_dof_vals_valid = self.get_snglfile_array_dict('chisq_dof')
        for ifo in self.ifos:
            chisq_vals = chisq_vals_valid[ifo][0]
            chisq_valid = chisq_vals_valid[ifo][1]
            chisq_dof_vals = chisq_dof_vals_valid[ifo][0]
            rchisq = chisq_vals / (2. * chisq_dof_vals - 2.)
            rchisq[np.logical_not(chisq_valid)] = -1.
            ofd.create_dataset(ifo + '_chisq', data=rchisq,
                               dtype=np.float64)

        # Single-detector fields
        for field in ['sg_chisq', 'end_time', 'sigmasq',
                      'psd_var_val']:
            logging.info(field)
            try:
                vals_valid = self.get_snglfile_array_dict(field)
            except KeyError:
                logging.info(field + " is not present in the "
                             "single-detector files")
            for ifo in self.ifos:
                vals = vals_valid[ifo][0]
                valid = vals_valid[ifo][1]
                vals[np.logical_not(valid)] = -1.
                ofd.create_dataset(ifo + '_'+field, data=vals,
                                   dtype=np.float32)

        snr_vals_valid = self.get_snglfile_array_dict('snr')
        network_snr_sq = np.zeros_like(snr_vals_valid[self.ifos[0]][0])
        for ifo in self.ifos:
            vals = snr_vals_valid[ifo][0]
            valid = snr_vals_valid[ifo][1]
            vals[np.logical_not(valid)] = -1.
            ofd.create_dataset(ifo + '_snr', data=vals,
                               dtype=np.float32)
            network_snr_sq[valid] += vals[valid] ** 2.0
        network_snr = np.sqrt(network_snr_sq)
        ofd.create_dataset('network_snr', data=network_snr, dtype=np.float32)

        logging.info("Triggered detectors")
        # Create a n_ifos by n_events matrix, with the ifo letter if the
        # event contains a trigger from the ifo, empty string if not
        triggered_matrix = [[ifo[0] if v else ''
                             for v in snr_vals_valid[ifo][1]]
                            for ifo in self.ifos]
        # Combine the ifo letters to make a single string per event
        triggered_detectors = [''.join(triggered).encode('ascii')
                               for triggered in zip(*triggered_matrix)]
        ofd.create_dataset('trig', data=triggered_detectors,
                           dtype='<S3')

        logging.info("active detectors")
        # Create a n_ifos by n_events matrix, with the ifo letter if the
        # ifo was active at the event time, empty string if not
        active_matrix = [[ifo[0] if t in self.active_segments[ifo]
                          else '' for t in time]
                         for ifo in self.ifos]
        # Combine the ifo letters to make a single string per event
        active_detectors = [''.join(active_at_time).encode('ascii')
                            for active_at_time in zip(*active_matrix)]
        ofd.create_dataset('obs', data=active_detectors,
                           dtype='<S3')

        ofd.close()


class ReadByTemplate(object):
    # Default assignment to {} is OK for a variable used only in __init__
    def __init__(self, filename, bank=None, segment_name=None, veto_files=None,
                 gating_veto_windows={}):
        self.filename = filename
        self.file = h5py.File(filename, 'r')
        self.ifo = tuple(self.file.keys())[0]
        self.valid = None
        self.bank = h5py.File(bank, 'r') if bank else {}

        # Determine the segments which define the boundaries of valid times
        # to use triggers
        key = '%s/search/' % self.ifo
        s, e = self.file[key + 'start_time'][:], self.file[key + 'end_time'][:]
        self.segs = veto.start_end_to_segments(s, e).coalesce()
        if segment_name is None:
            segment_name = []
        if veto_files is None:
            veto_files = []
        for vfile, name in zip(veto_files, segment_name):
            veto_segs = veto.select_segments_by_definer(vfile, ifo=self.ifo,
                                                        segment_name=name)
            self.segs = (self.segs - veto_segs).coalesce()
        if self.ifo in gating_veto_windows:
            gating_veto = gating_veto_windows[self.ifo].split(',')
            gveto_before = float(gating_veto[0])
            gveto_after = float(gating_veto[1])
            if gveto_before > 0 or gveto_after < 0:
                raise ValueError("Gating veto window values must be negative "
                                 "before gates and positive after gates.")
            if not (gveto_before == 0 and gveto_after == 0):
                autogate_times = np.unique(
                        self.file[self.ifo + '/gating/auto/time'][:])
                if self.ifo + '/gating/file' in self.file:
                    detgate_times = self.file[self.ifo + '/gating/file/time'][:]
                else:
                    detgate_times = []
                gate_times = np.concatenate((autogate_times, detgate_times))
                gating_veto_segs = veto.start_end_to_segments(
                        gate_times + gveto_before,
                        gate_times + gveto_after
                ).coalesce()
                self.segs = (self.segs - gating_veto_segs).coalesce()
        self.valid = veto.segments_to_start_end(self.segs)

    def get_data(self, col, num):
        """Get a column of data for template with id 'num'.

        Parameters
        ----------
        col: str
            Name of column to read
        num: int
            The template id to read triggers for

        Returns
        -------
        data: numpy.ndarray
            The requested column of data
        """
        ref = self.file['%s/%s_template' % (self.ifo, col)][num]
        return self.file['%s/%s' % (self.ifo, col)][ref]

    def set_template(self, num):
        """Set the active template to read from.

        Parameters
        ----------
        num: int
            The template id to read triggers for.

        Returns
        -------
        trigger_id: numpy.ndarray
            The indices of this templates triggers.
        """
        self.template_num = num
        times = self.get_data('end_time', num)

        # Determine which of these template's triggers are kept after
        # applying vetoes
        if self.valid:
            self.keep = veto.indices_within_times(times, self.valid[0],
                                                  self.valid[1])
#            logging.info('applying vetoes')
        else:
            self.keep = np.arange(0, len(times))

        if self.bank != {}:
            self.param = {}
            if 'parameters' in self.bank.attrs:
                for col in self.bank.attrs['parameters']:
                    self.param[col] = self.bank[col][self.template_num]
            else:
                for col in self.bank:
                    self.param[col] = self.bank[col][self.template_num]

        # Calculate the trigger id by adding the relative offset in self.keep
        # to the absolute beginning index of this templates triggers stored
        # in 'template_boundaries'
        trigger_id = self.keep + \
                         self.file['%s/template_boundaries' % self.ifo][num]
        return trigger_id

    def __getitem__(self, col):
        """ Return the column of data for current active template after
        applying vetoes

        Parameters
        ----------
        col: str
            Name of column to read

        Returns
        -------
        data: numpy.ndarray
            The requested column of data
        """
        if self.template_num is None:
            raise ValueError('You must call set_template to first pick the '
                             'template to read data from')
        data = self.get_data(col, self.template_num)
        data = data[self.keep] if self.valid else data
        return data


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
        err_msg = "Do not recognize --chisq-choice %s" % chisq_choice
        raise ValueError(err_msg)
    return chisq

def save_dict_to_hdf5(dic, filename):
    """
    Parameters
    ----------
    dic:
        python dictionary to be converted to hdf5 format
    filename:
        desired name of hdf5 file
    """
    with h5py.File(filename, 'w') as h5file:
        recursively_save_dict_contents_to_group(h5file, '/', dic)

def recursively_save_dict_contents_to_group(h5file, path, dic):
    """
    Parameters
    ----------
    h5file:
        h5py file to be written to
    path:
        path within h5py file to saved dictionary
    dic:
        python dictionary to be converted to hdf5 format
    """
    for key, item in dic.items():
        if isinstance(item, (np.ndarray, np.int64, np.float64, str, int, float,
                             bytes, tuple, list)):
            h5file[path + str(key)] = item
        elif isinstance(item, dict):
            recursively_save_dict_contents_to_group(h5file, path + key + '/', item)
        else:
            raise ValueError('Cannot save %s type' % type(item))

def load_hdf5_to_dict(h5file, path):
    """
    Parameters
    ----------
    h5file:
        h5py file to be loaded as a dictionary
    path:
        path within h5py file to load: '/' for the whole h5py file

    Returns
    -------
    dic:
        dictionary with hdf5 file group content
    """
    dic = {}
    for key, item in h5file[path].items():
        if isinstance(item, h5py.Dataset):
            dic[key] = item[()]
        elif isinstance(item, h5py.Group):
            dic[key] = load_hdf5_to_dict(h5file, path + key + '/')
        else:
            raise ValueError('Cannot load %s type' % type(item))
    return dic

def combine_and_copy(f, files, group):
    """ Combine the same column from multiple files and save to a third"""
    # ensure that the files input is stable for iteration order
    assert isinstance(files, (list, tuple))
    f[group] = np.concatenate([fi[group][:] if group in fi else \
                                   np.array([], dtype=np.uint32) for fi in files])

def name_all_datasets(files):
    assert isinstance(files, (list, tuple))
    datasets = []
    for fi in files:
        datasets += get_all_subkeys(fi, '/')
    return set(datasets)

def get_all_subkeys(grp, key):
    subkey_list = []
    subkey_start = key
    if key == '':
        grpk = grp
    else:
        grpk = grp[key]
    for sk in grpk.keys():
        path = subkey_start + '/' + sk
        if isinstance(grp[path], h5py.Dataset):
            subkey_list.append(path.lstrip('/'))
        else:
            subkey_list += get_all_subkeys(grp, path)
    # returns an empty list if there is no dataset or subgroup within the group
    return subkey_list

#
# =============================================================================
#
#                          Checkpointing utilities
#
# =============================================================================
#


def dump_state(state, fp, path=None, dsetname='state',
               protocol=pickle.HIGHEST_PROTOCOL):
    """Dumps the given state to an hdf5 file handler.

    The state is stored as a raw binary array to ``{path}/{dsetname}`` in the
    given hdf5 file handler. If a dataset with the same name and path is
    already in the file, the dataset will be resized and overwritten with the
    new state data.

    Parameters
    ----------
    state : any picklable object
        The sampler state to dump to file. Can be the object returned by
        any of the samplers' `.state` attribute (a dictionary of dictionaries),
        or any picklable object.
    fp : h5py.File
        An open hdf5 file handler. Must have write capability enabled.
    path : str, optional
        The path (group name) to store the state dataset to. Default (None)
        will result in the array being stored to the top level.
    dsetname : str, optional
        The name of the dataset to store the binary array to. Default is
        ``state``.
    protocol : int, optional
        The protocol version to use for pickling. See the :py:mod:`pickle`
        module for more details.
    """
    memfp = BytesIO()
    pickle.dump(state, memfp, protocol=protocol)
    dump_pickle_to_hdf(memfp, fp, path=path, dsetname=dsetname)


def dump_pickle_to_hdf(memfp, fp, path=None, dsetname='state'):
    """Dumps pickled data to an hdf5 file object.

    Parameters
    ----------
    memfp : file object
        Bytes stream of pickled data.
    fp : h5py.File
        An open hdf5 file handler. Must have write capability enabled.
    path : str, optional
        The path (group name) to store the state dataset to. Default (None)
        will result in the array being stored to the top level.
    dsetname : str, optional
        The name of the dataset to store the binary array to. Default is
        ``state``.
    """
    memfp.seek(0)
    bdata = np.frombuffer(memfp.read(), dtype='S1')
    if path is not None:
        dsetname = path + '/' + dsetname
    if dsetname not in fp:
        fp.create_dataset(dsetname, shape=bdata.shape, maxshape=(None,),
                          dtype=bdata.dtype)
    elif bdata.size != fp[dsetname].shape[0]:
        fp[dsetname].resize((bdata.size,))
    fp[dsetname][:] = bdata


def load_state(fp, path=None, dsetname='state'):
    """Loads a sampler state from the given hdf5 file object.

    The sampler state is expected to be stored as a raw bytes array which can
    be loaded by pickle.

    Parameters
    ----------
    fp : h5py.File
        An open hdf5 file handler.
    path : str, optional
        The path (group name) that the state data is stored to. Default (None)
        is to read from the top level.
    dsetname : str, optional
        The name of the dataset that the state data is stored to. Default is
        ``state``.
    """
    if path is not None:
        fp = fp[path]
    bdata = fp[dsetname][()].tobytes()
    return pickle.load(BytesIO(bdata))


__all__ = ('HFile', 'DictArray', 'StatmapData', 'MultiifoStatmapData',
           'FileData', 'DataFromFiles', 'SingleDetTriggers',
           'ForegroundTriggers', 'ReadByTemplate', 'chisq_choices',
           'get_chisq_from_file_choice', 'save_dict_to_hdf5',
           'recursively_save_dict_contents_to_group', 'load_hdf5_to_dict',
           'combine_and_copy', 'name_all_datasets', 'get_all_subkeys',
           'dump_state', 'dump_pickle_to_hdf', 'load_state')
