# convenience classes for accessing hdf5 trigger files
# the 'get_column()' method is implemented parallel to 
# the existing pylal.SnglInspiralUtils functions

import h5py
import numpy as np
import logging

from lal import LIGOTimeGPS

from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import ilwd
from glue.ligolw import utils as ligolw_utils
from glue.ligolw.utils import process as ligolw_process

from pycbc import version as pycbc_version
from pycbc.tmpltbank import return_search_summary
from pycbc.tmpltbank import return_empty_sngl
from pycbc import pnutils

class FileData(object):

    def __init__(self, fname, group=None, columnlist=None, filter_func=None):
        '''
        Parameters
        ----------
        group : string
            Name of group to be read from the file
        columnlist : list of strings
            Names of columns to be read; if None, use all existing columns 
        filter_func : string 
            String should evaluate to a Boolean expression using attributes
            of the class instance derived from columns: ex. 'self.snr < 6.5'
        '''
        if not fname: raise RuntimeError("Didn't get a file!")
        self.fname = fname
        self.h5file = h5py.File(fname, "r")
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
        '''
        Create a mask implementing the requested filter on the datasets

        Returns
        -------
        array of Boolean
            True for dataset indices to be returned by the get_column method
        '''
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
        '''
        Parameters
        ----------
        col : string
            Name of the dataset to be returned

        Returns
        -------
        numpy array
            Values from the dataset, filtered if requested
        '''
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
        '''
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
        '''
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

class ForegroundTriggers(object):
    # FIXME: A lot of this is hardcoded to expect two ifos
    def __init__(self, coinc_file, bank_file, sngl_files=None, n_loudest=None):
        self.coinc_file = FileData(coinc_file, group='foreground')
        self.sngl_files = {}
        if sngl_files is not None:
            for file in sngl_files:
                curr_dat = FileData(file)
                curr_ifo = curr_dat.group_key
                self.sngl_files[curr_ifo] = curr_dat
        self.bank_file = h5py.File(bank_file, "r")
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
        trigid1 = self.get_coincfile_array('trigger_id1')
        trigid2 = self.get_coincfile_array('trigger_id2')
        for ifo, file in self.sngl_files.items():
            try:
                ifo = file.h5file.keys()[0]
                fs = file.h5file
                if (fs[ifo]['template_id'][:][trigid1] ==\
                                                       self.template_id).all():
                    self._trig_ids[ifo] = trigid1
                else:
                    raise IndexError()
            except IndexError:
                # Trig1 doesn't fit, try trigid2
                try:
                    if (fs[ifo]['template_id'][:][trigid2] ==\
                                                       self.template_id).all():
                        self._trig_ids[ifo] = trigid2
                    else:
                        raise IndexError()
                except IndexError:
                     # Neither fit, fail
                    err_msg = "Cannot match id1 or id2 to the single triggers."
                    raise ValueError(err_msg)
        return self._trig_ids

    def get_coincfile_array(self, variable):
        return self.coinc_file.get_column(variable)[self.sort_arr]

    def get_bankfile_array(self, variable):
        return self.bank_file[variable][:][self.template_id]

    def get_snglfile_array_dict(self, variable):
        return_dict = {}
        for ifo in self.sngl_files.keys():
            curr = self.sngl_files[ifo].get_column(variable)[self.trig_id[ifo]]
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
        coinc_def_row.description = "sngl_inspiral-sngl_inspiral coincidences"
        coinc_def_row.coinc_def_id = coinc_def_id
        coinc_def_row.search_coinc_type = 0
        coinc_def_table.append(coinc_def_row)

        bank_col_names = ['mass1', 'mass2', 'spin1z', 'spin2z']
        bank_col_vals = {}
        for name in bank_col_names:
            bank_col_vals[name] = self.get_bankfile_array(name)

        coinc_event_names = ['ifar','time1','fap', 'stat']
        coinc_event_vals = {}
        for name in coinc_event_names:
            coinc_event_vals[name] = self.get_coincfile_array(name)
        
        sngl_col_names = ['snr', 'chisq', 'chisq_dof', 'bank_chisq',
                          'bank_chisq_dof', 'cont_chisq', 'cont_chisq_dof',
                          'end_time']
        sngl_col_vals = {}
        for name in sngl_col_names:
            sngl_col_vals[name] = self.get_snglfile_array_dict(name)

        for idx, coinc_idx in enumerate(self.sort_arr):
            # Set up IDs and mapping values
            curr_tmplt_id = self.template_id[idx]
            coinc_id = lsctables.CoincID(idx)

            # Set up sngls
            # FIXME: Need to do this to set up mapping from trigid1 to ifo
            self.trig_id            
            # FIXME: As two-ifo is hardcoded loop over all ifos
            sngl_combined_mchirp = 0
            sngl_combined_mtot = 0
            for ifo in ifos:
                sngl_id = self._trig_ids[ifo][idx]
                event_id = lsctables.SnglInspiralID(sngl_id)
                sngl = return_empty_sngl()
                curr_sngl_file = self.sngl_files[ifo].h5file[ifo]
                for name in sngl_col_names:
                    val = sngl_col_vals[name][ifo][idx]
                    setattr(sngl, name, val)
                for name in bank_col_names:
                    val = bank_col_vals[name][idx]
                    setattr(sngl, name, val)
                sngl.mtotal, sngl.eta = pnutils.mass1_mass2_to_mtotal_eta(
                        sngl.mass1, sngl.mass2)
                sngl.mchirp, junk = pnutils.mass1_mass2_to_mchirp_eta(
                        sngl.mass1, sngl.mass2)
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
