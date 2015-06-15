# convenience classes for accessing hdf5 trigger files
# the 'get_column()' method is implemented parallel to 
# the existing pylal.SnglInspiralUtils functions

import h5py
import numpy as np
import logging

class FileData(object):

    def __init__(self, fname, group, columnlist=None, filter_func=None):
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
        if not group: raise RuntimeError("Didn't get a group!")
        self.fname = fname
        self.h5file = h5py.File(fname, "r")
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
        try:
            vals = self.group[col]
            if self.filter_func:
                return vals[self.mask]
            else:
                return vals[:]
        except KeyError:  # if the column doesn't exist, as happens for zero triggers
            return np.array([])  # this may or may not have the right effect

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

