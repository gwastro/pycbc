# convenience classes for accessing hdf5 trigger files
# the 'get_column()' method is implemented parallel to 
# the existing pylal.SnglInspiralUtils functions

import h5py
import numpy as np


class FileData(object):

    def __init__(self, fname, group=None, columnlist=None, filter_func=None):
        '''
        group: names of group to be read from the file
        columns: list of names of columns to be read
        filter_func: Boolean expression using column attributes and math functions
                     eg. self.snr < 6.5
        '''
        self.fname = fname
        if not self.fname: raise RuntimeError("Didn't get a file!")
        if not columnlist: raise RuntimeError("Didn't get any columns!")
        if not group: raise RuntimeError("Didn't get a group!")
        self.filter_func = filter_func
        self.h5file = h5py.File(fname, "r")
        self.group = self.h5file[group]
        # restrict columns to those requested, otherwise use all columns
        self.columns = columnlist if columnlist is not None \
                       else self.group.keys()
        self._mask = None

    def close(self):
        self.h5file.close()

    @property
    def mask(self):
        if self.filter_func is None:
            raise RuntimeError("Can't get a mask without a filter function!")
        else:
            # only evaluate if there is no pre-calculated value
            if self._mask is None:
                # get the required columns into the namespace as numpy arrays
                for column in self.columns:
                    if column in self.filter_func:
                        setattr(self, column, self.group[column][:])
                self._mask = eval(self.filter_func)
            return self._mask

    def get_column(self, col):
        try:
            vals = self.group[col][:]
            if self.filter_func:
                return vals[self.mask]
            else:
                return vals
        except KeyError:  # if the column doesn't exist, as happens for zero triggers
            return np.array([])  # this may or may not have the right effect ..

class DataFromFiles(object):

    def __init__(self, filelist, group=None, columnlist=None, filter_func=None,
                 verbose=False):
        self.files = filelist
        self.group = group
        self.columns = columnlist
        self.filter_func = filter_func
        self.verbose = verbose

    def get_column(self, col):
        if self.verbose: print 'getting %s :' % col
        vals = []
        for f in self.files:
            d = FileData(f, group=self.group, columnlist=self.columns,
                         filter_func=self.filter_func)
            vals.append(d.get_column(col))
            # can't have more than ~1000 h5py.File objects open at once, so close tis
            d.close()
        if self.verbose: print '    got %i values' % sum(len(v) for v in vals)
        return np.concatenate(vals)

