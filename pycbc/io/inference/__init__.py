# Copyright (C) 2018  Collin Capano
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

"""I/O utilities for GWIn
"""

from __future__ import absolute_import

import os
import shutil
import logging
import h5py as _h5py

from .emcee import EmceeFile
from .txt import InferenceTXTFile

filetypes = {
    EmceeFile.name: EmceeFile,
}


def loadfile(path, mode=None, filetype=None, **kwargs):
    """Loads the given file using the appropriate InferenceFile class.

    If ``filetype`` is not provided, this will try to retreive the ``filetype``
    from the file's ``attrs``. If the file does not exist yet, an IOError will
    be raised if ``filetype`` is not provided.

    Parameters
    ----------
    path : str
        The filename to load.
    mode : str, optional
        What mode to load the file with, e.g., 'w' for write, 'r' for read,
        'a' for append. Default will default to h5py.File's mode, which is 'a'.
    filetype : str, optional
        Force the file to be loaded with the given class name. This must be
        provided if creating a new file.

    Returns
    -------
    filetype instance
        An open file handler to the file. The class used for IO with the file
        is determined by the ``filetype`` keyword (if provided) or the
        ``filetype`` stored in the file (if not provided).
    """
    if filetype is None:
        # try to read the file to get its filetype
        try:
            with _h5py.File(path, 'r') as fp:
                filetype = fp.attrs['filetype']
        except IOError:
            # file doesn't exist, filetype must be provided
            raise IOError("The file appears not to exist. In this case, "
                          "filetype must be provided.")
    return filetypes[filetype](path, mode=mode, **kwargs)

#
# =============================================================================
#
#                         HDF Utilities
#
# =============================================================================
#


def check_integrity(filename):
    """Checks the integrity of an InferenceFile.

    Checks done are:

        * can the file open?
        * do all of the datasets in the samples group have the same shape?
        * can the first and last sample in all of the datasets in the samples
          group be read?

    If any of these checks fail, an IOError is raised.

    Parameters
    ----------
    filename: str
        Name of an InferenceFile to check.

    Raises
    ------
    ValueError
        If the given file does not exist.
    KeyError
        If the samples group does not exist.
    IOError
        If any of the checks fail.
    """
    # check that the file exists
    if not os.path.exists(filename):
        raise ValueError("file {} does not exist".format(filename))
    # if the file is corrupted such that it cannot be opened, the next line
    # will raise an IOError
    with loadfile(filename, 'r') as fp:
        # check that all datasets in samples have the same shape
        parameters = fp[fp.samples_group].keys()
        group = fp.samples_group + '/{}'
        # use the first parameter as a reference shape
        ref_shape = fp[group.format(parameters[0])].shape
        if not all(fp[group.format(param)].shape == ref_shape
                   for param in parameters):
            raise IOError("not all datasets in the samples group have the "
                          "same shape")
        # check that we can read the first/last sample
        firstidx = tuple([0]*len(ref_shape))
        lastidx = tuple([-1]*len(ref_shape))
        for param in parameters:
            fp[group.format(param)][firstidx]
            fp[group.format(param)][lastidx]


def validate_checkpoint_files(checkpoint_file, backup_file):
    """Checks if the given checkpoint and/or backup files are valid.

    The checkpoint file is considered valid if:

        * it passes all tests run by ``check_integrity``;
        * it has at least one sample written to it (indicating at least one
          checkpoint has happened).

    The same applies to the backup file. The backup file must also have the
    same number of samples as the checkpoint file, otherwise, the backup is
    considered invalid.

    If the checkpoint (backup) file is found to be valid, but the backup
    (checkpoint) file is not valid, then the checkpoint (backup) is copied to
    the backup (checkpoint). Thus, this function ensures that checkpoint and
    backup files are either both valid or both invalid.

    Parameters
    ----------
    checkpoint_file : string
        Name of the checkpoint file.
    backup_file : string
        Name of the backup file.

    Returns
    -------
    checkpoint_valid : bool
        Whether or not the checkpoint (and backup) file may be used for loading
        samples.
    """
    # check if checkpoint file exists and is valid
    try:
        check_integrity(checkpoint_file)
        checkpoint_valid = True
    except (ValueError, KeyError, IOError):
        checkpoint_valid = False
    # backup file
    try:
        check_integrity(backup_file)
        backup_valid = True
    except (ValueError, KeyError, IOError):
        backup_valid = False
    # check if there are any samples in the file; if not, we'll just start from
    # scratch
    if checkpoint_valid:
        with loadfile(checkpoint_file, 'r') as fp:
            try:
                group = '{}/{}'.format(fp.samples_group, fp.variable_params[0])
                nsamples = fp[group].size
                checkpoint_valid = nsamples != 0
            except KeyError:
                checkpoint_valid = False
    # check if there are any samples in the backup file
    if backup_valid:
        with loadfile(backup_file, 'r') as fp:
            try:
                group = '{}/{}'.format(fp.samples_group, fp.variable_params[0])
                backup_nsamples = fp[group].size
                backup_valid = backup_nsamples != 0
            except KeyError:
                backup_valid = False
    # check that the checkpoint and backup have the same number of samples;
    # if not, assume the checkpoint has the correct number
    if checkpoint_valid and backup_valid:
        backup_valid = nsamples == backup_nsamples
    # decide what to do based on the files' statuses
    if checkpoint_valid and not backup_valid:
        # copy the checkpoint to the backup
        logging.info("Backup invalid; copying checkpoint file")
        shutil.copy(checkpoint_file, backup_file)
        backup_valid = True
    elif backup_valid and not checkpoint_valid:
        logging.info("Checkpoint invalid; copying backup file")
        # copy the backup to the checkpoint
        shutil.copy(backup_file, checkpoint_file)
        checkpoint_valid = True
    return checkpoint_valid
