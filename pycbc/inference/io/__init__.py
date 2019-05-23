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

"""I/O utilities for pycbc inference
"""

from __future__ import absolute_import
from __future__ import print_function

import os
import argparse
import shutil
import textwrap
import numpy
import logging
import h5py as _h5py
from pycbc.io.record import FieldArray, _numpy_function_lib
from pycbc import transforms as _transforms
from pycbc import waveform as _waveform

from pycbc.inference.option_utils import (ParseLabelArg, ParseParametersArg)
from .emcee import EmceeFile
from .emcee_pt import EmceePTFile
from .cpnest import CPNestFile
from .multinest import MultinestFile
from .posterior import PosteriorFile
from .txt import InferenceTXTFile

filetypes = {
    EmceeFile.name: EmceeFile,
    EmceePTFile.name: EmceePTFile,
    CPNestFile.name: CPNestFile,
    MultinestFile.name: MultinestFile,
    PosteriorFile.name: PosteriorFile
}


def get_file_type(filename):
    """ Returns I/O object to use for file.

    Parameters
    ----------
    filename : str
        Name of file.

    Returns
    -------
    file_type : {InferenceFile, InferenceTXTFile}
        The type of inference file object to use.
    """
    txt_extensions = [".txt", ".dat", ".csv"]
    hdf_extensions = [".hdf", ".h5", ".bkup", ".checkpoint"]
    for ext in hdf_extensions:
        if filename.endswith(ext):
            with _h5py.File(filename, 'r') as fp:
                filetype = fp.attrs['filetype']
            return filetypes[filetype]
    for ext in txt_extensions:
        if filename.endswith(ext):
            return InferenceTXTFile
    raise TypeError("Extension is not supported.")


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
            fileclass = get_file_type(path)
        except IOError:
            # file doesn't exist, filetype must be provided
            raise IOError("The file appears not to exist. In this case, "
                          "filetype must be provided.")
    else:
        fileclass = filetypes[filetype]
    return fileclass(path, mode=mode, **kwargs)


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
        # but only do the check if parameters have been written
        if len(parameters) > 0:
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
                _ = fp[group.format(param)][firstidx]
                _ = fp[group.format(param)][lastidx]


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


#
# =============================================================================
#
#                         Command-line Utilities
#
# =============================================================================
#
def get_common_parameters(input_files, collection=None):
    """Gets a list of variable params that are common across all input files.

    If no common parameters are found, a ``ValueError`` is raised.

    Parameters
    ----------
    input_files : list of str
        List of input files to load.
    collection : str, optional
        What group of parameters to load. Can be the name of a list of
        parameters stored in the files' attrs (e.g., "variable_params"), or
        "all". If "all", will load all of the parameters in the files'
        samples group. Default is to load all.

    Returns
    -------
    list :
        List of the parameter names.
    """
    if collection is None:
        collection = "all"
    parameters = []
    for fn in input_files:
        fp = loadfile(fn, 'r')
        if collection == 'all':
            ps = fp[fp.samples_group].keys()
        else:
            ps = fp.attrs[collection]
        parameters.append(set(ps))
        fp.close()
    parameters = list(set.intersection(*parameters))
    if parameters == []:
        raise ValueError("no common parameters found for collection {} in "
                         "files {}".format(collection, ', '.join(input_files)))
    return parameters


class NoInputFileError(Exception):
    """Raised in custom argparse Actions by arguments needing input-files when
    no file(s) were provided."""
    pass


class PrintFileParams(argparse.Action):
    """Argparse action that will load input files and print possible parameters
    to screen. Once this is done, the program is forced to exit immediately.

    The behvior is similar to --help, except that the input-file is read.

    .. note::
        The ``input_file`` attribute must be set in the parser namespace before
        this action is called. Otherwise, a ``NoInputFileError`` is raised.
    """
    def __init__(self, skip_args=None, nargs=0, **kwargs):
        if nargs != 0:
            raise ValueError("nargs for this action must be 0")
        super(PrintFileParams, self).__init__(nargs=nargs, **kwargs)
        self.skip_args = skip_args

    def __call__(self, parser, namespace, values, option_string=None):
        # get the input file(s)
        input_files = namespace.input_file
        if input_files is None:
            # see if we should raise an error
            try:
                raise_err = not parser.no_input_file_err
            except AttributeError:
                raise_err = True
            if raise_err:
                raise NoInputFileError("must provide at least one input file")
            else:
                # just return to stop further processing
                return
        filesbytype = {}
        fileparsers = {}
        for fn in input_files:
            fp = loadfile(fn, 'r')
            try:
                filesbytype[fp.name].append(fn)
            except KeyError:
                filesbytype[fp.name] = [fn]
                # get any extra options
                fileparsers[fp.name], _ = fp.extra_args_parser(
                    skip_args=self.skip_args, add_help=False)
            fp.close()
        # now print information about the intersection of all parameters
        parameters = get_common_parameters(input_files, collection='all')
        print("\n"+textwrap.fill("Parameters available with this (these) "
                                 "input file(s):"), end="\n\n")
        print(textwrap.fill(' '.join(sorted(parameters))),
              end="\n\n")
        # information about the pycbc functions
        pfuncs = sorted(FieldArray.functionlib.fget(FieldArray).keys())
        print(textwrap.fill("Available pycbc functions (see "
                            "http://pycbc.org/pycbc/latest/html for "
                            "more details):"), end="\n\n")
        print(textwrap.fill(', '.join(pfuncs)), end="\n\n")
        # numpy funcs
        npfuncs = sorted([name for (name, obj) in _numpy_function_lib.items()
                          if isinstance(obj, numpy.ufunc)])
        print(textwrap.fill("Available numpy functions:"),
              end="\n\n")
        print(textwrap.fill(', '.join(npfuncs)), end="\n\n")
        # misc
        consts = "e euler_gamma inf nan pi"
        print(textwrap.fill("Recognized constants:"),
              end="\n\n")
        print(consts, end="\n\n")
        print(textwrap.fill("Python arthimetic (+ - * / // ** %), "
                            "binary (&, |, etc.), and comparison (>, <, >=, "
                            "etc.) operators may also be used."), end="\n\n")
        # print out the extra arguments that may be used
        outstr = textwrap.fill("The following are additional command-line "
                               "options that may be provided, along with the "
                               "input files that understand them:")
        print("\n"+outstr, end="\n\n")
        for ftype, fparser in fileparsers.items():
            fnames = ', '.join(filesbytype[ftype])
            if fparser is None:
                outstr = textwrap.fill(
                    "File(s) {} use no additional options.".format(fnames))
                print(outstr, end="\n\n")
            else:
                fparser.usage = fnames
                fparser.print_help()
        parser.exit(0)


class ResultsArgumentParser(argparse.ArgumentParser):
    """Wraps argument parser, and preloads arguments needed for loading samples
    from a file.

    This parser class should be used by any program that wishes to use the
    standard arguments for loading samples. It provides functionality to parse
    file specific options. These file-specific arguments are not included in
    the standard ``--help`` (since they depend on what input files are given),
    but can be seen by running ``--file-help/-H``. The ``--file-help`` will
    also print off information about what parameters may be used given the
    input files.

    As with the standard ``ArgumentParser``, running this class's
    ``parse_args`` will result in an error if arguments are provided that are
    not recognized by the parser, nor by any of the file-specific arguments.
    For example, ``parse_args`` would work on the command
    ``--input-file results.hdf --walker 0`` if
    ``results.hdf`` was created by a sampler that recognizes a ``--walker``
    argument, but would raise an error if ``results.hdf`` was created by a
    sampler that does not recognize a ``--walker`` argument. The extra
    arguments that are recognized are determined by the sampler IO class's
    ``extra_args_parser``.

    Some arguments may be excluded from the parser using the ``skip_args``
    optional parameter.

    Parameters
    ----------
    skip_args : list of str, optional
        Do not add the given arguments to the parser. Arguments should be
        specified as the option string minus the leading '--'; e.g.,
        ``skip_args=['thin-start']`` would cause the ``thin-start`` argument
        to not be included. May also specify sampler-specific arguments. Note
        that ``input-file``, ``file-help``, and ``parameters`` are always
        added.
    defaultparams : {'variable_params', 'all'}, optional
        If no ``--parameters`` provided, which collection of parameters to
        load. If 'all' will load all parameters in the file's
        ``samples_group``. If 'variable_params' or None (the default) will load
        the variable parameters.
    autoparamlabels : bool, optional
        Passed to ``add_results_option_group``; see that function for details.
    \**kwargs :
        All other keyword arguments are passed to ``argparse.ArgumentParser``.
    """
    def __init__(self, skip_args=None, defaultparams=None,
                 autoparamlabels=True, **kwargs):
        super(ResultsArgumentParser, self).__init__(**kwargs)
        # add attribute to communicate to arguments what to do when there is
        # no input files
        self.no_input_file_err = False
        if skip_args is None:
            skip_args = []
        self.skip_args = skip_args
        if defaultparams is None:
            defaultparams = 'variable_params'
        self.defaultparams = defaultparams
        # add the results option grup
        self.add_results_option_group(autoparamlabels=autoparamlabels)

    @property
    def actions(self):
        """Exposes the actions this parser can do as a dictionary.

        The dictionary maps the ``dest`` to actions.
        """
        return {act.dest: act for act in self._actions}

    def _unset_required(self):
        """Convenience function to turn off required arguments for first parse.
        """
        self._required_args = [act for act in self._actions if act.required]
        for act in self._required_args:
            act.required = False

    def _reset_required(self):
        """Convenience function to turn required arguments back on.
        """
        for act in self._required_args:
            act.required = True

    def parse_known_args(self, args=None, namespace=None):
        """Parse args method to handle input-file dependent arguments."""
        # run parse args once to make sure the name space is populated with
        # the input files. We'll turn off raising NoInputFileErrors on this
        # pass
        self.no_input_file_err = True
        self._unset_required()
        opts, extra_opts = super(ResultsArgumentParser, self).parse_known_args(
            args, namespace)
        # now do it again
        self.no_input_file_err = False
        self._reset_required()
        opts, extra_opts = super(ResultsArgumentParser, self).parse_known_args(
            args, opts)
        # populate the parameters option if it wasn't specified
        if opts.parameters is None:
            parameters = get_common_parameters(opts.input_file,
                                               collection=self.defaultparams)
            # now call parse parameters action to populate the namespace
            self.actions['parameters'](self, opts, parameters)
        # parse the sampler-specific options and check for any unknowns
        unknown = []
        for fn in opts.input_file:
            fp = loadfile(fn, 'r')
            sampler_parser, _ = fp.extra_args_parser(skip_args=self.skip_args)
            if sampler_parser is not None:
                opts, still_unknown = sampler_parser.parse_known_args(
                    extra_opts, namespace=opts)
                unknown.append(set(still_unknown))
        # the intersection of the unknowns are options not understood by
        # any of the files
        if len(unknown) > 0:
            unknown = set.intersection(*unknown)
        return opts, list(unknown)

    def add_results_option_group(self, autoparamlabels=True):
        """Adds the options used to call pycbc.inference.io.results_from_cli
        function to the parser.

        These are options releated to loading the results from a run of
        pycbc_inference, for purposes of plotting and/or creating tables.

        Any argument strings included in the ``skip_args`` attribute will not
        be added.

        Parameters
        ----------
        autoparamlabels : bool, optional
            If True, the ``--parameters`` option will use labels from
            ``waveform.parameters`` if a parameter name is the same as a
            parameter there. Otherwise, will just use whatever label is
            provided. Default is True.
        """
        results_reading_group = self.add_argument_group(
            title="Arguments for loading results",
            description="Additional, file-specific arguments may also be "
            "provided, depending on what input-files are given. See "
            "--file-help for details.")
        results_reading_group.add_argument(
            "--input-file", type=str, required=True, nargs="+",
            action=ParseLabelArg, metavar='FILE[:LABEL]',
            help="Path to input HDF file(s). A label may be specified for "
                 "each input file to use for plots when multiple files are "
                 "specified.")
        # advanced help
        results_reading_group.add_argument(
            "-H", "--file-help",
            action=PrintFileParams, skip_args=self.skip_args,
            help="Based on the provided input-file(s), print all available "
                 "parameters that may be retrieved and all possible functions "
                 "on those parameters. Also print available additional "
                 "arguments that may be passed. This option is like an "
                 "advanced --help: if run, the program will just print the "
                 "information to screen, then exit.")
        if autoparamlabels:
            paramparser = ParseParametersArg
            lblhelp = (
                "If LABEL is the same as a parameter in "
                "pycbc.waveform.parameters, the label "
                "property of that parameter will be used (e.g., if LABEL "
                "were 'mchirp' then {} would be used). "
                .format(_waveform.parameters.mchirp.label))
        else:
            paramparser = ParseLabelArg
            lblhelp = ''
        results_reading_group.add_argument(
            "--parameters", type=str, nargs="+", metavar="PARAM[:LABEL]",
            action=paramparser,
            help="Name of parameters to load. If none provided will load all "
                 "of the model params in the input-file. If provided, the "
                 "parameters can be any of the model params or posterior "
                 "stats (loglikelihood, logprior, etc.) in the input file(s), "
                 "derived parameters from them, or any function of them. If "
                 "multiple files are provided, any parameter common to all "
                 "files may be used. Syntax for functions is python; any math "
                 "functions in the numpy libary may be used. Can optionally "
                 "also specify a LABEL for each parameter. If no LABEL is "
                 "provided, PARAM will used as the LABEL. {}"
                 "To see all possible parameters that may be used with the "
                 "given input file(s), as well as all avaiable functions, "
                 "run --file-help, along with one or more input files."
                 .format(lblhelp))
        return results_reading_group


def results_from_cli(opts, load_samples=True, **kwargs):
    """Loads an inference result file along with any labels associated with it
    from the command line options.

    Parameters
    ----------
    opts : ArgumentParser options
        The options from the command line.
    load_samples : bool, optional
        Load the samples from the file.

    Returns
    -------
    fp_all : (list of) BaseInferenceFile type
        The result file as an hdf file. If more than one input file,
        then it returns a list.
    parameters : list of str
        List of the parameters to use, parsed from the parameters option.
    labels : dict
        Dictionary of labels to associate with the parameters.
    samples_all : (list of) FieldArray(s) or None
        If load_samples, the samples as a FieldArray; otherwise, None.
        If more than one input file, then it returns a list.
    \**kwargs :
        Any other keyword arguments that are passed to read samples using
        samples_from_cli
    """

    # lists for files and samples from all input files
    fp_all = []
    samples_all = []

    input_files = opts.input_file
    if isinstance(input_files, str):
        input_files = [input_files]

    # loop over all input files
    for input_file in input_files:
        logging.info("Reading input file %s", input_file)

        # read input file
        fp = loadfile(input_file, "r")

        # load the samples
        if load_samples:
            logging.info("Loading samples")

            # check if need extra parameters for a non-sampling parameter
            file_parameters, ts = _transforms.get_common_cbc_transforms(
                opts.parameters, fp.variable_params)

            # read samples from file
            samples = fp.samples_from_cli(opts, parameters=file_parameters, **kwargs)

            logging.info("Using {} samples".format(samples.size))

            # add parameters not included in file
            samples = _transforms.apply_transforms(samples, ts)

        # else do not read samples
        else:
            samples = None

        # add results to lists from all input files
        if len(input_files) > 1:
            fp_all.append(fp)
            samples_all.append(samples)

        # else only one input file then do not return lists
        else:
            fp_all = fp
            samples_all = samples

    return fp_all, opts.parameters, opts.parameters_labels, samples_all


def injections_from_cli(opts):
    """Gets injection parameters from the inference file(s).

    Parameters
    ----------
    opts : argparser
        Argparser object that has the command-line objects to parse.

    Returns
    -------
    FieldArray
        Array of the injection parameters from all of the input files given
        by ``opts.input_file``.
    """
    input_files = opts.input_file
    if isinstance(input_files, str):
        input_files = [input_files]
    injections = None
    # loop over all input files getting the injection files
    for input_file in input_files:
        fp = loadfile(input_file, 'r')
        these_injs = fp.read_injections()
        if injections is None:
            injections = these_injs
        else:
            injections = injections.append(these_injs)
    return injections
