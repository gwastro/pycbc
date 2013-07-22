# Copyright (C) 2012--2013  Alex Nitz, Josh Willis
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

#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
"""
This module contains a few helper functions designed to make writing PyCBC
unit tests easier, while still allowing the tests to be run on CPU, CUDA,
and OpenCL.

All tests starting with 'test_' in the test subdirectory of pycbc are run
whenever the command 'python setup.py test' is given.  That command will
attempt to call each test, passing it the argument '-s <scheme>' where
scheme is each of 'cpu', 'cuda', and 'opencl' in turn. Unit tests
designed to validate code that should run under multiple schemes should accept
each of these options, rerunning the same tests under each successive scheme.

This will usually be done by putting something like:

   _scheme, _context = parse_args_all_schemes('MyFeature')

before the definition of the unit test class.  In the definition of 'setUp'
(which is a mandatory function that must be defined in the class, per the
Python unittest module design) one would then usually have:

   self.scheme = _scheme
   self.context = _context

and those properties of any instance can then be used by later tests defined
in the class; for example, beginning a block with 'with self.context:' to
ensure the appropriate context manager is used for that scheme.

Some other unit tests may be for features or sub-packages that are not GPU
capable, and cannot be meaningfully tested on the GPU.  Those tests, after
importing pycbc.test.utils, should instead just call:

   parse_args_cpu_only('MyFeature')

This call is needed because the tests must still be able to accept the arguments
specifying a GPU environment (since setup.py does not know which tests are GPU
capable and which are not) but when called with a GPU scheme will exit immediately.

Both functions take a single string as an argument.  That string is used to customize
the heading of all of the tests (according to feature and scheme) to make the output
of running all of the unit tests somewhat easier to parse when they are all run at
once.
"""

import pycbc
from sys import exit
from optparse import OptionParser, OptionValueError
from pycbc.scheme import CPUScheme, CUDAScheme, OpenCLScheme

def _check_scheme_all(option, opt_str, scheme, parser):
    if scheme=='cuda' and not pycbc.HAVE_CUDA:
        raise optparse.OptionValueError("CUDA not found")

    if scheme=='opencl' and not pycbc.HAVE_OPENCL:
        raise optparse.OptionValueError("OpenCL not found")
    setattr (parser.values, option.dest, scheme)


def parse_args_all_schemes(feature_str):
    _parser = OptionParser()
    _parser.add_option('--scheme','-s', action='callback', type = 'choice',
                       choices = ('cpu','cuda','opencl'),
                       default = 'cpu', dest = 'scheme', callback = _check_scheme_all,
                       help = 'specifies processing scheme, can be cpu [default], cuda, or opencl')
    _parser.add_option('--device-num','-d', action='store', type = 'int',
                       dest = 'devicenum', default=0,
                       help = 'specifies a GPU device to use for CUDA or OpenCL, 0 by default')
    (_opt_list, _args) = _parser.parse_args()

    # Changing the optvalues to a dict makes them easier to read
    _options = vars(_opt_list)

    _scheme = _options['scheme']

    if _scheme == 'cpu':
        _context = CPUScheme()
    if _scheme == 'cuda':
        _context = CUDAScheme(device_num=_options['devicenum'])
    if _scheme == 'opencl':
        _context = OpenCLScheme(device_num=_options['devicenum'])

    _scheme_dict = { 'cpu': 'CPU', 'cuda': 'CUDA', 'opencl' : 'OpenCL'}

    print "Running {0} unit tests for {1}:".format(_scheme_dict[_scheme],feature_str)

    return [_scheme,_context]

del _check_scheme_all

def _check_scheme_cpu(option, opt_str, scheme, parser):
    if scheme=='cuda':
        exit(0)

    if scheme=='opencl':
        exit(0)
    setattr (parser.values, option.dest, scheme)


def parse_args_cpu_only(feature_str):
    _parser = OptionParser()
    _parser.add_option('--scheme','-s', action='callback', type = 'choice',
                       choices = ('cpu','cuda','opencl'),
                       default = 'cpu', dest = 'scheme', callback = _check_scheme_cpu,
                       help = 'specifies processing scheme, can be cpu [default], cuda, or opencl')
    _parser.add_option('--device-num','-d', action='store', type = 'int',
                       dest = 'devicenum', default=0,
                       help = 'specifies a GPU device to use for CUDA or OpenCL, 0 by default')
    (_opt_list, _args) = _parser.parse_args()

    # In this case, the only reason we parsed the arguments was to exit if we were given
    # a GPU scheme.  So if we get here we're on the CPU, and should print out our message
    # and return.

    print "Running {0} unit tests for {1}:".format('CPU',feature_str)

    return

del _check_scheme_cpu

# Clean up our namespace

del exit, CPUScheme, CUDAScheme, OpenCLScheme, pycbc, OptionParser, OptionValueError
