# Copyright (C) 2017 Duncan Brown
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Generals
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""
This modules contains a function to provide an argparse action that reports
extremely verbose version information for PyCBC, lal, and lalsimulation.
"""

import os, sys
import argparse
import inspect
import subprocess

def print_link(library):
    err_msg = "Could not execute runtime linker to determine\n" + \
              "shared library paths for library:\n  " + library + "\n"
    FNULL = open(os.devnull, 'w')
    try:
        link = subprocess.check_output(['ldd', library],
                                       stderr=FNULL)
    except OSError:
        try:
            link = subprocess.check_output(['otool', '-L', library],
                                           stderr=FNULL)
        except:
            link = err_msg
    except:
        link = err_msg
    return link


class Version(argparse.Action):
    """ print the pycbc, lal and lalsimulation versions """
    def __init__(self, nargs=0, **kw):
        super(Version, self).__init__(nargs=nargs, **kw)


    def __call__(self, parser, namespace, values, option_string=None):

        import pycbc
        version_str="--- PyCBC Version --------------------------\n" + \
            pycbc.version.git_verbose_msg + \
            "\n\nImported from: " + inspect.getfile(pycbc)

        version_str += "\n\n--- LAL Version ----------------------------\n"
        try:
            import lal.git_version
            lal_module = inspect.getfile(lal)
            lal_library = os.path.join( os.path.dirname(lal_module),
                '_lal.so')
            version_str += lal.git_version.verbose_msg + \
            "\n\nImported from: " + lal_module + \
            "\n\nRuntime libraries:\n" + print_link(lal_library)
        except ImportError:
            version_str += "\nLAL not installed in environment\n"

        version_str += "\n\n--- LALSimulation Version-------------------\n"
        try:
            import lalsimulation.git_version
            lalsimulation_module = inspect.getfile(lalsimulation)
            lalsimulation_library = os.path.join( os.path.dirname(lalsimulation_module),
                '_lalsimulation.so')
            version_str += lalsimulation.git_version.verbose_msg + \
            "\n\nImported from: " + lalsimulation_module + \
            "\n\nRuntime libraries:\n" + print_link(lalsimulation_library)
        except ImportError:
            version_str += "\nLALSimulation not installed in environment\n"

        print(version_str)
        sys.exit(0)
