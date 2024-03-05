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

import os
import sys
import glob
import argparse
import inspect
import subprocess


def print_link(library):
    err_msg = (
        "Could not execute runtime linker to determine\n"
        f"shared library paths for library:\n  {library}\n"
    )
    try:
        # Linux
        link = subprocess.check_output(
            ['ldd', library],
            stderr=subprocess.DEVNULL,
            text=True
        )
    except OSError:
        try:
            # macOS
            link = subprocess.check_output(
                ['otool', '-L', library],
                stderr=subprocess.DEVNULL,
                text=True
            )
        except:
            link = err_msg
    except:
        link = err_msg
    return link


def get_lal_info(module, lib_glob):
    """Return a string reporting the version and runtime library information
    for a LAL Python import.
    """
    module_path = inspect.getfile(module)
    version_str = (
        module.git_version.verbose_msg +
        "\n\nImported from: " + module_path +
        "\n\nRuntime libraries:\n"
    )
    possible_lib_paths = glob.glob(
        os.path.join(os.path.dirname(module_path), lib_glob)
    )
    for lib_path in possible_lib_paths:
        version_str += print_link(lib_path)
    return version_str


class Version(argparse.Action):
    """Subclass of argparse.Action that prints version information for PyCBC,
    LAL and LALSimulation.
    """
    def __init__(self, nargs=0, **kw):
        super(Version, self).__init__(nargs=nargs, **kw)

    def __call__(self, parser, namespace, values, option_string=None):
        import pycbc

        version_str = (
            "--- PyCBC Version --------------------------\n" +
            pycbc.version.git_verbose_msg +
            "\n\nImported from: " + inspect.getfile(pycbc)
        )

        version_str += "\n\n--- LAL Version ----------------------------\n"
        try:
            import lal.git_version
        except ImportError:
            version_str += "\nLAL not installed in environment\n"
        else:
            version_str += get_lal_info(lal, '_lal*.so')

        version_str += "\n\n--- LALSimulation Version-------------------\n"
        try:
            import lalsimulation.git_version
        except ImportError:
            version_str += "\nLALSimulation not installed in environment\n"
        else:
            version_str += get_lal_info(lalsimulation, '_lalsimulation*.so')

        print(version_str)
        sys.exit(0)


__all__ = ['Version']
