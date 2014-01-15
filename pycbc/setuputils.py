# Copyright (C) 2012 Alex Nitz
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
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

""" Helper functions for configuring PyCBC
"""
import os, sys, commands

def pkg_config(pkg_libraries):
    """Use pkg-config to query for the location of libraries, library directories,
       and header directories

       Arguments:
           pkg_libries(list): A list of packages as strings

       Returns:
           libraries(list), library_dirs(list), include_dirs(list)
    """
    libraries=[]
    library_dirs=[] 
    include_dirs=[]

    # Check that we have the packages
    for pkg in pkg_libraries:
        if os.system('pkg-config --exists %s 2>/dev/null' % pkg) == 0:
            pass
        else:
            print "Could not find library {0}".format(pkg)
            sys.exit(1)

    # Get the pck-config flags
    if len(pkg_libraries)>0 :
        # PKG_CONFIG_ALLOW_SYSTEM_CFLAGS explicitly lists system paths.
        # On system-wide LAL installs, this is needed for swig to find lalswig.i
        for token in commands.getoutput("PKG_CONFIG_ALLOW_SYSTEM_CFLAGS=1 pkg-config --libs --cflags %s" % ' '.join(pkg_libraries)).split():
            if token.startswith("-l"):
                libraries.append(token[2:])
            elif token.startswith("-L"):
                library_dirs.append(token[2:])
            elif token.startswith("-I"):
                include_dirs.append(token[2:])

    return libraries, library_dirs, include_dirs

def pkg_config_header_strings(pkg_libraries):
    """ Returns a list of header strings that could be passed to a compiler
    """
    libs, lib_dirs, header_dirs = pkg_config(pkg_libraries)
  
    header_strings = []

    for header_dir in header_dirs:
        header_strings.append("-I" + header_dir)

    return header_strings        

