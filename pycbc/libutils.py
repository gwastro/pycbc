# Copyright (C) 2014 Josh Willis
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

"""
This module provides a simple interface for loading a shared library via ctypes,
allowing it to be specified in an OS-independent way and searched for preferentially
according to the paths that pkg-config specifies.
"""
import os, fnmatch, ctypes, commands, sys
from ctypes.util import find_library
from collections import deque

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

def pkg_config_check_exists(package):
    return (os.system('pkg-config --exists {0} 2>/dev/null'.format(package)) == 0)

def pkg_config_libdirs(packages):
    """
    Returns a list of all library paths that pkg-config says should be included when
    linking against the list of packages given as 'packages'. An empty return list means
    that the package may be found in the standard system locations, irrespective of
    pkg-config.
    """
    # First, check that we can call pkg-config on each package in the list
    for pkg in packages:
        if not pkg_config_check_exists(pkg):
            raise ValueError("Package {0} cannot be found on the pkg-config search path".format(pkg))

    libdirs = []
    for token in commands.getoutput("PKG_CONFIG_ALLOW_SYSTEM_LIBS=1 pkg-config --libs-only-L {0}".format(' '.join(packages))).split():
        if token.startswith("-L"):
            libdirs.append(token[2:])
    return libdirs

def get_libpath_from_dirlist(libname, dirs):
    """
    This function tries to find the architecture-independent library given by libname in the first
    available directory in the list dirs. 'Architecture-independent' means omitting any prefix such
    as 'lib' or suffix such as 'so' or 'dylib' or version number.  Within the first directory in which
    a matching pattern can be found, the lexicographically first such file is returned, as a string
    giving the full path name.  The only supported OSes at the moment are posix and mac, and this
    function does not attempt to determine which is being run.  So if for some reason your directory
    has both '.so' and '.dylib' libraries, who knows what will happen.  If the library cannot be found,
    None is returned.
    """
    dirqueue = deque(dirs)
    while (len(dirqueue) > 0):
        nextdir = dirqueue.popleft()
        possible = []
        # Our directory might be no good, so try/except
        try:
            for libfile in os.listdir(nextdir):
                if fnmatch.fnmatch(libfile,'lib'+libname+'.so*') or fnmatch.fnmatch(libfile,'lib'+libname+'.dylib*'):
                    possible.append(libfile)
        except OSError:
            pass
        # There might be more than one library found, we want the highest-numbered
        if (len(possible) > 0):
            possible.sort()
            return os.path.join(nextdir,possible[-1])
    # If we get here, we didn't find it...
    return None

def get_ctypes_library(libname, packages, mode=None):
    """
    This function takes a library name, specified in architecture-independent fashion (i.e.
    omitting any prefix such as 'lib' or suffix such as 'so' or 'dylib' or version number) and
    a list of packages that may provide that library, and according first to LD_LIBRARY_PATH,
    then the results of pkg-config, and falling back to the system search path, will try to
    return a CDLL ctypes object.  If 'mode' is given it will be used when loading the library.
    """
    libdirs = []
    # First try to get from LD_LIBRARY_PATH
    if "LD_LIBRARY_PATH" in os.environ:
        libdirs += os.environ["LD_LIBRARY_PATH"].split(":")
    # Next try to append via pkg_config
    try:
        libdirs += pkg_config_libdirs(packages)
    except ValueError:
        pass

    # Note that the function below can accept an empty list for libdirs, in which case
    # it will return None
    fullpath = get_libpath_from_dirlist(libname,libdirs)

    if fullpath is None:
        # This won't actually return a full-path, but it should be something
        # that can be found by CDLL
        fullpath = find_library(libname)

    if fullpath is None:
        # We got nothin'
        return None
    else:
        if mode is None:
            return ctypes.CDLL(fullpath)
        else:
            return ctypes.CDLL(fullpath,mode=mode)
