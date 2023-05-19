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

import importlib, inspect
import os, fnmatch, ctypes, sys, subprocess
from ctypes.util import find_library
from collections import deque
from subprocess import getoutput


# Be careful setting the mode for opening libraries! Some libraries (e.g.
# libgomp) seem to require the DEFAULT_MODE is used. Others (e.g. FFTW when
# MKL is also present) require that os.RTLD_DEEPBIND is used. If seeing
# segfaults around this code, play around with this!
DEFAULT_RTLD_MODE = ctypes.DEFAULT_MODE


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
            print("Could not find library {0}".format(pkg))
            sys.exit(1)

    # Get the pck-config flags
    if len(pkg_libraries)>0 :
        # PKG_CONFIG_ALLOW_SYSTEM_CFLAGS explicitly lists system paths.
        # On system-wide LAL installs, this is needed for swig to find lalswig.i
        for token in getoutput("PKG_CONFIG_ALLOW_SYSTEM_CFLAGS=1 pkg-config --libs --cflags %s" % ' '.join(pkg_libraries)).split():
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
    _, _, header_dirs = pkg_config(pkg_libraries)

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

    # don't try calling pkg-config if NO_PKGCONFIG is set in environment
    if os.environ.get("NO_PKGCONFIG", None):
        return []

    # if calling pkg-config failes, don't continue and don't try again.
    with open(os.devnull, "w") as FNULL:
        try:
            subprocess.check_call(["pkg-config", "--version"], stdout=FNULL)
        except:
            print(
                "PyCBC.libutils: pkg-config call failed, "
                "setting NO_PKGCONFIG=1",
                file=sys.stderr,
            )
            os.environ['NO_PKGCONFIG'] = "1"
            return []

    # First, check that we can call pkg-config on each package in the list
    for pkg in packages:
        if not pkg_config_check_exists(pkg):
            raise ValueError("Package {0} cannot be found on the pkg-config search path".format(pkg))

    libdirs = []
    for token in getoutput("PKG_CONFIG_ALLOW_SYSTEM_LIBS=1 pkg-config --libs-only-L {0}".format(' '.join(packages))).split():
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
                if fnmatch.fnmatch(libfile,'lib'+libname+'.so*') or \
                        fnmatch.fnmatch(libfile,'lib'+libname+'.dylib*') or \
                        fnmatch.fnmatch(libfile,'lib'+libname+'.*.dylib*') or \
                        fnmatch.fnmatch(libfile,libname+'.dll') or \
                        fnmatch.fnmatch(libfile,'cyg'+libname+'-*.dll'):
                    possible.append(libfile)
        except OSError:
            pass
        # There might be more than one library found, we want the highest-numbered
        if (len(possible) > 0):
            possible.sort()
            return os.path.join(nextdir,possible[-1])
    # If we get here, we didn't find it...
    return None

def get_ctypes_fullpath(libname, packages):
    libdirs = []
    # First try to get from LD_LIBRARY_PATH
    if "LD_LIBRARY_PATH" in os.environ:
        libdirs += os.environ["LD_LIBRARY_PATH"].split(":")
    # Next try to append via pkg_config
    try:
        libdirs += pkg_config_libdirs(packages)
    except ValueError:
        pass
    # We might be using conda/pip/virtualenv or some combination. This can
    # leave lib files in a directory that LD_LIBRARY_PATH or pkg_config
    # can miss.
    libdirs.append(os.path.join(sys.prefix, "lib"))

    # Note that the function below can accept an empty list for libdirs, in
    # which case it will return None
    fullpath = get_libpath_from_dirlist(libname, libdirs)

    if fullpath is None:
        # This won't actually return a full-path, but it should be something
        # that can be found by CDLL
        fullpath = find_library(libname)

    return fullpath

# Next, where possible we setup the capability to use dlmopen (which is present in
# GNU libc) to open libraries within a unique namespace, to avoid symbol collision
# when other packages may import the same library
#
# The following is based off of this github comment by user ihnorton:
#   https://github.com/pytorch/pytorch/issues/31300#issuecomment-567545126
#

_libc = ctypes.CDLL('')
if hasattr(_libc, 'dlmopen'):
    _dlmopen = _libc.dlmopen
    _dlmopen.restype = ctypes.c_void_p
else:
    _dlmopen = None

if hasattr(_libc, 'dlinfo'):
    _dlinfo = _libc.dlinfo
    _dlinfo.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_void_p]
else:
    _dlinfo = None

# From <dlfcn.h> in GNU libc
LM_ID_NEWLM = -1
RTLD_DI_LMID = 1

# Also in <dlfcn.h> is:
#     typedef long int Lmid_t
# so we will use that when calling dlmopen and dlinfo

HAVE_DLMOPEN = ((_dlinfo is not None ) and (_dlmopen is not None))

def get_ctypes_library_and_namespace(libpath, namespace=None):
    """
    This function takes a full library path (such as returned by
    libutils.get_ctypes_fullpath) and attempts to open it in a dedicated
    namespace. This is a GNU extension of the dynamic loader that is not
    available on all platforms.

    Parameters:
        libpath: The fully-qualified path to the library to open. May be
            of type 'str' or 'bytes'
        namespace: A ctypes.c_long value specifying the namespace, or
            'None'. Typically this function should be called the first
            time with 'None' for the namespace, and then the returned
            namespace should be reused on subsequent calls to libraries
            that should be opened in the same namespace

    Output:
        libobj: A ctypes.CDLL opened via 'dlmopen' in either a new, unique
            namespace (if 'namespace' was 'None' when called), or the
            specified namespace otherwise. If a failure occurred but no
            exception was raised, then this will also be 'None', and callers
            should check for that.
        namespace: The unique namespace in which the CDLL was opened. It
            will be new if 'None' was passed as 'namespace' on invocation,
            and the same as 'namespace' otherwise.
    """

    if not HAVE_DLMOPEN:
        return (None, None)

    # We insist on RTLD_DEEPBIND, because some versions of glibc segfault
    # on RTLD_GLOBAL. If you care about RTLD_GLOBAL you don't have a reason
    # to use namespaces. We also use RTLD_NOW for consistency with ctypes.
    flags = os.RTLD_DEEPBIND | os.RTLD_NOW

    if namespace is None:
        _dlmopen_ns = LM_ID_NEWLM
    else:
        _dlmopen_ns = namespace

    if not isinstance(libpath, bytes):
        truepath = libpath.encode('utf-8')
    else:
        trupath = libpath
        
    libhandle = _dlmopen(_dlmopen_ns, truepath, flags)
    if libhandle is None:
        return (None, None)
        
    libobj = ctypes.CDLL(truepath, handle=libhandle)
    if libobj is None:
        return (None, None)
    
    if namespace is None:
        namespace = ctypes.c_long(0)
        retval = _dlinfo(libhandle, RTLD_DI_LMID, ctypes.byref(namespace))
        if retval != 0:
            raise RuntimeError("Could not get value of new namespace for library {0}".format(libpath))

    return (libobj, namespace)

def get_ctypes_library(libname, packages, mode=DEFAULT_RTLD_MODE):
    """
    This function takes a library name, specified in architecture-independent fashion (i.e.
    omitting any prefix such as 'lib' or suffix such as 'so' or 'dylib' or version number) and
    a list of packages that may provide that library, and according first to LD_LIBRARY_PATH,
    then the results of pkg-config, and falling back to the system search path, will try to
    return a CDLL ctypes object.  If 'mode' is given it will be used when loading the library.
    """

    fullpath = get_ctypes_fullpath(libname, packages)
    
    if fullpath is None:
        # We got nothin'
        return None
    else:
        if mode is None:
            return ctypes.CDLL(fullpath)
        else:
            return ctypes.CDLL(fullpath, mode=mode)

def get_ctypes_library_optional_namespace(libname, packages,
                                          mode=DEFAULT_RTLD_MODE,
                                          use_namespace=False,
                                          namespace=None):
    """
    This function is a wrapper around both of 'get_ctypes_library' and
    'get_ctypes_library_namespace'. If 'use_namespace' is 'False', then
    it returns the tuple:
        (get_ctypes_library(libname, packages, mode), 'None')
    and the value of 'namespace' is ignored.

    If 'use_namespace' is True, then it returns the tuple:
        get_ctypes_library_and_namespace(libpath, namespace)
    and the value of 'mode' is ignored (the library is always opened
    with os.RTLD_DEEPBIND | os.RTLD_NOW). In this case, libpath is
    first determined by calling get_ctypes_fullpath(libname, packages).

    Thus, in most cases, this function can be called with 'use_namespace'
    set to 'libutils.HAVE_DLMOPEN', and expect a pair of return values,
    where the second only need be preserved if you expect to reuse a
    namespace. In all cases, if the first value of the return tuple is
    'None', then the call failed.

    In principle, if HAVE_DLMOPEN is 'True' but the call still fails, then
    this function can be called again with 'use_namespace = False', and
    the alternate path using 'get_ctypes_library' will be used on the
    second call.
    """

    if use_namespace:
        libpath = get_ctypes_fullpath(libname, packages)
        if libpath is None:
            return (None, None)
        else:
            return get_ctypes_library_and_namespace(libpath, namespace)
    else:
        libobj = get_ctypes_library(libname, packages, mode)
        return (libobj, None)
        
def import_optional(library_name):
    """ Try to import library but and return stub if not found

    Parameters
    ----------
    library_name: str
        The name of the python library to import

    Returns
    -------
    library: library or stub
        Either returns the library if importing is sucessful or it returns
        a stub which raises an import error and message when accessed.
    """
    try:
        return importlib.import_module(library_name)
    except ImportError:
        # module wasn't found so let's return a stub instead to inform
        # the user what has happened when they try to use related functions
        class no_module(object):
            def __init__(self, library):
                self.library = library

            def __getattribute__(self, attr):
                if attr == 'library':
                    return super().__getattribute__(attr)

                lib = self.library

                curframe = inspect.currentframe()
                calframe = inspect.getouterframes(curframe, 2)
                fun = calframe[1][3]
                msg =""" The function {} tried to access
                         '{}' of library '{}', however,
                        '{}' is not currently installed. To enable this
                        functionality install '{}' (e.g. through pip
                        / conda / system packages / source).
                      """.format(fun, attr, lib, lib, lib)
                raise ImportError(inspect.cleandoc(msg))
        return no_module(library_name)
