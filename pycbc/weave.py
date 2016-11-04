# Copyright (C) 2015 Larne Pekowsky, Alex Nitz
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

"""
This module provides methods for controlling scipy.weave
"""
import os.path, sys
import logging
import shutil, atexit, signal
import fcntl
import pycbc

## Blatently taken from weave to implement a crude file locking scheme
def pycbc_compile_function(code,arg_names,local_dict,global_dict,
                     module_dir,
                     compiler='',
                     verbose=1,
                     support_code=None,
                     headers=None,
                     customize=None,
                     type_converters=None,
                     auto_downcast=1,
                     **kw):
    """ Dummy wrapper around scipy weave compile to implement file locking
    """
    from scipy.weave.inline_tools import _compile_function
    headers = [] if headers is None else headers
    lockfile_dir = pycbc._cache_dir_path
    lockfile_name = os.path.join(lockfile_dir, 'code_lockfile')
    logging.info("attempting to aquire lock '%s' for "
                 "compiling code" % lockfile_name)
    if not os.path.exists(lockfile_dir):
        os.makedirs(lockfile_dir)
    lockfile = open(lockfile_name, 'w')
    fcntl.lockf(lockfile, fcntl.LOCK_EX)
    logging.info("we have aquired the lock")
    
    func = _compile_function(code,arg_names, local_dict, global_dict,
                     module_dir, compiler, verbose,
                     support_code, headers, customize,
                     type_converters,
                     auto_downcast, **kw)

    fcntl.lockf(lockfile, fcntl.LOCK_UN)
    logging.info("the lock has been released")

    return func


def insert_weave_option_group(parser):
    """
    Adds the options used to specify weave options.
    
    Parameters
    ----------
    parser : object
        OptionParser instance
    """
    optimization_group = parser.add_argument_group("Options for controlling "
                                   "scipy.weave")
    
    optimization_group.add_argument("--per-process-weave-cache",
                    action="store_true",
                    default=False,
                    help="""If given, each process will use a separate directory
                         for scipy.weave compilation.  This is slower, but safer if
                         several instances may be starting on the same machine at
                         the same time.""")

    optimization_group.add_argument("--clear-weave-cache-at-start", 
                    action="store_true",
                    default=False,
                    help="If given, delete the contents of the weave cache "
                         "when the process starts")

    optimization_group.add_argument("--clear-weave-cache-at-end", 
                    action="store_true",
                    default=False,
                    help="If given, delete the contents of the weave cache "
                         "when the process exits")

    optimization_group.add_argument("--fixed-weave-cache",
                    action="store_true",
                    default=False,
                    help="If given, use fixed directory PWD/pycbc_inspiral for "
                         " the weave cache")

def _clear_weave_cache():
    """Deletes the weave cache specified in os.environ['PYTHONCOMPILED']"""

    cache_dir = os.environ['PYTHONCOMPILED']
    if os.path.exists(cache_dir):
        shutil.rmtree(cache_dir)
    logging.info("Cleared weave cache %s", cache_dir)


def verify_weave_options(opt, parser):
    """Parses the CLI options, verifies that they are consistent and 
    reasonable, and acts on them if they are

    Parameters
    ----------
    opt : object
        Result of parsing the CLI with OptionParser, or any object with the
        required attributes 
    parser : object
        OptionParser instance.
    """

    # PYTHONCOMPILED is initially set in pycbc.__init__
    cache_dir = os.environ['PYTHONCOMPILED']

    # Check whether to use a fixed directory for scipy.weave
    if opt.fixed_weave_cache:
        if os.environ.get("FIXED_WEAVE_CACHE", None):
            cache_dir = os.environ["FIXED_WEAVE_CACHE"]
        elif getattr(sys, 'frozen', False):
            cache_dir = sys._MEIPASS
        else:
            cache_dir = os.path.join(os.getcwd(),"pycbc_inspiral")
        os.environ['PYTHONCOMPILED'] = cache_dir
        logging.debug("fixed_weave_cache: Setting weave cache to %s", cache_dir)
        sys.path = [cache_dir] + sys.path
        try: os.makedirs(cache_dir)
        except OSError: pass

    # Check whether to use a private directory for scipy.weave
    if opt.per_process_weave_cache:
        cache_dir = os.path.join(cache_dir, str(os.getpid()))
        os.environ['PYTHONCOMPILED'] = cache_dir
        logging.info("Setting weave cache to %s", cache_dir)

    if not os.path.exists(cache_dir):
        try:
            os.makedirs(cache_dir)
        except:
            logging.error("Unable to create weave cache %s", cache_dir)
            sys.exit(1)
        
    if opt.clear_weave_cache_at_start:
        _clear_weave_cache()
        os.makedirs(cache_dir)

    if opt.clear_weave_cache_at_end:
        atexit.register(_clear_weave_cache)
        signal.signal(signal.SIGTERM, _clear_weave_cache)



