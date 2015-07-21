# Copyright (C) 2015 Larne Pekowsky
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
import os, sys
import logging
import shutil, atexit, signal


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


def _clear_weave_cache():
    """Deletes the weave cache specified in os.environ['PYTHONCOMPILED']"""

    cache_dir = os.environ['PYTHONCOMPILED']
    if os.path.exists(cache_dir):
        shutil.rmtree(cache_dir)
    logging.info("Cleared weave cache %s" % cache_dir)


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

    # Check whether to use a private directory for scipy.weave
    if opt.per_process_weave_cache:
        cache_dir = os.path.join(cache_dir, str(os.getpid()))
        os.environ['PYTHONCOMPILED'] = cache_dir
        logging.info("Setting weave cache to %s" % cache_dir)

    if not os.path.exists(cache_dir):
        try:
            os.makedirs(cache_dir)
        except:
            logging.error("Unable to create weave cache %s" % cache_dir)
            sys.exit(1)
        
    if opt.clear_weave_cache_at_start:
        _clear_weave_cache()
        os.makedirs(cache_dir)

    if opt.clear_weave_cache_at_end:
        atexit.register(_clear_weave_cache)
        signal.signal(signal.SIGTERM, _clear_weave_cache)



