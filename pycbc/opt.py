# Copyright (C) 2015 Joshua Willis
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
This module defines optimization flags and determines hardware features that some
other modules and packages may use in addition to some optimized utilities.
"""
import os, sys
import logging
from collections import OrderedDict
from pycbc.libutils import get_lscpu_caches

logger = logging.getLogger('pycbc.opt')
caches_backend = 'lscpu'
# Work around different Python versions to get runtime
# info on hardware cache sizes
_USE_SUBPROCESS = False
# nix OS?
nix=False

if sys.platform=='darwin' or sys.platform=='linux':
    nix=True

#print(f"opt: OS is nix? {nix}")

if not nix:
    # If Windows, get L2 cache size from env
    # Ignore other vars
    LEVEL2_CACHE_SIZE = int(os.environ["LEVEL2_CACHE_SIZE"])
    logger.info("opt: using LEVEL2_CACHE_SIZE %d from environment",
                LEVEL2_CACHE_SIZE)
else:
    # darwin or linux
    # If python3, subprocess can be used
    try:
        import subprocess
        _USE_SUBPROCESS = True
    except ModuleNotFoundError:
        import commands
        _USE_SUBPROCESS = False

    #print(f"subprocess module found? {_USE_SUBPROCESS}")
    if _USE_SUBPROCESS:
        # Get cache sizes from lscpu, linesize from getconf
        # getconf does not return the correct cache size
        # on modern CPUs
        if sys.platform!='darwin':
            # Flag for later
            LEVEL2_CACHE_LINESIZE=None

            def getconf(confvar):
                """ getconf for cache line sizes """
                return int(subprocess.check_output(['getconf', confvar]))
        
            if caches_backend=='lscpu':
                """ lscpu for caches and their assoc """
                cache_info = get_lscpu_caches()

                def get_lscpu_val(confvar, caches_info=cache_info):
                    """ lscpu overload of getval """
                    return caches_info[confvar]
            
                getval = get_lscpu_val
                
            elif caches_backend=='getconf':
                # getconf overload of getval
                getval = getconf
                caches_info=None
            else:
                raise KeyError(f"Unknown cache backend {caches_backend}")
        
    else:
        def getconf(confvar):
            """ getconf overload, but with older commands module """
            retlist = commands.getstatusoutput('getconf ' + confvar)
            return int(retlist[1])
        
        getval=getconf
    
    if sys.platform!='darwin':
        LEVEL1_DCACHE_SIZE = getval('LEVEL1_DCACHE_SIZE')
        LEVEL2_CACHE_SIZE = getval('LEVEL2_CACHE_SIZE')
        LEVEL3_CACHE_SIZE = getval('LEVEL3_CACHE_SIZE')

        LEVEL1_DCACHE_ASSOC = getval('LEVEL1_DCACHE_ASSOC')
        LEVEL2_CACHE_ASSOC = getval('LEVEL2_CACHE_ASSOC')
        LEVEL3_CACHE_ASSOC = getval('LEVEL3_CACHE_ASSOC')

        # Can use getconf for cache line sizes
        # but it fails to fetch it for L3
        LEVEL1_DCACHE_LINESIZE = getval('LEVEL1_DCACHE_LINESIZE')
        LEVEL2_CACHE_LINESIZE = getval('LEVEL2_CACHE_LINESIZE')
        LEVEL3_CACHE_LINESIZE = getval('LEVEL3_CACHE_LINESIZE')
    else:
        # Get cache linesize from sysctl
        # On Apple M chips, different Lev cache linesizes can be different!
        # Also different cores (P vs E) can have different cache sizes!
        # Cache assocs are are not usually exposed! 
        # So get only sys reported lev2 size here instead
        LEVEL2_CACHE_LINESIZE=int(subprocess.check_output(['sysctl', '-n', 'hw.cachelinesize']))

    # Left here for testing. 
    #print("Cache sizes")
    #print(LEVEL1_DCACHE_SIZE,LEVEL2_CACHE_SIZE, LEVEL3_CACHE_SIZE)
    #print("Cache assoc")
    #print(LEVEL1_DCACHE_ASSOC, LEVEL2_CACHE_ASSOC, LEVEL3_CACHE_ASSOC)
    #print("Cache linesizes")
    #print(LEVEL1_DCACHE_LINESIZE, LEVEL2_CACHE_LINESIZE, LEVEL3_CACHE_LINESIZE)


def insert_optimization_option_group(parser):
    """
    Adds the options used to specify optimization-specific options.

    Parameters
    ----------
    parser : object
        OptionParser instance
    """
    optimization_group = parser.add_argument_group("Options for selecting "
                                   "optimization-specific settings")

    optimization_group.add_argument("--cpu-affinity", help="""
                    A set of CPUs on which to run, specified in a format suitable
                    to pass to taskset.""")
    optimization_group.add_argument("--cpu-affinity-from-env", help="""
                    The name of an enivornment variable containing a set
                    of CPUs on which to run,  specified in a format suitable
                    to pass to taskset.""")


def verify_optimization_options(opt, parser):
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

    # Pin to specified CPUs if requested
    requested_cpus = None

    if opt.cpu_affinity_from_env is not None:
        if opt.cpu_affinity is not None:
            logging.error("Both --cpu_affinity_from_env and --cpu_affinity specified")
            sys.exit(1)

        requested_cpus = os.environ.get(opt.cpu_affinity_from_env)

        if requested_cpus is None:
            logging.error("CPU affinity requested from environment variable %s "
                          "but this variable is not defined" % opt.cpu_affinity_from_env)
            sys.exit(1)

        if requested_cpus == '':
            logging.error("CPU affinity requested from environment variable %s "
                          "but this variable is empty" % opt.cpu_affinity_from_env)
            sys.exit(1)

    if requested_cpus is None:
        requested_cpus = opt.cpu_affinity

    if requested_cpus is not None:
        command = 'taskset -pc %s %d' % (requested_cpus, os.getpid())
        retcode = os.system(command)

        if retcode != 0:
            logging.error('taskset command <%s> failed with return code %d' % \
                          (command, retcode))
            sys.exit(1)

        logging.info("Pinned to CPUs %s " % requested_cpus)

class LimitedSizeDict(OrderedDict):
    """ Fixed sized dict for FIFO caching"""

    def __init__(self, *args, **kwds):
        self.size_limit = kwds.pop("size_limit", None)
        OrderedDict.__init__(self, *args, **kwds)
        self._check_size_limit()

    def __setitem__(self, key, value):
        OrderedDict.__setitem__(self, key, value)
        self._check_size_limit()

    def _check_size_limit(self):
        if self.size_limit is not None:
            while len(self) > self.size_limit:
                self.popitem(last=False)
