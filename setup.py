#!/usr/bin/python
# Copyright (C) 2011 Karsten Wiesner
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


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#

"""
setup.py file for PyCBC package
"""

import os

from distutils import log
from distutils.core import setup
from distutils.command import config
from distutils.command import build
from distutils.command.build import build
from distutils.command.clean import clean
from distutils import file_util
from distutils import spawn

ver = 0.1

class pycbc_build(build):
    
    def run(self):

        # run setup.py for all extensions
        
        current_dir = os.getcwd()
        
        os.chdir('/Users/kawies/dev/src/pycbc/pycbc/datavector/clayer/cpu')
        
        spawn.spawn(['python'] + ['setup.py'] + ['build'])
        # os.system('python hello.py')

        os.chdir(current_dir)
        
        build.run(self)
    
        # filecopy the _xyz.so libs to appropriate dirs in build


class pycbc_clean(clean):
    
    def run(self):
        
        self.all = 1       # clean everything in build/
        clean.run(self)
                           # clean specific files
        self.clean_files = []
                            
        for f in self.clean_files:
            print "removing {0}".format(f)
            try:
                os.unlink(f)
            except:
                log.warn("'%s' does not exist -- can't clean it" % f)


# do the actual work of building the package
setup (
    name = 'pycbc',
    version = ver,
    description = 'gravitational wave CBC analysis toolkit',
    author = 'Ligo Virgo Collaboration - pyCBC team',
    author_email = 'https://sugwg-git.phy.syr.edu/dokuwiki/doku.php?id=pycbc:home',
    #ext_modules = [datavectorcpu_ext, matchedfiltercpu_ext],
    #py_modules = ['datavecstim_opencl', 'datavecterm_cpu'],
    # library_dirs = [],
    # libraries = [[ 'pycbc', {        #### what is this?:
    #   'sources' : pycbc_c_sources,
    #   'include_dirs' : ['include'],
    #   'macros' : [] }]],
    # headers = ['include/pycbc.h'],
    package_dir = {'pycbc' : 'pycbc'},
    packages = ['pycbc',
                'pycbc.chisqveto',
                'pycbc.datavector',
                'pycbc.fft',
                'pycbc.matchedfilter',
                'pycbc.overwhiteningfilter',
                'pycbc.singledetectorevent',
                'pycbc.straindata',
                'pycbc.templatebank',], # list package directories
    cmdclass = {
            'build' : pycbc_build,
            'clean' : pycbc_clean 
       }
)
