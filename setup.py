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

# all supported extension directories (specifying also the order of the build!) 
main_package_dir = os.getcwd()
extension_dirs = [
            main_package_dir + '/pycbc/clayer/cpu',
            main_package_dir + '/pycbc/clayer/opencl',
            main_package_dir + '/pycbc/datavector/clayer/cpu',
            main_package_dir + '/pycbc/datavector/clayer/opencl',
            # todo setup.py main_package_dir + '/pycbc/straindata/clayer/cpu',
            main_package_dir + '/pycbc/matchedfilter/clayer/opencl']


class pycbc_build(build):
    
    def run(self):

        # run setup.py for all extensions
        current_dir = os.getcwd()
        for dir in extension_dirs:
            print 'changing dir to {0}'.format(dir)
            os.chdir(dir)
            spawn.spawn(['python'] + ['setup.py'] + ['build'])
        os.chdir(current_dir)

        # generic build of python modules
        build.run(self)
    
        # todo filecopy the _xyz.so libs to appropriate dirs in build


class pycbc_clean(clean):
    
    def run(self):
        
        # clean everything in build/
        self.all = 1
        clean.run(self)

        # clean all extensions
        current_dir = os.getcwd()
        for dir in extension_dirs:
            print 'changing dir to {0}'.format(dir) 
            os.chdir(dir)
            spawn.spawn(['python'] + ['setup.py'] + ['clean'])
        os.chdir(current_dir)

        # clean *.pyc files
        os.system('find . -name *.pyc -exec rm {} \;')

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
