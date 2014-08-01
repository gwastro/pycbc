#-----------------------------------------------------------------------------
# Copyright (c) 2013, PyInstaller Development Team.
#
# Distributed under the terms of the GNU General Public License with exception
# for distributing bootloader.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import os.path
from PyInstaller.hooks.hookutils import (collect_data_files, collect_submodules)
from distutils.sysconfig import get_python_inc


# IPython (tested with 0.13) requires the following files:
#   ./site-packages/IPython/config/profile/README_STARTUP
datas = collect_data_files('scipy.weave')
datas += collect_data_files('numpy.core')
datas += [(get_python_inc(), 'include')]

hiddenimports = collect_submodules('scipy.weave')
