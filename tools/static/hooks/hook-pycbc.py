#-----------------------------------------------------------------------------
# Copyright (c) 2013, PyInstaller Development Team.
#
# Distributed under the terms of the GNU General Public License with exception
# for distributing bootloader.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import os
from PyInstaller.hooks.hookutils import (collect_data_files, collect_submodules)

#Some of our libaries are not being picked up automatically (need to invest.)
#In the meantime we can pull them manually, in the same way we normally
#find them
from pycbc.libutils import get_libpath_from_dirlist, pkg_config_libdirs

def find_lib_path(libname, packages):
    libdirs = []
    if "LD_LIBRARY_PATH" in os.environ:
      libdirs += os.environ["LD_LIBRARY_PATH"].split(":")
    try:
        libdirs += pkg_config_libdirs(packages)
    except ValueError:
        pass

    path = get_libpath_from_dirlist(libname, libdirs)
    
    if path is not None:
        return [(path, '')]
    else:
        return []

# pull in the pycbc imports it can't find
hiddenimports = ['pycbc.fft.fft_cpu',
                 'pycbc.filter.matchedfilter_cpu',
                 'pycbc.vetoes.chisq_cpu',
                 'pycbc.waveform.spa_tmplt_cpu',
                 'pycbc.types.array_cpu',
                 'pycbc.fft.fftw',
                 'pycbc.fft.mkl',
                 'pycbc.fft.lalfft',
                 'pycbc.fft.npfft',
                 'pycbc.fft.__init__',
                 'pycbc.events.threshold_cpu',
                 ]
datas = []

# pull in all the mkl .so files
#datas += find_lib_path('mkl_rt', [])
#datas += find_lib_path('mkl_core', [])
#datas += find_lib_path('mkl_intel_thread', [])
#datas += find_lib_path('mkl_intel_lp64', [])
#datas += find_lib_path('mkl_avx2', [])
#datas += find_lib_path('mkl_def', [])

# try to pull in the openmp fftw files
datas += find_lib_path('fftw3', ['fftw3'])
datas += find_lib_path('fftw3f', ['fft3f'])
datas += find_lib_path('fftw3f_omp', ['fftw3f'])
datas += find_lib_path('fftw3_omp', ['fftw3'])
