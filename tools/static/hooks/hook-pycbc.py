#-----------------------------------------------------------------------------
# Copyright (c) 2013, PyInstaller Development Team.
#
# Distributed under the terms of the GNU General Public License with exception
# for distributing bootloader.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import os
try:
    from PyInstaller.utils.hooks import (collect_data_files, collect_submodules)
except ImportError:
    from PyInstaller.hooks.hookutils import (collect_data_files, collect_submodules)

# Executables that need MKL
needs_mkl = ['pycbc_inspiral','pycbc_single_template']

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
                 'pycbc.fft.backend_cpu',
                 'pycbc.fft.backend_mkl',
                 'pycbc.fft.fftw',
                 'pycbc.fft.mkl',
                 'pycbc.fft.npfft',
                 'pycbc.fft.__init__',
                 'pycbc.events.threshold_cpu',
                 'scipy.linalg.cython_blas',
                 'scipy.linalg.cython_lapack',
                 'scipy.special._ufuncs_cxx',
                 'h5py',
                 'h5py._conv',
                 'h5py._stub',
                 'mpld3'
                 ]

datas = []

# Add html assets to all executables
cwd     = os.environ.get("PYCBC_HOOKS_DIR", os.getcwd())
basedir = cwd.replace('tools/static','')
rootdir = basedir + 'pycbc/results'

for root, subdirs, files in os.walk(rootdir):
    for filename in files:
        if not filename.endswith('.py') and not filename.endswith('.pyc'):
            file_path  = os.path.join(root, filename)
            store_path = '/'.join(file_path.split('/')[:-1])
            store_path = store_path.replace(basedir, '')
            datas.append( (file_path, store_path) )

# Add em-bright data file
file_path = basedir + 'pycbc/tmpltbank/ns_sequences/equil_2H.dat'
store_path = '/'.join(file_path.split('/')[:-1])
store_path = store_path.replace(basedir, '')
datas.append( (file_path, store_path) )

if os.environ.get("NOW_BUILDING", None) in needs_mkl:
    # pull in all the mkl .so files
    datas += find_lib_path('mkl_rt', [])
    datas += find_lib_path('mkl_core', [])
    datas += find_lib_path('mkl_intel_thread', [])
    datas += find_lib_path('mkl_intel_lp64', [])
    datas += find_lib_path('mkl_avx2', [])
    datas += find_lib_path('mkl_def', [])
    datas += find_lib_path('iomp5', [])
    datas += find_lib_path('mkl_mc3', [])

# try to pull in the openmp fftw files
#datas += find_lib_path('fftw3', ['fftw3'])
#datas += find_lib_path('fftw3f', ['fft3f'])
#datas += find_lib_path('fftw3f_omp', ['fftw3f'])
#datas += find_lib_path('fftw3_omp', ['fftw3'])
