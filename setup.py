#!/usr/bin/python
# Copyright (C) 2012 Alex Nitz, Andrew Miller, Josh Willis
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
setup.py file for PyCBC package
"""
import os
import fnmatch
import sys
import subprocess
import shutil
from trace import fullmodname
from distutils.core import setup, Command, Extension
from distutils.command.clean import clean as _clean
from distutils.command.install import install as _install
from distutils.command.build import build as _build
from numpy import get_include as np_get_include
from pycbc.setuputils import pkg_config
from distutils.file_util import write_file

# Now use the above function to set up our extension library's needs:
ext_libraries, ext_library_dirs, ext_include_dirs = pkg_config(["lal", "lalsimulation", "lalinspiral"])

requires = ['lal.lal', 'lalinspiral.lalinspiral', 'lalsimulation.lalsimulation']
requires +=  ['numpy', 'scipy', 'glue', 'argparse']

def find_package_data(dirname):
    def find_paths(dirname):
        items = []
        for fname in os.listdir(dirname):
            path = os.path.join(dirname, fname)
            if os.path.isdir(path):
                items += find_paths(path)
            elif not path.endswith(".py") and not path.endswith(".pyc"):
                items.append(path)
        return items
    items = find_paths(dirname)
    return [os.path.relpath(path, dirname) for path in items]

def find_package_data(dirname):
    def find_paths(dirname):
        items = []
        for fname in os.listdir(dirname):
            path = os.path.join(dirname, fname)
            if os.path.isdir(path):
                items += find_paths(path)
            elif not path.endswith(".py") and not path.endswith(".pyc"):
                items.append(path)
        return items
    items = find_paths(dirname)
    return [os.path.relpath(path, dirname) for path in items]

# Add swig-generated files to the list of things to clean, so they
# get regenerated each time.
class clean(_clean):
    def finalize_options (self):
        _clean.finalize_options(self)
        self.clean_files = []
        self.clean_folders = ['docs/_build']
    def run(self):
        _clean.run(self)
        for f in self.clean_files:
            try:
                os.unlink(f)
                print 'removed {0}'.format(f)
            except:
                pass

        for fol in self.clean_folders:
            shutil.rmtree(fol, ignore_errors=True)
            print 'removed {0}'.format(fol)

class install(_install):
    def run(self):
        # Check for some of the required python packages
        for ppkg in requires:
            try:
                __import__(ppkg)
            except:
                raise RuntimeError('Failed to locate required package: %s  ' % ppkg)

        etcdirectory = os.path.join(self.install_data, 'etc')
        if not os.path.exists(etcdirectory):
            os.makedirs(etcdirectory)

        filename = os.path.join(etcdirectory, 'pycbc-user-env.sh')
        self.execute(write_file,
                     (filename, [self.extra_dirs]),
                     "creating %s" % filename)

        env_file = open(filename, 'w')
        print >> env_file, "# Source this file to access PyCBC"
        print >> env_file, "PATH=" + self.install_scripts + ":$PATH"
        print >> env_file, "PYTHONPATH=" + self.install_libbase + ":$PYTHONPATH"
        print >> env_file, "export PYTHONPATH"
        print >> env_file, "export PATH"
        env_file.close()
        _install.run(self)


# Override build order, so swig is handled first.
class build(_build):
    # override order of build sub-commands to work around
    # <http://bugs.python.org/issue7562>
    sub_commands = [ ('build_py', _build.has_pure_modules),
                     ('build_scripts', _build.has_scripts) ]

    def run(self):
        _build.run(self)

test_results = []
# Run all of the testing scripts
class TestBase(Command):
    user_options = []
    test_modules = []
    def initialize_options(self):
        self.scheme = None
        self.build_dir = None
    def finalize_options(self):
        #Populate the needed variables
        self.set_undefined_options('build',('build_lib', 'build_dir'))

    def find_test_modules(self,pattern):
       # Find all the unittests that match a given string pattern
        modules= []
        for path, dirs, files in os.walk("test"):
            for filename in fnmatch.filter(files, pattern):
                #add the test directories to the path
                sys.path.append(os.path.join(path))
                #save the module name for importing
                modules.append(fullmodname(filename))
        return modules

    def run(self):
        self.run_command('build')
        # Get the list of cpu test modules
        self.test_modules = self.find_test_modules("test*.py")

        # Run from the build directory
        if 'PYTHONPATH' in os.environ:
            os.environ['PYTHONPATH'] = self.build_dir + ":" + os.environ['PYTHONPATH']
        else:
            os.environ['PYTHONPATH'] = self.build_dir

        test_results.append("\n" + (self.scheme + " tests ").rjust(30))
        for test in self.test_modules:
            test_command = 'python ' + 'test/' + test + '.py -s ' + self.scheme
            a = subprocess.call(test_command,env=os.environ,shell=True)
            if a != 0:
                result_str = str(test).ljust(30) + ": Fail : " + str(a)
            else:
                result_str = str(test).ljust(30) + ": Pass"
            test_results.append(result_str)

        for test in test_results:
            print test

class test(Command):
    def has_cuda(self):
        import pycbc
        return pycbc.HAVE_CUDA
    def has_opencl(self):
        import pycbc
        return pycbc.HAVE_OPENCL

    sub_commands = [('test_cpu',None),('test_cuda',has_cuda),('test_opencl',has_opencl)]
    user_options = []
    description = "run the available tests for all compute schemes (cpu,cuda,opencl)"
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        for cmd_name in self.get_sub_commands():
            self.run_command(cmd_name)

class test_cpu(TestBase):
    description = "run all CPU tests"
    def initialize_options(self):
        TestBase.initialize_options(self)
        self.scheme = 'cpu'

class test_cuda(TestBase):
    description = "run CUDA tests"
    def initialize_options(self):
        TestBase.initialize_options(self)
        self.scheme = 'cuda'

class test_opencl(TestBase):
    description = "run OpenCL tests"
    def initialize_options(self):
        TestBase.initialize_options(self)
        self.scheme = 'opencl'

# write versioning info
def generate_version_info():
    """Get VCS info and write version info to version.py
    """
    from pycbc import _version_helper

    vcs_info = _version_helper.generate_git_version_info()

    with open('pycbc/version.py', 'w') as f:
        f.write("# Generated by setup.py for PyCBC on %s.\n\n"
                % vcs_info.build_date)

        # print general info
        f.write('version = \'%s\'\n' % vcs_info.version)
        f.write('date = \'%s\'\n' % vcs_info.date)
        f.write('release = %s\n' % vcs_info.release)

        # print git info
        f.write('\ngit_hash = \'%s\'\n' % vcs_info.hash)
        f.write('git_branch = \'%s\'\n' % vcs_info.branch)
        f.write('git_tag = \'%s\'\n' % vcs_info.tag)
        f.write('git_author = \'%s\'\n' % vcs_info.author)
        f.write('git_committer = \'%s\'\n' % vcs_info.committer)
        f.write('git_status = \'%s\'\n' % vcs_info.status)
        f.write('git_builder = \'%s\'\n' % vcs_info.builder)
        f.write('git_build_date = \'%s\'\n' % vcs_info.build_date)
        f.write('git_verbose_msg = """Branch: %s\nTag: %s\nId: %s\nBuilder: %s\nBuild date: %s\nRepository status is %s"""' %(vcs_info.branch,vcs_info.tag,vcs_info.hash,vcs_info.builder,vcs_info.build_date,vcs_info.status) )
    return vcs_info.version

class build_docs(Command):
    user_options = []
    description = "Build the documentation pages"
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        subprocess.check_call("cd docs; cp conf_std.py conf.py; sphinx-apidoc -o ./ -f -A 'PyCBC dev team' -V '0.1' ../pycbc && make html",
                            stderr=subprocess.STDOUT, shell=True)

class build_docs_test(Command):
    user_options = []
    description = "Build the documentation pages in testing mode"
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        subprocess.check_call("cd docs; cp conf_test.py conf.py; sphinx-apidoc -o ./ -f -A 'PyCBC dev team' -V '0.1' ../pycbc && make html",
                            stderr=subprocess.STDOUT, shell=True)                            
                           
# do the actual work of building the package
VERSION = generate_version_info()

setup (
    name = 'PyCBC',
    version = VERSION,
    description = 'Gravitational wave CBC analysis toolkit',
    author = 'Ligo Virgo Collaboration - PyCBC team',
    url = 'https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/',
    cmdclass = { 'test'  : test,
                 'build_docs' : build_docs,
                 'build_docs_test' : build_docs_test,
                 'install' : install,
                 'test_cpu':test_cpu,
                 'test_cuda':test_cuda,
                 'test_opencl':test_opencl,
                 'clean' : clean,
                 'build' : build},
    requires = requires,
    scripts  = [
               'bin/lalapps/lalapps_inspiral_ahope',
               'bin/lalapps/lalapps_tmpltbank_ahope',
               'bin/pycbc_banksim',
               'bin/pycbc_faithsim',
               'bin/pycbc_inspiral',
               'bin/pycbc_multi_inspiral',
               'bin/pycbc_make_banksim',
               'bin/pycbc_splitbank',
               'bin/pycbc_split_inspinj',
               'bin/pycbc_geom_aligned_2dstack',
               'bin/pycbc_geom_aligned_bank',
               'bin/pycbc_geom_nonspinbank',
               'bin/pycbc_aligned_bank_cat',
               'bin/pycbc_aligned_stoch_bank',
               'bin/pycbc_make_faithsim',
               'bin/pycbc_get_ffinal',
               'bin/pycbc_timeslides',
               'bin/pycbc_sqlite_simplify',
               'bin/pycbc_calculate_far',
               'bin/pycbc_compute_durations',
               'bin/pycbc_pipedown_plots',
               'bin/pycbc_tmpltbank_to_chi_params',
               'bin/pycbc_bank_verification',
               'bin/pycbc_run_sqlite',
               'bin/pycbc_inspinjfind',
               'bin/pycbc_write_results_page',
               'bin/gstlal/pycbc_calculate_likelihood',
               'bin/gstlal/pycbc_combine_likelihood',
               'bin/gstlal/pycbc_compute_far_from_snr_chisq_histograms',
               'bin/gstlal/pycbc_gen_ranking_data',
               'bin/gstlal/pycbc_pickle_horizon_distances',
               'bin/pycbc_make_coinc_workflow',
               'bin/pycbc_make_daily_workflow',
               'bin/pycbc_ligolw_find_playground',
               'bin/hdfcoinc/pycbc_make_hdf_coinc_workflow',
               'bin/pycbc_basic_pegasus_plan',
               'bin/hdfcoinc/pycbc_coinc_mergetrigs',
               'bin/hdfcoinc/pycbc_coinc_findtrigs',
               'bin/hdfcoinc/pycbc_coinc_bank2hdf',
               'bin/hdfcoinc/pycbc_coinc_trig2hdf',
               'bin/hdfcoinc/pycbc_coinc_statmap',
               'bin/hdfcoinc/pycbc_page_foreground',
               'bin/hdfcoinc/pycbc_page_foundmissed',
               'bin/hdfcoinc/pycbc_page_snrifar',
               'bin/hdfcoinc/pycbc_page_sensitivity',
               'bin/hdfcoinc/pycbc_page_banktriggerrate',
               'bin/hdfcoinc/pycbc_coinc_hdfinjfind',
               'bin/hdfcoinc/pycbc_page_snrchi',
               'bin/hdfcoinc/pycbc_results_server',
               'bin/hdfcoinc/pycbc_page_inspiralrange',
               'bin/hdfcoinc/pycbc_page_segments',
               'bin/pycbc_submit_dax',
               'bin/hdfcoinc/pycbc_page_coinc_snrchi',
               ],
    packages = [
               'pycbc',
               'pycbc.fft',
               'pycbc.types',
               'pycbc.filter',
               'pycbc.psd',
               'pycbc.waveform',
               'pycbc.events',
               'pycbc.noise',
               'pycbc.vetoes',
               'pycbc.tmpltbank',
               'pycbc.workflow',
               'pycbc.results',
               ],
     package_data = {'pycbc.workflow': find_package_data('pycbc/workflow'), 
	             'pycbc.results': find_package_data('pycbc/results')},
)


