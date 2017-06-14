#!/usr/bin/env python
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

from __future__ import print_function

import os, fnmatch, sys, subprocess, shutil

from trace import fullmodname

try:
    from setuptools.command.install import install as _install
    from setuptools.command.install_egg_info import install_egg_info as egg_info
    USE_SETUPTOOLS = True
except:
    from distutils.command.install import install as _install
    USE_SETUPTOOLS = False

from distutils.errors import DistutilsError
from distutils.core import setup, Command, Extension
from distutils.command.clean import clean as _clean
from distutils.file_util import write_file
from distutils.version import LooseVersion

try:
    import numpy.version
    if LooseVersion(numpy.version.version) < LooseVersion("1.6.4"):
        print(" Numpy >= 1.6.4 is required for pycbc dependencies. \n"
              " We found version %s already installed. Please update \n"
              " to a more recent version and then retry PyCBC  \n"
              " installation. \n"
              " \n"
              " Using pip: [pip install 'numpy>=1.6.4' --upgrade --user] \n"
              "" % numpy.version.version)
        exit(1)
except ImportError:
    pass

requires = ['lal.lal', 'lalsimulation.lalsimulation', 'glue']
setup_requires = []
install_requires =  setup_requires + ['Mako>=1.0.1',
                      'argparse>=1.3.0',
                      'decorator>=3.4.2',
                      'scipy>=0.13.0',
                      'weave>=0.16.0',
                      'unittest2',
                      'matplotlib>=1.3.1',
                      'numpy>=1.6.4',
                      'pillow',
                      'h5py>=2.5',
                      'jinja2',
                      'mpld3>=0.3',
                      'pyRXP>=2.1.0',
                      'pycbc-glue-obsolete==1.1.0',
                      'kombine==0.8.1',
                      'emcee>=2.2.0',
                      'corner>=2.0.1',
                      ]

#FIXME Remove me when we bump to h5py > 2.5
try:
    import h5py
except ImportError:
    setup_requires.append('cython')
else:
    import h5py.version
    if h5py.version.version < '2.5':
        setup_requires.append('cython')


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
                print('removed {0}'.format(f))
            except:
                pass

        for fol in self.clean_folders:
            shutil.rmtree(fol, ignore_errors=True)
            print('removed {0}'.format(fol))

class install(_install):
    def run(self):
        etcdirectory = os.path.join(self.install_data, 'etc')
        if not os.path.exists(etcdirectory):
            os.makedirs(etcdirectory)

        filename = os.path.join(etcdirectory, 'pycbc-user-env.sh')
        self.execute(write_file,
                     (filename, [self.extra_dirs]),
                     "creating %s" % filename)

        env_file = open(filename, 'w')
        print("# Source this file to access PyCBC", file=env_file)
        print("PATH=" + self.install_scripts + ":$PATH", file=env_file)
        print("PYTHONPATH=" + self.install_libbase + ":$PYTHONPATH",
              file=env_file)
        print("export PYTHONPATH", file=env_file)
        print("export PATH", file=env_file)
        env_file.close()

        _install.run(self)

def do_setup(*args):
    return True

_install._called_from_setup=do_setup

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
            test_command = [sys.executable,
                            'test/' + test + '.py',
                            '-s', self.scheme]
            a = subprocess.call(test_command, env=os.environ)
            if a != 0:
                result_str = str(test).ljust(30) + ": Fail : " + str(a)
            else:
                result_str = str(test).ljust(30) + ": Pass"
            test_results.append(result_str)

        for test in test_results:
            print(test)

class test(Command):
    def has_cuda(self):
        import pycbc
        return pycbc.HAVE_CUDA

    sub_commands = [('test_cpu',None),('test_cuda',has_cuda)]
    user_options = []
    description = "run the available tests for all compute schemes (cpu, cuda)"
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

# write versioning info
def get_version_info():
    """Get VCS info and write version info to version.py
    """
    from pycbc import _version_helper

    # If this is a pycbc git repo always populate versoin information using GIT
    try:
        vcs_info = _version_helper.generate_git_version_info()

        with open('pycbc/version.py', 'w') as f:
            f.write("# coding: utf-8\n")
            f.write("# Generated by setup.py for PyCBC on %s.\n\n"
                    % vcs_info.build_date)

            # print general info
            f.write('version = \'%s\'\n' % vcs_info.version)
            f.write('date = \'%s\'\n' % vcs_info.date)
            f.write('release = %s\n' % vcs_info.release)
            f.write('last_release = \'%s\'\n' % vcs_info.last_release)

            # print git info
            f.write('\ngit_hash = \'%s\'\n' % vcs_info.hash)
            f.write('git_branch = \'%s\'\n' % vcs_info.branch)
            f.write('git_tag = \'%s\'\n' % vcs_info.tag)
            f.write('git_author = \'%s\'\n' % vcs_info.author)
            f.write('git_committer = \'%s\'\n' % vcs_info.committer)
            f.write('git_status = \'%s\'\n' % vcs_info.status)
            f.write('git_builder = \'%s\'\n' % vcs_info.builder)
            f.write('git_build_date = \'%s\'\n' % vcs_info.build_date)
            f.write('git_verbose_msg = """Branch: %s\n'
                    'Tag: %s\n'
                    'Id: %s\n'
                    'Builder: %s\n'
                    'Build date: %s\n'
                    'Repository status is %s"""\n' %(vcs_info.branch,
                                                   vcs_info.tag,
                                                   vcs_info.hash,
                                                   vcs_info.builder,
                                                   vcs_info.build_date,
                                                   vcs_info.status))
            f.write('from pycbc._version import *\n')
            version = vcs_info.version

    # If this is a release or another kind of source distribution of PyCBC
    except:
        version = '1.7.2'
        release = 'True'
        date = hash = branch = tag = author = committer = status = builder = build_date = ''

        with open('pycbc/version.py', 'w') as f:
            f.write("# Generated by setup.py for PyCBC.\n\n")

            # print general infov
            f.write('version = \'%s\'\n' % version)
            f.write('date = \'%s\'\n' % date)
            f.write('release = %s\n' % release)

            # print git info
            f.write('\ngit_hash = \'%s\'\n' % hash)
            f.write('git_branch = \'%s\'\n' % branch)
            f.write('git_tag = \'%s\'\n' % tag)
            f.write('git_author = \'%s\'\n' % author)
            f.write('git_committer = \'%s\'\n' % committer)
            f.write('git_status = \'%s\'\n' % status)
            f.write('git_builder = \'%s\'\n' % builder)
            f.write('git_build_date = \'%s\'\n' % build_date)
            f.write('git_verbose_msg = """Version: %s Release: %s \n'
                    ' """\n' % (version, release))
            f.write('from pycbc._version import *\n')

    from pycbc import version
    version = version.version
    return version

class build_docs(Command):
    user_options = []
    description = "Build the documentation pages"
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        subprocess.check_call("cd docs; cp Makefile.std Makefile; cp conf_std.py conf.py; sphinx-apidoc "
                              " -o ./ -f -A 'PyCBC dev team' -V '0.1' ../pycbc && make html",
                            stderr=subprocess.STDOUT, shell=True)

class build_gh_pages(Command):
    user_options = []
    description = "Build the documentation pages for GitHub"
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        subprocess.check_call("mkdir -p _gh-pages/latest && touch _gh-pages/.nojekyll && "
                              "cd docs; cp Makefile.gh_pages Makefile; cp conf_std.py conf.py; sphinx-apidoc "
                              " -o ./ -f -A 'PyCBC dev team' -V '0.1' ../pycbc && make html",
                            stderr=subprocess.STDOUT, shell=True)

cmdclass = { 'test'  : test,
             'build_docs' : build_docs,
             'build_gh_pages' : build_gh_pages,
             'install' : install,
             'test_cpu':test_cpu,
             'test_cuda':test_cuda,
             'clean' : clean,
            }

extras_require = {'cuda': ['pycuda>=2015.1', 'scikit-cuda']}

# do the actual work of building the package
VERSION = get_version_info()

setup (
    name = 'PyCBC',
    version = VERSION,
    description = 'Analyze gravitational-wave data, find signals, and study their parameters.',
    author = 'Ligo Virgo Collaboration - PyCBC team',
    author_email = 'alex.nitz@ligo.org',
    url = 'https://ligo-cbc.github.io',
    download_url = 'https://github.com/ligo-cbc/pycbc/tarball/v%s' % VERSION,
    keywords = ['ligo', 'physics', 'gravity', 'signal processing', 'gravitational waves'],
    cmdclass = cmdclass,
    setup_requires = setup_requires,
    extras_require = extras_require,
    install_requires = install_requires,
    scripts  = [
               'bin/minifollowups/pycbc_injection_minifollowup',
               'bin/minifollowups/pycbc_foreground_minifollowup',
               'bin/minifollowups/pycbc_sngl_minifollowup',
               'bin/minifollowups/pycbc_single_template_plot',
               'bin/minifollowups/pycbc_plot_chigram',
               'bin/minifollowups/pycbc_page_coincinfo',
               'bin/minifollowups/pycbc_page_injinfo',
               'bin/minifollowups/pycbc_page_snglinfo',
               'bin/minifollowups/pycbc_plot_trigger_timeseries',
               'bin/pycbc_banksim',
               'bin/pycbc_banksim_combine_banks',
               'bin/pycbc_banksim_match_combine',
               'bin/pycbc_faithsim',
               'bin/pycbc_inspiral',
               'bin/pycbc_inspiral_skymax',
               'bin/pycbc_live',
               'bin/pycbc_live_nagios_monitor',
               'bin/pycbc_single_template',
               'bin/pycbc_multi_inspiral',
               'bin/pycbc_make_banksim',
               'bin/pycbc_splitbank',
               'bin/pycbc_hdf5_splitbank',
               'bin/pycbc_split_inspinj',
               'bin/bank/pycbc_brute_bank',
               'bin/bank/pycbc_geom_aligned_2dstack',
               'bin/bank/pycbc_geom_aligned_bank',
               'bin/bank/pycbc_geom_nonspinbank',
               'bin/bank/pycbc_aligned_bank_cat',
               'bin/bank/pycbc_aligned_stoch_bank',
               'bin/bank/pycbc_coinc_bank2hdf',
               'bin/bank/pycbc_tmpltbank_to_chi_params',
               'bin/bank/pycbc_bank_verification',
               'bin/pycbc_make_faithsim',
               'bin/pycbc_get_ffinal',
               'bin/pycbc_inj_cut',
               'bin/pycbc_upload_xml_to_gracedb',
               'bin/pycbc_dark_vs_bright_injections',
               'bin/pycbc_make_html_page',
               'bin/pycbc_optimal_snr',
               'bin/pycbc_fit_sngl_trigs',
               'bin/pycbc_randomize_inj_dist_by_optsnr',
               'bin/hdfcoinc/pycbc_calculate_psd',
               'bin/hdfcoinc/pycbc_average_psd',
               'bin/hdfcoinc/pycbc_coinc_mergetrigs',
               'bin/hdfcoinc/pycbc_coinc_findtrigs',
               'bin/hdfcoinc/pycbc_coinc_statmap',
               'bin/hdfcoinc/pycbc_coinc_statmap_inj',
               'bin/hdfcoinc/pycbc_page_foreground',
               'bin/hdfcoinc/pycbc_page_foundmissed',
               'bin/hdfcoinc/pycbc_page_ifar',
               'bin/hdfcoinc/pycbc_page_snrifar',
               'bin/hdfcoinc/pycbc_page_snrratehist',
               'bin/hdfcoinc/pycbc_page_sensitivity',
               'bin/hdfcoinc/pycbc_page_banktriggerrate',
               'bin/hdfcoinc/pycbc_coinc_hdfinjfind',
               'bin/hdfcoinc/pycbc_page_snrchi',
               'bin/hdfcoinc/pycbc_page_segments',
               'bin/hdfcoinc/pycbc_page_segtable',
               'bin/hdfcoinc/pycbc_page_segplot',
               'bin/hdfcoinc/pycbc_page_vetotable',
               'bin/hdfcoinc/pycbc_plot_psd_file',
               'bin/hdfcoinc/pycbc_plot_psd_timefreq',
               'bin/hdfcoinc/pycbc_plot_range',
               'bin/hdfcoinc/pycbc_foreground_censor',
               'bin/hdfcoinc/pycbc_plot_hist',
               'bin/hdfcoinc/pycbc_page_recovery',
               'bin/hdfcoinc/pycbc_page_injtable',
               'bin/hdfcoinc/pycbc_strip_injections',
               'bin/hdfcoinc/pycbc_page_coinc_snrchi',
               'bin/hdfcoinc/pycbc_distribute_background_bins',
               'bin/hdfcoinc/pycbc_combine_statmap',
               'bin/hdfcoinc/pycbc_stat_dtphase',
               'bin/hdfcoinc/pycbc_plot_singles_vs_params',
               'bin/hdfcoinc/pycbc_plot_singles_timefreq',
               'bin/hdfcoinc/pycbc_plot_throughput',
               'bin/hdfcoinc/pycbc_plot_background_coincs',
               'bin/hdfcoinc/pycbc_plot_bank_bins',
               'bin/hdfcoinc/pycbc_merge_psds',
               'bin/hdfcoinc/pycbc_plot_gating',
               'bin/hdfcoinc/pycbc_fit_sngls_by_template',
               'bin/hdfcoinc/pycbc_fit_sngls_over_param',
               'bin/hdfcoinc/pycbc_fit_sngls_binned',
               'bin/hwinj/pycbc_generate_hwinj',
               'bin/hwinj/pycbc_generate_hwinj_from_xml',
               'bin/hwinj/pycbc_plot_hwinj',
               'bin/hwinj/pycbc_insert_frame_hwinj',
               'bin/pycbc_submit_dax',
               'bin/mvsc/pycbc_mvsc_get_features',
               'bin/pycbc_coinc_time',
               'bin/pygrb/pycbc_make_offline_grb_workflow',
               'bin/pygrb/pycbc_make_grb_summary_page',
               'bin/pycbc_condition_strain',
               'bin/workflows/pycbc_make_inference_workflow',
               'bin/inference/pycbc_inference',
               'bin/inference/pycbc_inference_plot_acceptance_rate',
               'bin/inference/pycbc_inference_plot_acf',
               'bin/inference/pycbc_inference_plot_acl',
               'bin/inference/pycbc_inference_plot_geweke',
               'bin/inference/pycbc_inference_plot_gelman_rubin',
               'bin/inference/pycbc_inference_plot_movie',
               'bin/inference/pycbc_inference_plot_posterior',
               'bin/inference/pycbc_inference_plot_prior',
               'bin/inference/pycbc_inference_plot_samples',
               'bin/inference/pycbc_inference_table_summary',
               'bin/plotting/pycbc_plot_waveform',
               'bin/plotting/pycbc_banksim_plot_eff_fitting_factor',
               'bin/plotting/pycbc_banksim_table_point_injs',
               'bin/plotting/pycbc_banksim_plot_fitting_factors',
               'bin/workflows/pycbc_create_sbank_workflow',
               'bin/workflows/pycbc_create_uberbank_workflow',
               'bin/workflows/pycbc_make_coinc_search_workflow',
               'bin/workflows/pycbc_make_psd_estimation_workflow',
               'bin/workflows/pycbc_create_bank_verifier_workflow',
               'bin/pycbc_compress_bank',
               'bin/pycbc_ringinj',
               'tools/einsteinathome/pycbc_build_eah.sh'
               ],
    packages = [
               'pycbc',
               'pycbc.distributions',
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
               'pycbc.io',
               'pycbc.inference',
               'pycbc.inject',
               'pycbc.frame',
               ],
     package_data = {'pycbc.workflow': find_package_data('pycbc/workflow'),
	             'pycbc.results': find_package_data('pycbc/results'),
                     'pycbc.tmpltbank': find_package_data('pycbc/tmpltbank')},
)

