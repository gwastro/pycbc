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
    if numpy.version.version != '1.9.3': 
        print (" Numpy = 1.9.3 is required for pycbc dependencies. \n"
              " We found version %s already installed. Please update \n"
              " to this version and then retry PyCBC installation. \n"
              " \n"
              " Using pip: [pip install 'numpy==1.9.3 --upgrade --user] \n"
              "" % numpy.version.version)
        exit(1)
except ImportError:
    pass
                           
requires = ['lal.lal', 'lalsimulation.lalsimulation', 'glue', 'pylal']
setup_requires = []
install_requires =  setup_requires + [
                      'decorator==4.0.4',
                      'numpy==1.9.3',
                      'numpydoc==0.5',
                      'scipy==0.16.0',
                      'unittest2==1.1.0',
                      'h5py==2.5.0',
                      'alabaster==0.7.6',
                      'argparse==1.3.0',
                      'Babel==2.1.1',
                      'boto==2.5.2',
                      'Cython==0.23.2',
                      'docutils==0.12',
                      'Flask==0.10',
                      'Flask-Cache==0.13.1',
                      'Flask-SQLAlchemy==0.16',
                      'funcsigs==0.4',
                      'itsdangerous==0.21',
                      'Jinja2==2.7',
                      'linecache2==1.0.0',
                      'M2Crypto==0.22.3',
                      'Mako==1.0.2',
                      'MarkupSafe==0.18',
                      'matplotlib==1.4.3',
                      'mock==1.3.0',
                      'mpld3==0.3git',
                      'MySQL-python==1.2.5',
                      'nose==1.3.7',
                      'pam==0.1.4',
                      'pbr==1.8.0',
                      'pegasus-wms==4.5.2',
                      'Pillow==2.9.0',
                      'pip==7.1.2',
                      'psycopg2==2.6',
                      'pyRXP==2.1.0',
                      'pycbc-glue==0.9.6',
                      'pycbc-pylal==0.9.5',
                      'Pygments==2.0.2',
                      'pyOpenSSL==0.13',
                      'pyparsing==2.0.3',
                      'python-cjson==1.1.0',
                      'python-dateutil==2.4.2',
                      'pytz==2015.6',
                      'requests==1.2.3',
                      'setuptools==18.2',
                      'six==1.9.0',
                      'snowballstemmer==1.2.0',
                      'Sphinx==1.3.1',
                      'sphinx-rtd-theme==0.1.9',
                      'sphinxcontrib-programoutput==0.8',
                      'SQLAlchemy==0.8.0',
                      'traceback2==1.4.0',
                      'Werkzeug==0.9.3',
                      'WTForms==1.0.3']
links = ['https://github.com/ligo-cbc/mpld3/tarball/master#egg=mpld3-0.3git']

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
                print 'removed {0}'.format(f)
            except:
                pass

        for fol in self.clean_folders:
            shutil.rmtree(fol, ignore_errors=True)
            print 'removed {0}'.format(fol)

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
        print >> env_file, "# Source this file to access PyCBC"
        print >> env_file, "PATH=" + self.install_scripts + ":$PATH"
        print >> env_file, "PYTHONPATH=" + self.install_libbase + ":$PYTHONPATH"
        print >> env_file, "export PYTHONPATH"
        print >> env_file, "export PATH"
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
            f.write('git_verbose_msg = """Branch: %s\n'
                    'Tag: %s\n'
                    'Id: %s\n'
                    'Builder: %s\n'
                    'Build date: %s\n'
                    'Repository status is %s"""' %(vcs_info.branch, 
                                                   vcs_info.tag,
                                                   vcs_info.hash,
                                                   vcs_info.builder,
                                                   vcs_info.build_date,
                                                   vcs_info.status))
            version = vcs_info.version
            
    # If this is a release or another kind of source distribution of PyCBC
    except:
        version = '1.2.0'
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
                    ' """' % (version, release))
        
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
        subprocess.check_call("cd docs; cp Makefile.gh_pages Makefile; cp conf_std.py conf.py; sphinx-apidoc "
                              " -o ./ -f -A 'PyCBC dev team' -V '0.1' ../pycbc && make html",
                            stderr=subprocess.STDOUT, shell=True)

class build_docs_test(Command):
    user_options = []
    description = "Build the documentation pages in testing mode"
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        subprocess.check_call("cd docs; cp Makefile.std Makefile; cp conf_test.py conf.py; sphinx-apidoc "
                              " -o ./ -f -A 'PyCBC dev team' -V '0.1' ../pycbc && make html",
                            stderr=subprocess.STDOUT, shell=True)                            

cmdclass = { 'test'  : test,
             'build_docs' : build_docs,
             'build_gh_pages' : build_gh_pages,
             'build_docs_test' : build_docs_test,
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
    description = 'Gravitational wave CBC analysis toolkit',
    author = 'Ligo Virgo Collaboration - PyCBC team',
    author_email = 'alex.nitz@ligo.org',
    url = 'https://github.com/ligo-cbc/pycbc',
    download_url = 'https://github.com/ligo-cbc/pycbc/tarball/v%s' % VERSION,
    keywords = ['ligo', 'physics', 'gravity', 'signal processing'],
    cmdclass = cmdclass,
    setup_requires = setup_requires,
    extras_require = extras_require,
    install_requires = install_requires,
    dependency_links = links,
    scripts  = [
               'bin/minifollowups/pycbc_single_template_plot',
               'bin/minifollowups/pycbc_page_coincinfo',
               'bin/minifollowups/pycbc_plot_trigger_timeseries',
               'bin/lalapps/lalapps_inspiral_ahope',
               'bin/lalapps/lalapps_tmpltbank_ahope',
               'bin/pycbc_banksim',
               'bin/pycbc_faithsim',
               'bin/pycbc_inspiral',
               'bin/pycbc_single_template',
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
               'bin/pycbc_upload_xml_to_gracedb',
               'bin/pycbc_dark_vs_bright_injections',
               'bin/gstlal/pycbc_calculate_likelihood',
               'bin/gstlal/pycbc_combine_likelihood',
               'bin/gstlal/pycbc_compute_far_from_snr_chisq_histograms',
               'bin/gstlal/pycbc_gen_ranking_data',
               'bin/gstlal/pycbc_pickle_horizon_distances',
               'bin/pycbc_make_coinc_pipedown_workflow',
               'bin/pycbc_make_html_page',
               'bin/pycbc_ligolw_find_playground',
               'bin/hdfcoinc/pycbc_make_coinc_search_workflow',
               'bin/hdfcoinc/pycbc_make_psd_estimation_workflow',
               'bin/hdfcoinc/pycbc_calculate_psd',
               'bin/hdfcoinc/pycbc_average_psd',
               'bin/pycbc_optimal_snr',
               'bin/pycbc_fit_sngl_trigs',
               'bin/hdfcoinc/pycbc_coinc_mergetrigs',
               'bin/hdfcoinc/pycbc_coinc_findtrigs',
               'bin/hdfcoinc/pycbc_coinc_bank2hdf',
               'bin/hdfcoinc/pycbc_coinc_trig2hdf',
               'bin/hdfcoinc/pycbc_coinc_statmap',
               'bin/hdfcoinc/pycbc_coinc_statmap_inj',
               'bin/hdfcoinc/pycbc_page_foreground',
               'bin/hdfcoinc/pycbc_page_foundmissed',
               'bin/hdfcoinc/pycbc_page_ifar',
               'bin/hdfcoinc/pycbc_page_snrifar',
               'bin/hdfcoinc/pycbc_page_sensitivity',
               'bin/hdfcoinc/pycbc_page_banktriggerrate',
               'bin/hdfcoinc/pycbc_coinc_hdfinjfind',
               'bin/hdfcoinc/pycbc_page_snrchi',
               'bin/hdfcoinc/pycbc_page_segments',
               'bin/hdfcoinc/pycbc_page_segtable',
               'bin/hdfcoinc/pycbc_page_segplot',
               'bin/hdfcoinc/pycbc_page_vetotable',
               'bin/hdfcoinc/pycbc_plot_psd_file',
               'bin/hdfcoinc/pycbc_plot_range',
               'bin/hdfcoinc/pycbc_foreground_censor',
               'bin/hdfcoinc/pycbc_plot_hist',
               'bin/hdfcoinc/pycbc_page_recovery',
               'bin/hwinj/pycbc_generate_hwinj',
               'bin/hwinj/pycbc_generate_hwinj_from_xml',
               'bin/hwinj/pycbc_plot_hwinj',
               'bin/hwinj/pycbc_insert_frame_hwinj',
               'bin/hdfcoinc/pycbc_strip_injections',
               'bin/sngl/pycbc_ligolw_cluster',
               'bin/sngl/pycbc_plot_bank',
               'bin/sngl/pycbc_plot_glitchgram',
               'bin/sngl/pycbc_plot_histogram',
               'bin/sngl/pycbc_plot_params',
               'bin/sngl/pycbc_plot_rates',
               'bin/sngl/pycbc_plot_timeseries',
               'bin/hdfcoinc/pycbc_page_injtable',
               'bin/pycbc_submit_dax',
               'bin/pycbc_submit_dax_stampede',
               'bin/hdfcoinc/pycbc_page_coinc_snrchi',
               'bin/hdfcoinc/pycbc_distribute_background_bins',
               'bin/hdfcoinc/pycbc_combine_statmap',
               'bin/hdfcoinc/pycbc_stat_dtphase',
               'bin/hdfcoinc/pycbc_plot_singles_vs_params',
               'bin/hdfcoinc/pycbc_plot_singles_timefreq',
               'bin/mvsc/pycbc_mvsc_get_features',
               'bin/pycbc_coinc_time',
               'bin/hdfcoinc/pycbc_plot_background_coincs',
               'bin/hdfcoinc/pycbc_plot_bank_bins',
               'bin/pygrb/pycbc_make_offline_grb_workflow',
               'bin/pygrb/pycbc_make_grb_summary_page',
               'bin/hdfcoinc/pycbc_merge_psds',
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
               'pycbc.io',
               ],
     package_data = {'pycbc.workflow': find_package_data('pycbc/workflow'), 
	             'pycbc.results': find_package_data('pycbc/results'),
                     'pycbc.tmpltbank': find_package_data('pycbc/tmpltbank')},
)

