#!/usr/bin/python

# Copyright (C) 2015 Ian Harry
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Generals
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

import logging
import subprocess
import urllib.parse

import lal, lalframe
import pycbc.version

def get_library_version_info():
    """This will return a list of dictionaries containing versioning
    information about the various LIGO libraries that PyCBC will use in an
    analysis run."""
    library_list = []

    def add_info_new_version(info_dct, curr_module, extra_str):
        vcs_object = getattr(curr_module, extra_str +'VCSInfo')
        info_dct['ID'] =  vcs_object.vcsId
        info_dct['Status'] = vcs_object.vcsStatus
        info_dct['Version'] = vcs_object.version
        info_dct['Tag'] = vcs_object.vcsTag
        info_dct['Author'] = vcs_object.vcsAuthor
        info_dct['Branch'] = vcs_object.vcsBranch
        info_dct['Committer'] = vcs_object.vcsCommitter
        info_dct['Date'] = vcs_object.vcsDate

    lalinfo = {}
    lalinfo['Name'] = 'LAL'
    try:
        lalinfo['ID'] = lal.VCSId
        lalinfo['Status'] = lal.VCSStatus
        lalinfo['Version'] = lal.VCSVersion
        lalinfo['Tag'] = lal.VCSTag
        lalinfo['Author'] = lal.VCSAuthor
        lalinfo['Branch'] = lal.VCSBranch
        lalinfo['Committer'] = lal.VCSCommitter
        lalinfo['Date'] = lal.VCSDate
    except AttributeError:
        add_info_new_version(lalinfo, lal, '')
    library_list.append(lalinfo)

    lalframeinfo = {}
    try:
        lalframeinfo['Name'] = 'LALFrame'
        lalframeinfo['ID'] = lalframe.FrameVCSId
        lalframeinfo['Status'] = lalframe.FrameVCSStatus
        lalframeinfo['Version'] = lalframe.FrameVCSVersion
        lalframeinfo['Tag'] = lalframe.FrameVCSTag
        lalframeinfo['Author'] = lalframe.FrameVCSAuthor
        lalframeinfo['Branch'] = lalframe.FrameVCSBranch
        lalframeinfo['Committer'] = lalframe.FrameVCSCommitter
        lalframeinfo['Date'] = lalframe.FrameVCSDate
    except AttributeError:
        add_info_new_version(lalframeinfo, lalframe, 'Frame')
    library_list.append(lalframeinfo)

    lalsimulationinfo = {}
    lalsimulationinfo['Name'] = 'LALSimulation'
    try:
        import lalsimulation
        lalsimulationinfo['ID'] = lalsimulation.SimulationVCSId
        lalsimulationinfo['Status'] = lalsimulation.SimulationVCSStatus
        lalsimulationinfo['Version'] = lalsimulation.SimulationVCSVersion
        lalsimulationinfo['Tag'] = lalsimulation.SimulationVCSTag
        lalsimulationinfo['Author'] = lalsimulation.SimulationVCSAuthor
        lalsimulationinfo['Branch'] = lalsimulation.SimulationVCSBranch
        lalsimulationinfo['Committer'] = lalsimulation.SimulationVCSCommitter
        lalsimulationinfo['Date'] = lalsimulation.SimulationVCSDate
    except AttributeError:
        add_info_new_version(lalsimulationinfo, lalsimulation, 'Simulation')
    except ImportError:
        pass
    library_list.append(lalsimulationinfo)

    pycbcinfo = {}
    pycbcinfo['Name'] = 'PyCBC'
    pycbcinfo['ID'] = pycbc.version.version
    pycbcinfo['Status'] = pycbc.version.git_status
    pycbcinfo['Version'] = pycbc.version.release or ''
    pycbcinfo['Tag'] = pycbc.version.git_tag
    pycbcinfo['Author'] = pycbc.version.git_author
    pycbcinfo['Builder'] = pycbc.version.git_builder
    pycbcinfo['Branch'] = pycbc.version.git_branch
    pycbcinfo['Committer'] = pycbc.version.git_committer
    pycbcinfo['Date'] = pycbc.version.git_build_date
    library_list.append(pycbcinfo)

    return library_list

def get_code_version_numbers(executable_names, executable_files):
    """Will extract the version information from the executables listed in
    the executable section of the supplied ConfigParser object.

    Returns
    --------
    dict
        A dictionary keyed by the executable name with values giving the
        version string for each executable.
    """
    code_version_dict = {}
    for exe_name, value in zip(executable_names, executable_files):
        value = urllib.parse.urlparse(value)
        logging.info("Getting version info for %s", exe_name)
        version_string = None
        if value.scheme in ['gsiftp', 'http', 'https']:
            code_version_dict[exe_name] = "Using bundle downloaded from %s" % value
        elif value.scheme == 'singularity':
            txt = (
                "Executable run from a singularity image. See config file "
                "and site catalog for details of what image was used."
            )
            code_version_dict[exe_name] = txt
        else:
            try:
                version_string = subprocess.check_output(
                    [value.path, '--version'],
                    stderr=subprocess.STDOUT
                ).decode()
            except subprocess.CalledProcessError:
                version_string = "Executable fails on {} --version"
                version_string = version_string.format(value.path)
            except OSError:
                version_string = "Executable doesn't seem to exist(!)"
            code_version_dict[exe_name] = version_string
    return code_version_dict
