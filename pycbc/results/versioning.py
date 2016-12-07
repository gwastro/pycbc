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
import os
import subprocess
from pycbc.results import save_fig_with_metadata, html_escape

import lal, lalframe, lalmetaio, lalinspiral, lalsimulation
import glue.git_version, pycbc.version

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

    lalmetaioinfo = {}
    lalmetaioinfo['Name'] = 'LALMetaIO'
    try:
        lalmetaioinfo['ID'] = lalmetaio.MetaIOVCSId
        lalmetaioinfo['Status'] = lalmetaio.MetaIOVCSStatus
        lalmetaioinfo['Version'] = lalmetaio.MetaIOVCSVersion
        lalmetaioinfo['Tag'] = lalmetaio.MetaIOVCSTag
        lalmetaioinfo['Author'] = lalmetaio.MetaIOVCSAuthor
        lalmetaioinfo['Branch'] = lalmetaio.MetaIOVCSBranch
        lalmetaioinfo['Committer'] = lalmetaio.MetaIOVCSCommitter
        lalmetaioinfo['Date'] = lalmetaio.MetaIOVCSDate
    except AttributeError:
        add_info_new_version(lalmetaioinfo, lalmetaio, 'MetaIO')
    library_list.append(lalmetaioinfo)

    lalinspiralinfo = {}
    lalinspiralinfo['Name'] = 'LALInspiral'
    try:
        lalinspiralinfo['ID'] = lalinspiral.InspiralVCSId
        lalinspiralinfo['Status'] = lalinspiral.InspiralVCSStatus
        lalinspiralinfo['Version'] = lalinspiral.InspiralVCSVersion
        lalinspiralinfo['Tag'] = lalinspiral.InspiralVCSTag
        lalinspiralinfo['Author'] = lalinspiral.InspiralVCSAuthor
        lalinspiralinfo['Branch'] = lalinspiral.InspiralVCSBranch
        lalinspiralinfo['Committer'] = lalinspiral.InspiralVCSCommitter
        lalinspiralinfo['Date'] = lalinspiral.InspiralVCSDate
    except AttributeError:
        add_info_new_version(lalinspiralinfo, lalinspiral, 'Inspiral')
    library_list.append(lalinspiralinfo)

    lalsimulationinfo = {}
    lalsimulationinfo['Name'] = 'LALSimulation'
    try:
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
    library_list.append(lalsimulationinfo)

    glueinfo = {}
    glueinfo['Name'] = 'Glue'
    glueinfo['ID'] = glue.git_version.id
    glueinfo['Status'] = glue.git_version.status
    glueinfo['Version'] = glue.git_version.version
    glueinfo['Tag'] = glue.git_version.tag
    glueinfo['Author'] = glue.git_version.author
    glueinfo['Builder'] = glue.git_version.builder
    glueinfo['Branch'] = glue.git_version.branch
    glueinfo['Committer'] = glue.git_version.committer
    glueinfo['Date'] = glue.git_version.date
    library_list.append(glueinfo)

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

def write_library_information(path):
    library_list = get_library_version_info()
    for curr_lib in library_list:
        lib_name = curr_lib['Name']
        text = ''
        for key, value in curr_lib.items():
            text+='<li> %s : %s </li>\n' %(key,value)
        kwds = {'render-function' : 'render_text',
                'title' : '%s Version Information'%lib_name,
        }
        
        save_fig_with_metadata(html_escape(text), 
          os.path.join(path,'%s_version_information.html' %(lib_name)), **kwds)

def get_code_version_numbers(cp):
    """Will extract the version information from the executables listed in
    the executable section of the supplied ConfigParser object.

    Returns
    --------
    dict
        A dictionary keyed by the executable name with values giving the
        version string for each executable.
    """
    from pycbc.workflow.core import check_output_error_and_retcode
    code_version_dict = {}
    for item, value in cp.items('executables'):
        path, exe_name = os.path.split(value)
        version_string = None
        if value.startswith('gsiftp://') or value.startswith('http://'):
           code_version_dict[exe_name] = "Using bundle downloaded from %s" % value
        else:
            try:
                # FIXME: Replace with this version when python 2.7 is guaranteed
                # version_output = subprocess.check_output([value, '--version'],
                #                                         stderr=subprocess.STDOUT) 
                # Start of legacy block
                if value.startswith('file://'):
                    value = value[7:]
                output, error, retcode = \
                               check_output_error_and_retcode([value, '--version'])
                if not retcode == 0:
                    raise subprocess.CalledProcessError(retcode, '')
                # End of legacy block
                version_output = (output + error).replace('\n', ' ').split()
                # Look for a version
                if "Id:" in version_output:
                    index = version_output.index("Id:") + 1
                    version_string = 'Version is %s.' %(version_output[index],)
                elif "LALApps:" in version_output:
                    index = version_output.index("LALApps:") + 3
                    version_string = 'Version (lalapps) is %s.' %(version_output[index],)
                if version_string is None:
                    version_string = "Cannot identify version string in output."
            except subprocess.CalledProcessError:
                version_string = "Executable fails on %s --version" %(value)
            except OSError:
                version_string = "Executable doesn't seem to exist(!)"
            code_version_dict[exe_name] = version_string
    return code_version_dict

def write_code_versions(path, cp):
    code_version_dict = get_code_version_numbers(cp)
    html_text = ''
    for key,value in code_version_dict.items():
        html_text+= '<li><b>%s</b>: %s </li>\n' %(key,value.replace('@', '&#64;'))
    kwds = {'render-function' : 'render_text',
            'title' : 'Version Information from Executables',
    }
    save_fig_with_metadata(html_escape(html_text), 
        os.path.join(path,'version_information_from_executables.html'), **kwds)

def create_versioning_page(path, cp):
    logging.info("Entering versioning module")
    if not os.path.exists(path):
        os.mkdir(path)
    write_library_information(path)
    write_code_versions(path, cp)
    logging.info("Leaving versioning module")
