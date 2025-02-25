# Copyright (C) 2023 Gareth Cabourn Davies
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
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
Module to generate/manage the executable used for version information
in workflows
"""
import os
from pycbc.workflow.core import Executable


class VersioningExecutable(Executable):
    """
    Executable for getting version information
    """
    current_retention_level = Executable.FINAL_RESULT


def make_versioning_page(workflow, config_parser, out_dir, tags=None):
    """
    Make executable for versioning information
    """
    vers_exe = VersioningExecutable(
        workflow.cp,
        'page_versioning',
        out_dir=out_dir,
        ifos=workflow.ifos,
        tags=tags,
    )
    node = vers_exe.create_node()
    config_names = []
    exes = []
    for name, path in config_parser.items('executables'):
        exe_to_test = os.path.basename(path)
        if exe_to_test in exes:
            # executable is already part of the list,
            # find which index and add the name to the
            # one already stored
            path_idx = exes.index(exe_to_test)
            name_orig = config_names[path_idx]
            config_names[path_idx] = f"{name_orig},{name}"
        else:
            config_names.append(name)
            exes.append(exe_to_test)
    node.add_list_opt('--executables', exes)
    node.add_list_opt('--executables-names', config_names)
    node.new_output_file_opt(workflow.analysis_time, '.html', '--output-file')
    workflow.add_node(node)

    return node, node.output_files
