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
from urllib.request import pathname2url
from urllib.parse import urljoin
from pycbc.workflow.core import File, FileList, makedir, Executable, Node

class VersioningExecutable(Executable):
    current_retention_level = Executable.FINAL_RESULT

    def create_node(self):
        node=Node(self)
        return node

def make_versioning_page(workflow, cp, out_dir, tags=None):
    vers_exe = VersioningExecutable(
        workflow.cp,
        'page_versioning',
        out_dir=out_dir,
        ifos=workflow.ifos,
        tags=tags,
    )
    node = vers_exe.create_node()
    exe_names = []
    exe_paths = []
    for name, path in cp.items('executables'):
        if path in exe_paths: continue
        exe_names.append(name)
        file_url = urljoin('file:', pathname2url(path))
        exe_to_test = File(workflow.ifos, '',
                        workflow.analysis_time, file_url=file_url)
        exe_to_test.add_pfn(file_url, site='local')
        exe_paths.append(exe_to_test)
    node.add_input_list_opt('--executables-files', FileList(exe_paths))
    node.add_list_opt('--executables-names', exe_names)
    node.new_output_file_opt(workflow.analysis_time, '.html', '--output-file')
    workflow.add_node(node)

    return node
