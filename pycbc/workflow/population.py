# Copyright (C) 2019 Vaibhav Tiwari
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
This module is responsible for setting up plotting jobs.
https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/NOTYETCREATED.html
"""
import urlparse, urllib
from pycbc.workflow.core import File, FileList, makedir, Executable

class PopExecutable(Executable):
    """ Executable for Population and Rates
    """
    current_retention_level = Executable.FINAL_RESULT

    # plots and final results should get the highest priority
    # on the job queue
    def create_node(self):
        node = Executable.create_node(self)
        node.set_priority(1000)
        return node

def estimate_rate_posterior(workflow, trig_files, bank_file, injection_files, 
                            ratesamples_file, out_dir, exclude=None, require=None, tags=None):
    tags = [] if tags is None else tags
    makedir(out_dir)
    secs = requirestr(workflow.cp.get_subsections('calculate_rate_posterior'), require)
    secs = excludestr(secs, exclude)
    files = FileList([])
    for tag in secs:
            node.add_input_opt('--statmap-file', trig_file)
            node.add_input_opt('--bank-file', bank_file)
            node.add_input_opt('--sim-files', injection_files)
            node.add_input_opt('--prior-samples', trig_file)
            output_file = File(parent.ifo_list, self.exe_name, segment,
                           extension='.hdf5', store_file=self.retain_files,
                           directory=self.out_dir, tags=tags)
            node.add_output_opt('--output-file', output_file)
            workflow += node
            files += node.output_file
    return files