#!/usr/bin/python

# Copyright (C) 2015 Christopher M. Biwer
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

import os.path
import pycbc.results
from jinja2 import Environment, FileSystemLoader
from pycbc.results.metadata import save_html_with_metadata
from pycbc.results.render import get_embedded_config

def render_workflow_html_template(filename, subtemplate, filelists):
    """ Takes a template and list of tuples where the elements are File objects.
    """

    dir = os.path.dirname(filename)

    # render subtemplate
    subtemplate_dir = pycbc.results.__path__[0] + '/templates/wells'
    env = Environment(loader=FileSystemLoader(subtemplate_dir))
    env.globals.update(get_embedded_config=get_embedded_config)
    env.globals.update(len=len)
    subtemplate = env.get_template(subtemplate)
    context = {'filelists' : filelists,
               'dir' : dir}
    output = subtemplate.render(context)

    # save as html page
    kwds = {'render-function' : 'render_tmplt'}
    save_html_with_metadata(str(output), filename, None, kwds)
