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

import os.path, types

from ConfigParser import ConfigParser
from jinja2 import Environment, FileSystemLoader

import pycbc.results
from pycbc.workflow.segment import fromsegmentxml


def get_embedded_config(filename):
    """ Attempt to load config data attached to file
    """
    def check_option(self, section, name):
        return (self.has_section(section) and
               (self.has_option(section, name) or (name in self.defaults())))
    
    try:
        cp = pycbc.results.load_metadata_from_file(filename)
    except TypeError:
        cp = ConfigParser()
        
    cp.check_option = types.MethodType(check_option, cp)
    return cp

def setup_template_render(path, config_path):
    """ This function is the gateway for rendering a template for a file.
    """

    # initialization
    cp = get_embedded_config(path)
    output = ''

    # read configuration file
    if os.path.exists(config_path):
        cp.read(config_path)

        # render template
        filename = os.path.basename(path)
        if cp.has_option(filename, 'render-function'):
            render_function_name = cp.get(filename, 'render-function')
            render_function = eval(render_function_name)
            output = render_function(path, cp)
        else:
            output = render_default(path, cp)

    # if no configuration file is present
    # then render the default template
    else:
        output = render_default(path, cp)

    return output

def render_default(path, cp):
    """ This is the default function that will render a template to a string of HTML. The
    string will be for a drop-down tab that contains a link to the file.

    If the file extension requires information to be read, then that is passed to the
    content variable (eg. a segmentlistdict).
    """

    # define filename and slug from path
    filename = os.path.basename(path)
    slug = filename.replace('.', '_')

    # initializations
    content = None

    if path.endswith('.xml') or path.endswith('.xml.gz'):
        # segment or veto file return a segmentslistdict instance
        if 'SEG' in path or 'VETO' in path:
            with open(path, 'r') as xmlfile:
                content = fromsegmentxml(xmlfile, dict=True)
    elif path.endswith('.ini'):
        with open(path, 'rb') as f_handle:
            content = f_handle.read()

    # render template
    template_dir = pycbc.results.__path__[0] + '/templates/files'
    env = Environment(loader=FileSystemLoader(template_dir))
    env.globals.update(abs=abs)
    template = env.get_template('file_default.html')
    context = {'filename' : filename,
               'slug'     : slug,
               'cp'       : cp,
               'content'  : content}
    output = template.render(context)

    return output

def render_glitchgram(path, cp):
    """ Render a glitchgram file template.
    """

    # define filename and slug from path
    filename = os.path.basename(path)
    slug = filename.replace('.', '_')

    # render template
    template_dir = pycbc.results.__path__[0] + '/templates/files'
    env = Environment(loader=FileSystemLoader(template_dir))
    env.globals.update(abs=abs)
    template = env.get_template(cp.get(filename, 'template'))
    context = {'filename' : filename,
               'slug'     : slug,
               'cp'       : cp}
    output = template.render(context)

    return output
