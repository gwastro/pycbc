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

from ConfigParser import ConfigParser
from jinja2 import Environment, FileSystemLoader

import pycbc.results

def setup_template_render(path, config_path):

   # initialization
   cp = ConfigParser()
   output = ''

   # read configuration file
   if config_path:
       cp.read(config_path)

       # render template
       if cp.has_option(path, 'render-function'):
           render_function_name = cp.get_opt(path)
           render_function = eval(render_function_name)
           output = render_function(path, cp)
       else:
           output = render_default(path, cp)

   # if no configuration file is present
   # then render the default template
   else:
       output = render_deault(path, cp)

   return output

def render_default(path, cp):

   # define slug
   slug = path.split('/')[-1].split('.')[0]

   # render template
   template_dir = pycbc.results.__path__[0] + '/templates/files'
   env = Environment(loader=FileSystemLoader(template_dir))
   template = env.get_template('file_default.html')
   context = {'path' : path,
              'slug' : slug}
   output = template.render(context)

   return output

