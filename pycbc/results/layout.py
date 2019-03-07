# Copyright (C) 2015 Alexander Harvey Nitz
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
""" This module contains result page layout and numbering helper functions
"""
import os.path
from six.moves import zip_longest

def two_column_layout(path, cols, **kwargs):
    """ Make a well layout in a two column format

    Parameters
    ----------
    path: str
        Location to make the well html file
    cols: list of tuples
        The format of the items on the well result section. Each tuple
        contains the two files that are shown in the left and right hand
        side of a row in the well.html page.
    """
    path = os.path.join(os.getcwd(), path, 'well.html')
    from pycbc.results.render import render_workflow_html_template
    render_workflow_html_template(path, 'two_column.html', cols, **kwargs)

def single_layout(path, files, **kwargs):
    """ Make a well layout in  single column format

    path: str
        Location to make the well html file
    files: list of pycbc.workflow.core.Files
        This list of images to show in order within the well layout html file.
    """
    two_column_layout(path, [(f,) for f in files], **kwargs)

def grouper(iterable, n, fillvalue=None):
    """ Group items into chunks of n length
    """
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)

def group_layout(path, files, **kwargs):
    """ Make a well layout in chunks of two from a list of files

    path: str
        Location to make the well html file
    files: list of pycbc.workflow.core.Files
        This list of images to show in order within the well layout html file.
        Every two are placed on the same row.
    """
    if len(files) > 0:
        two_column_layout(path, list(grouper(files, 2)), **kwargs)

class SectionNumber(object):
    """ Class to help with numbering sections in an output page.
    """
    def __init__(self, base, secs):
        """ Create section numbering instance

        Parameters
        ----------
        base: path
            The path of the of output html results directory
        secs: list of strings
            List of the subsections of the output html page
        """
        self.base = base
        self.secs = secs
        self.name = {}
        self.count = {}
        self.num = {}

        for num, sec in enumerate(secs):
            self.name[sec] = '%s._%s' % (num + 1, sec)
            self.num[sec] = num
            self.count[sec] = 1

    def __getitem__ (self, path):
        """ Return the path to use for the given subsection with numbering
        included. The numbering increments for each new subsection request. If
        a section is re-requested, it gets the original numbering.
        """
        if path in self.name:
            name = self.name[path]
        else:
            sec, subsec = path.split('/')
            subnum = self.count[sec]
            num = self.num[sec]
            name = '%s/%s.%02d_%s' % (self.name[sec], num + 1, subnum, subsec)
            self.count[sec] += 1
            self.name[path] = name
        path = os.path.join(os.getcwd(), self.base, name)
        return path
