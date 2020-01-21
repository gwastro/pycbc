#!/usr/bin/env python

# Copyright (C) 2020 Collin Capano
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

"""Prints an RST table of data options for inference models.
"""
from __future__ import (print_function, absolute_import)
try:
    from StringIO import StringIO
except ModuleNotFoundError:  # python 3
    from io import StringIO

import re

from pycbc.inference.models import data_utils

from _dict_to_rst import (rst_dict_table, format_class)

# convenience class for storing row data
class Row(object):
    divider = '  '

    def __init__(self):
        self.groupmsg = ''
        self.option = ''
        self.metavar = ''
        self.helpmsg = ''

    @property
    def isgroup(self):
        return self.groupmsg != ''

    @property
    def grouplen(self):
        return max(map(len, self.groupmsg.split('\n')))

    @property
    def helplen(self):
        return max(map(len, self.helpmsg.split('\n')))

    def format(self, maxlen, optlen, metalen, helplen):
        if self.isgroup:
            out = ['{msg:^{width}}'.format(msg=msg, width=maxlen)
                   for msg in self.groupmsg.split('\n')]
        else:
            tmplt = '{msg:<{rpad}}'
            out = []
            for ii, helpmsg in enumerate(self.helpmsg.split('\n')):
                if ii == 0:
                    optstr = self.option
                    metastr = self.metavar
                else:
                    optstr = ''
                    metastr = ''
                optstr = tmplt.format(msg=optstr, rpad=optlen)
                metastr = tmplt.format(msg=metastr, rpad=metalen)
                helpstr = tmplt.format(msg=helpmsg, rpad=helplen)
                out.append(self.divider.join([optstr, metastr, helpstr]))
        return '\n'.join(out)

    def __len__(self):
        if self.isgroup:
            baselen = self.grouplen
        else:
            baselen = len(self.option) + len(self.metavar) + self.helplen
        return baselen + 2*len(self.divider)


# create a data parser that has the options
parser = data_utils.create_data_parser()

# dump the help message to a file buffer
fp = StringIO()
parser.print_help(file=fp)

# Regular expressions to interpret the help message:
# Lines with a "--option stuff" (i.e., a single space after the option) include
# metadata. Lines with "--option   msg" (i.e., multiple spaces after the
# option) contain no metadata, and just go straight to the help message.
regx_optmeta = re.compile(
    r'^\s+((-\S, )*)--(?P<option>\S+)\s(?P<metavar>\S.+)')
regx_optmsg = re.compile(r'^\s+((-\S, )*)--(?P<option>\S+)\s+(?P<msg>.+)')
# Note: optmsg will match optmeta lines, so need to test for optmeta first.

# Lines that start with whitespace but do not match optmeta or optmsg will be
# assumed to be the rest of the help message for either an option or a group.
regx_helpmsg = re.compile(r'^\s+(?P<msg>.+)')
# Note: this will match both optmeta and optmsg lines, so need to test for
# for those before this.

# Lines that do not start with whitespace will be considered to be the start of
# a new option group. This is mutually exclusive of all the previous regxs,
# since they all require whitespace at the start of a line.
regx_newgroup = re.compile(r'^(?P<msg>\S.+)')

# now format the string buffer into a rst table
fp.seek(0)
# we want to skip everything up to the "optional arguments:"
skip = True
while skip:
    line = fp.readline()
    m = regx_newgroup.match(line)
    if m is not None:
        skip = m.group('msg') != 'optional arguments:'

# advance past the 'optional arguments:' and the 'help' line
line = fp.readline()
line = fp.readline()

# now read through the rest of the lines, converting options into a list of
# tuples with order (option, meta data, help message), grouped by option groups
# add a header row
header = Row()
header.option = 'Name'
header.metavar = 'Syntax'
header.helpmsg = 'Description'
table = [header]

while line:
    # determine if the line is a new group
    newgroup = regx_newgroup.match(line)
    if newgroup:
        groupmsg = [newgroup.group('msg')]
        # continue reading until we get a blank line or an option
        line = fp.readline()
        while not (line == '' or line == '\n' or regx_optmsg.match(line)):
            groupmsg.append(regx_helpmsg.match(line).group('msg'))
            line = fp.readline()
        # compile the group message
        row = Row()
        row.groupmsg = '\n'.join(groupmsg)
        table.append(row)
    # check if the line contains an option
    m = regx_optmsg.match(line)
    if m:
        row = Row()
        helpmsg = []
        row.option = m.group('option')
        # check if the line is an option with metadata
        meta = regx_optmeta.match(line)
        if meta:
            row.metavar = meta.group('metavar')
        else:
            helpmsg.append(m.group('msg'))
        # continue reading to get the help message
        line = fp.readline()
        while not (line == '' or line == '\n' or regx_optmsg.match(line)):
            helpmsg.append(regx_helpmsg.match(line).group('msg'))
            line = fp.readline()
        # compile the help message
        row.helpmsg = '\n'.join(helpmsg)
        # add the row to the table
        table.append(row)
    # read the next line for the loop
    line = fp.readline()

# Now format the table
# get the maximum row length
maxlen = max(map(len, table))
# get the maximum length of each column
optlen = max([len(row.option) for row in table])
metalen = max([len(row.metavar) for row in table])
helplen = max([row.helplen for row in table])

# create row separators
# "major" will have == lines
majorbreak = Row()
majorbreak.option = '='*optlen
majorbreak.metavar = '='*metalen
majorbreak.helpmsg = '='*helplen
# "minor" will have -- lines
minorbreak = Row()
minorbreak.option = '-'*optlen
minorbreak.metavar = '-'*metalen
minorbreak.helpmsg = '-'*helplen

formatargs = [maxlen, optlen, metalen, helplen]

# Write the formatted table to file
filename = 'inference_data_opts-table.rst'
out = open(filename, 'w')
print(majorbreak.format(*formatargs), file=out)
# print the header
print(header.format(*formatargs), file=out)
print(majorbreak.format(*formatargs), file=out)
# print everything else in the table
for ii in range(1, len(table)):
    row = table[ii]
    if row.isgroup or ii == len(table)-1 or table[ii+1].isgroup:
        # use a majorbreak
        rbreak = majorbreak
    else:
        # use a minorbreak
        rbreak = minorbreak
    print(row.format(*formatargs), file=out)
    print(rbreak.format(*formatargs), file=out)
out.close()
