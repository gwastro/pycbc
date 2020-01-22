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
import textwrap

from pycbc.inference.models import data_utils

# wrapper for long metavars
metavar_txtwrap = textwrap.TextWrapper(width=34, break_long_words=False)

# convenience class for storing row data
class Row(object):

    def __init__(self, divider=None, lpad=0, wrap_option=True):
        if divider is None:
            divider = ' | '
        self.divider = divider
        self.lborder = ' '*lpad + divider[1:]
        self.rborder = divider[:-1]
        self._groupmsg = ''
        self.wrap_option = wrap_option
        self._option = ''
        self._metavar = ''
        self._helpmsg = ''

    @property
    def option(self):
        return self._option

    @option.setter
    def option(self, option):
        if self.wrap_option:
            # add `` around option string
            option = '``'+option+'``'
        self._option = option

    @property
    def metavar(self):
        return self._metavar

    @metavar.setter
    def metavar(self, metavar):
        # text wrap if metavar is more than 34 characters wide
        self._metavar = '\n'.join(metavar_txtwrap.wrap(metavar))
            
    @staticmethod
    def replace_doubledash(str_):
        """Replaces all instances of --arg with ``arg`` in a string."""
        pattern = r"--\w+-?\w*"
        for s in re.findall(pattern, str_):
            rep = '``' + s.replace('--', '') + '``'
            str_ = re.sub(s, rep, str_)
        return str_

    @property
    def helpmsg(self):
        return self._helpmsg

    @helpmsg.setter
    def helpmsg(self, msg):
        # replace all instances of --arg with ``arg``
        self._helpmsg = self.replace_doubledash(msg)

    @property
    def groupmsg(self):
        return self._groupmsg

    @groupmsg.setter
    def groupmsg(self, msg):
        # replace all instances of --arg with ``arg``
        self._groupmsg = self.replace_doubledash(msg)

    @property
    def isgroup(self):
        return self.groupmsg != ''

    @property
    def grouplen(self):
        return max(map(len, self.groupmsg.split('\n')))

    @property
    def metavarlen(self):
        return max(map(len, self.metavar.split('\n')))

    @property
    def helplen(self):
        return max(map(len, self.helpmsg.split('\n')))

    def format(self, maxlen, optlen, metalen, helplen):
        if self.isgroup:
            out = ['{lbdr}{msg:<{width}}{rbdr}'.format(lbdr=self.lborder,
                                                       msg=msg, width=maxlen,
                                                       rbdr=self.rborder)
                   for msg in self.groupmsg.split('\n')]
        else:
            tmplt = '{msg:<{rpad}}'
            out = []
            metavar = self.metavar.split('\n')
            helpmsg = self.helpmsg.split('\n')
            nlines = max(len(metavar), len(helpmsg))
            for ii in range(nlines):
                if ii == 0:
                    optstr = self.option
                else:
                    optstr = ''
                if ii < len(metavar):
                    metastr = metavar[ii]
                else:
                    metastr = ''
                if ii < len(helpmsg):
                    helpstr = helpmsg[ii]
                else:
                    helpstr = ''
                optstr = tmplt.format(msg=optstr, rpad=optlen)
                metastr = tmplt.format(msg=metastr, rpad=metalen)
                helpstr = tmplt.format(msg=helpstr, rpad=helplen)
                #rowstr = self.divider.join([optstr, helpstr, metastr])
                rowstr = self.divider.join([optstr, metastr, helpstr])
                # add borders
                rowstr = '{}{}{}'.format(self.lborder, rowstr, self.rborder)
                out.append(rowstr)
        return '\n'.join(out)

    def __len__(self):
        if self.isgroup:
            baselen = self.grouplen
        else:
            baselen = len(self.option) + self.metavarlen + self.helplen
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
header = Row(wrap_option=False)
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
        # remove the list of PSDs from the fake strain and psd model arguments,
        # referring instead to the psd table
        if (m.group('option') == 'fake-strain' or
            m.group('option') == 'psd-model'):
            # strip off everything after the "Choose from"
            helpmsg = row.helpmsg.replace('\n', ' ')
            idx = re.search(r"Choose +from", helpmsg).end()
            helpmsg = helpmsg[:idx] + " any available PSD model"
            # for fake strain, add zero Noise
            if m.group('option') == 'fake-strain':
                helpmsg += ', or ``zeroNoise``.'
            else:
                helpmsg += '.'
            row.helpmsg = textwrap.fill(helpmsg, width=54)
        # add the row to the table
        table.append(row)
    else:
        # read the next line for the loop
        line = fp.readline()

# Now format the table
# get the maximum length of each column
optlen = max([len(row.option) for row in table])
metalen = max([row.metavarlen for row in table])
helplen = max([row.helplen for row in table])
maxlen = optlen + metalen + helplen + 6  # the 6 is for the 2 dividers

# create row separators
# "major" will have == lines
majorline = Row(divider='=+=', wrap_option=False)
majorline.option = '='*optlen
majorline.metavar = '='*metalen
majorline.helpmsg = '='*helplen
# "minor" will have -- lines
minorline = Row(divider='-+-', wrap_option=False)
minorline.option = '-'*optlen
minorline.metavar = '-'*metalen
minorline.helpmsg = '-'*helplen

formatargs = [maxlen, optlen, metalen, helplen]

# Write the formatted table to file
filename = 'inference_data_opts-table.rst'
out = open(filename, 'w')
print(minorline.format(*formatargs), file=out)
# print the header
print(header.format(*formatargs), file=out)
print(majorline.format(*formatargs), file=out)
# print everything else in the table
for ii in range(1, len(table)):
    row = table[ii]
    print(row.format(*formatargs), file=out)
    print(minorline.format(*formatargs), file=out)
out.close()
