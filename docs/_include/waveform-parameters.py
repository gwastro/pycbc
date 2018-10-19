#!/usr/bin/env python

# Copyright (C) 2018 Duncan Macleod, Collin Capano
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

"""Prints the usable waveform parameters as a RST table.
"""

# NOTE: the manual call to OrdereDict can be removed in favour of
# `ParameterList.description_dict` when gwastro/pycbc#2125 is merged
# and released

from __future__ import (print_function, absolute_import)
from pycbc import waveform
from _dict_to_rst import rst_dict_table

allparams = (waveform.td_waveform_params +
             waveform.fd_waveform_params +
             waveform.location_params)

tbl = rst_dict_table(allparams.description_dict,
                     key_format='``\'{0}\'``'.format,
                     header=('Parameter', 'Description'),
                     sort=False)

filename = 'waveform-parameters.rst'
with open(filename, 'w') as fp:
    print(tbl, file=fp)
