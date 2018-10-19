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

"""Prints an RST table of available distributions from the distributions
module.
"""
from __future__ import (print_function, absolute_import)
from pycbc.distributions import distribs
from _dict_to_rst import (rst_dict_table, format_class)

tbl = rst_dict_table(distribs, key_format='``\'{0}\'``'.format,
                     header=('Name', 'Class'),
                     val_format=format_class)

filename = 'distributions-table.rst'
with open(filename, 'w') as fp:
    print(tbl, file=fp)
