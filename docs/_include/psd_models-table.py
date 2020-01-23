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

"""Prints an RST table of available psd models.
"""

from __future__ import (print_function, absolute_import)

from pycbc import psd

from _dict_to_rst import (rst_dict_table, format_function)

psds = {p: getattr(psd.analytical, p)
        for p in psd.analytical.get_psd_model_list()}

tbl = rst_dict_table(psds, key_format='``{0}``'.format,
                     header=('Name', 'Function'),
                     val_format=format_function)

filename = 'psd_models-table.rst'
with open(filename, 'w') as fp:
    print(tbl, file=fp)
