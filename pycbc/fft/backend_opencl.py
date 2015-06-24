#  Copyright (C) 2014 Josh Willis
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with with program; see the file COPYING. If not, write to the
#  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
#  MA  02111-1307  USA

import pycbc
from .core import _list_available

_backend_dict = {'pyfft' : 'cl_pyfft'}
_backend_list = ['pyfft']

_alist = []
_adict = {}

if pycbc.HAVE_OPENCL:
    _alist, _adict = _list_available(_backend_list, _backend_dict)
    opencl_backend = 'pyfft'
else:
    opencl_backend = None

def set_backend(backend_list):
    global opencl_backend
    for backend in backend_list:
        if backend in _alist:
            opencl_backend = backend
            break

def get_backend():
    return _adict[opencl_backend]

