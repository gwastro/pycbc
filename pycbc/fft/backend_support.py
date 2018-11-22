# Copyright (C) 2012  Josh Willis, Andrew Miller
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

#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
"""
This package provides a front-end to various fast Fourier transform
implementations within PyCBC.
"""

import pycbc
import pycbc.scheme


# These are global variables, that are modified by the various scheme-
# dependent submodules, to maintain a list of all possible backends
# for all possible schemes that are available at runtime.  This list
# and dict are then used when parsing command-line options.

_all_backends_list = []
_all_backends_dict = {}

# The following is the function called by each scheme's setup to add whatever new
# backends may have been found to the global list.  Since some backends may be
# shared, we must first check to make sure that the item in the list is not already
# in the global list, and we assume that the keys to the dict are in one-to-one
# correspondence with the items in the list.

def _update_global_available(new_list, new_dict, global_list, global_dict):
    for item in new_list:
        if item not in global_list:
            global_list.append(item)
            global_dict.update({item:new_dict[item]})

def get_backend_modules():
    return _all_backends_dict.values()

def get_backend_names():
    return _all_backends_dict.keys()

BACKEND_PREFIX="pycbc.fft.backend_"

@pycbc.scheme.schemed(BACKEND_PREFIX)
def set_backend(backend_list):
    err_msg = "This function is a stub that should be overridden using "
    err_msg += "the scheme. You shouldn't be seeing this error!"
    raise ValueError(err_msg)

@pycbc.scheme.schemed(BACKEND_PREFIX)
def get_backend():
    err_msg = "This function is a stub that should be overridden using "
    err_msg += "the scheme. You shouldn't be seeing this error!"
    raise ValueError(err_msg)

# Import all scheme-dependent backends, to get _all_backends accurate:

for scheme_name in ["cpu", "mkl", "cuda"]:
    try:
        mod = __import__('pycbc.fft.backend_' + scheme_name, fromlist = ['_alist', '_adict'])
        _alist = getattr(mod, "_alist")
        _adict = getattr(mod, "_adict")
        _update_global_available(_alist, _adict, _all_backends_list,
                                 _all_backends_dict)
    except ImportError:
        pass
