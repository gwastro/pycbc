# Copyright (C) 2020 Alexander Nitz
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
"""
This modules provides classes for evaluating Gaussian distributions.
"""

import importlib
import numpy.random
from pycbc import VARARGS_DELIM

class External(object):
    """ Distribution defined by external cdfinv and logpdf functions
 
    To add to an inference configuration file:
 
    .. code-block:: ini

        [prior-param1+param2]
        module = custom_mod
        logpdf = custom_function_name
        cdfinv = custom_function_name2
   
    Parameters
    ----------
    params : list
        list of parameter names 
    logpdf : function
        function which returns the logpdf
    cdfinv : function
        function which applies the invcdf
        
    Examples
    --------
    To instantate by hand and example of function format
    
    >>> import numpy
    >>> params = 
    >>> def logpdf(x=None, y=None):
    ...     p = numpy.ones(len(x))
    ...     res = {}
    ...     res['x'] = p
    ...     res['y'] = p
    ...     return res
    >>>
    >>> def cdfinv(**kwds):
    ...     return kwds
    >>> e = External(['x', 'y'], logpdf, cdfinv)
    >>> e.rvs(size=10) 
    """
    name = "external"

    def __init__(self, params, logpdf, cdfinv):
        self.params = params
        self.logpdf = logpdf
        self.cdfinv = cdfinv

    def rvs(self, size=1):
        draw = {}
        for param in self.params:
            draw[param] = numpy.random.uniform(0, 1, size=size)
        return self.cdfinv(**draw)

    @classmethod
    def from_config(cls, cp, section, variable_args):
        tag = variable_args
        params = variable_args.split(VARARGS_DELIM)
        modulestr = cp.get_opt_tag(section, 'module', tag)
        logpdfstr = cp.get_opt_tag(section, 'logpdf', tag)
        cdfinvstr = cp.get_opt_tag(section, 'cdfinv', tag)
        
        mod = importlib.import_module(modulestr)
        logpdf = getattr(mod, logpdfstr)
        cdfinv = getattr(mod, cdfinvstr)
        return cls(params, logpdf, cdfinv)

__all__ = ['external']
