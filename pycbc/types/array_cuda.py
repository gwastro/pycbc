# Copyright (C) 2012  Alex Nitz
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

from pycuda.elementwise import ElementwiseKernel
from pycuda.reduction import ReductionKernel
from pycuda.tools import register_dtype
from pycuda.tools import context_dependent_memoize
from pycuda.tools import dtype_to_ctype
from pytools import match_precision
import numpy as np

@context_dependent_memoize
def get_norm_kernel(dtype_x, dtype_out):
    return ElementwiseKernel(
            "%(tp_x)s *x, %(tp_z)s *z" % {
                "tp_x": dtype_to_ctype(dtype_x),
                "tp_z": dtype_to_ctype(dtype_out),
                },
            "z[i] = norm(x[i])",
            "norm")

def squared_norm(a):
    dtype_out = match_precision(np.dtype('float64'),a.dtype)
    out = a._new_like_me(dtype=dtype_out)
    krnl = get_norm_kernel(a.dtype,dtype_out)
    krnl(a,out)
    return out     
 
# Define PYCUDA MAXLOC for both single and double precission ################## 
       
maxloc_preamble = """

    struct MAXLOCN{
        TTYPE max;
        LTYPE   loc;
        
        __device__
        MAXLOCN(){}
        
        __device__
        MAXLOCN(MAXLOCN const &src): max(src.max), loc(src.loc){}
        __device__
        MAXLOCN(MAXLOCN const volatile &src): max(src.max), loc(src.loc){}
        
        __device__
        MAXLOCN volatile &operator=( MAXLOCN const &src) volatile{
            max = src.max;
            loc = src.loc;
            return *this;
        }
    };
    
    __device__
    MAXLOCN maxloc_red(MAXLOCN a, MAXLOCN b){
        if (a.max > b.max)
            return a;
        else 
            return b;  
    }
    
    __device__
    MAXLOCN maxloc_start(){
        MAXLOCN t;
        t.max=0;
        t.loc=0;
        return t;
    }
    
    __device__
    MAXLOCN maxloc_map(TTYPE val, LTYPE loc){
        MAXLOCN t;
        t.max = val;
        t.loc = loc;
        return t;
    }

    """
    
maxloc_preamble_single = """
    #define MAXLOCN maxlocs
    #define TTYPE float
    #define LTYPE int
""" + maxloc_preamble

maxloc_preamble_double = """
    #define MAXLOCN maxlocd
    #define TTYPE double
    #define LTYPE long
""" + maxloc_preamble
    
maxloc_dtype_double = np.dtype([("max", np.float64), ("loc", np.int64)])
maxloc_dtype_single = np.dtype([("max", np.float32), ("loc", np.int32)])

register_dtype(maxloc_dtype_single, "maxlocs")
register_dtype(maxloc_dtype_double, "maxlocd")

mls = ReductionKernel(maxloc_dtype_single, neutral = "maxloc_start()",
        reduce_expr="maxloc_red(a, b)", map_expr="maxloc_map(x[i], i)",
        arguments="float *x", preamble=maxloc_preamble_single)

mld = ReductionKernel(maxloc_dtype_double, neutral = "maxloc_start()",
        reduce_expr="maxloc_red(a, b)", map_expr="maxloc_map(x[i], i)",
        arguments="double *x", preamble=maxloc_preamble_double)
        
max_loc = {'single':mls,'double':mld}




