# Copyright (C) 2012  Alex Nitz
#
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
"""This file contains backported functionality from future versions of librarays.
Mostly done with monkey-patching.
"""

# Add in missing numpy functionality
import numpy

def unique(ar, return_index=False, return_inverse=False):
    """
        KILL ME!!!!
    """
    import numpy as np
    try:
        ar = ar.flatten()
    except AttributeError:
        if not return_inverse and not return_index:
            return np.sort(list(set(ar)))
        else:
            ar = np.asanyarray(ar).flatten()

    if ar.size == 0:
        if return_inverse and return_index:
            return ar, np.empty(0, np.bool), np.empty(0, np.bool)
        elif return_inverse or return_index:
            return ar, np.empty(0, np.bool)
        else:
            return ar

    if return_inverse or return_index:
        if return_index:
            perm = ar.argsort(kind='mergesort')
        else:
            perm = ar.argsort()
        aux = ar[perm]
        flag = np.concatenate(([True], aux[1:] != aux[:-1]))
        if return_inverse:
            iflag = np.cumsum(flag) - 1
            iperm = perm.argsort()
            if return_index:
                return aux[flag], perm[flag], iflag[iperm]
            else:
                return aux[flag], iflag[iperm]
        else:
            return aux[flag], perm[flag]

    else:
        ar.sort()
        flag = np.concatenate(([True], ar[1:] != ar[:-1]))
        return ar[flag]

numpy.unique = unique

def in1d(ar1, ar2, assume_unique=False, invert=False):
    """Stolen from numpy, please kill me and upgrade numpy....
    """
    import numpy as np
    # Ravel both arrays, behavior for the first array could be different
    ar1 = np.asarray(ar1).ravel()
    ar2 = np.asarray(ar2).ravel()

    # This code is significantly faster when the condition is satisfied.
    if len(ar2) < 10 * len(ar1) ** 0.145:
        if invert:
            mask = np.ones(len(ar1), dtype=np.bool)
            for a in ar2:
                mask &= (ar1 != a)
        else:
            mask = np.zeros(len(ar1), dtype=np.bool)
            for a in ar2:
                mask |= (ar1 == a)
        return mask

    # Otherwise use sorting
    if not assume_unique:
        ar1, rev_idx = np.unique(ar1, return_inverse=True)
        ar2 = np.unique(ar2)

    ar = np.concatenate( (ar1, ar2) )
    # We need this to be a stable sort, so always use 'mergesort'
    # here. The values from the first array should always come before
    # the values from the second array.
    order = ar.argsort(kind='mergesort')
    sar = ar[order]
    if invert:
        bool_ar = (sar[1:] != sar[:-1])
    else:
        bool_ar = (sar[1:] == sar[:-1])
    flag = np.concatenate( (bool_ar, [invert]) )
    indx = order.argsort(kind='mergesort')[:len( ar1 )]

    if assume_unique:
        return flag[indx]
    else:
        return flag[indx][rev_idx]

numpy.in1d = in1d

def block_diag(*arrs):
    import numpy as np
    if arrs == ():
        arrs = ([],)
    arrs = [np.atleast_2d(a) for a in arrs]

    bad_args = [k for k in range(len(arrs)) if arrs[k].ndim > 2]
    if bad_args:
        raise ValueError("arguments in the following positions have dimension "
                            "greater than 2: %s" % bad_args)

    shapes = np.array([a.shape for a in arrs])
    out = np.zeros(np.sum(shapes, axis=0), dtype=arrs[0].dtype)

    r, c = 0, 0
    for i, (rr, cc) in enumerate(shapes):
        out[r:r + rr, c:c + cc] = arrs[i]
        r += rr
        c += cc
    return out
