# Copyright (C) 2013 Ian W. Harry
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

from __future__ import division
import copy
import numpy
import lal

def generate_hexagonal_lattice(maxv1, minv1, maxv2, minv2, mindist):
    """
    This function generates a 2-dimensional lattice of points using a hexagonal
    lattice.

    Parameters
    -----------
    maxv1 : float
        Largest value in the 1st dimension to cover
    minv1 : float
        Smallest value in the 1st dimension to cover
    maxv2 : float
        Largest value in the 2nd dimension to cover
    minv2 : float
        Smallest value in the 2nd dimension to cover
    mindist : float
        Maximum allowed mismatch between a point in the parameter space and the
        generated bank of points.

    Returns
    --------
    v1s : numpy.array
        Array of positions in the first dimension
    v2s : numpy.array
        Array of positions in the second dimension
    """
    if minv1 > maxv1:
      raise ValueError("Invalid input to function.")
    if minv2 > maxv2:
      raise ValueError("Invalid input to function.")
    # Place first point
    v1s = [minv1]
    v2s = [minv2]
    initPoint = [minv1,minv2]
    # Place first line
    initLine = [initPoint]
    tmpv1 = minv1
    while (tmpv1 < maxv1):
        tmpv1 = tmpv1 + (3 * mindist)**(0.5)
        initLine.append([tmpv1,minv2])
        v1s.append(tmpv1)
        v2s.append(minv2)
    initLine = numpy.array(initLine)
    initLine2 = copy.deepcopy(initLine)
    initLine2[:,0] += 0.5 * (3*mindist)**0.5
    initLine2[:,1] += 1.5 * (mindist)**0.5
    for i in xrange(len(initLine2)):
        v1s.append(initLine2[i,0])
        v2s.append(initLine2[i,1])
    tmpv2_1 = initLine[0,1]
    tmpv2_2 = initLine2[0,1]
    while tmpv2_1 < maxv2 and tmpv2_2 < maxv2:
        tmpv2_1 = tmpv2_1 + 3.0 * (mindist)**0.5
        tmpv2_2 = tmpv2_2 + 3.0 * (mindist)**0.5 
        initLine[:,1] = tmpv2_1
        initLine2[:,1] = tmpv2_2
        for i in xrange(len(initLine)):
            v1s.append(initLine[i,0])
            v2s.append(initLine[i,1])
        for i in xrange(len(initLine2)):
            v1s.append(initLine2[i,0])
            v2s.append(initLine2[i,1])
    v1s = numpy.array(v1s)
    v2s = numpy.array(v2s)
    return v1s, v2s

def generate_anstar_3d_lattice(maxv1, minv1, maxv2, minv2, maxv3, minv3, \
                               mindist):
    """
    This function calls into LAL routines to generate a 3-dimensional array
    of points using the An^* lattice.

    Parameters
    -----------
    maxv1 : float
        Largest value in the 1st dimension to cover
    minv1 : float
        Smallest value in the 1st dimension to cover
    maxv2 : float
        Largest value in the 2nd dimension to cover
    minv2 : float
        Smallest value in the 2nd dimension to cover
    maxv3 : float
        Largest value in the 3rd dimension to cover
    minv3 : float
        Smallest value in the 3rd dimension to cover
    mindist : float
        Maximum allowed mismatch between a point in the parameter space and the
        generated bank of points.

    Returns
    --------
    v1s : numpy.array
        Array of positions in the first dimension
    v2s : numpy.array
        Array of positions in the second dimension
    v3s : numpy.array
        Array of positions in the second dimension
    """
    # Lalpulsar not a requirement for the rest of pycbc, so check if we have it
    # here in this function.
    try:
        import lalpulsar
    except:
        raise ImportError("A SWIG-wrapped install of lalpulsar is needed to use the anstar tiling functionality.")

    tiling = lalpulsar.CreateLatticeTiling(3)
    lalpulsar.SetLatticeTilingConstantBound(tiling, 0, minv1, maxv1)
    lalpulsar.SetLatticeTilingConstantBound(tiling, 1, minv2, maxv2)
    lalpulsar.SetLatticeTilingConstantBound(tiling, 2, minv3, maxv3)
    # Make a 3x3 Euclidean lattice
    a = lal.gsl_matrix(3,3)
    a.data[0,0] = 1
    a.data[1,1] = 1
    a.data[2,2] = 1
    try:
        # old versions of lalpulsar used an enumeration
        lattice = lalpulsar.TILING_LATTICE_ANSTAR
    except AttributeError:
        # newer versions of lalpulsar use a string
        lattice = 'An-star'
    lalpulsar.SetTilingLatticeAndMetric(tiling, lattice, a, mindist)
    try:
        iterator = lalpulsar.CreateLatticeTilingIterator(tiling, 3)
    except TypeError:
        # old versions of lalpulsar required the flags argument
        # (set to 0 for defaults)
        iterator = lalpulsar.CreateLatticeTilingIterator(tiling, 3, 0)

    vs1 = []
    vs2 = []
    vs3 = []
    count = 0
    curr_point = lal.gsl_vector(3)
    while (lalpulsar.NextLatticeTilingPoint(iterator, curr_point) > 0):
        vs1.append(curr_point.data[0])
        vs2.append(curr_point.data[1])
        vs3.append(curr_point.data[2])
    return vs1, vs2, vs3

