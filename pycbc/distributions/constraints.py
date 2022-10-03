# Copyright (C) 2017 Christopher M. Biwer
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
This modules provides classes for evaluating multi-dimensional constraints.
"""

import re
import scipy.spatial
import numpy
import h5py
from pycbc import transforms
from pycbc.io import record


class Constraint(object):
    """Creates a constraint that evaluates to True if parameters obey
    the constraint and False if they do not.
    """
    name = "custom"

    def __init__(self, constraint_arg, static_args=None, transforms=None,
                 **kwargs):
        static_args = (
            {} if static_args is None
            else dict(sorted(
                static_args.items(), key=lambda x: len(x[0]), reverse=True))
            )
        for arg, val in static_args.items():
            swp = f"'{val}'" if isinstance(val, str) else str(val)
            # Substitute static arg name for value if it appears in the
            # constraint_arg string at the beginning of a word and is not
            # followed by an underscore or equals sign.
            # This ensures that static_args that are also kwargs in function calls are
            # handled correctly, i.e., the kwarg is not touched while its value is replaced
            # with the static_arg value.
            constraint_arg = re.sub(
                r'\b{}(?!\_|\=)'.format(arg), swp, constraint_arg)
        self.constraint_arg = constraint_arg
        self.transforms = transforms
        for kwarg in kwargs.keys():
            setattr(self, kwarg, kwargs[kwarg])

    def __call__(self, params):
        """Evaluates constraint.
        """
        # cast to FieldArray
        if isinstance(params, dict):
            params = record.FieldArray.from_kwargs(**params)
        elif not isinstance(params, record.FieldArray):
           raise ValueError("params must be dict or FieldArray instance")

        # try to evaluate; this will assume that all of the needed parameters
        # for the constraint exists in params
        try:
            out = self._constraint(params)
        except NameError:
            # one or more needed parameters don't exist; try applying the
            # transforms
            params = transforms.apply_transforms(params, self.transforms) \
                     if self.transforms else params
            out = self._constraint(params)
        if isinstance(out, record.FieldArray):
            out = out.item() if params.size == 1 else out
        return out

    def _constraint(self, params):
        """ Evaluates constraint function.
        """
        return params[self.constraint_arg]


class SupernovaeConvexHull(Constraint):
    """Pre defined constraint for core-collapse waveforms that checks
    whether a given set of coefficients lie within the convex hull of
    the coefficients of the principal component basis vectors.
    """
    name = "supernovae_convex_hull"
    required_parameters = ["coeff_0", "coeff_1"]

    def __init__(self, constraint_arg, transforms=None, **kwargs):
        super(SupernovaeConvexHull,
              self).__init__(constraint_arg, transforms=transforms, **kwargs)

        if 'principal_components_file' in kwargs:
            pc_filename = kwargs['principal_components_file']
            hull_dimention = numpy.array(kwargs['hull_dimention'])
            self.hull_dimention = int(hull_dimention)
            pc_file = h5py.File(pc_filename, 'r')
            pc_coefficients = numpy.array(pc_file.get('coefficients'))
            pc_file.close()
            hull_points = []
            for dim in range(self.hull_dimention):
                hull_points.append(pc_coefficients[:, dim])
            hull_points = numpy.array(hull_points).T
            pc_coeffs_hull = scipy.spatial.Delaunay(hull_points)
            self._hull = pc_coeffs_hull

    def _constraint(self, params):

        output_array = []
        points = numpy.array([params["coeff_0"],
                              params["coeff_1"],
                              params["coeff_2"]])
        for coeff_index in range(len(params["coeff_0"])):
            point = points[:, coeff_index][:self.hull_dimention]
            output_array.append(self._hull.find_simplex(point) >= 0)
        return numpy.array(output_array)


# list of all constraints
constraints = {
    Constraint.name : Constraint,
    SupernovaeConvexHull.name : SupernovaeConvexHull,
}
