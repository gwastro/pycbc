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

from pycbc import transforms
from pycbc.io import record

class Constraint(object):
    """Creates a constraint that evaluates to True if parameters obey
    the constraint and False if they do not.
    """
    name = "custom"

    def __init__(self, constraint_arg, transforms=None, **kwargs):
        self.constraint_arg = constraint_arg
        self.transforms = transforms
        self._scratch = None
        for kwarg in kwargs.keys():
            setattr(self, kwarg, kwargs[kwarg])

    def _getscratch(self, params):
        """Gets FieldArray scratch space.

        FieldArray creation can be slow. This tries to speed up successive
        calls by creating a scratch FieldArray space based on the given
        parameters.
        """
        if self._scratch is None or (set(params.keys()) != 
                                     set(self._scratch.fieldnames)):
            self._scratch = record.FieldArray.from_kwargs(**params)
        else:
            for p in params:
                try:
                    self._scratch[p][:] = params[p]
                except ValueError:
                    # can happen if there is a size mismatch, just
                    # create a new one
                    self._scratch = record.FieldArray.from_kwargs(**params)
                    break
        return self._scratch

    def __call__(self, params):
        """Evaluates constraint.
        """
        # cast to FieldArray
        if isinstance(params, dict):
            params = self._getscratch(params)
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


# list of all constraints
constraints = {
    Constraint.name : Constraint,
}
