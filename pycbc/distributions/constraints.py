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
    """ Creates a constraint that evaluates to True if parameters obey
    the constraint and False if they do not.
    """
    name = "custom"
    required_parameters = []
    def __init__(self, constraint_arg, transforms=None, **kwargs):
        self.constraint_arg = constraint_arg
        self.transforms = transforms
        for kwarg in kwargs.keys():
            setattr(self, kwarg, kwargs[kwarg])

    def __call__(self, params):
        """ Evaluates constraint.
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

class MtotalLT(Constraint):
    """ Pre-defined constraint that check if total mass is less than a value.
    """
    name = "mtotal_lt"
    required_parameters = ["mass1", "mass2"]

    def _constraint(self, params):
        """ Evaluates constraint function.
        """
        return params["mass1"] + params["mass2"] < self.mtotal

class CartesianSpinSpace(Constraint):
    """ Pre-defined constraint that check if Cartesian parameters
    are within acceptable values.
    """
    name = "cartesian_spin_space"
    required_parameters = ["mass1", "mass2", "spin1x", "spin1y", "spin1z",
                           "spin2x", "spin2y", "spin2z"]

    def _constraint(self, params):
        """ Evaluates constraint function.
        """
        if (params["spin1x"]**2 + params["spin1y"]**2 +
                params["spin1z"]**2)**2 > 1:
            return False
        elif (params["spin2x"]**2 + params["spin2y"]**2 +
                  params["spin2z"]**2)**2 > 1:
            return False
        else:
            return True

class EffectiveSpinSpace(Constraint):
    """ Pre-defined constraint that check if effective spin parameters
    are within acceptable values.
    """
    name = "effective_spin_space"
    required_parameters = ["mass1", "mass2", "q", "xi1", "xi2",
                           "chi_eff", "chi_a"]

    def _constraint(self, params):
        """ Evaluates constraint function.
        """

        # ensure that mass1 > mass2
        if params["mass1"] < params["mass2"]:
            return False

        # constraint for secondary mass
        a = ((4.0 * params["q"]**2 + 3.0 * params["q"])
                 / (4.0 + 3.0 * params["q"]) * params["xi2"])**2
        b = ((1.0 + params["q"]**2) / 4.0
                 * (params["chi_eff"] + params["chi_a"])**2)
        if a + b > 1:
            return False

        # constraint for primary mass
        a = params["xi1"]**2
        b = ((1.0 + params["q"])**2 / (4.0 * params["q"]**2)
                 * (params["chi_eff"] - params["chi_a"])**2)
        if a + b > 1:
            return False

        return True

# list of all constraints
constraints = {
    Constraint.name : Constraint,
    MtotalLT.name : MtotalLT,
    CartesianSpinSpace.name : CartesianSpinSpace,
    EffectiveSpinSpace.name : EffectiveSpinSpace,
}
