# Copyright (C) 2022  Collin Capano
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

"""Hierarchical model definitions."""

from .base import BaseModel



#
# =============================================================================
#
#                       Hierarhical model definition
#
# =============================================================================
#


class HierarhcicalModel(BaseModel):
    r"""Model that is a combination of other models.

    Sub-models are treated as being independent of each other, although
    they can share parameters. In other words, the hiearchical likelihood is:

    .. math::

        p(\bold{D}|\bold{\vartheta}, \bold{H}) =
            \prod_{I}^{K} p(\bold{d}_I|\bold{\vartheta}, H_{I})

    Submodels are provided as a dictionary upon initialization with a unique
    label assigned to each model, e.g.,
    `{'event1' -> model1, 'event2' -> model2}`. Variable and static parameters
    that are specific to each submodel should be prepended with `{label}_`,
    where `{label}__` is the label associated with the given submodel. Shared
    parameters across multiple models have no labels prepended. Parameters with
    labels prepended will override parameters without labels for the event
    they are assigned.

    Parameters
    ----------
    variable_params: (tuple of) string(s)
        A tuple of parameter names that will be varied.
    submodels: dict
        Dictionary of model labels -> model instances of all the submodels.
    """
    name = 'hierarchical'

    def __init__(self, variable_params, submodels, **kwargs):
        self.submodels = submodels
        super().__init__(variable_params, **kwargs)
        # store a map of model names -> parameters for quick look up later
        self.param_map = map_params(self.variable_params)

    @BaseModel.variable_params.setter
    def variable_params(self, variable_params):
        # overrides BaseModel's variable params to use HierarchicalParam
        # instances for the variable parameters
        if isinstance(variable_params, str):
            variable_params = [variable_params]
        self._variable_params = tuple(HierarchicalParam(p, self.submodels)
                                      for p in variable_params)

    @BaseModel.static_params.setter
    def static_params(self, static_params):
        if static_params is None:
            static_params = {}
        self._static_params = {HierarchicalParam(p, self.submodels): val
                               for p, val in static_params.items()}

    def update(self, **params):
        """Updates the current parameter positions, resets stats, and updates
        the sub-models.

        Model labels are stripped from parameters before being passed to
        the submodels.
        """
        # update the hierarchical parameters
        super().update(**params)
        # update each of the submodels
        for lbl, model in self.submodels.items():
            model.update(**{p.subname: params[p.fullname]
                            for p in self.param_map[lbl]})

    def _loglikelihood(self):
        # takes the sum of the constitutent models' loglikelihoods
        logl = 0.
        for lbl, model in self.submodels.items():
            logl += model.loglikelihood
        return logl


class HierarchicalParam(str):
    """Sub-class of str for hierarchical parameter names.

    This adds attributes that keep track of the model label(s) the parameter
    is associated with, along with the name that is passed to the models.
    """
    delim = '__'
    model_delim = '_'

    def __new__(cls, fullname, possible_models):
        fullname = str(fullname)
        obj = str.__new__(cls, fullname)
        obj.fullname = fullname
        models, subp = HierarchicalParam.parse(fullname, possible_models)
        obj.models = models
        obj.subname = subp
        return obj

    @classmethod
    def parse(cls, fullname, possible_models):
        """Parses the full parameter name into the models the parameter is
        associated with and the parameter name that is passed to the models.

        Parameters
        ----------
        fullname : str
            The full name of the parameter, which includes both the model
            label(s) and the parameter name.
        possible_models : set
            Set of model labels the parameter can be associated with.

        Returns
        -------
        models : list
            List of the model labels the parameter is associated with.
        subp : str
            Parameter name that is passed to the models. This is the parameter
            name with the model label(s) stripped from it.
        """
        # make sure possible models is a set
        possible_models = set(possible_models)
        p = fullname.split(cls.delim, 1)
        if len(p) == 1:
            # is a global fullname, associate with all
            subp = fullname
            models = possible_models.copy()
        else:
            models, subp = p
            # convert into set of model label(s)
            models = set(models.split(cls.model_delim))
            # make sure the given labels are in the list of possible models
            unknown = models - possible_models
            if any(unknown):
                raise ValueError('unrecognized model label(s) {} present in '
                                 'parameter {}'.format(', '.join(unknown),
                                                       fullname))
        return models, subp


def map_params(params):
    """Creates a map of models -> parameters.

    Parameters
    ----------
    params : list of HierarchicalParam instances
        The list of hierarchical parameter names to prase.

    Returns
    -------
    dict :
        Dictionary of model labels -> associated parameters.
    """
    param_map = {}
    for p in params:
        for lbl in p.models:
            try:
                param_map[lbl].append(p)
            except KeyError:
                param_map[lbl] = [p]
    return param_map
