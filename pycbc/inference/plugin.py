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

"""Utilities for plugin model discovery."""


def retrieve_model_plugins(model_dict):
    """Retrieves and processes external model plugins.
    """
    import pkg_resources
    # Check for fd waveforms
    for plugin in pkg_resources.iter_entry_points('pycbc.inference.models'):
        add_custom_model(plugin.resolve(), model_dict)


def add_custom_model(model, models, force=False):
    """Makes a custom model available to PyCBC.

    The provided model will be added to the dictionary of models that PyCBC
    knows about, using the model's ``name`` attribute. If the ``name`` is the
    same as a model that already exists in PyCBC, a :py:exc:`RuntimeError` will
    be raised unless the ``force`` option is set to ``True``.

    Parameters
    ----------
    model : pycbc.inference.models.base.BaseModel
        The model to use. The model should be a sub-class of
        :py:class:`BaseModel <pycbc.inference.models.base.BaseModel>` to ensure
        it has the correct API for use within ``pycbc_inference``.
    force : bool, optional
        Add the model even if its ``name`` attribute is the same as a model
        that is already in :py:data:`pycbc.inference.models.models`. Otherwise,
        a :py:exc:`RuntimeError` will be raised. Default is ``False``.
    """
    #from pycbc.inference.models import models
    if model.name in models and not force:
        raise RuntimeError("Cannot load plugin model {}; the name is already "
                           "in use.".format(model.name))
    models[model.name] = model 
