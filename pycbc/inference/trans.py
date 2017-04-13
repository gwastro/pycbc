# Copyright (C) 2017  Christopher M. Biwer
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
This modules provides classes and functions for convert sampling parameters
to a set of base parameters.
"""

import numpy
from pycbc import conversions
from pycbc import coordinates
from pycbc import cosmology
from pycbc.io import record
from pycbc.waveform import parameters

class BaseConversion(object):
    """ A base class for converting between two sets of parameters.
    """
    _inputs = set([])
    _outputs = set([])

    @property
    def inputs(self):
        """ Returns a set of input parameters.
        """
        return set(self._inputs)

    @property
    def outputs(self):
        """ Returns a set of output parameters.
        """
        return set(self._outputs)

    @classmethod
    def _convert(cls, maps):
        """ This function converts from inputs to outputs.
        """
        raise NotImplementedError("Not added.")

    @classmethod
    def _convert_inverse(cls, maps):
        """ The inverse conversions of _convert. This function converts from
        outputs to inputs.
        """
        raise NotImplementedError("Not added.")

    @classmethod
    def convert(cls, old_maps):
        """ Convert inputs to outputs. This function accepts either
        a WaveformArray or dict. It will return the same output type as the
        input mapping object. Internally it calls _convert.

        Parameters
        ----------
        old_maps : {WaveformArray, dict}
            Mapping object to add new keys.

        Returns
        -------
        {WaveformArray, dict}
            Mapping object with new keys.
        """
        new_maps = cls._convert(old_maps)
        return cls.format_output(old_maps, new_maps)

    @staticmethod
    def format_output(old_maps, new_maps):
        """ This function takes the returned dict from _convert and converts
        it to the same datatype as the input to _convert.

        Parameters
        ----------
        old_maps : {WaveformArray, dict}
            The mapping object to add new maps to.
        new_maps : dict
            A dict with key as parameter name and value is numpy.array.

        Returns
        -------
        {WaveformArray, dict}
            The old_maps object with new keys from new_maps.
        """

        # if input is WaveformArray then return WaveformArray
        if isinstance(old_maps, record.WaveformArray):
            keys = new_maps.keys()
            values = [new_maps[key] for key in keys]
            old_maps = old_maps.add_fields(values, keys)
            return old_maps

        # if input is dict then return dict
        elif isinstance(old_maps, dict):
            return old_maps.update(new_maps)

        # else error
        else:
            raise TypeError("Input type must be WaveformArray or dict.")

    def inverse(self):
        """ Inverts the conversions being done. Inputs become outputs and
        vice versa. The function convert will now call the inverse
        transformation.
        """
        self._inputs, self._outputs = self._outputs, self._inputs
        self._convert, self._convert_inverse = \
                                      self._convert_inverse, self._convert

class MchirpQToMass1Mass2(BaseConversion):
    """ Converts chirp mass and mass ratio to component masses.
    """
    _inputs = [parameters.mchirp, parameters.q]
    _outputs = [parameters.mass1, parameters.mass2]

    @classmethod
    def _convert(cls, maps):
        out = {}
        out[parameters.mass1] = conversions.mass1_from_mchirp_q(
                                                maps[parameters.mchirp],
                                                maps[parameters.q])
        out[parameters.mass2] = conversions.mass2_from_mchirp_q(
                                                maps[parameters.mchirp],
                                                maps[parameters.q])
        return out

    @classmethod
    def _convert_inverse(cls, maps):
        out = {}
        out[parameters.mchirp] = \
                 conversions.mchirp_from_mchirp_q(maps[parameters.mchirp],
                                                  maps[parameters.q])
        m_p = conversions.primary_mass(maps[parameters.mass1],
                                       maps[parameters.mass2])
        m_s = conversions.secondary_mass(maps[parameters.mass1],
                                         maps[parameters.mass2])
        out[parmeters.q] = m_p / m_s
        return out

class SphericalSpin1ToCartesianSpin1(BaseConversion):
    """ Converts spherical spin parameters (magnitude and two angles) to
    catesian spin parameters. This class only converts spsins for the first
    component mass.
    """
    _inputs = [parameters.spin1_a, parameters.spin1_azimuthal,
               parameters.spin1_polar]
    _outputs = [parameters.spin2x, parameters.spin2y, parameters.spin2z]

    @classmethod
    def _convert(cls, maps):
        a, az, po = cls.ordered_inputs
        data = coordinates.spherical_to_cartesian(maps[a], maps[az], maps[po])
        out = {param : val for param, val in zip(ordered_outputs, data)}
        return out

class SphericalSpin2ToCartesianSpin2(SphericalSpin1ToCartesianSpin1):
    """ Converts spherical spin parameters (magnitude and two angles) to
    catesian spin parameters. This class only converts spsins for the second
    component mass.
    """
    _inputs = [parameters.spin2_a, parameters.spin2_azimuthal,
               parameters.spin2_polar]
    _outputs = [parameters.spin2x, parameters.spin2y, parameters.spin2z]

class MassSpinToCartesianSpin(BaseConversion):
    """ Converts component masses, and mass-weighted spin parameters
    (eg. effective spin) to cartesian spin coordinates.
    """
    # mass-spin parameters not in pycbc.waveform.parameters yet
    _inputs = ["mass1", "mass2", "chi_eff", "chi_a", "xi1", "xi2",
               "phi_a", "phi_s"]
    _outputs = [parameters.spin1x, parameters.spin1y, parameters.spin1z,
                parameters.spin2x, parameters.spin2y, parameters.spin2z]

    @classmethod
    def _convert(cls, maps):
        out = {}
        out[parameters.spin1x] = conversions.spin1x_from_xi1_phi_a_phi_s(
                               maps["xi1"], maps["phi_a"], maps["phi_s"])
        out[parameters.spin1y] = conversions.spin1y_from_xi1_phi_a_phi_s(
                               maps["xi1"], maps["phi_a"], maps["phi_s"])
        out[parameters.spin1z] = \
                         conversions.spin1z_from_mass1_mass2_chi_eff_chi_a(
                               maps[parameters.mass1], maps[parameters.mass2],
                               maps["chi_eff"], maps["chi_a"])
        out[parameters.spin2x] = \
                         conversions.spin2x_from_mass1_mass2_xi2_phi_a_phi_s(
                               maps[parameters.mass1], maps[parameters.mass2],
                               maps["xi2"], maps["phi_a"], maps["phi_s"])
        out[parameters.spin2y] = \
                         conversions.spin2y_from_mass1_mass2_xi2_phi_a_phi_s(
                               maps[parameters.mass1], maps[parameters.mass2],
                               maps["xi2"], maps["phi_a"], maps["phi_s"])
        out[parameters.spin2z] = \
                         conversions.spin2z_from_mass1_mass2_chi_eff_chi_a(
                               maps[parameters.mass1], maps[parameters.mass2],
                               maps["chi_eff"], maps["chi_a"])
        return out

class DistanceToRedshift(BaseConversion):
    """ Converts distance to redshift.
    """
    _inputs = [parameters.distance]
    _outputs = [parameters.redshift]

    @classmethod
    def _convert(cls, maps):
        out = {parameters.redshift : cosmology.redshift(maps["distance"])}
        return out

# list of all Conversions
converts = [MchirpQToMass1Mass2, SphericalSpin1ToCartesianSpin1,
            SphericalSpin2ToCartesianSpin2, MassSpinToCartesianSpin,
            DistanceToRedshift]

def add_base_parameters(sampling_params):
    """ Adds a standard set of base parameters to a mapping object
    (ie. WaveformArray or dict). The standard set of base parameters includes
    mass1, mass2, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z, and redshift.

    This function loop over all Conversion classes in this module and sees if
    the conversions is required. If it is required, then the new keys are
    added.

    Parameters
    ----------
    sampling_params : {WaveformArray, dict}
        Mapping object to add new keys.

    Returns
    -------
    sampling_params : {WaveformArray, dict}
       Mapping object with new fields.
    """
    converters = [converter() for converter in converts]
    for converter in converters:
        current_params = set(sampling_params.fieldnames)
        if (converter.inputs.issubset(current_params) and
                not converter.outputs.issubset(current_params)):
            sampling_params = converter.convert(sampling_params)
    return sampling_params
