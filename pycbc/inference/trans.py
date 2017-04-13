# Copyright (C) 2016  Christopher M. Biwer
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
    """ Base class.
    """
    inputs = set([])
    outputs = set([])

    @classmethod
    def convert(cls, maps):
        raise NotImplementedError("Not added.")

    @classmethod
    def convert_inverse(cls, maps):
        raise NotImplementedError("Not added.")

    @staticmethod
    def format_output(old_maps, new_maps):
        """ Convert inverse.
        """

        # if input is WaveformArray then return WaveformArray
        if isinstance(old_maps, record.WaveformArray):
            keys = maps.keys()
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
        """ Inverse conversion of the class name.
        """
        self.inputs, self.outputs = self.outputs, self.inputs
        self.convert, self.convert_inverse = self.convert_inverse, self.convert

class MchirpQToMass1Mass2(BaseConversion):
    """ Converts mchirp and q to mass1 and mass2.
    """
    inputs = set([parameters.mchirp, parameters.q])
    outputs = set([parameters.mass1, parameters.mass2])

    @classmethod
    def convert(cls, maps):
        out = {}
        out[parameters.mass1] = conversions.mass1_from_mchirp_q(
                                                maps[parameters.mchirp],
                                                maps[parameters.q])
        out[parameters.mass2] = conversions.mass2_from_mchirp_q(
                                                maps[parameters.mchirp],
                                                maps[parameters.q])
        return cls.format_output(maps, out)

    @classmethod
    def convert_inverse(cls, maps):
        out = {}
        out[parameters.mchirp] = \
                 conversions.mchirp_from_mchirp_q(maps[parameters.mchirp],
                                                  maps[parameters.q])
        m_p = conversions.primary_mass(maps[parameters.mass1],
                                       maps[parameters.mass2])
        m_s = conversions.secondary_mass(maps[parameters.mass1],
                                         maps[parameters.mass2])
        out[parmeters.q] = m_p / m_s
        return cls.format_output(maps, out)

class SphericalSpin1ToCartesianSpin1(BaseConversion):
    """ Converts spin1x, spin1y, and spin1z to spin1_a, spin1_azimuthal,
    and spin1_polar.
    """
    ordered_inputs = [parameters.spin1_a, parameters.spin1_azimuthal,
                      parameters.spin1_polar]
    ordered_outputs = [parameters.spin2x, parameters.spin2y, parameters.spin2z]
    inputs = set(ordered_inputs)
    outputs = set(ordered_outputs)

    @classmethod
    def convert(cls, maps):
        a, az, po = cls.ordered_inputs
        data = coordinates.spherical_to_cartesian(maps[a], maps[az], maps[po])
        out = {param : val for param, val in zip(ordered_outputs, data)}
        return cls.format_output(maps, out)

class SphericalSpin2ToCartesianSpin2(SphericalSpin1ToCartesianSpin1):
    """ Converts spin2x, spin2y, and spin2z to spin2_a, spin2_azimuthal,
    and spin2_polar.
    """
    ordered_inputs = [parameters.spin2_a, parameters.spin2_azimuthal,
                      parameters.spin2_polar]
    ordered_outputs = [parameters.spin2x, parameters.spin2y, parameters.spin2z]
    inputs = set(ordered_inputs)
    outputs = set(ordered_outputs)

class MassSpinToCartesianSpin(BaseConversion):
    """ Converts mass1, mass2, chi_eff, chi_a, xi1, xi2, phi_a, and phi_s to
    spin1x, spin1y, spin1z, spin2x, spin2y, and spin2z.
    """
    # mass-spin parameters not in pycbc.waveform.parameters yet
    inputs = set(["mass1", "mass2", "chi_eff", "chi_a", "xi1", "xi2",
                  "phi_a", "phi_s"])
    outputs = set([parameters.spin1x, parameters.spin1y, parameters.spin1z,
                   parameters.spin2x, parameters.spin2y, parameters.spin2z])

    @classmethod
    def convert(cls, maps):
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
        return cls.format_output(maps, out)

class DistanceToRedshift(BaseConversion):
    """ Converts distance to redshift.
    """
    inputs = set([parameters.distance])
    outputs = set([parameters.redshift])

    @classmethod
    def convert(cls, maps):
        out = {parameters.redshift : cosmology.redshift(maps["distance"])}
        return cls.format_output(maps, out)

# list of all Conversions
converts = [MchirpQToMass1Mass2, SphericalSpin1ToCartesianSpin1,
            SphericalSpin2ToCartesianSpin2,
            MassSpinToCartesianSpin]

def add_base_parameters(sampling_params):
    """ Adds a standard set of base parameters to the WaveformArray for
    plotting. Standard set of base parameters includes mass1, mass2,
    spin1x, spin1y, spin1z, spin2x, spin2y, spin2z, and redshift.

    Parameters
    ----------
    sampling_params : WaveformArray
        WaveformArray to add new fields.

    Returns
    -------
    WaveformArray
       WaveformArray with new fields.
    """

    # convert mchirp and q sampling to base parameters
    current_params = set(sampling_params.fieldnames)
    if (MchirpQToMass1Mass2.inputs.issubset(current_params) and
            not MchirpQToMass1Mass2.outputs.issubset(current_params)):
        sampling_params = MchirpQToMass1Mass2.convert(sampling_params)

    # convert spherical spins sampling to base parameters
    for converter in [SphericalSpin1ToCartesianSpin1,
                      SphericalSpin2ToCartesianSpin2]:
        current_params = set(sampling_params.fieldnames)
        if (converter.inputs.issubset(current_params) and
                not converter.outputs.issubset(current_params)):
            sampling_params = converter.convert(sampling_params)

    # convert mass-spin sampling to base parameters
    current_params = set(sampling_params.fieldnames)
    if (MassSpinToCartesianSpin.inputs.issubset(current_params) and
            not MassSpinToCartesianSpin.outputs.issubset(current_params)):
        sampling_params = MassSpinToCartesianSpin.convert(sampling_params)

    # convert distance sampling to base parameters
    current_params = set(sampling_params.fieldnames)
    if (DistanceToRedshift.inputs.issubset(current_params) and
            not DistanceToRedshift.outputs.issubset(current_params)):
        sampling_params = DistanceToRedshift.convert(sampling_params)

    return sampling_params
