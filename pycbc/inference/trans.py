""" Convert sampling parameters to a set of base parameters.
"""

import numpy
from pycbc import conversions
from pycbc import coordinates
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
    def format_output(inputs, maps):
        """ Convert inverse.
        """

        # if input is WaveformArray then return WaveformArray
        if isinstance(inputs, record.WaveformArray):
            keys = maps.keys()
            values = [maps[key] for key in maps]
            inputs = inputs.add_fields(values, keys)
            return inputs

        # if input is dict then return dict
        elif isinstance(inputs, dict):
            return inputs.update(maps)

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
        mass1 = conversions.mass1_from_mchirp_q(maps[parameters.mchirp],
                                                maps[parameters.q])
        mass2 = conversions.mass2_from_mchirp_q(maps[parameters.mchirp],
                                                maps[parameters.q])
        return cls.format_output(
                    maps, {parameters.mass1 : mass1, parameters.mass2 : mass2})

    @classmethod
    def convert_inverse(cls, maps):
        mchirp = conversions.mchirp_from_mchirp_q(maps[parameters.mchirp],
                                                  maps[parameters.q])
        m_p = conversions.primary_mass(maps[parameters.mass1],
                                       maps[parameters.mass2])
        m_s = conversions.secondary_mass(maps[parameters.mass1],
                                         maps[parameters.mass2])
        q = m_p / m_s
        return cls.format_output(
                    maps, {parameters.mchirp : mchirp, parameters.q : q})

class SphericalSpin1ToCartesianSpin1(BaseConversion):
    """ Converts spin1x, spin1y, and spin1z to spin1_a, spin1_azimuthal,
    and spin1_polar.
    """
    ordered_inputs = [parameters.spin1_a, parameters.spin1_azimuthal,
                      parameters.spin1_polar]
    inputs = set(ordered_inputs)
    outputs = set([parameters.spin1x, parameters.spin1y, parameters.spin1z])

    @classmethod
    def convert(cls, maps):
        out = {}
        a, az, po = cls.ordered_inputs
        a_val, az_val, po_val = coordinates.spherical_to_cartesian(
                                                   maps[a], maps[az], maps[po])
        return cls.format_output(maps, {a : a_val, az : az_val, po : po_val})

class SphericalSpin2ToCartesianSpin2(BaseConversion):
    """ Converts spin1x, spin1y, and spin1z to spin1_a, spin1_azimuthal,
    and spin1_polar.
    """
    ordered_inputs = [parameters.spin2_a, parameters.spin2_azimuthal,
                  parameters.spin2_polar]
    inputs = set(ordered_inputs)
    outputs = set([parameters.spin2x, parameters.spin2y, parameters.spin2z])

class MassSpinToCartesianSpin(BaseConversion):
    """ Converts mass1, mass2, chi_eff, chi_a, xi1, xi2, phi_a, and phi_s to
    spin1x, spin1y, spin1z, spin2x, spin2y, and spin2z.
    """
    inputs = set([])
    outputs = set([])

    @classmethod
    def convert(cls, maps):
        spin1x = conversions.spin1x_from_xi1_phi_a_phi_s(
                               maps["xi1"], maps["phi_a"], maps["phi_s"])
        spin1y = conversions.spin1y_from_xi1_phi_a_phi_s(
                               maps["xi1"], maps["phi_a"], maps["phi_s"])
        spin1z = conversions.spin1z_from_mass1_mass2_chi_eff_chi_a(
                               maps["mass1"], maps["mass2"],
                               maps["chi_eff"], maps["chi_a"])
        spin2x = conversions.spin2x_from_mass1_mass2_xi2_phi_a_phi_s(
                               maps["mass1"], maps["mass2"], maps["xi2"],
                               maps["phi_a"], maps["phi_s"])
        spin2y = conversions.spin2y_from_mass1_mass2_xi2_phi_a_phi_s(
                               maps["mass1"], maps["mass2"], maps["xi2"],
                               maps["phi_a"], maps["phi_s"])
        spin2z = conversions.spin2z_from_mass1_mass2_chi_eff_chi_a(
                               maps["mass1"], maps["mass2"],
                               maps["chi_eff"], maps["chi_a"])
        return cls.format_output(
                    maps, {"spin1x" : spin1x, "spin1y" : spin1y,
                           "spin1z" : spin1z, "spin2x" : spin2x,
                           "spin2y" : spin2y, "spin2z" : spin2z})

# list of all Conversions
converts = [MchirpQToMass1Mass2, SphericalSpin1ToCartesianSpin1,
            SphericalSpin2ToCartesianSpin2,
            MassSpinToCartesianSpin]

def add_base_parameters(sampling_params):
    """ Adds a standard set of base parameters to the WaveformArray for
    plotting. Standard set of base parameters includes mass1, mass2,
    spin1x, spin1y, spin1z, spin2x, spin2y, and spin2z.

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

    return sampling_params
