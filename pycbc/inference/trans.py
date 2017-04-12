""" Convert sampling parameters to a set of base parameters used for plotting.
"""

import numpy
from pycbc import conversions
from pycbc import coordinates
from pycbc.waveform import parameters

class BaseConversion(object):

    inputs = set([])
    outputs = set([])

    def convert(self, **kwargs):
        raise NotImplementedError("Not added.")

    def convert_inverse(self, **kwargs):
        raise NotImplementedError("Not added.")

    def inverse(self):
        self.inputs, self.outputs = self.outputs, self.inputs
        self.convert, self.convert_inverse = self.convert_inverse, self.convert

class MchirpQToMass1Mass2(object):
    """ Converts mchirp and q to mass1 and mass2.
    """

    inputs = set([parameters.mchirp, parameters.q])
    outputs = set([parameters.mass1, parameters.mass2])

    @staticmethod
    def convert(**kwargs):
        mass1 = conversions.mass1_from_mchirp_q(kwargs["mchirp"], kwargs["q"])
        mass2 = conversions.mass2_from_mchirp_q(kwargs["mchirp"], kwargs["q"])
        return {"mass1" : mass1, "mass2" : mass2}

    @staticmethod
    def convert_inverse(**kwargs):
        mchirp = conversions.mchirp_from_mchirp_q(kwargs["mchirp"], kwargs["q"])
        m_p = conversions.primary_mass(kwargs["mass1"], kwargs["mass2"])
        m_s = conversions.secondary_mass(kwargs["mass1"], kwargs["mass2"])
        q = m_p / m_s
        return {"mchirp" : mchirp, "q" : q}

class SphericalSpinToCartesianSpin(object):
    """ Converts spin1x, spin1y, spin1z, spin2x, spin2y, and spin2z to
    a_1, spin1_azimuthal, spin1_polar, a_2, spin2_azimuthal, and spin2_polar.
    """

    inputs = set([])
    outputs = set([])

    @staticmethod
    def convert(**kwargs):
        out = {}
        parameters = [["a_1", "spin1_azimuthal", "spin1_polar"],
                      ["a_2", "spin2_azimuthal", "spin2_polar"]]
        for a, az, po in parameters:
            if set([a, az, po]).issubset(set(kwargs.keys())):
                a_val, az_val, po_val = \
                     coordinates.spherical_to_cartesian(kwargs[a],
                                                        kwargs[az],
                                                        kwargs[po])
                out.update({a : a_val, az : az_val, po : po_val})
        return out

    @staticmethod
    def convert_inverse(**kwargs):
        raise NotImplementedError("Not added.")

class MassSpinToCartesianSpin(object):
    """ Converts mass1, mass2, chi_eff, chi_a, xi1, xi2, phi_a, and phi_s to
    spin1x, spin1y, spin1z, spin2x, spin2y, and spin2z.
    """

    inputs = set([])
    outputs = set([])

    @staticmethod
    def convert(**kwargs):
        spin1x = conversions.spin1x_from_xi1_phi_a_phi_s(
                               kwargs["xi1"], kwargs["phi_a"], kwargs["phi_s"])
        spin1y = conversions.spin1y_from_xi1_phi_a_phi_s(
                               kwargs["xi1"], kwargs["phi_a"], kwargs["phi_s"])
        spin1z = conversions.spin1z_from_mass1_mass2_chi_eff_chi_a(
                               kwargs["mass1"], kwargs["mass2"],
                               kwargs["chi_eff"], kwargs["chi_a"])
        spin2x = conversions.spin2x_from_mass1_mass2_xi2_phi_a_phi_s(
                               kwargs["mass1"], kwargs["mass2"], kwargs["xi2"],
                               kwargs["phi_a"], kwargs["phi_s"])
        spin2y = conversions.spin2y_from_mass1_mass2_xi2_phi_a_phi_s(
                               kwargs["mass1"], kwargs["mass2"], kwargs["xi2"],
                               kwargs["phi_a"], kwargs["phi_s"])
        spin2z = conversions.spin2z_from_mass1_mass2_chi_eff_chi_a(
                               kwargs["mass1"], kwargs["mass2"],
                               kwargs["chi_eff"], kwargs["chi_a"])
        return {"spin1x" : spin1x, "spin1y" : spin1y, "spin1z" : spin1z,
                "spin2x" : spin2x, "spin2y" : spin2y, "spin2z" : spin2z}

    @staticmethod
    def convert_inverse(**kwargs):
        raise NotImplementedError("Not added.")

converts = [MchirpQToMass1Mass2, SphericalSpinToCartesianSpin,
            MassSpinToCartesianSpin]

def convert(sampling_params):
    """ Converts from sampling parameters to base parameters for plotting.
    """

    # use dict
    params_dict = {param : sampling_params[param]
                   for param in sampling_params.fieldnames}

    # convert mchirp and q sampling to base parameters
    current_params = set(sampling_params.fieldnames)
    test_params = set(["mchirp", "q"])
    if test_params.issubset(current_params):
        new_fields = MchirpQToMass1Mass2.convert(**params_dict)
        keys = new_fields.keys()
        values = [new_fields[key] for key in new_fields]
        sampling_params = sampling_params.add_fields(values, keys)
        params_dict.update(new_fields)

    # convert spherical spins sampling to base parameters
    current_params = set(sampling_params.fieldnames)
    test_params1 = set(["a_1", "spin1_azimuthal", "spin1_polar"])
    test_params2 = set(["a_2", "spin2_azimuthal", "spin2_polar"])
    if (test_params1.issubset(current_params) or
            test_params2.issubset(current_params)):
        new_fields = SphericalSpinToCartesianSpin.convert(**params_dict)
        keys = new_fields.keys()
        values = [new_fields[key] for key in new_fields]
        sampling_params = sampling_params.add_fields(values, keys)
        params_dict.update(new_fields)

    # convert mass-spin sampling to base parameters
    current_params = set(sampling_params.fieldnames)
    test_params = set(["chi_eff", "chi_a", "xi1", "xi2", "phi_s", "phi_a"])
    if test_params.issubset(current_params):
        new_fields = MassSpinToCartesianSpin.convert(**params_dict)
        keys = new_fields.keys()
        values = [new_fields[key] for key in new_fields]
        sampling_params = sampling_params.add_fields(values, keys)
        params_dict.update(new_fields)

    return sampling_params


