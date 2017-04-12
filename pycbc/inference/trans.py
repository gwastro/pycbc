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
        mass1 = conversions.mass1_from_mchirp_q(kwargs[parameters.mchirp],
                                                kwargs[parameters.q])
        mass2 = conversions.mass2_from_mchirp_q(kwargs[parameters.mchirp],
                                                kwargs[parameters.q])
        return {"mass1" : mass1, "mass2" : mass2}

    @staticmethod
    def convert_inverse(**kwargs):
        mchirp = conversions.mchirp_from_mchirp_q(kwargs[parameters.mchirp],
                                                  kwargs[parameters.q])
        m_p = conversions.primary_mass(kwargs[parameters.mass1],
                                       kwargs[parameters.mass2])
        m_s = conversions.secondary_mass(kwargs[parameters.mass1],
                                         kwargs[parameters.mass2])
        q = m_p / m_s
        return {"mchirp" : mchirp, "q" : q}

class SphericalSpin1ToCartesianSpin1(object):
    """ Converts spin1x, spin1y, spin1z, spin2x, spin2y, and spin2z to
    a_1, spin1_azimuthal, spin1_polar, a_2, spin2_azimuthal, and spin2_polar.
    """

    ordered_inputs = [parameters.spin1_a, parameters.spin1_azimuthal,
                  parameters.spin1_polar]
    inputs = set(ordered_inputs)
    outputs = set([parameters.spin1x, parameters.spin1y, parameters.spin1z])

    @classmethod
    def convert(cls, **kwargs):
        out = {}
        a, az, po = cls.ordered_inputs
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

class SphericalSpin2ToCartesianSpin2(object):
    """ Converts spin1x, spin1y, spin1z, spin2x, spin2y, and spin2z to
    a_1, spin1_azimuthal, spin1_polar, a_2, spin2_azimuthal, and spin2_polar.
    """

    ordered_inputs = [parameters.spin2_a, parameters.spin2_azimuthal,
                  parameters.spin2_polar]
    inputs = set(ordered_inputs)
    outputs = set([parameters.spin2x, parameters.spin2y, parameters.spin2z])

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

    # use dict
    params_dict = {param : sampling_params[param]
                   for param in sampling_params.fieldnames}

    # convert mchirp and q sampling to base parameters
    current_params = set(sampling_params.fieldnames)
    if (MchirpQToMass1Mass2.inputs.issubset(current_params) and
            not MchirpQToMass1Mass2.outputs.issubset(current_params)):
        new_fields = MchirpQToMass1Mass2.convert(**params_dict)
        keys = new_fields.keys()
        values = [new_fields[key] for key in new_fields]
        sampling_params = sampling_params.add_fields(values, keys)
        params_dict.update(new_fields)

    # convert spherical spins sampling to base parameters
    current_params = set(sampling_params.fieldnames)
    test_params = set(["spin1_a", "spin1_azimuthal", "spin1_polar"])
    if (SphericalSpin1ToCartesianSpin1.inputs.issubset(current_params) and
            not SphericalSpin1ToCartesianSpin1.outputs.issubset(current_params)):
        new_fields = SphericalSpin1ToCartesianSpin1.convert(**params_dict)
        keys = new_fields.keys()
        values = [new_fields[key] for key in new_fields]
        sampling_params = sampling_params.add_fields(values, keys)
        params_dict.update(new_fields)

    # FIXME: reused
    # convert spherical spins sampling to base parameters
    current_params = set(sampling_params.fieldnames)
    test_params = set(["spin2_a", "spin2_azimuthal", "spin2_polar"])
    if (SphericalSpin2ToCartesianSpin2.inputs.issubset(current_params) and
            not SphericalSpin2ToCartesianSpin2.outputs.issubset(current_params)):
        new_fields = SphericalSpin2ToCartesianSpin2.convert(**params_dict)
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


