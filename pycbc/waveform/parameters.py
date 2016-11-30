# Copyright (C) 2016 Collin Capano
#
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
"""Classes to define common parameters used for waveform generation.
"""

from UserList import UserList

#
# =============================================================================
#
#                                  Base definitions
#
# =============================================================================
#

class Parameter(str):
    """A class that stores information about a parameter. This is done by
    sub-classing string, adding additional attributes.
    """

    def __new__(cls, name, dtype=None, default=None, label=None,
                 description="No description."):
        obj = str.__new__(cls, name)
        obj.name = name
        obj.dtype = dtype
        obj.default = default
        obj.label = label
        obj.description = description
        return obj

    def docstr(self, prefix='', include_label=True):
        """Returns a string summarizing the parameter. Format is:
        <prefix>``name`` : {``default``, ``dtype``}
        <prefix>   ``description`` Label: ``label``.
        """
        outstr = "%s%s : {%s, %s}\n" %(prefix, self.name, str(self.default),
            str(self.dtype).replace("<type '", '').replace("'>", '')) + \
            "%s    %s" %(prefix, self.description)
        if include_label:
            outstr += " Label: %s" %(self.label)
        return outstr


class ParameterList(UserList):
    """A list of parameters. Each element in the list is expected to be a
    Parameter instance.
    """

    @property
    def names(self):
        """Returns a list of the names of each parameter."""
        return [x.name for x in self]

    @property
    def aslist(self):
        """Cast to basic list."""
        return list(self)

    @property
    def asdict(self):
        """Returns a dictionary of the parameters keyed by the parameters."""
        return dict([[x, x] for x in self])

    def defaults(self, include_nulls=False):
        """Returns a list of the name and default value of each parameter,
        as tuples. If include_nulls is False, only parameters that do not
        have None as a default are returnted. Otherwise, all parameters
        are returned.
        """
        return [(x, x.default) for x in self
                    if include_nulls or x.default is not None]

    def default_dict(self, include_nulls=False):
        """Returns a dictionary of the name and default value of each
        parameter.
        """
        return dict(self.defaults(include_nulls=include_nulls))

    @property
    def nodefaults(self):
        """Returns a ParameterList of the parameters that have None for
        defaults.
        """
        return ParameterList([x for x in self if x.default is None])

    @property
    def dtypes(self):
        """Returns a list of the name and dtype of each parameter,
        as tuples.
        """
        return [(x, x.dtype) for x in self]

    @property
    def dtype_dict(self):
        """Returns a dictionary of the name and dtype of each parameter."""
        return dict(self.dtypes)

    @property
    def descriptions(self):
        """Returns a list of the name and description of each parameter,
        as tuples.
        """
        return [(x, x.description) for x in self]

    @property
    def description_dict(self):
        """Return a dictionary of the name and description of each parameter.
        """
        return dict(self.descriptions)

    @property
    def labels(self):
        """Returns a list of each parameter and its label, as tuples."""
        return [(x, x.label) for x in self]

    @property
    def label_dict(self):
        """Return a dictionary of the name and label of each parameter.
        """
        return dict(self.labels)

    def docstr(self, prefix='', include_label=True):
        """Returns the ``docstr`` of each parameter joined together."""
        return '\n'.join([x.docstr(prefix, include_label) for x in self])

#
# =============================================================================
#
#                        Parameter definitions
#
# =============================================================================
#

#
#   CBC intrinsic parameters
#
mass1 = Parameter("mass1",
                dtype=float, default=None, label=r"$m_1~(\mathrm{M}_\odot)$",
                description="The mass of the first component object in the "
                          "binary (in solar masses).")
mass2 = Parameter("mass2",
                dtype=float, default=None, label=r"$m_2~(\mathrm{M}_\odot)$",
                description="The mass of the second component object in the "
                            "binary (in solar masses).")
spin1x = Parameter("spin1x",
                dtype=float, default=0., label=r"$\chi_{1x}$",
                description="The x component of the first binary component's "
                            "dimensionless spin.")
spin1y = Parameter("spin1y",
                dtype=float, default=0., label=r"$\chi_{1y}$",
                description="The y component of the first binary component's "
                            "dimensionless spin.")
spin1z = Parameter("spin1z",
                dtype=float, default=0., label=r"$\chi_{1z}$",
                description="The z component of the first binary component's "
                            "dimensionless spin.")
spin2x = Parameter("spin2x",
                dtype=float, default=0., label=r"$\chi_{2x}$",
                description="The x component of the second binary component's "
                            "dimensionless spin.")
spin2y = Parameter("spin2y",
                dtype=float, default=0., label=r"$\chi_{2y}$",
                description="The y component of the second binary component's "
                            "dimensionless spin.")
spin2z = Parameter("spin2z",
                dtype=float, default=0., label=r"$\chi_{2z}$",
                description="The z component of the second binary component's "
                            "dimensionless spin.")
eccentricity = Parameter("eccentricity",
                dtype=float, default=0., label=r"$e$",
                description="Eccentricity.")

# derived parameters (these are not used for waveform generation) for masses
mchirp = Parameter("mchirp",
                dtype=float, label=r"$\mathcal{M}~(\mathrm{M}_\odot)$",
                description="The chirp mass of the binary (in solar masses).")
eta = Parameter("eta",
                dtype=float, label=r"$\eta$",
                description="The symmetric mass ratio of the binary.")
mtotal = Parameter("mtotal",
                dtype=float, label=r"$M~(\mathrm{M}_\odot)$",
                description="The total mass of the binary (in solar masses).")
q = Parameter("q",
                dtype=float, label=r"$q$",
                description="The mass ratio, m1/m2, where m1 >= m2.")
m_p = Parameter("m_p",
                dtype=float, label=r"$m_{\mathrm{pr}}$",
                description="Mass of the primary object (in solar masses).")
m_s = Parameter("m_s",
                dtype=float, label=r"$m_{\mathrm{sc}}$",
                description="Mass of the secondary object (in solar masses).")

# derived parameters for component spins
chi_eff = Parameter("chi_eff",
                dtype=float, label=r"$\chi_\mathrm{eff}$",
                description="Effective spin of the binary.")
spin_px = Parameter("spin_px",
                dtype=float, label=r"$\chi_{\mathrm{pr}\,x}$",
                description="The x component of the dimensionless spin of the "
                            "primary object.")
spin_py = Parameter("spin_py",
                dtype=float, label=r"$\chi_{\mathrm{pr}\,y}$",
                description="The y component of the dimensionless spin of the "
                            "primary object.")
spin_pz = Parameter("spin_pz",
                dtype=float, label=r"$\chi_{\mathrm{pr}\,z}$",
                description="The z component of the dimensionless spin of the "
                            "primary object.")
spin_sx = Parameter("spin_sx",
                dtype=float, label=r"$\chi_{\mathrm{sc}\,x}$",
                description="The x component of the dimensionless spin of the "
                            "secondary object.")
spin_sy = Parameter("spin_sy",
                dtype=float, label=r"$\chi_{\mathrm{sc}\,y}$",
                description="The y component of the dimensionless spin of the "
                            "secondary object.")
spin_sz = Parameter("spin_sz",
                dtype=float, label=r"$\chi_{\mathrm{sc}\,z}$",
                description="The z component of the dimensionless spin of the "
                            "secondary object.")
lambda1 = Parameter("lambda1",
                dtype=float, default=0., label=r"$\lambda_1$",
                description="The tidal deformability parameter of object 1.")
lambda2 = Parameter("lambda2",
                dtype=float, default=0., label=r"$\lambda_2$",
                description="The tidal deformability parameter of object 2.")
dquad_mon1 = Parameter("dquad_mon1",
                dtype=float, default=0., label=r"$qm_1$",
                description="Quadrupole-monopole parameter / m_1^5 -1.")
dquad_mon2 = Parameter("dquad_mon2",
                dtype=float, default=0., label=r"$qm_2$",
                description="Quadrupole-monopole parameter / m_2^5 -1.")

# derived parameters for component spin magnitude and angles
spin1_a = Parameter("spin1_a",
                    dtype=float, label=r"$\a_{1}$",
                    description="The dimensionless spin magnitude "
                                "$|\vec{s}/m_{1}^2|$.")
spin2_a = Parameter("spin2_a",
                    dtype=float, label=r"$\a_{2}$",
                    description="The dimensionless spin magnitude "
                                "$|\vec{s}/m_{2}^2|$.")
spin1_azimuthal = Parameter(
                      "spin1_azimuthal",
                      dtype=float, label=r"$\theta_1^\mathrm{azimuthal}$",
                      description="The azimuthal spin angle for mass 1.")
spin2_azimuthal = Parameter(
                      "spin2_azimuthal",
                      dtype=float, label=r"$\theta_2^\mathrm{azimuthal}$",
                      description="The azimuthal spin angle for mass 2.")
spin1_polar = Parameter("spin1_polar",
                        dtype=float, label=r"$\theta_1^\mathrm{polar}$",
                        description="The polar spin angle for mass 1.")
spin2_polar = Parameter("spin2_polar",
                        dtype=float, label=r"$\theta_2^\mathrm{polar}$",
                        description="The polar spin angle for mass 2.")


#
#   Parameters needed for CBC waveform generation
#
f_lower = Parameter("f_lower",
                dtype=float, default=None, label=r"$f_0$ (Hz)",
                description="The starting frequency of the waveform (in Hz).")
f_final = Parameter("f_final",
                dtype=float, default=0, label=r"$f_{\mathrm{final}}$ (Hz)",
                description="The ending frequency of the waveform. The "
                            "default (0) indicates that the choice is made by "
                            "the respective approximant.")
f_ref = Parameter("f_ref",
                dtype=float, default=0, label=r"$f_{\mathrm{ref}}$ (Hz)",
                description="The reference frequency.")
delta_f = Parameter("delta_f",
                dtype=float, default=None, label=r"$\Delta f$ (Hz)",
                description="The frequency step used to generate the waveform "
                            "(in Hz).")
delta_t = Parameter("delta_t",
                dtype=float, default=None, label=r"$\Delta t$ (s)",
                description="The time step used to generate the waveform "
                            "(in s).")
sample_points = Parameter("sample_points",
                dtype="Array", default=None, label=None,
                description="An array of the frequencies (in Hz) at which to "
                            "generate the waveform.")
approximant = Parameter("approximant",
                dtype=str, default=None, label=None,
                description="A string that indicates the chosen approximant.")
phase_order = Parameter("phase_order",
                dtype=int, default=-1, label=None,
                description="The pN order of the orbital phase. The default "
                            "of -1 indicates that all implemented orders are "
                            "used.")
spin_order = Parameter("spin_order",
                dtype=int, default=-1, label=None,
                description="The pN order of the spin corrections. The "
                            "default of -1 indicates that all implemented "
                            "orders are used.")
tidal_order = Parameter("tidal_order",
                dtype=int, default=-1, label=None,
                description="The pN order of the tidal corrections. The "
                            "default of -1 indicates that all implemented "
                            "orders are used.")
amplitude_order = Parameter("amplitude_order",
                dtype=int, default=-1, label=None,
                description="The pN order of the amplitude. The default of -1 "
                            "indicates that all implemented orders are used.")
eccentricity_order = Parameter("eccentricity_order",
                dtype=int, default=-1, label=None,
                description="The pN order of the eccentricity corrections."
                        "The default of -1 indicates that all implemented orders are used.")
numrel_data = Parameter("numrel_data",
                dtype=str, default="", label=None,
                description="Sets the NR flags; only needed for NR waveforms.")

#
#   General location parameters
#
distance = Parameter("distance",
                dtype=float, default=1., label=r"$d_L$ (Mpc)",
                description="Luminosity distance to the binary (in Mpc).")
coa_phase = Parameter("coa_phase",
                dtype=float, default=0., label=r"$\phi_c$",
                description="Coalesence phase of the binary (in rad).")
inclination = Parameter("inclination",
                dtype=float, default=0., label=r"$\iota$",
                description="Inclination (rad), defined as the angle between "
                            "the total angular momentum J and the "
                            "line-of-sight.")
long_asc_nodes = Parameter("long_asc_nodes",
                dtype=float, default=0., label=r"$\Omega$",
                description="Longitude of ascending nodes axis (rad).")
mean_per_ano = Parameter("mean_per_ano",
                dtype=float, default=0., label=r"$\delta$",
                description="Mean anomaly of the periastron (rad).")
tc = Parameter("tc",
                dtype=float, default=None, label=r"$t_c$ (s)",
                description="Coalescence time (s).")
ra = Parameter("ra",
                dtype=float, default=None, label=r"$\alpha$",
                description="Right ascension (rad).")
dec = Parameter("dec",
                dtype=float, default=None, label=r"$\delta$",
                description="Declination (rad).")
polarization = Parameter("polarization",
                dtype=float, default=None, label=r"$\psi$",
                description="Polarization (rad).")

#
#   Non mandatory flags with default values
#
frame_axis = Parameter("frame_axis",
                dtype=int, default=0,
                description="Allow to choose among orbital_l, view and total_j")
modes_choice = Parameter("modes_choice",
                dtype=int, default=0,
                description="Allow to turn on  among orbital_l, view and total_j")
side_bands = Parameter("side_bands",
                dtype=int, default=0,
                description="Flag for generating sidebands")

#
# =============================================================================
#
#                        Parameter list definitions
#
# =============================================================================
#

# parameters describing the location of a binary w.r.t. the Earth. Note: we
# do not include distance here. This is because these parameters are not
# passed to the waveform generators in lalsimulation, but are instead applied
# after a waveform is generated. Distance, however, is a parameter used by
# the waveform generators.
location_params = ParameterList([tc, ra, dec, polarization])

# parameters describing the orientation of a binary w.r.t. the radiation
# frame. Note: we include distance here, as it is typically used for generating
# waveforms.
orientation_params = ParameterList([distance, coa_phase, inclination, long_asc_nodes, mean_per_ano])

# the extrinsic parameters of a waveform
extrinsic_params = orientation_params + location_params 

# intrinsic parameters of a CBC waveform
cbc_intrinsic_params = ParameterList([
    mass1, mass2, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z,
    eccentricity, lambda1, lambda2, dquad_mon1, dquad_mon2])

# the parameters of a cbc in the radiation frame
cbc_rframe_params = cbc_intrinsic_params + orientation_params

# common generation parameters are parameters needed to generate either
# a TD, FD, or frequency sequence waveform
common_generation_params = ParameterList([
    approximant, f_ref, phase_order, spin_order, tidal_order, amplitude_order, eccentricity_order])

# Flags having discrete values, optional to generate either
# a TD, FD, or frequency sequence waveform
flags_generation_params = ParameterList([frame_axis, modes_choice, side_bands])

# the following are parameters needed to generate an FD or TD waveform that
# is equally sampled
common_gen_equal_sampled_params = ParameterList([f_lower]) + \
    common_generation_params + flags_generation_params

# the following are parameters needed to generate an FD waveform
fd_waveform_params = cbc_rframe_params + ParameterList([delta_f]) + \
    common_gen_equal_sampled_params + ParameterList([f_final])

# the following are parameters needed to generate a TD waveform
td_waveform_params = cbc_rframe_params + ParameterList([delta_t]) + \
    common_gen_equal_sampled_params + ParameterList([numrel_data]) + \
    flags_generation_params

# the following are parameters needed to generate a frequency series waveform
fd_waveform_sequence_params = cbc_rframe_params + \
    ParameterList([sample_points]) + common_generation_params + \
    flags_generation_params
