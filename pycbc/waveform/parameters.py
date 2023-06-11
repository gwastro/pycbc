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

from collections import OrderedDict
try:
    from collections import UserList
except ImportError:
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
        dtype_str = str(self.dtype).replace("<type '", '').replace("'>", '')
        dtype_str = dtype_str.replace("<class '", '')
        outstr = "%s%s : {%s, %s}\n%s    %s" % (
                prefix, self.name, str(self.default), dtype_str, prefix,
                self.description)
        if include_label:
            outstr += " Label: %s" % (self.label)
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

    def defaults(self):
        """Returns a list of the name and default value of each parameter,
        as tuples.
        """
        return [(x, x.default) for x in self]

    def default_dict(self):
        """Returns a dictionary of the name and default value of each
        parameter.
        """
        return OrderedDict(self.defaults())

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
        return OrderedDict(self.dtypes)

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
        return OrderedDict(self.descriptions)

    @property
    def labels(self):
        """Returns a list of each parameter and its label, as tuples."""
        return [(x, x.label) for x in self]

    @property
    def label_dict(self):
        """Return a dictionary of the name and label of each parameter.
        """
        return OrderedDict(self.labels)

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

# derived parameters (these are not used for waveform generation)
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
srcmass1 = Parameter("srcmass1", dtype=float,
                     label=r"$m_1^{\rm{src}}~(\mathrm{M}_\odot)$",
                     description="The mass of the first component object in "
                                 "the source frame (in solar masses).")
srcmass2 = Parameter("srcmass1", dtype=float,
                     label=r"$m_2^{\rm{src}}~(\mathrm{M}_\odot)$",
                     description="The mass of the second component object in "
                                 "the source frame (in solar masses).")
srcmchirp = Parameter("srcmchirp", dtype=float,
                      label=r"$\mathcal{M}^{\rm{src}}~(\mathrm{M}_\odot)$",
                      description="The chirp mass of the binary in the "
                                  "source frame (in solar masses).")
srcmtotal = Parameter("mtotal", dtype=float,
                      label=r"$M^{\rm{src}}~(\mathrm{M}_\odot)$",
                      description="The total mass of the binary in the "
                                  "source frame (in solar masses).")
primary_mass = Parameter("primary_mass",
                dtype=float, label=r"$m_{1}$",
                description="Mass of the primary object (in solar masses).")
secondary_mass = Parameter("secondary_mass",
                dtype=float, label=r"$m_{2}$",
                description="Mass of the secondary object (in solar masses).")

# derived parameters for component spins
chi_eff = Parameter("chi_eff",
                dtype=float, label=r"$\chi_\mathrm{eff}$",
                description="Effective spin of the binary.")
chi_p = Parameter("chi_p",
                dtype=float, label=r"$\chi_p$",
                description="Effective precessing spin of the binary.")
spin_px = Parameter("spin_px",
                dtype=float, label=r"$\chi_{1x}$",
                description="The x component of the dimensionless spin of the "
                            "primary object.")
spin_py = Parameter("spin_py",
                dtype=float, label=r"$\chi_{1y}$",
                description="The y component of the dimensionless spin of the "
                            "primary object.")
spin_pz = Parameter("spin_pz",
                dtype=float, label=r"$\chi_{1z}$",
                description="The z component of the dimensionless spin of the "
                            "primary object.")
spin_sx = Parameter("spin_sx",
                dtype=float, label=r"$\chi_{2x}$",
                description="The x component of the dimensionless spin of the "
                            "secondary object.")
spin_sy = Parameter("spin_sy",
                dtype=float, label=r"$\chi_{2y}$",
                description="The y component of the dimensionless spin of the "
                            "secondary object.")
spin_sz = Parameter("spin_sz",
                dtype=float, label=r"$\chi_{2z}$",
                description="The z component of the dimensionless spin of the "
                            "secondary object.")
lambda1 = Parameter("lambda1",
                dtype=float, default=None, label=r"$\Lambda_1$",
                description="The dimensionless tidal deformability parameter of object 1.")
lambda2 = Parameter("lambda2",
                dtype=float, default=None, label=r"$\Lambda_2$",
                description="The dimensionless tidal deformability parameter of object 2.")
dquad_mon1 = Parameter("dquad_mon1",
                dtype=float, default=None, label=r"$qm_1$",
                description="Quadrupole-monopole parameter / m_1^5 -1.")
dquad_mon2 = Parameter("dquad_mon2",
                dtype=float, default=None, label=r"$qm_2$",
                description="Quadrupole-monopole parameter / m_2^5 -1.")
lambda_octu1 = Parameter("lambda_octu1",
                dtype=float, default=None, label=r"$\Lambda_3^{(1)}$",
                description="The octupolar tidal deformability parameter of "
                            "object 1.")
lambda_octu2 = Parameter("lambda_octu2",
                dtype=float, default=None, label=r"$\Lambda_3^{(2)}$",
                description="The octupolar tidal deformability parameter of "
                            "object 2.")
quadfmode1 = Parameter("quadfmode1",
                dtype=float, default=None, label=r"$m_1 \omega_{02}^{(1)}$",
                description="The quadrupolar f-mode angular frequency of "
                            "object 1.")
quadfmode2 = Parameter("quadfmode2",
                dtype=float, default=None, label=r"$m_ \omega_{02}^{(2)}$",
                description="The quadrupolar f-mode angular frequency of "
                            "object 2.")
octufmode1 = Parameter("octufmode1",
                dtype=float, default=None, label=r"$m_1 \omega_{03}^{(1)}$",
                description="The octupolar f-mode angular frequency of "
                            "object 1.")
octufmode2 = Parameter("octufmode2",
                dtype=float, default=None, label=r"$m_ \omega_{03}^{(2)}$",
                description="The octupolar f-mode angular frequency of "
                            "object 2.")

# derived parameters for component spin magnitude and angles
spin1_a = Parameter("spin1_a",
                    dtype=float, label=r"$a_{1}$",
                    description="The dimensionless spin magnitude "
                                r"$|\vec{s}/m_{1}^2|$.")
spin2_a = Parameter("spin2_a",
                    dtype=float, label=r"$a_{2}$",
                    description="The dimensionless spin magnitude "
                                r"$|\vec{s}/m_{2}^2|$.")
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
f_final_func = Parameter("f_final_func",
                dtype=str, default="", label=None,
                description="Use the given frequency function to compute f_final "
                            "based on the parameters of the waveform.")
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
remnant_mass = Parameter("remnant_mass",
                        dtype=float, label=r"$m_{\mathrm{rem}}$",
                        description="Remnant mass of NS-BH merger. See "
                        "conversions.remnant_mass_"
                        "from_mass1_mass2_spin1x_spin1y_spin1z_eos")

#
#   General location parameters
#
distance = Parameter("distance",
                dtype=float, default=1., label=r"$d_L$ (Mpc)",
                description="Luminosity distance to the binary (in Mpc).")
chirp_distance = Parameter("chirp_distance",
                dtype=float, default=1., label=r"$d_c$ (Mpc)",
                description="Chirp distance to the binary (in Mpc).")
coa_phase = Parameter("coa_phase",
                dtype=float, default=0., label=r"$\phi_c$",
                description="Coalesence phase of the binary (in rad).")
inclination = Parameter("inclination",
                dtype=float, default=0., label=r"$\iota$",
                description="Inclination (rad), defined as the angle between "
                            "the orbital angular momentum L and the "
                            "line-of-sight at the reference frequency.")
thetajn = Parameter("thetajn",
                    dtype=float, default=0., label=r"$\theta_{JN}$",
                    description="The angle between the total angular momentum "
                                "J and the line-of-sight.")
long_asc_nodes = Parameter("long_asc_nodes",
                dtype=float, default=0., label=r"$\Omega$",
                description="Longitude of ascending nodes axis (rad).")
mean_per_ano = Parameter("mean_per_ano",
                dtype=float, default=0., label=r"$\delta$",
                description="Mean anomaly of the periastron (rad).")
tc = Parameter("tc",
                dtype=float, default=None, label=r"$t_c$ (s)",
                description="Coalescence time (s).")
delta_tc = Parameter("delta_tc", dtype=float,
                     label=r"$\Delta t_c~(\rm{s})$",
                     description="Coalesence time offset.")
ra = Parameter("ra",
                dtype=float, default=0., label=r"$\alpha$",
                description="Right ascension (rad).")
dec = Parameter("dec",
                dtype=float, default=0., label=r"$\delta$",
                description="Declination (rad).")
polarization = Parameter("polarization",
                dtype=float, default=0., label=r"$\psi$",
                description="Polarization (rad).")
redshift = Parameter("redshift",
                dtype=float, default=None, label=r"$z$",
                description="Redshift.")
comoving_volume = Parameter("comoving_volume", dtype=float,
                            label=r"$V_C~(\rm{Mpc}^3)$",
                            description="Comoving volume (in cubic Mpc).")
eclipticlatitude = Parameter("eclipticlatitude",
                dtype=float, default=0., label=r"$\beta$",
                description="eclipticlatitude wrt SSB coords.")
eclipticlongitude = Parameter("eclipticlongitude",
                dtype=float, default=0., label=r"$\lambda$",
                description="eclipticlongitude wrt SSB coords.")

#
#   Calibration parameters
#
delta_fs = Parameter("calib_delta_fs",
                     dtype=float,
                     description="Change in optical spring freq (Hz).")
delta_fc = Parameter("calib_delta_fc",
                     dtype=float,
                     description="Change in cavity pole freq (Hz).")
delta_qinv = Parameter("calib_delta_qinv",
                       dtype=float,
                       description="Change in inverse quality factor.")
kappa_c = Parameter("calib_kappa_c",
                    dtype=float)
kappa_tst_re = Parameter("calib_kappa_tst_re",
                         dtype=float)
kappa_tst_im = Parameter("calib_kappa_tst_im",
                         dtype=float)
kappa_pu_re = Parameter("calib_kappa_pu_re",
                        dtype=float)
kappa_pu_im = Parameter("calib_kappa_pu_im",
                        dtype=float)

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
mode_array = Parameter("mode_array",
                dtype=list, default=None,
                description="Choose which (l,m) modes to include when "
                            "generating a waveform. "
                            "Only if approximant supports this feature."
                            "By default pass None and let lalsimulation "
                            "use it's default behaviour."
                            "Example: mode_array = [ [2,2], [2,-2] ]")

#
#   Parametrized testing general relativity parameters
#
dchi0 = Parameter("dchi0",
                dtype=float, default=0., label=r"$d\chi_0$",
                description="0PN testingGR parameter.")
dchi1 = Parameter("dchi1",
                dtype=float, default=0., label=r"$d\chi_1$",
                description="0.5PN testingGR parameter.")
dchi2 = Parameter("dchi2",
                dtype=float, default=0., label=r"$d\chi_2$",
                description="1PN testingGR parameter.")
dchi3 = Parameter("dchi3",
                dtype=float, default=0., label=r"$d\chi_3$",
                description="1.5PN testingGR parameter.")
dchi4 = Parameter("dchi4",
                dtype=float, default=0., label=r"$d\chi_4$",
                description="2PN testingGR parameter.")
dchi5 = Parameter("dchi5",
                dtype=float, default=0., label=r"$d\chi_5$",
                description="2.5PN testingGR parameter.")
dchi5l = Parameter("dchi5l",
                dtype=float, default=0., label=r"$d\chi_5{l}$",
                description="2.5PN logrithm testingGR parameter.")
dchi6 = Parameter("dchi6",
                dtype=float, default=0., label=r"$d\chi_6$",
                description="3PN testingGR parameter.")
dchi6l = Parameter("dchi6l",
                dtype=float, default=0., label=r"$d\chi_{6l}$",
                description="3PN logrithm testingGR parameter.")
dchi7 = Parameter("dchi7",
                dtype=float, default=0., label=r"$d\chi_7$",
                description="3.5PN testingGR parameter.")
dalpha1 = Parameter("dalpha1",
                dtype=float, default=0., label=r"$d\alpha_1$",
                description="Merger-ringdown testingGR parameter.")
dalpha2 = Parameter("dalpha2",
                dtype=float, default=0., label=r"$d\alpha_2$",
                description="Merger-ringdown testingGR parameter.")
dalpha3 = Parameter("dalpha3",
                dtype=float, default=0., label=r"$d\alpha_3$",
                description="Merger-ringdown testingGR parameter.")
dalpha4 = Parameter("dalpha4",
                dtype=float, default=0., label=r"$d\alpha_4$",
                description="Merger-ringdown testingGR parameter.")
dalpha5 = Parameter("dalpha5",
                dtype=float, default=0., label=r"$d\alpha_5$",
                description="Merger-ringdown testingGR parameter.")
dbeta1 = Parameter("dbeta1",
                dtype=float, default=0., label=r"$d\beta_1$",
                description="Intermediate testingGR parameter.")
dbeta2 = Parameter("dbeta2",
                dtype=float, default=0., label=r"$d\beta_2$",
                description="Intermediate testingGR parameter.")
dbeta3 = Parameter("dbeta3",
                dtype=float, default=0., label=r"$d\beta_3$",
                description="Intermediate testingGR parameter.")
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
location_params = ParameterList([tc, ra, dec, polarization,
                                eclipticlatitude, eclipticlongitude])

# parameters describing the orientation of a binary w.r.t. the radiation
# frame. Note: we include distance here, as it is typically used for generating
# waveforms.
orientation_params = ParameterList\
    ([distance, coa_phase, inclination, long_asc_nodes, mean_per_ano])

# the extrinsic parameters of a waveform
extrinsic_params = orientation_params + location_params


# testing GR parameters
testingGR_params = ParameterList\
    ([dchi0, dchi1, dchi2, dchi3, dchi4, dchi5, dchi5l, dchi6, dchi6l,
      dchi7, dalpha1, dalpha2, dalpha3, dalpha4, dalpha5,
      dbeta1, dbeta2, dbeta3])

# intrinsic parameters of a CBC waveform. Some of these are not recognized
# by every waveform model
cbc_intrinsic_params = ParameterList\
    ([mass1, mass2, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z,
      eccentricity, lambda1, lambda2, dquad_mon1, dquad_mon2, lambda_octu1,
      lambda_octu2, quadfmode1, quadfmode2, octufmode1, octufmode2]) + \
    testingGR_params

# the parameters of a cbc in the radiation frame
cbc_rframe_params = cbc_intrinsic_params + orientation_params

# calibration parameters
calibration_params = ParameterList([
    delta_fc, delta_fs, delta_qinv, kappa_c, kappa_tst_re, kappa_tst_im,
    kappa_pu_re, kappa_pu_im])

# common generation parameters are parameters needed to generate either
# a TD, FD, or frequency sequence waveform
common_generation_params = ParameterList([
    approximant, f_ref, phase_order, spin_order, tidal_order, amplitude_order, eccentricity_order])

# Flags having discrete values, optional to generate either
# a TD, FD, or frequency sequence waveform
flags_generation_params = ParameterList([frame_axis, modes_choice, side_bands, mode_array])

# the following are parameters needed to generate an FD or TD waveform that
# is equally sampled
common_gen_equal_sampled_params = ParameterList([f_lower]) + \
    common_generation_params + flags_generation_params

# the following are parameters that can be used to generate an FD waveform
fd_waveform_params = cbc_rframe_params + ParameterList([delta_f]) + \
    common_gen_equal_sampled_params + ParameterList([f_final, f_final_func])

# the following are parameters that can be used to generate a TD waveform
td_waveform_params = cbc_rframe_params + ParameterList([delta_t]) + \
    common_gen_equal_sampled_params + ParameterList([numrel_data]) + \
    flags_generation_params

# The following are the minimum set of parameters that are required to
# generate a FD or TD waveform. All other parameters have some default value as
# defined above. Defaults of None simply mean that the value is not passed into
# the lal_dict structure and the waveform generator will take whatever default
# behaviour
td_required = ParameterList([f_lower, delta_t, approximant])
fd_required = ParameterList([f_lower, delta_f, approximant])

####
cbc_td_required = ParameterList([mass1, mass2, f_lower, delta_t, approximant])
cbc_fd_required = ParameterList([mass1, mass2, f_lower, delta_f, approximant])

# the following are parameters that can be used to generate a
# frequency series waveform
fd_waveform_sequence_params = cbc_rframe_params + \
    ParameterList([sample_points]) + common_generation_params + \
    flags_generation_params
