# Copyright (C) 2013 Ian W. Harry
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

import argparse
import logging
import os
import textwrap

import numpy

from pycbc import pnutils
from pycbc.neutron_stars import load_ns_sequence
from pycbc.tmpltbank.lambda_mapping import pycbc_valid_orders_help_descriptions
from pycbc.types import nonnegative_float, positive_float

logger = logging.getLogger(__name__)


class IndentedHelpFormatterWithNL(argparse.ArgumentDefaultsHelpFormatter):
    """
    This class taken from
    https://groups.google.com/forum/#!topic/comp.lang.python/bfbmtUGhW8I
    and is used to format the argparse help messages to deal with line breaking
    nicer. Specfically the pn-order help is large and looks crappy without this.
    This function is (C) Tim Chase
    """

    def format_description(self, description):
        """
        No documentation
        """
        if not description:
            return ""
        desc_width = self.width - self.current_indent
        indent = " " * self.current_indent
        # the above is still the same
        bits = description.split("\n")
        formatted_bits = [
            textwrap.fill(
                bit, desc_width, initial_indent=indent, subsequent_indent=indent
            )
            for bit in bits
        ]
        return "\n".join(formatted_bits) + "\n"

    def format_option(self, option):
        """
        No documentation
        """
        # The help for each option consists of two parts:
        #   * the opt strings and metavars
        #   eg. ("-x", or "-fFILENAME, --file=FILENAME")
        #   * the user-supplied help string
        #   eg. ("turn on expert mode", "read data from FILENAME")
        #
        # If possible, we write both of these on the same line:
        #   -x    turn on expert mode
        #
        # But if the opt string list is too long, we put the help
        # string on a second line, indented to the same column it would
        # start in if it fit on the first line.
        #   -fFILENAME, --file=FILENAME
        #       read data from FILENAME
        result = []
        opts = self.option_strings[option]
        opt_width = self.help_position - self.current_indent - 2
        if len(opts) > opt_width:
            opts = f"{'':{self.current_indent}}{opts}\n"
            indent_first = self.help_position
        else:  # start help on same line as opts
            opts = f"{'':{self.current_indent}}{opts:<{opt_width}}  "
            indent_first = 0
        result.append(opts)
        if option.help:
            help_text = self.expand_default(option)
            # Everything is the same up through here
            help_lines = []
            for para in help_text.split("\n"):
                help_lines.extend(textwrap.wrap(para, self.help_width))
            # Everything is the same after here
            result.append(f"{'':{indent_first}}{help_lines[0]}\n")
            result.extend(
                [f"{'':{self.help_position}}{line}\n" for line in help_lines[1:]]
            )
        elif opts[-1] != "\n":
            result.append("\n")
        return "".join(result)


def get_options_from_group(option_group):
    """
    Take an option group and return all the options that are defined in that
    group.
    """
    option_list = option_group._group_actions
    command_lines = []
    for option in option_list:
        option_strings = option.option_strings
        for string in option_strings:
            if string.startswith("--"):
                command_lines.append(string)
    return command_lines


def insert_base_bank_options(parser, match_req=True):
    """
    Adds essential common options for template bank generation to an
    ArgumentParser instance.
    """

    def match_type(s):
        err_msg = f"must be a number between 0 and 1 excluded, not {s!r}"
        try:
            value = float(s)
        except ValueError:
            raise argparse.ArgumentTypeError(err_msg)
        if value <= 0 or value >= 1:
            raise argparse.ArgumentTypeError(err_msg)
        return value

    parser.add_argument(
        "-m",
        "--min-match",
        type=match_type,
        required=match_req,
        help="Generate bank with specified minimum match. Required.",
    )
    parser.add_argument(
        "-O", "--output-file", required=True, help="Output file name. Required."
    )
    parser.add_argument(
        "--f-low-column",
        type=str,
        metavar="NAME",
        help="If given, store the lower frequency cutoff into "
        "column NAME of the single-inspiral table. "
        "(Requires an output file ending in .xml)",
    )
    parser.add_argument(
        "--output-f-final",
        action="store_true",
        default=False,
        help="Include 'f_final' in the output hdf file.",
    )


def insert_metric_calculation_options(parser):
    """
    Adds the options used to obtain a metric in the bank generation codes to an
    argparser as an OptionGroup. This should be used if you want to use these
    options in your code.
    """
    metric_opts = parser.add_argument_group(
        "Options related to calculating the parameter space metric"
    )
    metric_opts.add_argument(
        "--pn-order",
        action="store",
        type=str,
        required=True,
        help="Determines the PN order to use.  For a bank of "
        "non-spinning templates, spin-related terms in the "
        "metric will be zero.  REQUIRED. "
        f"Choices: {pycbc_valid_orders_help_descriptions}",
    )
    metric_opts.add_argument(
        "--f0",
        action="store",
        type=positive_float,
        default=70.0,
        help="f0 is used as a dynamic scaling factor when "
        "calculating integrals used in metric construction.  "
        "I.e. instead of integrating F(f) we integrate F(f/f0) "
        "then rescale by powers of f0.  The default value 70Hz "
        "should be fine for most applications.  OPTIONAL. "
        "UNITS=Hz",
    )
    metric_opts.add_argument(
        "--f-low",
        action="store",
        type=positive_float,
        required=True,
        help="Lower frequency cutoff used in computing the "
        "parameter space metric.  REQUIRED. UNITS=Hz",
    )
    metric_opts.add_argument(
        "--f-upper",
        action="store",
        type=positive_float,
        required=True,
        help="Upper frequency cutoff used in computing the "
        "parameter space metric.  REQUIRED. UNITS=Hz",
    )
    metric_opts.add_argument(
        "--delta-f",
        action="store",
        type=positive_float,
        required=True,
        help="Frequency spacing used in computing the parameter "
        r"space metric:  integrals of the form \int F(f) df "
        r"are approximated as \sum F(f) delta_f.  REQUIRED. "
        "UNITS=Hz",
    )
    metric_opts.add_argument(
        "--write-metric",
        action="store_true",
        default=False,
        help="If given write the metric components to disk as they are calculated.",
    )
    return metric_opts


def verify_metric_calculation_options(opts, parser):
    """
    Parses the metric calculation options given and verifies that they are
    correct.

    Parameters
    ----------
    opts : argparse.Values instance
        Result of parsing the input options with OptionParser
    parser : object
        The OptionParser instance.

    """
    if not opts.pn_order:
        parser.error("Must supply --pn-order")


class MetricParameters:
    """
    This class holds all of the options that are parsed in the function
    insert_metric_calculation_options
    and all products produced using these options. It can also be initialized
    from the __init__ function, providing directly the options normally
    provided on the command line.
    """

    _psd = None
    _metric = None
    _evals = None
    _evecs = None
    _evecs_cv = None

    def __init__(self, pn_order, f_low, f_upper, delta_f, f0=70, write_metric=False):
        """
        Initialize an instance of the MetricParameters by providing all
        options directly. See the help message associated with any code
        that uses the metric options for more details of how to set each of
        these, e.g. pycbc_aligned_stoch_bank --help
        """
        self.pn_order = pn_order
        self.f_low = f_low
        self.f_upper = f_upper
        self.delta_f = delta_f
        self.f0 = f0
        self._moments = None
        self.write_metric = write_metric

    @classmethod
    def from_argparse(cls, opts):
        """
        Initialize an instance of the MetricParameters class from an
        argparse.OptionParser instance. This assumes that
        insert_metric_calculation_options
        and
        verify_metric_calculation_options
        have already been called before initializing the class.
        """
        return cls(
            opts.pn_order,
            opts.f_low,
            opts.f_upper,
            opts.delta_f,
            f0=opts.f0,
            write_metric=opts.write_metric,
        )

    @property
    def psd(self):
        """
        A pyCBC FrequencySeries holding the appropriate PSD.
        Return the PSD used in the metric calculation.
        """
        if not self._psd:
            err_msg = "The PSD has not been set in the MetricParameters "
            err_msg += "instance."
            raise ValueError(err_msg)
        return self._psd

    @psd.setter
    def psd(self, in_psd):
        self._psd = in_psd

    @property
    def moments(self):
        """
        Moments structure
        This contains the result of all the integrals used in computing the
        metrics above. It can be used for additional metric calculations.
        This is composed of two compound
        dictionaries. The first entry indicates which moment is being
        calculated and the second entry indicates the upper frequency cutoff
        that was used.

        In all cases x = f/f0.

        For the first entries the options are:

        moments['J%d' %(i)][f_cutoff]
        This stores the integral of
        x**((-i)/3.) * delta X / PSD(x)

        moments['log%d' %(i)][f_cutoff]
        This stores the integral of
        (numpy.log(x**(1./3.))) x**((-i)/3.) * delta X / PSD(x)

        moments['loglog%d' %(i)][f_cutoff]
        This stores the integral of
        (numpy.log(x**(1./3.)))**2 x**((-i)/3.) * delta X / PSD(x)

        moments['loglog%d' %(i)][f_cutoff]
        This stores the integral of
        (numpy.log(x**(1./3.)))**3 x**((-i)/3.) * delta X / PSD(x)

        moments['loglog%d' %(i)][f_cutoff]
        This stores the integral of
        (numpy.log(x**(1./3.)))**4 x**((-i)/3.) * delta X / PSD(x)

        The second entry stores the frequency cutoff that was used when
        computing the integral.
        """
        return self._moments

    @moments.setter
    def moments(self, in_moments):
        self._moments = in_moments

    @property
    def evals(self):
        """
        The eigenvalues of the parameter space.
        This is a Dictionary of numpy.array
        Each entry in the dictionary corresponds to the different frequency
        ranges described in vary_fmax. If vary_fmax = False, the only entry
        will be f_upper, this corresponds to integrals in [f_low,f_upper). This
        entry is always present. Each other entry will use floats as keys to
        the dictionary. These floats give the upper frequency cutoff when it is
        varying.
        Each numpy.array contains the eigenvalues which, with the eigenvectors
        in evecs, are needed to rotate the
        coordinate system to one in which the metric is the identity matrix.
        """
        if self._evals is None:
            err_msg = "The metric eigenvalues have not been set in the "
            err_msg += "MetricParameters instance."
            raise ValueError(err_msg)
        return self._evals

    @evals.setter
    def evals(self, in_evals):
        if self.write_metric:
            for frequency in in_evals.keys():
                numpy.savetxt(f"metric_evals_{int(frequency)}.dat", in_evals[frequency])
        self._evals = in_evals

    @property
    def evecs(self):
        """
        The eigenvectors of the parameter space.
        This is a Dictionary of numpy.matrix
        Each entry in the dictionary is as described under evals.
        Each numpy.matrix contains the eigenvectors which, with the eigenvalues
        in evals, are needed to rotate the
        coordinate system to one in which the metric is the identity matrix.
        """
        if self._evecs is None:
            err_msg = "The metric eigenvectors have not been set in the "
            err_msg += "MetricParameters instance."
            raise ValueError(err_msg)
        return self._evecs

    @evecs.setter
    def evecs(self, in_evecs):
        if self.write_metric:
            for frequency in in_evecs.keys():
                numpy.savetxt(f"metric_evecs_{int(frequency)}.dat", in_evecs[frequency])
        self._evecs = in_evecs

    @property
    def metric(self):
        """
        The metric of the parameter space.
        This is a Dictionary of numpy.matrix
        Each entry in the dictionary is as described under evals.
        Each numpy.matrix contains the metric of the parameter space in the
        Lambda_i coordinate system.
        """
        if self._metric is None:
            err_msg = "The metric eigenvectors have not been set in the "
            err_msg += "MetricParameters instance."
            raise ValueError(err_msg)
        return self._metric

    @metric.setter
    def metric(self, in_metric):
        if self.write_metric:
            for frequency in in_metric.keys():
                numpy.savetxt(
                    f"metric_components_{int(frequency)}.dat", in_metric[frequency]
                )
        self._metric = in_metric

    @property
    def time_unprojected_metric(self):
        """
        The metric of the parameter space with the time dimension unprojected.
        This is a Dictionary of numpy.matrix
        Each entry in the dictionary is as described under evals.
        Each numpy.matrix contains the metric of the parameter space in the
        Lambda_i, t coordinate system. The time components are always in the
        last [-1] position in the matrix.
        """
        if self._time_unprojected_metric is None:
            err_msg = "The time unprojected metric has not been set in the "
            err_msg += "MetricParameters instance."
            raise ValueError(err_msg)
        return self._time_unprojected_metric

    @time_unprojected_metric.setter
    def time_unprojected_metric(self, in_metric):
        if self.write_metric:
            for frequency in in_metric.keys():
                numpy.savetxt(
                    f"metric_timeunprojected_{int(frequency)}.dat", in_metric[frequency]
                )
        self._time_unprojected_metric = in_metric

    @property
    def evecs_cv(self):
        """
        The eigenvectors of the principal directions of the mu space.
        This is a Dictionary of numpy.matrix
        Each entry in the dictionary is as described under evals.
        Each numpy.matrix contains the eigenvectors which, with the eigenvalues
        in evals, are needed to rotate the
        coordinate system to one in which the metric is the identity matrix.
        """
        if self._evecs_cv is None:
            err_msg = "The covariance eigenvectors have not been set in the "
            err_msg += "MetricParameters instance."
            raise ValueError(err_msg)
        return self._evecs_cv

    @evecs_cv.setter
    def evecs_cv(self, in_evecs):
        if self.write_metric:
            for frequency in in_evecs.keys():
                numpy.savetxt(
                    f"covariance_evecs_{int(frequency)}.dat", in_evecs[frequency]
                )
        self._evecs_cv = in_evecs


def insert_mass_range_option_group(parser, non_spin=False):
    """
    Adds the options used to specify mass ranges in the bank generation codes
    to an argparser as an OptionGroup. This should be used if you
    want to use these options in your code.

    Parameters
    ----------
    parser : object
        OptionParser instance.
    non_spin : boolean, optional (default=False)
        If this is provided the spin-related options will not be added.

    """
    mass_opts = parser.add_argument_group(
        "Options related to mass and spin limits for bank generation"
    )
    mass_opts.add_argument(
        "--min-mass1",
        action="store",
        type=positive_float,
        required=True,
        help="Minimum mass1: must be >= min-mass2. REQUIRED. UNITS=Solar mass",
    )
    mass_opts.add_argument(
        "--max-mass1",
        action="store",
        type=positive_float,
        required=True,
        help="Maximum mass1: must be >= max-mass2. REQUIRED. UNITS=Solar mass",
    )
    mass_opts.add_argument(
        "--min-mass2",
        action="store",
        type=positive_float,
        required=True,
        help="Minimum mass2. REQUIRED. UNITS=Solar mass",
    )
    mass_opts.add_argument(
        "--max-mass2",
        action="store",
        type=positive_float,
        required=True,
        help="Maximum mass2. REQUIRED. UNITS=Solar mass",
    )
    mass_opts.add_argument(
        "--max-total-mass",
        action="store",
        type=positive_float,
        default=None,
        help="Maximum total mass. OPTIONAL, if not provided "
        "the max total mass is determined by the component "
        "masses. UNITS=Solar mass",
    )
    mass_opts.add_argument(
        "--min-total-mass",
        action="store",
        type=positive_float,
        default=None,
        help="Minimum total mass. OPTIONAL, if not provided the "
        "min total mass is determined by the component masses."
        " UNITS=Solar mass",
    )
    mass_opts.add_argument(
        "--max-chirp-mass",
        action="store",
        type=positive_float,
        default=None,
        help="Maximum chirp mass. OPTIONAL, if not provided the "
        "max chirp mass is determined by the component masses."
        " UNITS=Solar mass",
    )
    mass_opts.add_argument(
        "--min-chirp-mass",
        action="store",
        type=positive_float,
        default=None,
        help="Minimum total mass. OPTIONAL, if not provided the "
        "min chirp mass is determined by the component masses."
        " UNITS=Solar mass",
    )
    mass_opts.add_argument(
        "--max-eta",
        action="store",
        type=positive_float,
        default=0.25,
        help="Maximum symmetric mass ratio. OPTIONAL, no upper bound"
        " on eta will be imposed if not provided. "
        "UNITS=Solar mass.",
    )
    mass_opts.add_argument(
        "--min-eta",
        action="store",
        type=nonnegative_float,
        default=0.0,
        help="Minimum symmetric mass ratio. OPTIONAL, no lower bound"
        " on eta will be imposed if not provided. "
        "UNITS=Solar mass.",
    )
    mass_opts.add_argument(
        "--ns-eos",
        action="store",
        default=None,
        help="Select the EOS to be used for the NS when calculating "
        "the remnant disk mass. Only 2H is currently supported. "
        "OPTIONAL",
    )
    mass_opts.add_argument(
        "--remnant-mass-threshold",
        action="store",
        type=nonnegative_float,
        default=None,
        help="Setting this filters EM dim NS-BH binaries: if the "
        "remnant disk mass does not exceed this value, the NS-BH "
        "binary is dropped from the target parameter space. "
        "When it is set to None (default value) the EM dim "
        "filter is not activated. OPTIONAL",
    )
    mass_opts.add_argument(
        "--use-eos-max-ns-mass",
        action="store_true",
        default=False,
        help="Cut the mass range of the smaller object to the maximum "
        "mass allowed by EOS. "
        "OPTIONAL",
    )
    mass_opts.add_argument(
        "--delta-bh-spin",
        action="store",
        type=positive_float,
        default=None,
        help="Grid spacing used for the BH spin z component when "
        "generating the surface of the minumum minimum symmetric "
        "mass ratio as a function of BH spin and NS mass required "
        "to produce a remnant disk mass that exceeds the threshold "
        "specificed in --remnant-mass-threshold. "
        "OPTIONAL (0.1 by default) ",
    )
    mass_opts.add_argument(
        "--delta-ns-mass",
        action="store",
        type=positive_float,
        default=None,
        help="Grid spacing used for the NS mass when generating the "
        "surface of the minumum minimum symmetric mass ratio as "
        "a function of BH spin and NS mass required to produce "
        "a remnant disk mass that exceeds the thrsehold specified "
        "in --remnant-mass-threshold. "
        "OPTIONAL (0.1 by default) ",
    )
    if non_spin:
        parser.add_argument_group(mass_opts)
        return mass_opts

    mass_opts.add_argument(
        "--max-ns-spin-mag",
        action="store",
        type=nonnegative_float,
        default=None,
        help="Maximum neutron star spin magnitude.  Neutron stars "
        "are defined as components lighter than the NS-BH "
        "boundary (3 Msun by default). REQUIRED if min-mass2 "
        "< ns-bh-boundary-mass",
    )
    mass_opts.add_argument(
        "--max-bh-spin-mag",
        action="store",
        type=nonnegative_float,
        default=None,
        help="Maximum black hole spin magnitude.  Black holes are "
        "defined as components at or above the NS-BH boundary "
        "(3 Msun by default). REQUIRED if max-mass1 >= "
        "ns-bh-boundary-mass",
    )
    # Mutually exclusive group prevents both options being set on command line
    # If --nsbh-flag is True then spinning bank generation must ignore the
    # default value of ns-bh-boundary-mass.
    action = mass_opts.add_mutually_exclusive_group(required=False)
    action.add_argument(
        "--ns-bh-boundary-mass",
        action="store",
        type=positive_float,
        help="Mass boundary between neutron stars and black holes. "
        "Components below this mass are considered neutron "
        "stars and are subject to the neutron star spin limits. "
        "Components at/above are subject to the black hole spin "
        "limits.  OPTIONAL, "
        f"default={MassRangeParameters.default_nsbh_boundary_mass:f}.  "
        "UNITS=Solar mass",
    )
    action.add_argument(
        "--nsbh-flag",
        action="store_true",
        default=False,
        help="Set this flag if generating a bank that contains only "
        "systems with 1 black hole and 1 neutron star. With "
        "this flag set the heavier body will always be subject "
        "to the black hole spin restriction and the lighter "
        "to the neutron star spin restriction, regardless of "
        "mass.  OPTIONAL.  If set, the value of "
        "--ns-bh-boundary-mass will be ignored.",
    )
    return mass_opts


def verify_mass_range_options(opts, parser, non_spin=False):
    """
    Parses the metric calculation options given and verifies that they are
    correct.

    Parameters
    ----------
    opts : argparse.Values instance
        Result of parsing the input options with OptionParser
    parser : object
        The OptionParser instance.
    non_spin : boolean, optional (default=False)
        If this is provided the spin-related options will not be checked.

    """
    # Mass1 must be the heavier!
    if opts.min_mass1 < opts.min_mass2:
        parser.error("min-mass1 cannot be less than min-mass2!")
    if opts.max_mass1 < opts.max_mass2:
        parser.error("max-mass1 cannot be less than max-mass2!")
    # If given are min/max total mass/chirp mass possible?
    if opts.min_total_mass and (opts.min_total_mass > opts.max_mass1 + opts.max_mass2):
        err_msg = f"Supplied minimum total mass {opts.min_total_mass:f} "
        err_msg += "greater than the sum of the two max component masses "
        err_msg += f" {opts.max_mass1:f} and {opts.max_mass2:f}."
        parser.error(err_msg)
    if opts.max_total_mass and (opts.max_total_mass < opts.min_mass1 + opts.min_mass2):
        err_msg = f"Supplied maximum total mass {opts.max_total_mass:f} "
        err_msg += "smaller than the sum of the two min component masses "
        err_msg += f" {opts.min_mass1:f} and {opts.min_mass2:f}."
        parser.error(err_msg)
    if (
        opts.max_total_mass
        and opts.min_total_mass
        and (opts.max_total_mass < opts.min_total_mass)
    ):
        parser.error("Min total mass must be larger than max total mass.")
    # Warn the user that his/her setup is such that EM dim NS-BH binaries
    # will not be targeted by the template bank that is being built.  Also
    # inform him/her about the caveats involved in this.
    if (
        hasattr(opts, "remnant_mass_threshold")
        and opts.remnant_mass_threshold is not None
    ):
        logger.info("""You have asked to exclude EM dim NS-BH systems from the
                       target parameter space. The script will assume that m1
                       is the BH and m2 is the NS: make sure that your settings
                       respect this convention. The script will also treat the
                       NS as non-spinning: use NS spins in the template bank
                       at your own risk!""")
        if opts.use_eos_max_ns_mass:
            logger.info("""You have asked to take into account the maximum NS
                        mass value for the EOS in use.""")
        # Find out if the EM constraint surface data already exists or not
        # and inform user whether this will be read from file or generated.
        # This is the minumum eta as a function of BH spin and NS mass
        # required to produce an EM counterpart
        if os.path.isfile("constraint_em_bright.npz"):
            logger.info("""The constraint surface for EM bright binaries
                        will be read in from constraint_em_bright.npz.""")

    # Assign min/max total mass from mass1, mass2 if not specified
    if (not opts.min_total_mass) or (
        (opts.min_mass1 + opts.min_mass2) > opts.min_total_mass
    ):
        opts.min_total_mass = opts.min_mass1 + opts.min_mass2
    if (not opts.max_total_mass) or (
        (opts.max_mass1 + opts.max_mass2) < opts.max_total_mass
    ):
        opts.max_total_mass = opts.max_mass1 + opts.max_mass2

    # It is vital that min and max total mass be set correctly.
    # This is becasue the heavily-used function get_random_mass will place
    # points first in total mass (to some power), and then in eta. If the total
    # mass limits are not well known ahead of time it will place unphysical
    # points and fail.
    # This test is a bit convoluted as we identify the maximum and minimum
    # possible total mass from chirp mass and/or eta restrictions.
    if opts.min_chirp_mass is not None:
        # Need to get the smallest possible min_tot_mass from this chirp mass
        # There are 4 possibilities for where the min_tot_mass is found on the
        # line of min_chirp_mass that interacts with the component mass limits.
        # Either it is found at max_m2, or at min_m1, or it starts on the equal
        # mass line within the parameter space, or it doesn't intersect
        # at all.
        # First let's get the masses at both of these possible points
        m1_at_max_m2 = pnutils.mchirp_mass1_to_mass2(
            opts.min_chirp_mass, opts.max_mass2
        )
        if m1_at_max_m2 < opts.max_mass2:
            # Unphysical, remove
            m1_at_max_m2 = -1
        m2_at_min_m1 = pnutils.mchirp_mass1_to_mass2(
            opts.min_chirp_mass, opts.min_mass1
        )
        if m2_at_min_m1 > opts.min_mass1:
            # Unphysical, remove
            m2_at_min_m1 = -1
        # Get the values on the equal mass line
        m1_at_equal_mass, m2_at_equal_mass = pnutils.mchirp_eta_to_mass1_mass2(
            opts.min_chirp_mass, 0.25
        )

        # Are any of these possible?
        if m1_at_max_m2 <= opts.max_mass1 and m1_at_max_m2 >= opts.min_mass1:
            min_tot_mass = opts.max_mass2 + m1_at_max_m2
        elif m2_at_min_m1 <= opts.max_mass2 and m2_at_min_m1 >= opts.min_mass2:
            min_tot_mass = opts.min_mass1 + m2_at_min_m1
        elif (
            m1_at_equal_mass <= opts.max_mass1
            and m1_at_equal_mass >= opts.min_mass1
            and m2_at_equal_mass <= opts.max_mass2
            and m2_at_equal_mass >= opts.min_mass2
        ):
            min_tot_mass = m1_at_equal_mass + m2_at_equal_mass
        # So either the restriction is low enough to be redundant, or is
        # removing all the parameter space
        elif m2_at_min_m1 < opts.min_mass2:
            # This is the redundant case, ignore
            min_tot_mass = opts.min_total_mass
        else:
            # And this is the bad case
            err_msg = "The minimum chirp mass provided is not possible given "
            err_msg += "restrictions on component masses."
            raise ValueError(err_msg)
        # Is there also an eta restriction?
        if opts.max_eta:
            # Get the value of m1,m2 at max_eta, min_chirp_mass
            max_eta_m1, max_eta_m2 = pnutils.mchirp_eta_to_mass1_mass2(
                opts.min_chirp_mass, opts.max_eta
            )
            max_eta_min_tot_mass = max_eta_m1 + max_eta_m2
            if max_eta_min_tot_mass > min_tot_mass:
                # Okay, eta does restrict this further. Still physical?
                min_tot_mass = max_eta_min_tot_mass
                if max_eta_m1 > opts.max_mass1:
                    err_msg = "The combination of component mass, chirp "
                    err_msg += "mass, eta and (possibly) total mass limits "
                    err_msg += "have precluded all systems."
                    raise ValueError(err_msg)
        # Update min_tot_mass if needed
        if min_tot_mass > opts.min_total_mass:
            opts.min_total_mass = float(min_tot_mass)

    # Then need to do max_chirp_mass and min_eta
    if opts.max_chirp_mass is not None:
        # Need to get the largest possible maxn_tot_mass from this chirp mass
        # There are 3 possibilities for where the max_tot_mass is found on the
        # line of max_chirp_mass that interacts with the component mass limits.
        # Either it is found at min_m2, or at max_m1, or it doesn't intersect
        # at all.
        # First let's get the masses at both of these possible points
        m1_at_min_m2 = pnutils.mchirp_mass1_to_mass2(
            opts.max_chirp_mass, opts.min_mass2
        )
        m2_at_max_m1 = pnutils.mchirp_mass1_to_mass2(
            opts.max_chirp_mass, opts.max_mass1
        )
        # Are either of these possible?
        if m1_at_min_m2 <= opts.max_mass1 and m1_at_min_m2 >= opts.min_mass1:
            max_tot_mass = opts.min_mass2 + m1_at_min_m2
        elif m2_at_max_m1 <= opts.max_mass2 and m2_at_max_m1 >= opts.min_mass2:
            max_tot_mass = opts.max_mass1 + m2_at_max_m1
        # So either the restriction is low enough to be redundant, or is
        # removing all the paramter space
        elif m2_at_max_m1 > opts.max_mass2:
            # This is the redundant case, ignore
            max_tot_mass = opts.max_total_mass
        else:
            # And this is the bad case
            err_msg = "The maximum chirp mass provided is not possible given "
            err_msg += "restrictions on component masses."
            raise ValueError(err_msg)
        # Is there also an eta restriction?
        if opts.min_eta:
            # Get the value of m1,m2 at max_eta, min_chirp_mass
            min_eta_m1, min_eta_m2 = pnutils.mchirp_eta_to_mass1_mass2(
                opts.max_chirp_mass, opts.min_eta
            )
            min_eta_max_tot_mass = min_eta_m1 + min_eta_m2
            if min_eta_max_tot_mass < max_tot_mass:
                # Okay, eta does restrict this further. Still physical?
                max_tot_mass = min_eta_max_tot_mass
                if min_eta_m1 < opts.min_mass1:
                    err_msg = "The combination of component mass, chirp "
                    err_msg += "mass, eta and (possibly) total mass limits "
                    err_msg += "have precluded all systems."
                    raise ValueError(err_msg)
        # Update min_tot_mass if needed
        if max_tot_mass < opts.max_total_mass:
            opts.max_total_mass = float(max_tot_mass)

    # Need to check max_eta alone for minimum and maximum mass
    if opts.max_eta:
        # Similar to above except this can affect both the minimum and maximum
        # total mass. Need to identify where the line of max_eta intersects
        # the parameter space, and if it affects mass restrictions.
        m1_at_min_m2 = pnutils.eta_mass1_to_mass2(
            opts.max_eta, opts.min_mass2, return_mass_heavier=True
        )
        m2_at_min_m1 = pnutils.eta_mass1_to_mass2(
            opts.max_eta, opts.min_mass1, return_mass_heavier=False
        )
        m1_at_max_m2 = pnutils.eta_mass1_to_mass2(
            opts.max_eta, opts.max_mass2, return_mass_heavier=True
        )
        m2_at_max_m1 = pnutils.eta_mass1_to_mass2(
            opts.max_eta, opts.max_mass1, return_mass_heavier=False
        )
        # Check for restrictions on the minimum total mass
        # Are either of these possible?
        if m1_at_min_m2 <= opts.max_mass1 and m1_at_min_m2 >= opts.min_mass1:
            min_tot_mass = opts.min_mass2 + m1_at_min_m2
        elif m2_at_min_m1 <= opts.max_mass2 and m2_at_min_m1 >= opts.min_mass2:
            # This case doesn't change the minimal total mass
            min_tot_mass = opts.min_total_mass
        # So either the restriction is low enough to be redundant, or is
        # removing all the paramter space
        elif m2_at_min_m1 > opts.max_mass2:
            # This is the redundant case, ignore
            min_tot_mass = opts.min_total_mass
        elif opts.max_eta == 0.25 and (
            m1_at_min_m2 < opts.min_mass2 or m2_at_min_m1 > opts.min_mass1
        ):
            # This just catches potential roundoff issues in the case that
            # max-eta is not used
            min_tot_mass = opts.min_total_mass
        else:
            # And this is the bad case
            err_msg = "The maximum eta provided is not possible given "
            err_msg += "restrictions on component masses."
            print(m1_at_min_m2, m2_at_min_m1, m1_at_max_m2, m2_at_max_m1)
            print(opts.min_mass1, opts.max_mass1, opts.min_mass2, opts.max_mass2)
            raise ValueError(err_msg)
        # Update min_tot_mass if needed
        if min_tot_mass > opts.min_total_mass:
            opts.min_total_mass = float(min_tot_mass)

        # Check for restrictions on the maximum total mass
        # Are either of these possible?
        if m2_at_max_m1 <= opts.max_mass2 and m2_at_max_m1 >= opts.min_mass2:
            max_tot_mass = opts.max_mass1 + m2_at_max_m1
        elif m1_at_max_m2 <= opts.max_mass1 and m1_at_max_m2 >= opts.min_mass1:
            # This case doesn't change the maximal total mass
            max_tot_mass = opts.max_total_mass
        # So either the restriction is low enough to be redundant, or is
        # removing all the paramter space, the latter case is already tested
        else:
            # This is the redundant case, ignore
            max_tot_mass = opts.max_total_mass
        if max_tot_mass < opts.max_total_mass:
            opts.max_total_mass = float(max_tot_mass)

    # Need to check min_eta alone for maximum and minimum total mass
    if opts.min_eta:
        # Same as max_eta.
        # Need to identify where the line of max_eta intersects
        # the parameter space, and if it affects mass restrictions.
        m1_at_min_m2 = pnutils.eta_mass1_to_mass2(
            opts.min_eta, opts.min_mass2, return_mass_heavier=True
        )
        m2_at_min_m1 = pnutils.eta_mass1_to_mass2(
            opts.min_eta, opts.min_mass1, return_mass_heavier=False
        )
        m1_at_max_m2 = pnutils.eta_mass1_to_mass2(
            opts.min_eta, opts.max_mass2, return_mass_heavier=True
        )
        m2_at_max_m1 = pnutils.eta_mass1_to_mass2(
            opts.min_eta, opts.max_mass1, return_mass_heavier=False
        )

        # Check for restrictions on the maximum total mass
        # Are either of these possible?
        if m1_at_max_m2 <= opts.max_mass1 and m1_at_max_m2 >= opts.min_mass1:
            max_tot_mass = opts.max_mass2 + m1_at_max_m2

        elif m2_at_max_m1 <= opts.max_mass2 and m2_at_max_m1 >= opts.min_mass2:
            # This case doesn't affect the maximum total mass
            max_tot_mass = opts.max_total_mass
        # So either the restriction is low enough to be redundant, or is
        # removing all the paramter space
        elif m2_at_max_m1 < opts.min_mass2:
            # This is the redundant case, ignore
            max_tot_mass = opts.max_total_mass
        else:
            # And this is the bad case
            err_msg = "The minimum eta provided is not possible given "
            err_msg += "restrictions on component masses."
            raise ValueError(err_msg)
        # Update min_tot_mass if needed
        if max_tot_mass < opts.max_total_mass:
            opts.max_total_mass = float(max_tot_mass)

        # Check for restrictions on the minimum total mass
        # Are either of these possible?
        if m2_at_min_m1 <= opts.max_mass2 and m2_at_min_m1 >= opts.min_mass2:
            min_tot_mass = opts.min_mass1 + m2_at_min_m1
        elif m1_at_min_m2 <= opts.max_mass1 and m1_at_min_m2 >= opts.min_mass1:
            # This case doesn't change the maximal total mass
            min_tot_mass = opts.min_total_mass
        # So either the restriction is low enough to be redundant, or is
        # removing all the paramter space, which is tested above
        else:
            # This is the redundant case, ignore
            min_tot_mass = opts.min_total_mass
        if min_tot_mass > opts.min_total_mass:
            opts.min_total_mass = float(min_tot_mass)

    if opts.max_total_mass < opts.min_total_mass:
        err_msg = "After including restrictions on chirp mass, component mass, "
        err_msg += "eta and total mass, no physical systems are possible."
        raise ValueError(err_msg)

    if opts.max_eta and opts.min_eta and (opts.max_eta < opts.min_eta):
        parser.error("--max-eta must be larger than --min-eta.")
    if non_spin:
        return

    if opts.max_ns_spin_mag is None:
        if opts.nsbh_flag:
            parser.error("Must supply --max_ns_spin_mag with --nsbh-flag")
        # Can ignore this if no NSs will be generated
        elif opts.min_mass2 < (
            opts.ns_bh_boundary_mass or MassRangeParameters.default_nsbh_boundary_mass
        ):
            parser.error(
                "Must supply --max-ns-spin-mag for the chosen value of --min_mass2"
            )
        else:
            opts.max_ns_spin_mag = opts.max_bh_spin_mag
    if opts.max_bh_spin_mag is None:
        if opts.nsbh_flag:
            parser.error("Must supply --max_bh_spin_mag with --nsbh-flag")
        # Can ignore this if no BHs will be generated
        if opts.max_mass1 >= (
            opts.ns_bh_boundary_mass or MassRangeParameters.default_nsbh_boundary_mass
        ):
            parser.error(
                "Must supply --max-bh-spin-mag for the chosen value of --max_mass1"
            )
        else:
            opts.max_bh_spin_mag = opts.max_ns_spin_mag


class MassRangeParameters:
    """
    This class holds all of the options that are parsed in the function
    insert_mass_range_option_group
    and all products produced using these options. It can also be initialized
    from the __init__ function providing directly the options normally
    provided on the command line
    """

    default_nsbh_boundary_mass = 3.0
    default_ns_eos = "2H"
    default_delta_bh_spin = 0.1
    default_delta_ns_mass = 0.1

    def __init__(
        self,
        min_mass1,
        max_mass1,
        min_mass2,
        max_mass2,
        max_ns_spin_mag=0,
        max_bh_spin_mag=0,
        max_tot_mass=None,
        min_tot_mass=None,
        max_eta=None,
        min_eta=0,
        max_chirp_mass=None,
        min_chirp_mass=None,
        ns_bh_boundary_mass=None,
        nsbh_flag=False,
        remnant_mass_threshold=None,
        ns_eos=None,
        use_eos_max_ns_mass=False,
        delta_bh_spin=None,
        delta_ns_mass=None,
    ):
        """
        Initialize an instance of the MassRangeParameters by providing all
        options directly. See the help message associated with any code
        that uses the metric options for more details of how to set each of
        these. For e.g. pycbc_aligned_stoch_bank --help
        """
        self.min_mass1 = min_mass1
        self.max_mass1 = max_mass1
        self.min_mass2 = min_mass2
        self.max_mass2 = max_mass2
        self.max_ns_spin_mag = max_ns_spin_mag
        self.max_bh_spin_mag = max_bh_spin_mag
        self.min_tot_mass = min_mass1 + min_mass2
        if min_tot_mass and (min_tot_mass > self.min_tot_mass):
            self.min_tot_mass = min_tot_mass
        self.max_tot_mass = max_mass1 + max_mass2
        if max_tot_mass and (max_tot_mass < self.max_tot_mass):
            self.max_tot_mass = max_tot_mass
        self.max_tot_mass = max_tot_mass
        self.min_tot_mass = min_tot_mass
        if max_eta:
            self.max_eta = max_eta
        else:
            self.max_eta = 0.25
        self.max_chirp_mass = max_chirp_mass
        self.min_chirp_mass = min_chirp_mass
        self.min_eta = min_eta
        self.ns_bh_boundary_mass = (
            ns_bh_boundary_mass or self.default_nsbh_boundary_mass
        )
        self.nsbh_flag = nsbh_flag
        self.remnant_mass_threshold = remnant_mass_threshold
        self.ns_eos = ns_eos or self.default_ns_eos
        self.delta_bh_spin = delta_bh_spin or self.default_delta_bh_spin
        self.delta_ns_mass = delta_ns_mass or self.default_delta_ns_mass
        self.use_eos_max_ns_mass = use_eos_max_ns_mass
        if self.remnant_mass_threshold is not None:
            if self.ns_eos != "2H":
                err_msg = """
                         By setting a value for --remnant-mass-threshold
                         you have asked to filter out EM dim NS-BH templates.
                         The EOS you chose is not supported currently: please
                         remove the --ns-eos option from your command line or
                         set it to '2H'.
                         """
                raise ValueError(err_msg)
            if use_eos_max_ns_mass:
                _, max_ns_g_mass = load_ns_sequence(self.ns_eos)
                if self.max_mass2 > max_ns_g_mass:
                    err_msg = f"""
                             The maximum NS mass supported by this EOS is
                             {max_ns_g_mass - 0.0000000001}. Please set --max-mass2 to this value or run
                             without the --use-eos-max-ns-mass flag.
                             """
                    raise ValueError(err_msg)
            self.delta_bh_spin = delta_bh_spin or self.default_delta_bh_spin
            self.delta_ns_mass = delta_ns_mass or self.default_delta_ns_mass

        # FIXME: This may be inaccurate if Eta limits are given
        # This will not cause any problems, but maybe could be fixed.
        self.min_comp_mass = self.min_mass2
        self.max_comp_mass = self.max_mass1

        # WARNING: We expect mass1 > mass2 ALWAYS
        # Check input:
        if (min_mass2 > min_mass1) or (max_mass2 > max_mass1):
            err_msg = "Mass1 must be larger than mass2. Check input options."
            raise ValueError(err_msg)

        if (min_mass2 > max_mass2) or (min_mass1 > max_mass1):
            err_msg = "Minimum masses cannot be larger than maximum masses."
            err_msg += "Check input options."
            raise ValueError(err_msg)

    @classmethod
    def from_argparse(cls, opts, non_spin=False):
        """
        Initialize an instance of the MassRangeParameters class from an
        argparse.OptionParser instance. This assumes that
        insert_mass_range_option_group
        and
        verify_mass_range_options
        have already been called before initializing the class.
        """
        if non_spin:
            return cls(
                opts.min_mass1,
                opts.max_mass1,
                opts.min_mass2,
                opts.max_mass2,
                max_tot_mass=opts.max_total_mass,
                min_tot_mass=opts.min_total_mass,
                max_eta=opts.max_eta,
                min_eta=opts.min_eta,
                max_chirp_mass=opts.max_chirp_mass,
                min_chirp_mass=opts.min_chirp_mass,
                remnant_mass_threshold=opts.remnant_mass_threshold,
                ns_eos=opts.ns_eos,
                use_eos_max_ns_mass=opts.use_eos_max_ns_mass,
                delta_bh_spin=opts.delta_bh_spin,
                delta_ns_mass=opts.delta_ns_mass,
            )
        return cls(
            opts.min_mass1,
            opts.max_mass1,
            opts.min_mass2,
            opts.max_mass2,
            max_tot_mass=opts.max_total_mass,
            min_tot_mass=opts.min_total_mass,
            max_eta=opts.max_eta,
            min_eta=opts.min_eta,
            max_ns_spin_mag=opts.max_ns_spin_mag,
            max_bh_spin_mag=opts.max_bh_spin_mag,
            nsbh_flag=opts.nsbh_flag,
            max_chirp_mass=opts.max_chirp_mass,
            min_chirp_mass=opts.min_chirp_mass,
            ns_bh_boundary_mass=opts.ns_bh_boundary_mass,
            remnant_mass_threshold=opts.remnant_mass_threshold,
            ns_eos=opts.ns_eos,
            use_eos_max_ns_mass=opts.use_eos_max_ns_mass,
            delta_bh_spin=opts.delta_bh_spin,
            delta_ns_mass=opts.delta_ns_mass,
        )

    def is_outside_range(self, mass1, mass2, spin1z, spin2z):
        """
        Test if a given location in mass1, mass2, spin1z, spin2z is within the
        range of parameters allowed by the massParams object.
        """
        # Mass1 test
        if mass1 * 1.001 < self.min_mass1:
            return 1
        if mass1 > self.max_mass1 * 1.001:
            return 1
        # Mass2 test
        if mass2 * 1.001 < self.min_mass2:
            return 1
        if mass2 > self.max_mass2 * 1.001:
            return 1
        # Spin1 test
        if self.nsbh_flag:
            if abs(spin1z) > self.max_bh_spin_mag * 1.001:
                return 1
        else:
            spin1z_mag = abs(spin1z)
            if not (
                (
                    mass1 * 1.001 > self.ns_bh_boundary_mass
                    and spin1z_mag <= self.max_bh_spin_mag * 1.001
                )
                or (
                    mass1 < self.ns_bh_boundary_mass * 1.001
                    and spin1z_mag <= self.max_ns_spin_mag * 1.001
                )
            ):
                return 1
        # Spin2 test
        if self.nsbh_flag:
            if abs(spin2z) > self.max_ns_spin_mag * 1.001:
                return 1
        else:
            spin2z_mag = abs(spin2z)
            if not (
                (
                    mass2 * 1.001 > self.ns_bh_boundary_mass
                    and spin2z_mag <= self.max_bh_spin_mag * 1.001
                )
                or (
                    mass2 < self.ns_bh_boundary_mass * 1.001
                    and spin2z_mag <= self.max_ns_spin_mag * 1.001
                )
            ):
                return 1
        # Total mass test
        m_tot = mass1 + mass2
        if m_tot > self.max_tot_mass * 1.001:
            return 1
        if m_tot * 1.001 < self.min_tot_mass:
            return 1

        # Eta test
        eta = mass1 * mass2 / (m_tot * m_tot)
        if eta > self.max_eta * 1.001:
            return 1
        if eta * 1.001 < self.min_eta:
            return 1

        # Chirp mass test
        chirp_mass = m_tot * eta ** (3.0 / 5.0)
        if self.min_chirp_mass is not None and chirp_mass * 1.001 < self.min_chirp_mass:
            return 1
        if self.max_chirp_mass is not None and chirp_mass > self.max_chirp_mass * 1.001:
            return 1

        return 0
