#!/usr/bin/env python

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

import optparse
import textwrap
from pycbc.tmpltbank.lambda_mapping import pycbcValidOrdersHelpDescriptions

class IndentedHelpFormatterWithNL(optparse.IndentedHelpFormatter):
    """
    This class taken from 
    https://groups.google.com/forum/#!topic/comp.lang.python/bfbmtUGhW8I
    and is used to format the optparse help messages to deal with line breaking
    nicer. Specfically the pn-order help is large and looks crappy without this.
    This function is (C) Tim Chase
    """
    def format_description(self, description):
        """
        No documentation
        """
        if not description: return ""
        desc_width = self.width - self.current_indent
        indent = " "*self.current_indent
        # the above is still the same
        bits = description.split('\n')
        formatted_bits = [
            textwrap.fill(bit,
                desc_width,
                initial_indent=indent,
                subsequent_indent=indent)
            for bit in bits]
        result = "\n".join(formatted_bits) + "\n"
        return result

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
            opts = "%*s%s\n" % (self.current_indent, "", opts)
            indent_first = self.help_position
        else: # start help on same line as opts
            opts = "%*s%-*s  " % (self.current_indent, "", opt_width, opts)
            indent_first = 0
        result.append(opts)
        if option.help:
            help_text = self.expand_default(option)
            # Everything is the same up through here
            help_lines = []
            for para in help_text.split("\n"):
                help_lines.extend(textwrap.wrap(para, self.help_width))
            # Everything is the same after here
            result.append("%*s%s\n" % (
                indent_first, "", help_lines[0]))
            result.extend(["%*s%s\n" % (self.help_position, "", line)
                for line in help_lines[1:]])
        elif opts[-1] != "\n":
            result.append("\n")
        return "".join(result)

def insert_metric_calculation_options(parser):
    """
    Adds the options used to obtain a metric in the bank generation codes to an
    optparser as an OptionGroup. This should be used if you want to use these
    options in your code.
    """
    metricOpts = optparse.OptionGroup(parser,
                   "Options related to calculating the parameter space "+\
                   "metric.")
    metricOpts.add_option("", "--pn-order", action="store", type="string",\
                   default=None,\
                   help="Determines the PN order to use. REQUIRED ARGUMENT: "+\
                        "choices are: %s" %(pycbcValidOrdersHelpDescriptions))
    metricOpts.add_option("", "--f0", action="store", type="float",\
                  default=70., help="f0 for use in metric calculation." +\
                                    "If using ethinca calculation this "+\
                                    "must be equal to f-low."
                                    "OPTIONAL: default= %default")
    metricOpts.add_option("", "--f-low", action="store", type="float",\
                  default=15., help="f_low for use in metric calculation," +\
                                    "OPTIONAL: default= %default")
    metricOpts.add_option("", "--f-upper", action="store", type="float",\
                  default=2000., help="Upper cutoff frequency when "+\
                  "generating the metric. OPTIONAL: default= %default")
    metricOpts.add_option("", "--delta-f", action="store", type="float",\
                  default=0.001, help="delta_f for use in metric calculation,"+\
                                      "linear interpolation used to get this,"+\
                                      "OPTIONAL: default= %default")
    parser.add_option_group(metricOpts)

def verify_metric_calculation_options(opts, parser):
    """
    Parses the metric calculation options given and verifies that they are
    correct.

       Parameters
    ----------
    opts : optparse.Values instance
        Result of parsing the input options with OptionParser,
    parser : object
        The OptionParser instance.
    """
    if not opts.pn_order:
        parser.error("Must supply --pn-order")
    try:
        if opts.calculate_ethinca_metric:
            if opts.f_low != opts.f0:
                parser.error("If calculating ethinca --f0 must be equal to "+\
                             "--f-low.")
    except AttributeError:
        pass


def insert_mass_range_option_group(parser,nonSpin=False):
    """
    Adds the options used to specify mass ranges in the bank generation codes
    to an optparser as an OptionGroup. This should be used if you
    want to use these options in your code.
 
    Parameters
    -----------
    parser : object
        OptionParser instance.
    nonSpin : boolean, optional (default=False)
        If this is provided the spin-related options will not be added.
    """
    massOpts = optparse.OptionGroup(parser,
                   "Options related to mass and spin limits.")
    massOpts.add_option("", "--min-mass1", action="store", type="float",\
                  default=None, help="Minimum mass1 to generate bank with"+\
                                     ", mass1 *must* be larger than mass2." +\
                                      "REQUIRED ARGUMENT")
    massOpts.add_option("", "--max-mass1", action="store", type="float",\
                  default=None, help="Maximum mass1 to generate bank with."+\
                                      "REQUIRED ARGUMENT.")
    massOpts.add_option("", "--min-mass2", action="store", type="float",\
                  default=None, help="Minimum mass2 to generate bank with."+\
                                      "REQUIRED ARGUMENT.")
    massOpts.add_option("", "--max-mass2", action="store", type="float",\
                  default=None, help="Maximum mass2 to generate bank with."+\
                                      "REQUIRED ARGUMENT.")
    massOpts.add_option("", "--max-total-mass", action="store", type="float",\
                  default=None, help="Minimum total mass to generate bank "+\
                                      "with. OPTIONAL, no limit on max total "+\
                                      "mass if not provided.")
    massOpts.add_option("", "--min-total-mass", action="store", type="float",\
                  default=None, help="Maximum total mass to generate bank "+\
                                      "with. OPTIONAL, no limit on min total "+\
                                      "mass if not provided.")
    massOpts.add_option("", "--max-eta", action="store", type="float",\
                  default=None, help="Minimum symmetric mass ratio to "+\
                                      "generate bank "+\
                                      "with. OPTIONAL, no limit on maximum "+\
                                      "symmetric mass ratio if not provided.")
    massOpts.add_option("", "--min-eta", action="store", type="float",\
                  default=None, help="Maximum symmetric mass ratio to "+\
                                      "generate bank "+\
                                      "with. OPTIONAL, no limit on minimum "+\
                                      "symmetric mass ratio if not provided.")

    if nonSpin:
        parser.add_option_group(massOpts)
        return

    parser.add_option("", "--max-ns-spin-mag", action="store", type="float",\
                  default=None, help="Maximum neutron star spin magnitude."+\
                                     "Neutron stars are defined as objects "+\
                                     "with mass less than 3 solar masses."+\
                                     "REQUIRED ARGUMENT.")
    parser.add_option("", "--max-bh-spin-mag", action="store", type="float",\
                  default=None, help="Maximum black hole spin magnitude."+\
                                     "Neutron stars are defined as objects "+\
                                     "with mass less than 3 solar masses."+\
                                     "REQUIRED ARGUMENT.")
    parser.add_option("", "--nsbh-flag", action="store_true", default=False,\
                    help="Set this flag if you are generating a bank that "+\
                         "contains only systems with 1 black hole and 1 "+\
                         "neutron star. This prevents you having a bunch "+\
                         "of 3,3 systems where the 'neutron star' takes "+\
                         "the black hole's value of spin."+\
                         "OPTIONAL default= %default")
    parser.add_option_group(massOpts)

def verify_mass_range_options(opts, parser, nonSpin=False):
    """
    Parses the metric calculation options given and verifies that they are
    correct.

       Parameters
    ----------
    opts : optparse.Values instance
        Result of parsing the input options with OptionParser,
    parser : object
        The OptionParser instance.
    nonSpin : boolean, optional (default=False)
        If this is provided the spin-related options will not be checked.
    """
    if not opts.min_mass1:
        parser.error("Must supply --min-mass1")
    if not opts.min_mass2:
        parser.error("Must supply --min-mass2")
    if not opts.max_mass1:
        parser.error("Must supply --max-mass1")
    if not opts.max_mass2:
        parser.error("Must supply --max-mass2")
    # Mass1 must be the heavier!
    if opts.min_mass1 < opts.min_mass2:
        parser.error("Mass1 must be heavier the mass2.")
    if opts.max_mass1 < opts.max_mass2:
        parser.error("Mass1 must be heavier the mass2.")
    if (not opts.min_total_mass) or \
            ((opts.min_mass1 + opts.min_mass2) > opts.min_total_mass):
        opts.min_total_mass = opts.min_mass1 + opts.min_mass2
    if (not opts.max_total_mass) or \
            ((opts.max_mass1 + opts.max_mass2) < opts.max_total_mass):
        opts.max_total_mass = opts.max_mass1 + opts.max_mass2
    if opts.max_eta and opts.min_eta:
        if opts.max_eta < opts.min_eta:
            parser.error("If given --max-eta must be larger than --min-eta.")
    if nonSpin:
        return
    if not opts.max_ns_spin_mag:
        # Can ignore this if no NSs will be generated
        if opts.min_mass2 < 3:
            parser.error("Must supply --max-ns-spin-mag")
        else:
            opts.max_ns_spin_mag = opts.max_bh_spin_mag
    if not opts.max_bh_spin_mag:
        # Can ignore this if no BHs will be generated
        if opts.max_mass1 > 3:
            parser.error("Must supply --max-bh-spin-mag")
        else:
            opts.max_bh_spin_mag = opts.max_ns_spin_mag


def insert_ethinca_calculation_option_group(parser):
    """
    Adds the options used to obtain the ethinca metric. This should be used
    if you want ethinca calculation done in your code.
 
    Parameters
    -----------
    parser : object
        OptionParser instance.
    """
    ethincaGroup = optparse.OptionGroup(parser, "Ethinca metric options",
                    "These options are used in the calculation of the gamma "
                    "components for ethinca calculations if needed.")
    ethincaGroup.add_option("","--calculate-ethinca-metric", \
                        action="store_true",\
                        default=False, help="If given the ethinca metric "+\
                        "will be calculated and stored in the Gamma entried "+\
                        "in the sngl_inspiral table. "+\
                        "OPTIONAL: default= %default")
    ethincaGroup.add_option("","--ethinca-calc-density", action="store",\
                        default=10, help="The ethinca metric is calculated "+\
                        "using a given value for f_max. Tmpltbank uses the "+\
                        "ISCO frequency for every template and recomputes "+\
                        "the metric components every time it does this. "+\
                        "This code generates the metric for discrete values "+\
                        "f_max and uses the closest one to each template's "+\
                        "ISCO frequency. This value sets the spacing between "+\
                        "these discrete values of frequency cutoff. "+\
                        "OPTIONAL: default= %default")
    parser.add_option_group(ethincaGroup)

def verify_ethinca_calculation_options(opts, parser):
    """
    Parses the ethinca calculation options given and verifies that they are
    correct.

       Parameters
    ----------
    opts : optparse.Values instance
        Result of parsing the input options with OptionParser,
    parser : object
        The OptionParser instance.
    nonSpin : boolean, optional (default=False)
        If this is provided the spin-related options will not be checked.
    """
    # Currently no checks here.
    pass


# FIXME: This will be replaced by something in pycbc.data (or whatever), when
# Alex writes this.
def insert_data_reading_options(parser):
    """
    This function is only here while we wait for Alex to write it properly in
    pycbc.data, or whatever the module is called.
    """
    dataReadingGroup=optparse.OptionGroup(parser, "Options for obtaining h(t)",
                  "These options are used for generating h(t) either by "
                  "reading from file or by generating it. This is only needed "
                  "if the PSD is to be estimated from the data, ie. if the "
                  "--psd-estimation option is given.")

    dataReadingGroup.add_option("--gps-start-time", \
                            help="The gps start time of the data", type=int)
    dataReadingGroup.add_option("--gps-end-time", \
                            help="The gps end time of the data", type=int)
    dataReadingGroup.add_option("--strain-high-pass", type=float, \
                            help="High pass frequency")
    dataReadingGroup.add_option("--pad-data", \
              help="Extra padding to remove highpass corruption (s)", type=int)
    dataReadingGroup.add_option("--sample-rate", type=int, \
                            help="The sample rate to use for h(t) generation.")
    dataReadingGroup.add_option("--frame-cache", type=str, \
                            help="Cache file containing the frame locations.")
    dataReadingGroup.add_option("--channel-name", type=str, \
                   help="The channel containing the gravitational strain data")
    parser.add_option_group(dataReadingGroup)

def verify_data_reading_options(opts, parser):
    """
    This function is only here while we wait for Alex to write it properly in
    pycbc.data, or whatever the module is called.
    """
    # This is a placeholder, does nothing yet.
    pass

