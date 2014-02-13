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
from pycbc.tmpltbank.lambda_mapping import *

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
                   help="Determines the PN order to use. Note that if you "+\
                        "placing a bank of non-spinning templates, any "+\
                        "spin-related terms in the metric will always "+\
                        "be zero. REQUIRED ARGUMENT: "+\
                        "choices are: %s" %(pycbcValidOrdersHelpDescriptions))
    metricOpts.add_option("", "--f0", action="store", type="float",\
                  default=70.,\
                  help="f0 is used as a dynamic rescaling factor when "+\
                       "calculating the integrals used in metric "+\
                       "construction. IE. instead of integrating F(f) we "+\
                       "integrate F(f/f0) and then remove f0 after the fact."+\
                       "The default option should be fine here for most "+\
                       "applications. OPTIONAL: default= %default."+\
                       "WARNING: If using ethinca calculation this must be "+\
                       "equal to f-low. UNITS=Hz")
    metricOpts.add_option("", "--f-low", action="store", type="float",\
                  default=None, help="The lower frequency cutoff to use "+\
                       "when computing the components of the parameter "+\
                       "space metric. REQUIRED ARGUMENT. UNITS=Hz")
    metricOpts.add_option("", "--f-upper", action="store", type="float",\
                  default=None, help="The upper frequency cutoff to use "+\
                       "when computing the components of the parameter "+\
                       "space metric. REQUIRED ARGUMENT. UNITS=Hz")
    metricOpts.add_option("", "--delta-f", action="store", type="float",\
                  default=None, help="The frequency spacing to use when "+\
                       "computing the components of the parameter space "+\
                       "metric. IE. the various integrals should be "+\
                       "\int F(f) df, however we only approximate this as "+\
                       "\sum F(f) delta_f. This sets delta_f. "+\
                       "REQUIRED ARGUMENT. UNITS=Hz")
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
    if not opts.f_low:
        parser.error("Must supply --f-low")
    if not opts.f_upper:
        parser.error("Must supply --f-upper")
    if not opts.delta_f:
        parser.error("Must supply --delta-f")
    try:
        if opts.calculate_ethinca_metric:
            if opts.f_low != opts.f0:
                parser.error("If calculating ethinca --f0 must be equal to "+\
                             "--f-low.")
    except AttributeError:
        pass


class metricParameters:
    """
    This class holds all of the options that are parsed in the function
    insert_metric_calculation_options
    and all products produced using these options. It can also be initialized
    from the __init__ function providing directly the options normally
    provided on the command line
    """
    _psd = None
    _metric = None
    _evals = None
    _evecs = None
    _evecsCV = None
    def __init__(self, pnOrder, fLow, fUpper, deltaF, f0=70):
        """
        Initialize an instance of the metricParameters by providing all
        options directly. See the help message associated with any code
        that uses the metric options for more details of how to set each of
        these. For e.g. pycbc_aligned_stoch_bank --help
        """
        self.pnOrder=pnOrder
        self.fLow=fLow
        self.fUpper=fUpper
        self.deltaF=deltaF
        self.f0=f0
        self._moments=None

    @classmethod
    def from_optparse(cls, opts):
        """
        Initialize an instance of the metricParameters class from an
        optparse.OptionParser instance. This assumes that
        insert_metric_calculation_options
        and
        verify_metric_calculation_options
        have already been called before initializing the class.
        """
        return cls(opts.pn_order, opts.f_low, opts.f_upper, opts.delta_f,\
                   f0=opts.f0)

    @property
    def psd(self):
        """
        A pyCBC FrequencySeries holding the appropriate PSD.
        Return the PSD used in the metric calculation.
        """
        if not self._psd:
            errMsg = "The PSD has not been set in the metricParameters "
            errMsg += "instance."
            raise ValueError(errMsg)
        return self._psd

    @psd.setter
    def psd(self, inPsd):
        self._psd = inPsd

    @property
    def moments(self):
        """
        Moments structure
        This contains the result of all the integrals used in computing the
        metrics above. It can be used for the ethinca components calculation,
        or other similar calculations. This is composed of two compound
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
    def moments(self, inMoments):
        self._moments=inMoments

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
        if not self._evals:
            errMsg = "The metric eigenvalues have not been set in the "
            errMsg += "metricParameters instance."
            raise ValueError(errMsg)
        return self._evals

    @evals.setter
    def evals(self, inEvals):
        self._evals = inEvals

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
        if not self._evecs:
            errMsg = "The metric eigenvectors have not been set in the "
            errMsg += "metricParameters instance."
            raise ValueError(errMsg)
        return self._evecs

    @evecs.setter
    def evecs(self, inEvecs):
        self._evecs = inEvecs

    @property
    def metric(self):
        """
        The eigenvectors of the parameter space.
        This is a Dictionary of numpy.matrix
        Each entry in the dictionary is as described under evals.
        Each numpy.matrix contains the metric of the parameter space in the
        Lambda_i coordinate system.
        """
        if not self._metric:
            errMsg = "The metric eigenvectors have not been set in the "
            errMsg += "metricParameters instance."
            raise ValueError(errMsg)
        return self._metric

    @metric.setter
    def metric(self, inMetric):
        self._metric = inMetric

    @property
    def evecsCV(self):
        """
        The eigenvectors of the principal directions of the mu space.
        This is a Dictionary of numpy.matrix
        Each entry in the dictionary is as described under evals.
        Each numpy.matrix contains the eigenvectors which, with the eigenvalues
        in evals, are needed to rotate the
        coordinate system to one in which the metric is the identity matrix.
        """
        if not self._evecsCV:
            errMsg = "The covariance eigenvectors have not been set in the "
            errMsg += "metricParameters instance."
            raise ValueError(errMsg)
        return self._evecsCV

    @evecsCV.setter
    def evecsCV(self, inEvecs):
        self._evecsCV = inEvecs



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
                                      "REQUIRED ARGUMENT. UNITS=Solar mass.")
    massOpts.add_option("", "--max-mass1", action="store", type="float",\
                  default=None, help="Maximum mass1 to generate bank with."+\
                                      "REQUIRED ARGUMENT. UNITS=Solar mass.")
    massOpts.add_option("", "--min-mass2", action="store", type="float",\
                  default=None, help="Minimum mass2 to generate bank with."+\
                                      "REQUIRED ARGUMENT. UNITS=Solar mass.")
    massOpts.add_option("", "--max-mass2", action="store", type="float",\
                  default=None, help="Maximum mass2 to generate bank with."+\
                                      "REQUIRED ARGUMENT. UNITS=Solar mass.")
    massOpts.add_option("", "--max-total-mass", action="store", type="float",\
                  default=None, help="Minimum total mass to generate bank "+\
                                      "with. OPTIONAL, no limit on max total "+\
                                      "mass if not provided. UNITS=Solar mass.")
    massOpts.add_option("", "--min-total-mass", action="store", type="float",\
                  default=None, help="Maximum total mass to generate bank "+\
                                      "with. OPTIONAL, no limit on min total "+\
                                      "mass if not provided. UNITS=Solar mass.")
    massOpts.add_option("", "--max-eta", action="store", type="float",\
                  default=0.25, help="Minimum symmetric mass ratio to "+\
                                      "generate bank "+\
                                      "with. OPTIONAL, no limit on maximum "+\
                                      "symmetric mass ratio if not provided."+\
                                      "UNITS=Solar mass.")
    massOpts.add_option("", "--min-eta", action="store", type="float",\
                  default=0., help="Maximum symmetric mass ratio to "+\
                                      "generate bank "+\
                                      "with. OPTIONAL, no limit on minimum "+\
                                      "symmetric mass ratio if not provided."+\
                                      "UNITS=Solar mass.")

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
    if opts.max_ns_spin_mag is None:
        # Can ignore this if no NSs will be generated
        if opts.min_mass2 < 3:
            parser.error("Must supply --max-ns-spin-mag")
        else:
            opts.max_ns_spin_mag = opts.max_bh_spin_mag
    if opts.max_bh_spin_mag is None:
        # Can ignore this if no BHs will be generated
        if opts.max_mass1 > 3:
            parser.error("Must supply --max-bh-spin-mag")
        else:
            opts.max_bh_spin_mag = opts.max_ns_spin_mag

class massRangeParameters(object):
    """
    This class holds all of the options that are parsed in the function
    insert_mass_range_option_group
    and all products produced using these options. It can also be initialized
    from the __init__ function providing directly the options normally
    provided on the command line
    """
    def __init__(self, minMass1, maxMass1, minMass2, maxMass2,
                 maxNSSpinMag=0, maxBHSpinMag=0, maxTotMass=None,
                 minTotMass=None, maxEta=None, minEta=0, nsbhFlag=False):
        """
        Initialize an instance of the massRangeParameters by providing all
        options directly. See the help message associated with any code
        that uses the metric options for more details of how to set each of
        these. For e.g. pycbc_aligned_stoch_bank --help
        """
        self.minMass1=minMass1
        self.maxMass1=maxMass1
        self.minMass2=minMass2
        self.maxMass2=maxMass2
        self.maxNSSpinMag=maxNSSpinMag
        self.maxBHSpinMag=maxBHSpinMag
        self.minTotMass = minMass1 + minMass2
        if minTotMass and (minTotMass > self.minTotMass):
            self.minTotMass = minTotMass
        self.maxTotMass = maxMass1 + maxMass2
        if maxTotMass and (maxTotMass < self.maxTotMass):
            self.maxTotMass = maxTotMass
        self.maxTotMass=maxTotMass
        self.minTotMass=minTotMass
        if maxEta:
            self.maxEta=maxEta
        else:
            self.maxEta=0.25
        self.minEta=minEta
        # FIXME: Check you have NSBHs if this is set.
        # In fact, why not set automatically if you *have* only NSBHs?
        self.nsbhFlag=nsbhFlag

        # FIXME: This may be inaccurate if Eta limits are given
        # This will not cause any problems, but maybe could be fixed.
        self.minCompMass = self.minMass2
        self.maxCompMass = self.maxMass1

        # WARNING: We expect mass1 > mass2 ALWAYS
        # Check input:
        if (minMass2 > minMass1) or (maxMass2 > maxMass1):
            errMsg = "Mass1 must be larger than mass2. Check input options."
            raise ValueError(errMsg)

        if (minMass2 > maxMass2) or (minMass1 > maxMass1):
            errMsg = "Minimum masses cannot be larger than maximum masses."
            errMsg += "Check input options."
            raise ValueError(errMsg)


    @classmethod
    def from_optparse(cls, opts, nonSpin=False):
        """
        Initialize an instance of the massRangeParameters class from an
        optparse.OptionParser instance. This assumes that
        insert_mass_range_option_group
        and
        verify_mass_range_options
        have already been called before initializing the class.
        """
        if nonSpin:
            return cls(opts.min_mass1, opts.max_mass1, opts.min_mass2,\
                       opts.max_mass2, maxTotMass=opts.max_total_mass,\
                       minTotMass=opts.min_total_mass, maxEta=opts.max_eta,\
                       minEta=opts.min_eta)
        else:
            return cls(opts.min_mass1, opts.max_mass1, opts.min_mass2,\
                       opts.max_mass2, maxTotMass=opts.max_total_mass,\
                       minTotMass=opts.min_total_mass, maxEta=opts.max_eta,\
                       minEta=opts.min_eta, maxNSSpinMag=opts.max_ns_spin_mag,\
                       maxBHSpinMag=opts.max_bh_spin_mag, \
                       nsbhFlag=opts.nsbh_flag)

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
                        "OPTIONAL: default= %default, UNITS=Hz.")
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
