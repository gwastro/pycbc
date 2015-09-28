# Copyright (C) 2013  Ian Harry
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
"""
This module provides a wrapper to the ConfigParser utilities for pycbc
workflow construction. This module is described in the page here:
https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/ahope/initialization_inifile.html
"""

import os
import re
import distutils.spawn
import ConfigParser
import glue.pipeline
import pycbc.workflow

def add_workflow_command_line_group(parser):
    """
    The standard way of initializing a ConfigParser object in workflow will be
    to do it from the command line. This is done by giving a

    --local-config-files filea.ini fileb.ini filec.ini

    command. You can also set config file override commands on the command
    line. This will be most useful when setting (for example) start and
    end times, or active ifos. This is done by

    --config-overrides section1:option1:value1 section2:option2:value2 ...

    This can also be given as

    --config-overrides section1:option1

    where the value will be left as ''.

    This function returns an argparse OptionGroup to ensure these options are
    parsed correctly and can then be sent directly to initialize an
    WorkflowConfigParser.

    Parameters
    -----------
    parser : argparse.ArgumentParser instance
        The initialized argparse instance to add the workflow option group to.
    """
    workflowArgs = parser.add_argument_group('workflow',
                                          'Options needed for workflow setup.')
    workflowArgs.add_argument("--config-files", nargs="+", action='store',
                           metavar="CONFIGFILE",
                           help="List of config files to be used in "
                                "analysis.")
    workflowArgs.add_argument("--config-overrides", nargs="*", action='store',
                           metavar="SECTION:OPTION:VALUE",
                           help="List of section,option,value combinations to add into the configuration file. Normally the gps start and end times might be provided this way, and user specific locations (ie. output directories). This can also be provided as SECTION:OPTION or SECTION:OPTION: both of which indicate that the corresponding value is left blank.")


class WorkflowConfigParser(glue.pipeline.DeepCopyableConfigParser):
    """
    This is a sub-class of glue.pipeline.DeepCopyableConfigParser, which lets
    us add a few additional helper features that are useful in workflows.
    """
    def __init__(self, configFiles=[], overrideTuples=[], parsedFilePath=None):
        """
        Initialize an WorkflowConfigParser. This reads the input configuration
        files, overrides values if necessary and performs the interpolation.
        See https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/ahope/initialization_inifile.html
        
        Parameters
        -----------
        configFiles : Path to .ini file, or list of paths
            The file(s) to be read in and parsed.
        overrideTuples : List of (section, option, value) tuples
            Add the (section, option, value) triplets provided
            in this list to the provided .ini file(s). If the section, option
            pair is already present, it will be overwritten.
        parsedFilePath : Path, optional (default=None)
            If given, write the parsed .ini file back to disk at this location.

        Returns
        --------
        WorkflowConfigParser
            Initialized WorkflowConfigParser instance.
        """
        glue.pipeline.DeepCopyableConfigParser.__init__(self)
        
        # Enable case sensitive options
        self.optionxform = str

        # The full path needs to be specified to avoid circular imports between
        # this file and pycbc.workflow.core
        configFiles = [pycbc.workflow.core.resolve_url(cFile) for cFile in configFiles]

        self.read_ini_file(configFiles)

        # Do overrides first
        for override in overrideTuples:
            if len(override) not in [2,3]:
                print override
                errMsg = "Overrides must be tuples of length 2 or 3."
                errMsg = "Got %s." %(str(override))
                raise ValueError(errMsg)
            section = override[0]
            option = override[1]
            value = ''
            if len(override) == 3:
                value = override[2]
            # Check for section existence, create if needed
            if not self.has_section(section):
                self.add_section(section)
            self.set(section, option, value)

        # Replace exe macros with full paths
        self.perform_exe_expansion()

        # Split sections like [inspiral&tmplt] into [inspiral] and [tmplt]
        self.split_multi_sections()

        # Check for any substitutions that can be made
        # FIXME: The python 3 version of ConfigParser can do this automatically
        # move over to that if it can be backported to python2.X.
        # We use the same formatting as the new configparser module when doing
        # ExtendedInterpolation
        # This is described at
        # http://docs.python.org/3.4/library/configparser.html
        self.perform_extended_interpolation()

        # Check for duplicate options in sub-sections
        self.sanity_check_subsections()

        # Dump parsed .ini file if needed
        if parsedFilePath:
            fp = open(parsedFilePath,'w')
            cp.write(fp)
            fp.close()


    @classmethod
    def from_args(cls, args):
        """
        Initialize a WorkflowConfigParser instance using the command line values
        parsed in args. args must contain the values provided by the
        workflow_command_line_group() function. If you are not using the standard
        workflow command line interface, you should probably initialize directly
        using __init__()

        Parameters
        -----------
        args : argparse.ArgumentParser
            The command line arguments parsed by argparse
        """
        # Identify the config files
        confFiles = []

        # files and URLs to resolve
        if args.config_files:
            confFiles += args.config_files
        
        # Identify the overrides
        confOverrides = args.config_overrides or []
        # and parse them
        parsedOverrides = []
        for override in confOverrides:
            splitOverride = override.split(":")
            if len(splitOverride) == 3:
                parsedOverrides.append(tuple(splitOverride))
            elif len(splitOverride) == 2:
                parsedOverrides.append(tuple(splitOverride + [""]))
            elif len(splitOverride) > 3:
                # Cannot have colons in either section name or variable name
                # but the value may contain colons
                rec_value = ':'.join(splitOverride[2:])
                parsedOverrides.append(tuple(splitOverride[:2] + [rec_value]))
            else:
                raise ValueError(
                    "Overrides must be of format section:option:value "
                    "or section:option. Cannot parse %s." % override)

        return cls(confFiles, parsedOverrides) 


    def read_ini_file(self, cpFile):
        """
        Read a .ini file and return it as a ConfigParser class.
        This function does none of the parsing/combining of sections. It simply
        reads the file and returns it unedited

        Stub awaiting more functionality - see configparser_test.py

        Parameters
        ----------
        cpFile : Path to .ini file, or list of paths
            The path(s) to a .ini file to be read in

        Returns
        -------
        cp : ConfigParser
            The ConfigParser class containing the read in .ini file
        """
        # Read the file
        self.read(cpFile)


    def perform_exe_expansion(self):
        """
        This function will look through the executables section of the
        ConfigParser object and replace any values using macros with full paths.

        For any values that look like 

        ${which:lalapps_tmpltbank}

        will be replaced with the equivalent of which(lalapps_tmpltbank)

        Otherwise values will be unchanged.
        """
        # Only works on executables section
        if self.has_section('executables'):
            for option, value in self.items('executables'):
                # Check the value
                newStr = self.interpolate_exe(value)
                if newStr != value:
                    self.set('executables', option, newStr)


    def interpolate_exe(self, testString):
        """
        Replace testString with a path to an executable based on the format.

        If this looks like

        ${which:lalapps_tmpltbank}
 
        it will return the equivalent of which(lalapps_tmpltbank)

        Otherwise it will return an unchanged string.

        Parameters
        -----------
        testString : string
            The input string

        Returns
        --------
        newString : string
            The output string.
        """
        # First check if any interpolation is needed and abort if not
        testString = testString.strip()
        if not (testString.startswith('${') and testString.endswith('}')):
            return testString

        # This may not be an exe interpolation, so even if it has ${XXX} form
        # I may not have to do anything
        newString = testString

        # Strip the ${ and }
        testString = testString[2:-1]

        testList = testString.split(':')

        # Maybe we can add a few different possibilities for substitution
        if len(testList) == 2:
            if testList[0] == 'which':
                newString = distutils.spawn.find_executable(testList[1])
                if not newString:
                    errMsg = "Cannot find exe %s in your path " %(testList[1])
                    errMsg += "and you specified ${which:%s}." %(testList[1])
                    raise ValueError(errMsg)

        return newString

    def get_subsections(self, section_name):
        """ Return a list of subsections for the given section name
        """
        subsections = [sec for sec in  self.sections() if sec.startswith(section_name + '-')]   

        for sec in subsections:
            sp = sec.split('-')
            # This is unusual, but a format [section-subsection-tag] is okay. Just
            # check that [section-subsection] section exists. If not it is possible
            # the user is trying to use an subsection name with '-' in it
            if (len(sp) > 1) and not self.has_section('%s-%s' % (sp[0], sp[1])):
                 raise ValueError( "Workflow uses the '-' as a delimiter so this is"
                     "interpreted as section-subsection-tag. "
                     "While checking section %s, no section with "
                     "name %s-%s was found. " 
                     "If you did not intend to use tags in an "
                     "'advanced user' manner, or do not understand what "
                     "this means, don't use dashes in section "
                     "names. So [injection-nsbhinj] is good. "
                     "[injection-nsbh-inj] is not." % (sec, sp[0], sp[1]))
        
        if len(subsections) > 0:
            return [sec.split('-')[1] for sec in subsections]
        elif self.has_section(section_name):
            return ['']
        else:
            return []

    def perform_extended_interpolation(self):
        """
        Filter through an ini file and replace all examples of 
        ExtendedInterpolation formatting with the exact value. For values like
        ${example} this is replaced with the value that corresponds to the
        option called example ***in the same section***

        For values like ${common|example} this is replaced with the value that
        corresponds to the option example in the section [common]. Note that
        in the python3 config parser this is ${common:example} but python2.7
        interprets the : the same as a = and this breaks things

        Nested interpolation is not supported here.
        """

        # Do not allow any interpolation of the section names
        for section in self.sections():
            for option,value in self.items(section):
                # Check the option name
                newStr = self.interpolate_string(option, section)
                if newStr != option:
                    self.set(section,newStr,value)
                    self.remove_option(section,option)
                # Check the value
                newStr = self.interpolate_string(value, section)
                if newStr != value:
                    self.set(section,option,newStr)


    def interpolate_string(self, testString, section):
        """
        Take a string and replace all example of ExtendedInterpolation
        formatting within the string with the exact value. 

        For values like ${example} this is replaced with the value that
        corresponds to the option called example ***in the same section***

        For values like ${common|example} this is replaced with the value that
        corresponds to the option example in the section [common]. Note that
        in the python3 config parser this is ${common:example} but python2.7
        interprets the : the same as a = and this breaks things

        Nested interpolation is not supported here.

        Parameters
        ----------
        testString : String
            The string to parse and interpolate
        section : String
            The current section of the ConfigParser object

        Returns
        ----------
        testString : String
            Interpolated string
        """

        # First check if any interpolation is needed and abort if not
        reObj = re.search(r"\$\{.*?\}", testString)
        while reObj:
            # Not really sure how this works, but this will obtain the first
            # instance of a string contained within ${....}
            repString = (reObj).group(0)[2:-1]
            # Need to test which of the two formats we have
            splitString = repString.split('|')
            if len(splitString) == 1:
                try:
                    testString = testString.replace('${'+repString+'}',\
                                            self.get(section,splitString[0]))
                except ConfigParser.NoOptionError:
                    print "Substitution failed"
                    raise
            if len(splitString) == 2:
                try:
                    testString = testString.replace('${'+repString+'}',\
                                       self.get(splitString[0],splitString[1]))
                except ConfigParser.NoOptionError:
                    print "Substitution failed"
                    raise
            reObj = re.search(r"\$\{.*?\}", testString)

        return testString


    def split_multi_sections(self):
        """
        Parse through the WorkflowConfigParser instance and splits any sections
        labelled with an "&" sign (for e.g. [inspiral&tmpltbank]) into
        [inspiral] and [tmpltbank] sections. If these individual sections
        already exist they  will be appended to. If an option exists in both the
        [inspiral] and [inspiral&tmpltbank] sections an error will be thrown
        """
        # Begin by looping over all sections
        for section in self.sections():
            # Only continue if section needs splitting
            if '&' not in section:
                continue
            # Get list of section names to add these options to
            splitSections = section.split('&')
            for newSec in splitSections:
                # Add sections if they don't already exist
                if not self.has_section(newSec):
                    self.add_section(newSec)
                self.add_options_to_section(newSec, self.items(section))
            self.remove_section(section)


    def add_options_to_section(self ,section, items, overwrite_options=False):
        """
        Add a set of options and values to a section of a ConfigParser object.
        Will throw an error if any of the options being added already exist,
        this behaviour can be overridden if desired

        Parameters
        ----------
        section : string
            The name of the section to add options+values to
        items : list of tuples
            Each tuple contains (at [0]) the option and (at [1]) the value to
            add to the section of the ini file
        overwrite_options : Boolean, optional
            By default this function will throw a ValueError if an option exists
            in both the original section in the ConfigParser *and* in the
            provided items.
            This will override so that the options+values given in items
            will replace the original values if the value is set to True.
            Default = True 
        """
        # Sanity checking
        if not self.has_section(section):
            raise ValueError('Section %s not present in ConfigParser.' \
                             %(section,))

        # Check for duplicate options first
        for option,value in items:
            if not overwrite_options:
                if option in self.options(section):
                    raise ValueError('Option exists in both original ' + \
                                  'ConfigParser section [%s] and ' %(section,) + \
                                  'input list: %s' %(option,))
            self.set(section,option,value)


    def sanity_check_subsections(self):
        """
        This function goes through the ConfigParset and checks that any options
        given in the [SECTION_NAME] section are not also given in any 
        [SECTION_NAME-SUBSECTION] sections.

        """
        # Loop over the sections in the ini file
        for section in self.sections():
            # Loop over the sections again
            for section2 in self.sections():
                # Check if any are subsections of section
                if section2.startswith(section + '-'):
                    # Check for duplicate options whenever this exists
                    self.check_duplicate_options(section, section2,
                                                 raise_error=True)


    def check_duplicate_options(self, section1, section2, raise_error=False):
        """
        Check for duplicate options in two sections, section1 and section2.
        Will return a list of the duplicate options.
    
        Parameters
        ----------
        section1 : string
            The name of the first section to compare
        section2 : string
            The name of the second section to compare
        raise_error : Boolean, optional (default=False)
            If True, raise an error if duplicates are present.

        Returns
        ----------
        duplicates : List
            List of duplicate options
        """
        # Sanity checking
        if not self.has_section(section1):
            raise ValueError('Section %s not present in ConfigParser.'\
                             %(section1,) )
        if not self.has_section(section2):
            raise ValueError('Section %s not present in ConfigParser.'\
                             %(section2,) )

        items1 = self.options(section1)
        items2 = self.options(section2)

        # The list comprehension here creates a list of all duplicate items
        duplicates = [x for x in items1 if x in items2]

        if duplicates and raise_error:
            raise ValueError('The following options appear in both section ' +\
                             '%s and %s: %s' \
                             %(section1,section2,' '.join(duplicates)))

        return duplicates


    def get_opt_tag(self, section, option, tag):
        """
        Convenience function accessing get_opt_tags() for a single tag: see 
        documentation for that function.
        NB calling get_opt_tags() directly is preferred for simplicity.

        Parameters
        -----------
        self : ConfigParser object
            The ConfigParser object (automatically passed when this is appended
            to the ConfigParser class)
        section : string
            The section of the ConfigParser object to read
        option : string
            The ConfigParser option to look for
        tag : string
            The name of the subsection to look in, if not found in [section]
 
        Returns
        --------
        string
            The value of the options being searched for
        """
        return self.get_opt_tags(section, option, [tag])


    def get_opt_tags(self, section, option, tags):
        """
        Supplement to ConfigParser.ConfigParser.get(). This will search for an
        option in [section] and if it doesn't find it will also try in
        [section-tag] for every value of tag in tags.
        Will raise a ConfigParser.Error if it cannot find a value.

        Parameters
        -----------
        self : ConfigParser object
            The ConfigParser object (automatically passed when this is appended
            to the ConfigParser class)
        section : string
            The section of the ConfigParser object to read
        option : string
            The ConfigParser option to look for
        tags : list of strings
            The name of subsections to look in, if not found in [section]
 
        Returns
        --------
        string
            The value of the options being searched for
        """
        # Need lower case tag name; also exclude cases with tag=None
        if tags:
            tags = [tag.lower() for tag in tags if tag is not None]

        try:
            return self.get(section, option)
        except ConfigParser.Error:
            errString = "No option '%s' in section [%s] " %(option,section)
            if not tags:
                raise ConfigParser.Error(errString + ".")
            returnVals = []
            sectionList = ["%s-%s" %(section, tag) for tag in tags]
            for tag in tags:
                if self.has_section('%s-%s' %(section, tag)):
                    if self.has_option('%s-%s' %(section, tag), option):
                        returnVals.append(self.get('%s-%s' %(section, tag),
                                                    option))
            if not returnVals:
                errString += "or in sections [%s]." %("] [".join(sectionList))
                raise ConfigParser.Error(errString)
            if len(returnVals) > 1:
                errString += "and multiple entries found in sections [%s]."\
                              %("] [".join(sectionList))
                raise ConfigParser.Error(errString)
            return returnVals[0]


    def has_option_tag(self, section, option, tag):
        """
        Convenience function accessing has_option_tags() for a single tag: see
        documentation for that function.
        NB calling has_option_tags() directly is preferred for simplicity.

        Parameters
        -----------
        self : ConfigParser object
            The ConfigParser object (automatically passed when this is appended
            to the ConfigParser class)
        section : string
            The section of the ConfigParser object to read
        option : string
            The ConfigParser option to look for
        tag : string
            The name of the subsection to look in, if not found in [section]
 
        Returns
        --------
        Boolean
            Is the option in the section or [section-tag]
        """
        return self.has_option_tags(section, option, [tag])


    def has_option_tags(self, section, option, tags):
        """
        Supplement to ConfigParser.ConfigParser.has_option().
        This will search for an option in [section] and if it doesn't find it
        will also try in [section-tag] for each value in tags.
        Returns True if the option is found and false if not.

        Parameters
        -----------
        self : ConfigParser object
            The ConfigParser object (automatically passed when this is appended
            to the ConfigParser class)
        section : string
            The section of the ConfigParser object to read
        option : string
            The ConfigParser option to look for
        tags : list of strings
            The names of the subsection to look in, if not found in [section]
 
        Returns
        --------
        Boolean
            Is the option in the section or [section-tag] (for tag in tags)
        """
        try:
            self.get_opt_tags(section, option, tags)
            return True
        except ConfigParser.Error:
            return False
