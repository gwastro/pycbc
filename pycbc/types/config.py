# Copyright (C) 2013,2017,2021 Ian Harry, Duncan Brown, Alex Nitz
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
This module provides a wrapper to the ConfigParser utilities for pycbc.
This module is described in the page here:
"""
import re
import itertools
import logging
from io import StringIO
import configparser as ConfigParser


class DeepCopyableConfigParser(ConfigParser.ConfigParser):
    """
    The standard SafeConfigParser no longer supports deepcopy() as of python
    2.7 (see http://bugs.python.org/issue16058). This subclass restores that
    functionality.
    """

    def __deepcopy__(self, memo):
        # http://stackoverflow.com/questions/23416370
        # /manually-building-a-deep-copy-of-a-configparser-in-python-2-7
        config_string = StringIO()
        self.write(config_string)
        config_string.seek(0)
        new_config = self.__class__()
        new_config.readfp(config_string)
        return new_config


class InterpolatingConfigParser(DeepCopyableConfigParser):
    """
    This is a sub-class of DeepCopyableConfigParser, which lets
    us add a few additional helper features that are useful in workflows.
    """

    def __init__(
        self,
        configFiles=None,
        overrideTuples=None,
        parsedFilePath=None,
        deleteTuples=None,
        skip_extended=False,
        sanitize_newline=True,
    ):
        """
         Initialize an InterpolatingConfigParser. This reads the input configuration
         files, overrides values if necessary and performs the interpolation.

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
         deleteTuples : List of (section, option) tuples
             Delete the (section, option) pairs provided
             in this list from provided .ini file(s). If the section only
             is provided, the entire section will be deleted.

        Returns
         --------
         InterpolatingConfigParser
             Initialized InterpolatingConfigParser instance.
        """
        if configFiles is None:
            configFiles = []
        if overrideTuples is None:
            overrideTuples = []
        if deleteTuples is None:
            deleteTuples = []
        DeepCopyableConfigParser.__init__(self)

        # Enable case sensitive options
        self.optionxform = str

        self.read_ini_file(configFiles)

        # Split sections like [inspiral&tmplt] into [inspiral] and [tmplt]
        self.split_multi_sections()

        # Populate shared options from the [sharedoptions] section
        self.populate_shared_sections()

        # Do deletes from command line
        for delete in deleteTuples:
            if len(delete) == 1:
                if self.remove_section(delete[0]) is False:
                    raise ValueError(
                        "Cannot delete section %s, "
                        "no such section in configuration." % delete
                    )

                logging.info(
                    "Deleting section %s from configuration", delete[0]
                )
            elif len(delete) == 2:
                if self.remove_option(delete[0], delete[1]) is False:
                    raise ValueError(
                        "Cannot delete option %s from section %s,"
                        " no such option in configuration." % delete
                    )

                logging.info(
                    "Deleting option %s from section %s in " "configuration",
                    delete[1],
                    delete[0],
                )
            else:
                raise ValueError(
                    "Deletes must be tuples of length 1 or 2. "
                    "Got %s." % str(delete)
                )

        # Do overrides from command line
        for override in overrideTuples:
            if len(override) not in [2, 3]:
                errmsg = "Overrides must be tuples of length 2 or 3."
                errmsg = "Got %s." % (str(override))
                raise ValueError(errmsg)
            section = override[0]
            option = override[1]
            value = ""
            if len(override) == 3:
                value = override[2]
            # Check for section existence, create if needed
            if not self.has_section(section):
                self.add_section(section)
            self.set(section, option, value)
            logging.info(
                "Overriding section %s option %s with value %s "
                "in configuration.",
                section,
                option,
                value,
            )

        # Check for any substitutions that can be made
        if not skip_extended:
            self.perform_extended_interpolation()

        # replace newlines in input with spaces
        # this enables command line conversion compatibility
        if sanitize_newline:
            self.sanitize_newline()

        # Check for duplicate options in sub-sections
        self.sanity_check_subsections()

        # Dump parsed .ini file if needed
        if parsedFilePath:
            fp = open(parsedFilePath, "w")
            self.write(fp)
            fp.close()

    @classmethod
    def from_cli(cls, opts):
        """Initialize the config parser using options parsed from the command
        line.

        The parsed options ``opts`` must include options provided by
        :py:func:`add_workflow_command_line_group`.

        Parameters
        -----------
        opts : argparse.ArgumentParser
            The command line arguments parsed by argparse
        """
        # read configuration file
        logging.info("Reading configuration file")
        if opts.config_overrides is not None:
            overrides = [
                tuple(override.split(":", 2))
                for override in opts.config_overrides
            ]
        else:
            overrides = None
        if opts.config_delete is not None:
            deletes = [
                tuple(delete.split(":")) for delete in opts.config_delete
            ]
        else:
            deletes = None
        return cls(opts.config_files, overrides, deleteTuples=deletes)

    def read_ini_file(self, fpath):
        """
        Read a .ini file and return it as a ConfigParser class.
        This function does none of the parsing/combining of sections. It simply
        reads the file and returns it unedited

        Stub awaiting more functionality - see configparser_test.py

        Parameters
        ----------
        fpath : Path to .ini file, or list of paths
            The path(s) to a .ini file to be read in

        Returns
        -------
        cp : ConfigParser
            The ConfigParser class containing the read in .ini file
        """
        # Read the file
        self.read(fpath)

    def get_subsections(self, section_name):
        """Return a list of subsections for the given section name"""
        # Keep only subsection names
        subsections = [
            sec[len(section_name) + 1:]
            for sec in self.sections()
            if sec.startswith(section_name + "-")
        ]

        for sec in subsections:
            sp = sec.split("-")
            # This is unusual, but a format [section-subsection-tag] is okay. Just
            # check that [section-subsection] section exists. If not it is possible
            # the user is trying to use an subsection name with '-' in it
            if (len(sp) > 1) and not self.has_section(
                "%s-%s" % (section_name, sp[0])
            ):
                raise ValueError(
                    "Workflow uses the '-' as a delimiter so "
                    "this is interpreted as section-subsection-tag. "
                    "While checking section %s, no section with "
                    "name %s-%s was found. "
                    "If you did not intend to use tags in an "
                    "'advanced user' manner, or do not understand what "
                    "this means, don't use dashes in section "
                    "names. So [injection-nsbhinj] is good. "
                    "[injection-nsbh-inj] is not." % (sec, sp[0], sp[1])
                )

        if len(subsections) > 0:
            return [sec.split("-")[0] for sec in subsections]
        elif self.has_section(section_name):
            return [""]
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
            for option, value in self.items(section):
                # Check the option name
                new_str = self.interpolate_string(option, section)
                if new_str != option:
                    self.set(section, new_str, value)
                    self.remove_option(section, option)
                # Check the value
                new_str = self.interpolate_string(value, section)
                if new_str != value:
                    self.set(section, option, new_str)

    def sanitize_newline(self):
        """
        Filter through an ini file and replace all examples of
        newlines with spaces. This is useful for command line conversion
        and allow multiline configparser inputs without added backslashes
        """

        # Do not allow any interpolation of the section names
        for section in self.sections():
            for option, value in self.items(section):
                new_value = value.replace('\n', ' ').replace('\r', ' ')
                self.set(section, option, new_value)

    def interpolate_string(self, test_string, section):
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
        test_string : String
            The string to parse and interpolate
        section : String
            The current section of the ConfigParser object

        Returns
        ----------
        test_string : String
            Interpolated string
        """

        # First check if any interpolation is needed and abort if not
        re_obj = re.search(r"\$\{.*?\}", test_string)
        while re_obj:
            # Not really sure how this works, but this will obtain the first
            # instance of a string contained within ${....}
            rep_string = (re_obj).group(0)[2:-1]
            # Need to test which of the two formats we have
            split_string = rep_string.split("|")
            if len(split_string) == 1:
                try:
                    test_string = test_string.replace(
                        "${" + rep_string + "}",
                        self.get(section, split_string[0]),
                    )
                except ConfigParser.NoOptionError:
                    print("Substitution failed")
                    raise
            if len(split_string) == 2:
                try:
                    test_string = test_string.replace(
                        "${" + rep_string + "}",
                        self.get(split_string[0], split_string[1]),
                    )
                except ConfigParser.NoOptionError:
                    print("Substitution failed")
                    raise
            re_obj = re.search(r"\$\{.*?\}", test_string)

        return test_string

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
            if "&" not in section:
                continue
            # Get list of section names to add these options to
            split_sections = section.split("&")
            for new_sec in split_sections:
                # Add sections if they don't already exist
                if not self.has_section(new_sec):
                    self.add_section(new_sec)
                self.add_options_to_section(new_sec, self.items(section))
            self.remove_section(section)

    def populate_shared_sections(self):
        """Parse the [sharedoptions] section of the ini file.

        That section should contain entries according to:

          * massparams = inspiral, tmpltbank
          * dataparams = tmpltbank

        This will result in all options in [sharedoptions-massparams] being
        copied into the [inspiral] and [tmpltbank] sections and the options
        in [sharedoptions-dataparams] being copited into [tmpltbank].
        In the case of duplicates an error will be raised.
        """
        if not self.has_section("sharedoptions"):
            # No sharedoptions, exit
            return
        for key, value in self.items("sharedoptions"):
            assert self.has_section("sharedoptions-%s" % (key))
            # Comma separated
            values = value.split(",")
            common_options = self.items("sharedoptions-%s" % (key))
            for section in values:
                if not self.has_section(section):
                    self.add_section(section)
                for arg, val in common_options:
                    if arg in self.options(section):
                        raise ValueError(
                            "Option exists in both original "
                            + "ConfigParser section [%s] and " % (section,)
                            + "sharedoptions section: %s %s"
                            % (arg, "sharedoptions-%s" % (key))
                        )
                    self.set(section, arg, val)
            self.remove_section("sharedoptions-%s" % (key))
        self.remove_section("sharedoptions")

    def add_options_to_section(self, section, items, overwrite_options=False):
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
            raise ValueError(
                "Section %s not present in ConfigParser." % (section,)
            )

        # Check for duplicate options first
        for option, value in items:
            if not overwrite_options:
                if option in self.options(section):
                    raise ValueError(
                        "Option exists in both original "
                        + "ConfigParser section [%s] and " % (section,)
                        + "input list: %s" % (option,)
                    )
            self.set(section, option, value)

    def sanity_check_subsections(self):
        """
        This function goes through the ConfigParset and checks that any options
        given in the [SECTION_NAME] section are not also given in any
        [SECTION_NAME-SUBSECTION] sections.

        """
        # Loop over the sections in the ini file
        for section in self.sections():
            # [pegasus_profile] specially is allowed to be overriden by
            # sub-sections
            if section == "pegasus_profile":
                continue

            # Loop over the sections again
            for section2 in self.sections():
                # Check if any are subsections of section
                if section2.startswith(section + "-"):
                    # Check for duplicate options whenever this exists
                    self.check_duplicate_options(
                        section, section2, raise_error=True
                    )

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
            raise ValueError(
                "Section %s not present in ConfigParser." % (section1,)
            )
        if not self.has_section(section2):
            raise ValueError(
                "Section %s not present in ConfigParser." % (section2,)
            )

        items1 = self.options(section1)
        items2 = self.options(section2)

        # The list comprehension here creates a list of all duplicate items
        duplicates = [x for x in items1 if x in items2]

        if duplicates and raise_error:
            raise ValueError(
                "The following options appear in both section "
                + "%s and %s: %s" % (section1, section2, " ".join(duplicates))
            )

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
            err_string = "No option '%s' in section [%s] " % (option, section)
            if not tags:
                raise ConfigParser.Error(err_string + ".")
            return_vals = []
            sub_section_list = []
            for sec_len in range(1, len(tags) + 1):
                for tag_permutation in itertools.permutations(tags, sec_len):
                    joined_name = "-".join(tag_permutation)
                    sub_section_list.append(joined_name)
            section_list = ["%s-%s" % (section, sb) for sb in sub_section_list]
            err_section_list = []
            for sub in sub_section_list:
                if self.has_section("%s-%s" % (section, sub)):
                    if self.has_option("%s-%s" % (section, sub), option):
                        err_section_list.append("%s-%s" % (section, sub))
                        return_vals.append(
                            self.get("%s-%s" % (section, sub), option)
                        )

            # We also want to recursively go into sections

            if not return_vals:
                err_string += "or in sections [%s]." % (
                    "] [".join(section_list)
                )
                raise ConfigParser.Error(err_string)
            if len(return_vals) > 1:
                err_string += (
                    "and multiple entries found in sections [%s]."
                    % ("] [".join(err_section_list))
                )
                raise ConfigParser.Error(err_string)
            return return_vals[0]

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
