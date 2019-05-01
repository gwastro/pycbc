# Copyright (C) 2013,2017 Ian Harry, Duncan Brown
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
import stat
import string
import shutil
import time
import logging
import requests
import distutils.spawn
import itertools
import six
from six.moves import configparser as ConfigParser
from six.moves.urllib.parse import urlparse
from six.moves import http_cookiejar as cookielib
from six.moves.http_cookiejar import (_warn_unhandled_exception,
                                      LoadError, Cookie)
import glue.pipeline
from bs4 import BeautifulSoup

def _really_load(self, f, filename, ignore_discard, ignore_expires):
    """
    This function is required to monkey patch MozillaCookieJar's _really_load
    function which does not understand the curl format cookie file created
    by ecp-cookie-init. It patches the code so that #HttpOnly_ get loaded.

    https://bugs.python.org/issue2190
    https://bugs.python.org/file37625/httponly.patch
    """
    now = time.time()

    magic = f.readline()
    if not re.search(self.magic_re, magic):
        f.close()
        raise LoadError(
            "%r does not look like a Netscape format cookies file" %
            filename)

    try:
        while 1:
            line = f.readline()
            if line == "": break

            # last field may be absent, so keep any trailing tab
            if line.endswith("\n"): line = line[:-1]

            sline = line.strip()
            # support HttpOnly cookies (as stored by curl or old Firefox).
            if sline.startswith("#HttpOnly_"):
                line = sline[10:]
            # skip comments and blank lines XXX what is $ for?
            elif (sline.startswith(("#", "$")) or sline == ""):
                continue

            domain, domain_specified, path, secure, expires, name, value = \
                    line.split("\t")
            secure = (secure == "TRUE")
            domain_specified = (domain_specified == "TRUE")
            if name == "":
                # cookies.txt regards 'Set-Cookie: foo' as a cookie
                # with no name, whereas cookielib regards it as a
                # cookie with no value.
                name = value
                value = None

            initial_dot = domain.startswith(".")
            assert domain_specified == initial_dot

            discard = False
            if expires == "":
                expires = None
                discard = True

            # assume path_specified is false
            c = Cookie(0, name, value,
                       None, False,
                       domain, domain_specified, initial_dot,
                       path, False,
                       secure,
                       expires,
                       discard,
                       None,
                       None,
                       {})
            if not ignore_discard and c.discard:
                continue
            if not ignore_expires and c.is_expired(now):
                continue
            self.set_cookie(c)

    except IOError:
        raise
    except Exception:
        _warn_unhandled_exception()
        raise LoadError("invalid Netscape format cookies file %r: %r" %
                        (filename, line))

# Now monkey patch the code
cookielib.MozillaCookieJar._really_load = _really_load # noqa

ecp_cookie_error = """The attempt to download the file at

{}

was redirected to the git.ligo.org sign-in page. This means that you likely
forgot to initialize your ECP cookie or that your LIGO.ORG credentials are
otherwise invalid. Create a valid ECP cookie for git.ligo.org by running

ecp-cookie-init LIGO.ORG https://git.ligo.org/users/auth/shibboleth/callback albert.einstein

before attempting to download files from git.ligo.org.
"""

def istext(s, text_characters=None, threshold=0.3):
    """
    Determines if the string is a set of binary data or a text file.
    This is done by checking if a large proportion of characters are > 0X7E
    (0x7F is <DEL> and unprintable) or low bit control codes. In other words
    things that you wouldn't see (often) in a text file. (ASCII past 0x7F
    might appear, but rarely).

    Code modified from
    https://www.safaribooksonline.com/library/view/python-cookbook-2nd/0596007973/ch01s12.html
    """
    # if s contains any null, it's not text:
    if six.PY2 and "\0" in s:
        return False
    # an "empty" string is "text" (arbitrary but reasonable choice):
    if not s:
        return True

    text_characters = "".join(map(chr, range(32, 127))) + "\n\r\t\b"
    if six.PY2:
        _null_trans = string.maketrans("", "")
        # Get the substring of s made up of non-text characters
        t = s.translate(_null_trans, text_characters)
    else:
        # Not yet sure how to deal with this in python3. Will need example.
        return True

        # trans = str.maketrans('', '', text_characters)
        # t = s.translate(trans)

    # s is 'text' if less than 30% of its characters are non-text ones:
    return len(t)/float(len(s)) <= threshold

def resolve_url(url, directory=None, permissions=None):
    """
    Resolves a URL to a local file, and returns the path to
    that file.
    """

    u = urlparse(url)

    # create the name of the destination file
    if directory is None:
        directory = os.getcwd()
    filename = os.path.join(directory,os.path.basename(u.path))

    if u.scheme == '' or u.scheme == 'file':
        # for regular files, make a direct copy
        if os.path.isfile(u.path):
            if os.path.isfile(filename):
                # check to see if src and dest are the same file
                src_inode = os.stat(u.path)[stat.ST_INO]
                dst_inode = os.stat(filename)[stat.ST_INO]
                if src_inode != dst_inode:
                    shutil.copy(u.path, filename)
            else:
                shutil.copy(u.path, filename)
        else:
            errmsg  = "Cannot open file %s from URL %s" % (u.path, url)
            raise ValueError(errmsg)

    elif u.scheme == 'http' or u.scheme == 'https':
        s = requests.Session()
        s.mount(str(u.scheme)+'://',
            requests.adapters.HTTPAdapter(max_retries=5))

        # look for an ecp cookie file and load the cookies
        cookie_dict = {}
        ecp_file = '/tmp/ecpcookie.u%d' % os.getuid()
        if os.path.isfile(ecp_file):
            cj = cookielib.MozillaCookieJar()
            cj.load(ecp_file, ignore_discard=True, ignore_expires=True)
        else:
            cj = []

        for c in cj:
            if c.domain == u.netloc:
                # load cookies for this server
                cookie_dict[c.name] = c.value
            elif u.netloc == "code.pycbc.phy.syr.edu" and \
              c.domain == "git.ligo.org":
                # handle the redirect for code.pycbc to git.ligo.org
                cookie_dict[c.name] = c.value

        r = s.get(url, cookies=cookie_dict, allow_redirects=True)
        if r.status_code != 200:
            errmsg = "Unable to download %s\nError code = %d" % (url,
                r.status_code)
            raise ValueError(errmsg)

        # if we are downloading from git.ligo.org, check that we
        # did not get redirected to the sign-in page
        if u.netloc == 'git.ligo.org' or u.netloc == 'code.pycbc.phy.syr.edu':
            # Check if we have downloaded a binary file.
            if istext(r.content):
                soup = BeautifulSoup(r.content, 'html.parser')
                desc = soup.findAll(attrs={"property":"og:url"})
                if len(desc) and \
                  desc[0]['content'] == 'https://git.ligo.org/users/sign_in':
                    raise ValueError(ecp_cookie_error.format(url))

        output_fp = open(filename, 'wb')
        output_fp.write(r.content)
        output_fp.close()

    else:
        # TODO: We could support other schemes such as gsiftp by
        # calling out to globus-url-copy
        errmsg  = "Unknown URL scheme: %s\n" % (u.scheme)
        errmsg += "Currently supported are: file, http, and https."
        raise ValueError(errmsg)

    if not os.path.isfile(filename):
        errmsg = "Error trying to create file %s from %s" % (filename,url)
        raise ValueError(errmsg)

    if permissions:
        if os.access(filename, os.W_OK):
            os.chmod(filename, permissions)
        else:
            # check that the file has at least the permissions requested
            s = os.stat(filename)[stat.ST_MODE]
            if (s & permissions) != permissions:
                errmsg = "Could not change permissions on %s (read-only)" % url
                raise ValueError(errmsg)

    return filename


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

    To remove a configuration option, use the command line argument

    --config-delete section1:option1

    which will delete option1 from [section1] or

    --config-delete section1

    to delete all of the options in [section1]

    Deletes are implemented before overrides.

    This function returns an argparse OptionGroup to ensure these options are
    parsed correctly and can then be sent directly to initialize an
    WorkflowConfigParser.

    Parameters
    -----------
    parser : argparse.ArgumentParser instance
        The initialized argparse instance to add the workflow option group to.
    """
    workflowArgs = parser.add_argument_group('Configuration',
                                             'Options needed for parsing '
                                             'config file(s).')
    workflowArgs.add_argument("--config-files", nargs="+", action='store',
                           metavar="CONFIGFILE",
                           help="List of config files to be used in "
                                "analysis.")
    workflowArgs.add_argument("--config-overrides", nargs="*", action='store',
                           metavar="SECTION:OPTION:VALUE",
                           help="List of section,option,value combinations to "
                           "add into the configuration file. Normally the gps "
                           "start and end times might be provided this way, "
                           "and user specific locations (ie. output directories). "
                           "This can also be provided as SECTION:OPTION or "
                           "SECTION:OPTION: both of which indicate that the "
                           "corresponding value is left blank.")
    workflowArgs.add_argument("--config-delete", nargs="*", action='store',
                           metavar="SECTION:OPTION",
                           help="List of section,option combinations to delete "
                           "from the configuration file. This can also be "
                           "provided as SECTION which deletes the enture section"
                           " from the configuration file or SECTION:OPTION "
                           "which deletes a specific option from a given "
                           "section.")


class WorkflowConfigParser(glue.pipeline.DeepCopyableConfigParser):
    """
    This is a sub-class of glue.pipeline.DeepCopyableConfigParser, which lets
    us add a few additional helper features that are useful in workflows.
    """
    def __init__(self, configFiles=None, overrideTuples=None, parsedFilePath=None, deleteTuples=None):
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
        deleteTuples : List of (section, option) tuples
            Delete the (section, option) pairs provided
            in this list from provided .ini file(s). If the section only
            is provided, the entire section will be deleted.

        Returns
        --------
        WorkflowConfigParser
            Initialized WorkflowConfigParser instance.
        """
        if configFiles is None:
            configFiles = []
        if overrideTuples is None:
            overrideTuples = []
        if deleteTuples is None:
            deleteTuples = []
        glue.pipeline.DeepCopyableConfigParser.__init__(self)

        # Enable case sensitive options
        self.optionxform = str

        configFiles = [resolve_url(cFile) for cFile in configFiles]

        self.read_ini_file(configFiles)

        # Replace exe macros with full paths
        self.perform_exe_expansion()

        # Split sections like [inspiral&tmplt] into [inspiral] and [tmplt]
        self.split_multi_sections()

        # Populate shared options from the [sharedoptions] section
        self.populate_shared_sections()

        # Do deletes from command line
        for delete in deleteTuples:
            if len(delete) == 1:
                if self.remove_section(delete[0]) is False:
                    raise ValueError("Cannot delete section %s, "
                        "no such section in configuration." % delete )
                else:
                    logging.info("Deleting section %s from configuration",
                                 delete[0])
            elif len(delete) == 2:
                if self.remove_option(delete[0],delete[1]) is False:
                    raise ValueError("Cannot delete option %s from section %s,"
                        " no such option in configuration." % delete )
                else:
                    logging.info("Deleting option %s from section %s in "
                                 "configuration", delete[1], delete[0])
            else:
                raise ValueError("Deletes must be tuples of length 1 or 2. "
                    "Got %s." % str(delete) )

        # Do overrides from command line
        for override in overrideTuples:
            if len(override) not in [2,3]:
                errmsg = "Overrides must be tuples of length 2 or 3."
                errmsg = "Got %s." % (str(override) )
                raise ValueError(errmsg)
            section = override[0]
            option = override[1]
            value = ''
            if len(override) == 3:
                value = override[2]
            # Check for section existence, create if needed
            if not self.has_section(section):
                self.add_section(section)
            self.set(section, option, value)
            logging.info("Overriding section %s option %s with value %s "
                "in configuration.", section, option, value )


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
            self.write(fp)
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

        # Identify the deletes
        confDeletes = args.config_delete or []
        # and parse them
        parsedDeletes = []
        for delete in confDeletes:
            splitDelete = delete.split(":")
            if len(splitDelete) > 2:
                raise ValueError(
                    "Deletes must be of format section:option "
                    "or section. Cannot parse %s." % str(delete))
            else:
                parsedDeletes.append(tuple(splitDelete))

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
                    "or section:option. Cannot parse %s." % str(override))

        return cls(confFiles, parsedOverrides, None, parsedDeletes)


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
                    errmsg = "Cannot find exe %s in your path " %(testList[1])
                    errmsg += "and you specified ${which:%s}." %(testList[1])
                    raise ValueError(errmsg)

        return newString

    def get_subsections(self, section_name):
        """ Return a list of subsections for the given section name
        """
        # Keep only subsection names
        subsections = [sec[len(section_name)+1:] for sec in self.sections()\
                       if sec.startswith(section_name + '-')]

        for sec in subsections:
            sp = sec.split('-')
            # This is unusual, but a format [section-subsection-tag] is okay. Just
            # check that [section-subsection] section exists. If not it is possible
            # the user is trying to use an subsection name with '-' in it
            if (len(sp) > 1) and not self.has_section('%s-%s' % (section_name,
                                                                 sp[0])):
                raise ValueError( "Workflow uses the '-' as a delimiter so "
                    "this is interpreted as section-subsection-tag. "
                    "While checking section %s, no section with "
                    "name %s-%s was found. "
                    "If you did not intend to use tags in an "
                    "'advanced user' manner, or do not understand what "
                    "this means, don't use dashes in section "
                    "names. So [injection-nsbhinj] is good. "
                    "[injection-nsbh-inj] is not." % (sec, sp[0], sp[1]))

        if len(subsections) > 0:
            return [sec.split('-')[0] for sec in subsections]
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
                    print("Substitution failed")
                    raise
            if len(splitString) == 2:
                try:
                    testString = testString.replace('${'+repString+'}',\
                                       self.get(splitString[0],splitString[1]))
                except ConfigParser.NoOptionError:
                    print("Substitution failed")
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
        if not self.has_section('sharedoptions'):
            # No sharedoptions, exit
            return
        for key, value in self.items('sharedoptions'):
            assert(self.has_section('sharedoptions-%s' %(key)))
            # Comma separated
            values = value.split(',')
            common_options = self.items('sharedoptions-%s' %(key))
            for section in values:
                if not self.has_section(section):
                    self.add_section(section)
                for arg, val in common_options:
                    if arg in self.options(section):
                        raise ValueError('Option exists in both original ' + \
                               'ConfigParser section [%s] and ' %(section,) + \
                               'sharedoptions section: %s %s' \
                               %(arg,'sharedoptions-%s' %(key)))
                    self.set(section, arg, val)
            self.remove_section('sharedoptions-%s' %(key))
        self.remove_section('sharedoptions')

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
            # [pegasus_profile] specially is allowed to be overriden by
            # sub-sections
            if section == 'pegasus_profile':
                continue

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
            err_string = "No option '%s' in section [%s] " %(option,section)
            if not tags:
                raise ConfigParser.Error(err_string + ".")
            return_vals = []
            sub_section_list = []
            for sec_len in range(1, len(tags)+1):
                for tag_permutation in itertools.permutations(tags, sec_len):
                    joined_name = '-'.join(tag_permutation)
                    sub_section_list.append(joined_name)
            section_list = ["%s-%s" %(section, sb) for sb in sub_section_list]
            err_section_list = []
            for sub in sub_section_list:
                if self.has_section('%s-%s' %(section, sub)):
                    if self.has_option('%s-%s' %(section, sub), option):
                        err_section_list.append("%s-%s" %(section, sub))
                        return_vals.append(self.get('%s-%s' %(section, sub),
                                                    option))

            # We also want to recursively go into sections

            if not return_vals:
                err_string += "or in sections [%s]." \
                               %("] [".join(section_list))
                raise ConfigParser.Error(err_string)
            if len(return_vals) > 1:
                err_string += "and multiple entries found in sections [%s]."\
                              %("] [".join(err_section_list))
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


    @staticmethod
    def add_config_opts_to_parser(parser):
        """Adds options for configuration files to the given parser."""
        parser.add_argument("--config-files", type=str, nargs="+",
                            required=True,
                            help="A file parsable by "
                                 "pycbc.workflow.WorkflowConfigParser.")
        parser.add_argument("--config-overrides", type=str, nargs="+",
                            default=None, metavar="SECTION:OPTION:VALUE",
                            help="List of section:option:value combinations "
                                 "to add into the configuration file.")


    @classmethod
    def from_cli(cls, opts):
        """Loads a config file from the given options, with overrides and
        deletes applied.
        """
        # read configuration file
        logging.info("Reading configuration file")
        if opts.config_overrides is not None:
            overrides = [override.split(":")
                         for override in opts.config_overrides]
        else:
            overrides = None
        if opts.config_delete is not None:
            deletes = [delete.split(":") for delete in opts.config_delete]
        else:
            deletes = None
        return cls(opts.config_files, overrides, deleteTuples=deletes)
