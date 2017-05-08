import re
import copy
try:
    import configparser as ConfigParser
except ImportError:
    import ConfigParser

def parse_workflow_ini_file(cpFile,parsed_filepath=None):
    """Read a .ini file in, parse it as described in the documentation linked
    to above, and return the parsed ini file.

    Parameters
    ----------
    cpFile : The path to a .ini file to be read in
    parsed_filepath: Boolean, optional
        If provided, the .ini file, after parsing, will be written to this
        location   

    Returns
    -------
    cp: The parsed ConfigParser class containing the read in .ini file
    """    
    # First read the .ini file
    cp = read_ini_file(cpFile)
    print(cp.sections())

    # Check for any substitutions that can be made
    # FIXME: The python 3 version of ConfigParser can do this automatically
    # move over to that if it can be backported to python2.X.
    # We use the same formatting as the new configparser module when doing
    # ExtendedInterpolation
    # This is described at http://docs.python.org/3.4/library/configparser.html
    #cp = perform_extended_interpolation(cp)

    # Split sections like [inspiral&tmplt] into [inspiral] and [tmplt]
    cp = split_multi_sections(cp)
    print(cp.sections())

    # Check for duplicate options in sub-sections
    sanity_check_subsections(cp)
    print(cp.sections())

    # Dump parsed .ini file if needed
    if parsed_filepath:
        fp = open(parsed_filepath,'w')
        cp.write(fp)
        fp.close()

    return cp


def read_ini_file(cpFile):
    """Read a .ini file and return it as a ConfigParser class.
    This function does none of the parsing/combining of sections. It simply
    reads the file and returns it unedited

    Parameters
    ----------
    cpFile : The path to a .ini file to be read in

    Returns
    -------
    cp: The ConfigParser class containing the read in .ini file
    """    

    # Initialise ConfigParser class
    cp = ConfigParser.ConfigParser(\
        interpolation=ConfigParser.ExtendedInterpolation())
    # Read the file
    fp = open(cpFile,'r')
    cp.read_file(fp)
    fp.close()
    return cp

def perform_extended_interpolation(cp,preserve_orig_file=False):
    """Filter through an ini file and replace all examples of 
    ExtendedInterpolation formatting with the exact value. For values like
    ${example} this is replaced with the value that corresponds to the option
    called example ***in the same section***

    For values like ${common|example} this is replaced with the value that
    corresponds to the option example in the section [common]. Note that
    in the python3 config parser this is ${common:example} but python2.7
    interprets the : the same as a = and this breaks things

    Nested interpolation is not supported here.

    Parameters
    ----------
    cp: ConfigParser object
    preserve_orig_file: Boolean, optional
        By default the input ConfigParser object will be modified in place. If
        this is set deepcopy will be used and the input will be preserved.
        Default = False

    Returns
    -------
    cp: parsed ConfigParser object
    """
    # Deepcopy the cp object if needed
    if preserve_orig_file:
        cp = copy.deepcopy(cp)

    # Do not allow any interpolation of the section names
    for section in cp.sections():
         for option,value in cp.items(section):
             # Check the option name
             newStr = interpolate_string(option,cp,section)
             if newStr != option:
                 cp.set(section,newStr,value)
                 cp.remove_option(section,option)
             # Check the value
             newStr = interpolate_string(value,cp,section)
             if newStr != value:
                 cp.set(section,option,newStr)

    return cp
        
def interpolate_string(testString,cp,section):
    """Take a string and replace all example of ExtendedInterpolation formatting
    within the string with the exact value. 

    For values like ${example} this is replaced with the value that corresponds
    to the option called example ***in the same section***

    For values like ${common|example} this is replaced with the value that
    corresponds to the option example in the section [common]. Note that
    in the python3 config parser this is ${common:example} but python2.7
    interprets the : the same as a = and this breaks things

    Nested interpolation is not supported here.

    Parameters
    ----------
    testString: String
        The string to parse and interpolate
    cp: ConfigParser
        The ConfigParser object to look for the interpolation strings within
    section: String
        The current section of the ConfigParser object

    Returns
    ----------
    testString: String
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
            testString = testString.replace('${'+repString+'}',\
                                    cp.get(section,splitString[0]))
        if len(splitString) == 2:
            testString = testString.replace('${'+repString+'}',\
                                    cp.get(splitString[0],splitString[1]))
        reObj = re.search(r"\$\{.*?\}", testString)

    return testString

def split_multi_sections(cp,preserve_orig_file=False):
    """Parse through a supplied ConfigParser object and splits any sections
    labelled with an "&" sign (for e.g. [inspiral&tmpltbank]) into [inspiral]
    and [tmpltbank] sections. If these individual sections already exist they
    will be appended to. If an option exists in both the [inspiral] and 
    [inspiral&tmpltbank] sections an error will be thrown
   
    Parameters
    ----------
    cp: The ConfigParser class
    preserve_orig_file: Boolean, optional
        By default the input ConfigParser object will be modified in place. If
        this is set deepcopy will be used and the input will be preserved.
        Default = False

    Returns
    ----------
    cp: The ConfigParser class 
    """
    # Deepcopy the cp object if needed
    if preserve_orig_file:
        cp = copy.deepcopy(cp)

    # Begin by looping over all sections
    for section in cp.sections():
        # Only continue if section needs splitting
        if '&' not in section:
            continue
        # Get list of section names to add these options to
        splitSections = section.split('&')
        for newSec in splitSections:
            # Add sections if they don't already exist
            if not cp.has_section(newSec):
                cp.add_section(newSec)
            add_options_to_section(cp,newSec,cp.items(section))
        cp.remove_section(section)
    return cp

def sanity_check_subsections(cp):
    """This function goes through the ConfigParset and checks that any options
    given in the [SECTION_NAME] section are not also given in any 
    [SECTION_NAME-SUBSECTION] sections.

    Parameters
    ----------
    cp: The ConfigParser class

    Returns
    ----------
    None
    """
    # Loop over the sections in the ini file
    for section in cp.sections():
        # Loop over the sections again
        for section2 in cp.sections():
            # Check if any are subsections of section
            if section2.startswith(section + '-'):
                # Check for duplicate options whenever this exists
                check_duplicate_options(cp,section,section2,raise_error=True)

def add_options_to_section(cp,section,items,preserve_orig_file=False,\
                           overwrite_options=False):
    """Add a set of options and values to a section of a ConfigParser object.
    Will throw an error if any of the options being added already exist, this
    behaviour can be overridden if desired

    Parameters
    ----------
    cp: The ConfigParser class
    section: string
        The name of the section to add options+values to
    items: list of tuples
        Each tuple contains (at [0]) the option and (at [1]) the value to add
        to the section of the ini file
    preserve_orig_file: Boolean, optional
        By default the input ConfigParser object will be modified in place. If
        this is set deepcopy will be used and the input will be preserved.
        Default = False
    overwrite_options: Boolean, optional
        By default this function will throw a ValueError if an option exists in
        both the original section in the ConfigParser *and* in the provided
        items. This will override so that the options+values given in items
        will replace the original values if the value is set to True.
        Default = True 
    Returns
    ----------
    cp: The ConfigParser class 
    """
    # Sanity checking
    if not cp.has_section(section):
        raise ValueError('Section %s not present in ConfigParser.' %(section,))

    # Deepcopy the cp object if needed
    if preserve_orig_file:
        cp = copy.deepcopy(cp)

    # Check for duplicate options first
    for option,value in items:
        if not overwrite_options:
            if option in cp.options(section):
                raise ValueError('Option %s exists in both original' + \
                                 'ConfigParser and input list' %(option,))
        cp.set(section,option,value)

    return cp

def check_duplicate_options(cp,section1,section2,raise_error=False):
    """Check for duplicate options in two sections, section1 and section2.
    Will return True if there are duplicate options and False if not
    

    Parameters
    ----------
    cp: The ConfigParser class
    section1: string
        The name of the first section to compare
    section2: string
        The name of the second section to compare
    raise_error: Boolean, optional
        If True, raise an error if duplicates are present.
        Default = False

    Returns
    ----------
    duplicate: List
        List of duplicate options
    """
    # Sanity checking
    if not cp.has_section(section1):
        raise ValueError('Section %s not present in ConfigParser.'%(section1,))
    if not cp.has_section(section2):
        raise ValueError('Section %s not present in ConfigParser.'%(section2,))

    items1 = cp.items(section1)
    items2 = cp.items(section2)

    # The list comprehension here creates a list of all duplicate items
    duplicates = [x for x in items1 if x in items2]

    if duplicates and raise_error:
        raise ValueError('The following options appear in both section ' +\
                         '%s and %s: %s' \
                         %(section1,section2,duplicates.join(' ')))
    
    return duplicates
