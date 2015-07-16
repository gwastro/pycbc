# Copyright (C) 2015 Ian Harry
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Generals
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""
This modules contains extensions for use with argparse
"""
import copy
import argparse
from collections import defaultdict

class DictWithDefaultReturn(defaultdict):
    default_set = False
    ifo_set = False

class MultiDetOptionAction(argparse.Action):
    # Initialise the same as the standard 'append' action
    def __init__(self,
                 option_strings,
                 dest,
                 nargs=None,
                 const=None,
                 default=None,
                 type=None,
                 choices=None,
                 required=False,
                 help=None,
                 metavar=None):
        if type is not None:
            self.internal_type=type
        else:
            self.internal_type=str
        new_default = DictWithDefaultReturn(lambda: default)
        #new_default.default_value=default
        if nargs == 0:
            raise ValueError('nargs for append actions must be > 0; if arg '
                             'strings are not supplying the value to append, '
                             'the append const action may be more appropriate')
        if const is not None and nargs != OPTIONAL:
            raise ValueError('nargs must be %r to supply const' % OPTIONAL)
        super(MultiDetOptionAction, self).__init__(
            option_strings=option_strings,
            dest=dest,
            nargs=nargs,
            const=const,
            default=new_default,
            type=str,
            choices=choices,
            required=required,
            help=help,
            metavar=metavar)

    def __call__(self, parser, namespace, values, option_string=None):
        # Again this is modified from the standard argparse 'append' action
        err_msg = "Issue with option: %s \n" %(self.dest,)
        err_msg += "Received value: %s \n" %(' '.join(values),)
        if getattr(namespace, self.dest, None) is None:
            setattr(namespace, self.dest, DictWithDefaultReturn())
        items = getattr(namespace, self.dest)
        items = copy.copy(items)
        for value in values:
            value = value.split(':')
            if len(value) == 2:
                # "Normal" case, all ifos supplied independetly as "H1,VALUE"
                if items.default_set:
                    err_msg += "If you are supplying a value for all ifos, you "
                    err_msg += "cannot also supply values for specific ifos."
                    raise ValueError(err_msg)
                items[value[0]] = self.internal_type(value[1])
                items.ifo_set = True
            elif len(value) == 1:
                # OR supply only one value and use this for all ifos
                if items.default_set:
                    err_msg += "If you are supplying a value for all ifos, you "
                    err_msg += "must only supply one value."
                    raise ValueError(err_msg)
                # Can't use a global and ifo specific options
                if items.ifo_set:
                    err_msg += "If you are supplying a value for all ifos, you "
                    err_msg += "cannot also supply values for specific ifos."
                    raise ValueError(err_msg)
                #items.default_value = self.internal_type(value[0])
                new_default = self.internal_type(value[0])
                items.default_factory = lambda: new_default
                items.default_set = True
            else:
                err_msg += "The character ':' is used to deliminate the "
                err_msg += "ifo and the value. Please do not use it more than "
                err_msg += "once."
                raise ValueError(err_msg)
        setattr(namespace, self.dest, items)

class MultiDetOptionActionSpecial(MultiDetOptionAction):
    """
    This class in an extension of the MultiDetOptionAction class to handle
    cases where the : is already a special character. For example the channel
    name is something like H1:CHANNEL_NAME. Here the channel name *must*
    be provided uniquely for each ifo. The dictionary key is set to H1 and the
    value to H1:CHANNEL_NAME for this example.
    """
    def __call__(self, parser, namespace, values, option_string=None):
        # Again this is modified from the standard argparse 'append' action
        err_msg = "Issue with option: %s \n" %(self.dest,)
        err_msg += "Received value: %s \n" %(' '.join(values),)
        if getattr(namespace, self.dest, None) is None:
            setattr(namespace, self.dest, {})
        items = getattr(namespace, self.dest)
        items = copy.copy(items)
        for value in values:
            value_split = value.split(':')
            if len(value_split) == 2:
                # "Normal" case, all ifos supplied independetly as "H1:VALUE"
                if items.has_key(value_split[0]):
                    err_msg += "Multiple values supplied for ifo %s.\n" \
                               %(value_split[0],)
                    err_msg += "Already have %s." %(items[value_split[0]])
                    raise ValueError(err_msg)
                else:
                    items[value_split[0]] = value
            elif len(value_split) == 3:
                # This is an unadvertised feature. It is used for cases where I
                # want to pretend H1 data is actually L1 (or similar). So if I
                # supply --channel-name H1:L1:LDAS-STRAIN I can use L1 data and
                # pretend it is H1 internally.
                if items.has_key(value_split[0]):
                    err_msg += "Multiple values supplied for ifo %s.\n" \
                               %(value_split[0],)
                    err_msg += "Already have %s." %(items[value_split[0]])
                    raise ValueError(err_msg)
                else:
                    items[value_split[0]] = ':'.join(value_split[1:3])
            else:
                err_msg += "The character ':' is used to deliminate the "
                err_msg += "ifo and the value. It must appear exactly "
                err_msg += "once."
                raise ValueError(err_msg)
        setattr(namespace, self.dest, items)


class MultiDetOptionAppendAction(MultiDetOptionAction):
    def __call__(self, parser, namespace, values, option_string=None):
        # Again this is modified from the standard argparse 'append' action
        if getattr(namespace, self.dest, None) is None:
            setattr(namespace, self.dest, {})
        items = getattr(namespace, self.dest)
        items = copy.copy(items)
        for value in values:
            value = value.split(':')
            if len(value) == 2:
                # "Normal" case, all ifos supplied independetly as "H1:VALUE"
                if items.has_key(value[0]):
                    items[value[0]].append(self.internal_type(value[1]))
                else:
                    items[value[0]] = [self.internal_type(value[1])]
            else:
                err_msg = "Issue with option: %s \n" %(self.dest,)
                err_msg += "Received value: %s \n" %(' '.join(values),)
                err_msg += "The character ':' is used to distinguish the "
                err_msg += "ifo and the value. It must be given exactly once "
                err_msg += "for all entries"
                raise ValueError(err_msg)
        setattr(namespace, self.dest, items)

def required_opts(opt, parser, opt_list, required_by=None):
    """Check that all the opts are defined 
    
    Parameters
    ----------
    opt : object
        Result of option parsing
    parser : object
        OptionParser instance.
    opt_list : list of strings
    required_by : string, optional
        the option that requires these options (if applicable)
    """
    for name in opt_list:
        attr = name[2:].replace('-', '_')
        if not hasattr(opt, attr) or (getattr(opt, attr) is None):
            err_str = "%s is missing " % name
            if required_by is not None:
                err_str += ", required by %s" % required_by
            parser.error(err_str)

def required_opts_multi_ifo(opt, parser, ifo, opt_list, required_by=None):
    """Check that all the opts are defined 
    
    Parameters
    ----------
    opt : object
        Result of option parsing
    parser : object
        OptionParser instance.
    ifo : string
    opt_list : list of strings
    required_by : string, optional
        the option that requires these options (if applicable)
    """
    for name in opt_list:
        attr = name[2:].replace('-', '_')
        try:
            if getattr(opt, attr)[ifo] is None:
                raise KeyError
        except KeyError:
            err_str = "%s is missing " % name
            if required_by is not None:
                err_str += ", required by %s" % required_by
            parser.error(err_str)

def ensure_one_opt(opt, parser, opt_list):
    """  Check that one and only one in the opt_list is defined in opt
    
    Parameters
    ----------
    opt : object
        Result of option parsing
    parser : object
        OptionParser instance.
    opt_list : list of strings
    """

    the_one = None
    for name in opt_list:
        attr = name[2:].replace('-', '_')
        if hasattr(opt, attr) and (getattr(opt, attr) is not None):
            if the_one is None:
                the_one = name
            else:
                parser.error("%s and %s are mutually exculsive" \
                              % (the_one, name))

    if the_one is None:
        parser.error("you must supply one of the following %s" \
                      % (', '.join(opt_list)))

def ensure_one_opt_multi_ifo(opt, parser, ifo, opt_list):
    """  Check that one and only one in the opt_list is defined in opt
    
    Parameters
    ----------
    opt : object
        Result of option parsing
    parser : object
        OptionParser instance.
    opt_list : list of strings
    """

    the_one = None
    for name in opt_list:
        attr = name[2:].replace('-', '_')
        try:
            if getattr(opt, attr)[ifo] is None:
                raise KeyError
        except KeyError:
            pass
        else:
            if the_one is None:
                the_one = name
            else:
                parser.error("%s and %s are mutually exculsive" \
                              % (the_one, name))

    if the_one is None:
        print opt, ifo
        parser.error("you must supply one of the following %s" \
                      % (', '.join(opt_list)))

def copy_opts_for_single_ifo(opt, ifo):
    """
    Takes the namespace object (opt) from the multi-detector interface and
    returns a namespace object for a single ifo that can be used with
    functions expecting output from the single-detector interface.
    """
    opt = copy.deepcopy(opt)
    for arg, val in vars(opt).items():
        if isinstance(val, DictWithDefaultReturn):
            setattr(opt, arg, getattr(opt, arg)[ifo])
    return opt

def convert_to_process_params_dict(opt):
    """
    Takes the namespace object (opt) from the multi-detector interface and
    returns a dictionary of command line options that will be handled correctly
    by the register_to_process_params ligolw function.
    """
    opt = copy.deepcopy(opt)
    for arg, val in vars(opt).items():
        if isinstance(val, DictWithDefaultReturn):
            new_val = []
            for key in val.keys():
                if isinstance(val[key], list):
                    for item in val[key]:
                        if item is not None:
                            new_val.append(':'.join([key, str(item)]))
                else:
                    if val[key] is not None:
                        new_val.append(':'.join([key, str(val[key])]))
            setattr(opt, arg, new_val)
    return vars(opt)

def positive_float(s):
    """
    Ensure argument is a positive real number and return it as float.

    To be used as type in argparse arguments.
    """
    err_msg = "must be a positive number, not %r" % s
    try:
        value = float(s)
    except ValueError, e:
        raise argparse.ArgumentTypeError(err_msg)
    if value <= 0:
        raise argparse.ArgumentTypeError(err_msg)
    return value

def nonnegative_float(s):
    """
    Ensure argument is a positive real number or zero and return it as float.

    To be used as type in argparse arguments.
    """
    err_msg = "must be either positive or zero, not %r" % s
    try:
        value = float(s)
    except ValueError, e:
        raise argparse.ArgumentTypeError(err_msg)
    if value < 0:
        raise argparse.ArgumentTypeError(err_msg)
    return value