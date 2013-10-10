#!/usr/bin/python
# Copyright (C) 2012 Alex Nitz, Andrew Miller, Tito Dal Canton
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
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

from optparse import OptionParser, OptionGroup

from pycbc.psd.read import *
from pycbc.psd.analytical import *
from pycbc.psd.estimate import *

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
            

def from_cli(opt, length, delta_f, low_frequency_cutoff, 
             strain=None, dyn_range_factor=1):
    """Parses the CLI options related to the noise PSD and returns a
    FrequencySeries with the corresponding PSD. If necessary, the PSD is
    linearly interpolated to achieve the resolution specified in the CLI.

    Parameters
    ----------
    opt : object
        Result of parsing the CLI with OptionParser, or any object with the
        required attributes (psd_model, psd_file, asd_file, psd_estimation,
        psd_segment_length, psd_segment_stride, psd_inverse_length, psd_output).
    length : int
        The length in samples of the output PSD.
    delta_f : float
        The frequency step of the output PSD.
    low_frequency_cutoff: float
        The low frequncy cutoff to use when calculating the PSD.
    strain : {None, TimeSeries}
        Time series containing the data from which the PSD should be measured,
        when psd_estimation is in use.
    dyn_range_factor : {1, float}
        For PSDs taken from models or text files, if `dyn_range_factor` is
        not None, then the PSD is multiplied by `dyn_range_factor` ** 2.

    Returns
    -------
    psd : FrequencySeries
        The frequency series containing the PSD.
    """
    f_low = low_frequency_cutoff
    sample_rate = int((length -1) * 2 * delta_f)

    if (opt.psd_model or opt.psd_file) and not opt.psd_estimation:
        # PSD from lalsimulation or file
        if opt.psd_model:
            psd = from_string(opt.psd_model, length, delta_f, f_low)           
        elif opt.psd_file:
            psd = from_txt(opt.psd_file, length, 
                           delta_f, f_low, is_asd_file=False)                          
        elif opt.asd_file:
            psd = from_txt(opt.psd_file, length, 
                           delta_f, f_low, is_asd_file=True)
        
        psd *= dyn_range_factor ** 2
           
    elif opt.psd_estimation and not (opt.psd_model or opt.psd_file):
        # estimate PSD from data
        psd = welch(strain, avg_method=opt.psd_estimation,
                    seg_len=int(opt.psd_segment_length * sample_rate),
                    seg_stride=int(opt.psd_segment_stride * sample_rate))
                    
        if delta_f != psd.delta_f:
            psd = interpolate(psd, delta_f) 
                    
    if opt.psd_inverse_length:
        psd = inverse_spectrum_truncation(psd, 
            int(opt.psd_inverse_length * sample_rate),
            low_frequency_cutoff=f_low)
            
    if opt.psd_output:
        (psd / (dyn_range_factor ** 2)).save(opt.psd_output)

    return psd

def insert_psd_option_group(parser):
    """
    Adds the options used to call the pycbc.psd.from_cli function to an
    optparser OptionGroup, which is the returned. This should be used if you
    want to use these options in your code.

    Returns
    --------
    psd_options : optparser.OptionGroup
        The optparser.OptionGroup containing all the necessary psd options.
    """
    psd_options = OptionGroup(parser,
                   "Options related to the noise PSD generation")
    psd_options.add_option("--psd-model",
                           help="Get PSD from given analytical model")
    psd_options.add_option("--psd-file",
                           help="Get PSD using given PSD ASCII file")
    psd_options.add_option("--asd-file",
                           help="Get PSD using given ASD ASCII file")
    psd_options.add_option("--psd-estimation",
                 help="Measure PSD from the data, using given average method",
                 choices=["mean", "median", "median-mean"])
    psd_options.add_option("--psd-segment-length", type=float, 
                          help="The segment length for PSD estimation (s)")
    psd_options.add_option("--psd-segment-stride", type=float, 
                          help="The segment stride for PSD estimation (s)")
    psd_options.add_option("--psd-inverse-length", type=float, 
                          help="The maximum length of the impulse response " +\
                          "of the overwhitening filter (s)")
    psd_options.add_option("--psd-output", help="Write PSD to specified file")
    parser.add_option_group(psd_options)
    
def verify_psd_options(opt, parser):
    """Parses the CLI options and verifies that they are consistent and 
    reasonable.
    
    Parameters
    ----------
    opt : object
        Result of parsing the CLI with OptionParser, or any object with the
        required attributes (psd_model, psd_file, asd_file, psd_estimation,
        psd_segment_length, psd_segment_stride, psd_inverse_length, psd_output).
    parser : object
        OptionParser instance.
    """
    ensure_one_opt(opt, parser, ['--psd-file', '--psd-model', 
                                '--psd-estimation', '--asd-file'])
                                
    if opt.psd_estimation:
        required_opts(opt, parser, 
                      ['--psd-segment-stride', '--psd-segment-length' ],
                      required_by = "--psd-estimation")
        
    
