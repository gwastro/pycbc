# Copyright (C) 2013 Alex Nitz
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
This modules contains functions reading, generating, and segmenting strain data
"""
import logging

from optparse import OptionGroup
from pycbc import psd, DYN_RANGE_FAC
from pycbc.types import float32
from pycbc.frame import read_frame
from pycbc.inject import InjectionSet
from pycbc.filter import resample_to_delta_t, highpass

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
            
def from_cli(opt):
    """Parses the CLI options related to strain data reading and conditioning.

    Parameters
    ----------
    opt : object
        Result of parsing the CLI with OptionParser, or any object with the
        required attributes  (gps-start-time, gps-end-time, strain-high-pass, 
        pad-data, sample-rate, frame-cache, channel-name, fake-strain, 
        fake-strain-seed).

    Returns
    -------
    strain : TimeSeries
        The time series containing the conditioned strain data.
    """
    if opt.frame_cache:
        logging.info("Reading Frames")
        strain = read_frame(opt.frame_cache, opt.channel_name, 
                            start_time=opt.gps_start_time-opt.pad_data, 
                            end_time=opt.gps_end_time+opt.pad_data)

        if opt.injection_file:
            logging.info("Applying injections")
            injections = InjectionSet(opt.injection_file)
            injections.apply(strain, opt.channel_name[0:2])

        logging.info("Highpass Filtering")
        strain = highpass(strain, frequency=opt.strain_high_pass)

        logging.info("Converting to float32")
        strain = (strain * DYN_RANGE_FAC).astype(float32)

        logging.info("Resampling data")
        strain = resample_to_delta_t(strain, 1.0/opt.sample_rate)

        logging.info("Remove Padding")
        start = opt.pad_data*opt.sample_rate
        end = len(strain)-opt.sample_rate*opt.pad_data
        strain = strain[start:end]
        
    if opt.fake_strain:
        logging.info("Generating Fake Strain")
        duration = opt.gps_end_time - opt.gps_start_time
        tlen = duration * opt.sample_rate
        pdf = 1.0/128
        plen = int(opt.sample_rate / pdf) / 2 + 1

        logging.info("Making PSD for strain")
        strain_psd = pycbc.psd.from_string(opt.fake_strain, plen, 
                                           pdf, opt.low_frequency_cutoff)
        
        logging.info("Making colored noise")
        strain = pycbc.noise.noise_from_psd(tlen, 1.0/opt.sample_rate, 
                                            strain_psd, 
                                            seed=opt.fake_strain_seed)

        if opt.injection_file:
            logging.info("Applying injections")
            injections = InjectionSet(opt.injection_file)
            injections.apply(strain, opt.channel_name[0:2])
        
        logging.info("Converting to float32")
        strain = (DYN_RANGE_FAC * strain).astype(float32)    
    return strain

def insert_strain_option_group(parser):
    """
    Adds the options used to call the pycbc.strain.from_cli function to an
    optparser as an OptionGroup. This should be used if you
    want to use these options in your code.
 
    Parameters
    -----------
    parser : object
        OptionParser instance.
    """
    
    data_reading_group = OptionGroup(parser, "Options for obtaining h(t)",
                  "These options are used for generating h(t) either by "
                  "reading from a file or by generating it. This is only needed "
                  "if the PSD is to be estimated from the data, ie. if the "
                  "--psd-estimation option is given.")

    # Required options
    data_reading_group.add_option("--gps-start-time", 
                            help="The gps start time of the data", type=int)
    data_reading_group.add_option("--gps-end-time", 
                            help="The gps end time of the data", type=int)
    data_reading_group.add_option("--strain-high-pass", type=float, 
                            help="High pass frequency")
    data_reading_group.add_option("--pad-data", 
              help="Extra padding to remove highpass corruption (s)", type=int)
    data_reading_group.add_option("--sample-rate", type=int, 
                            help="The sample rate to use for h(t) generation.")
    data_reading_group.add_option("--channel-name", type=str, 
                   help="The channel containing the gravitational strain data")
                   
    #Read from cache file              
    data_reading_group.add_option("--frame-cache", type=str, 
                            help="Cache file containing the frame locations.")
    
    #Generate gaussian noise with given psd           
    data_reading_group.add_option("--fake-strain", 
                help="Name of model PSD for generating fake gaussian noise." +\
                     " Choices are " + str(psd.get_list()) , 
                     choices=psd.get_list())
    data_reading_group.add_option("--fake-strain-seed", type=int, default=0,
                help="Seed value for the generation of fake colored" + \
                     " gaussian noise")
                                
    #optional       
    data_reading_group.add_option("--injection-file", type=str, 
                      help="(optional) Injection file used to add " + \
                           "waveforms into the strain")                 
                   
    parser.add_option_group(data_reading_group)

def verify_strain_options(opts, parser):
    """Parses the strain data CLI options and verifies that they are consistent 
    and reasonable.
    
    Parameters
    ----------
    opt : object
        Result of parsing the CLI with OptionParser, or any object with the
        required attributes (gps-start-time, gps-end-time, strain-high-pass, 
        pad-data, sample-rate, frame-cache, channel-name, fake-strain, 
        fake-strain-seed).
    parser : object
        OptionParser instance.
    """
    ensure_one_opt(opts, parser, ['--frame-cache', '--fake-strain'])
    
    required_opts(opts, parser, 
                  ['--gps-start-time', '--gps-end-time', '--strain-high-pass',
                   '--pad-data', '--sample-rate', '--channel-name',
                   ])               
                   
def StrainData(object):
    pass        
