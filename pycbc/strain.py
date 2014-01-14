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
import logging, numpy, lal
import pycbc.noise
import pycbc.psd
from pycbc import psd
from pycbc.types import float32
from pycbc.frame import read_frame
from pycbc.inject import InjectionSet
from pycbc.filter import resample_to_delta_t, highpass, make_frequency_series

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
            
def from_cli(opt, dyn_range_fac=1):
    """Parses the CLI options related to strain data reading and conditioning.

    Parameters
    ----------
    opt : object
        Result of parsing the CLI with OptionParser, or any object with the
        required attributes  (gps-start-time, gps-end-time, strain-high-pass, 
        pad-data, sample-rate, (frame-cache or frame-files), channel-name, 
        fake-strain, fake-strain-seed).
        
    dyn_range_fac: {float, 1}, optional
        A large constant to reduce the dynamic range of the strain.

    Returns
    -------
    strain : TimeSeries
        The time series containing the conditioned strain data.
    """
    if opt.frame_cache or opt.frame_files:
        if opt.frame_cache:
            frame_source = opt.frame_cache 
        if opt.frame_files:
            frame_source = opt.frame_files 

        logging.info("Reading Frames")
        strain = read_frame(frame_source, opt.channel_name, 
                            start_time=opt.gps_start_time-opt.pad_data, 
                            end_time=opt.gps_end_time+opt.pad_data)

        if opt.injection_file:
            logging.info("Applying injections")
            injections = InjectionSet(opt.injection_file)
            injections.apply(strain, opt.channel_name[0:2])

        logging.info("Highpass Filtering")
        strain = highpass(strain, frequency=opt.strain_high_pass)

        logging.info("Converting to float32")
        strain = (strain * dyn_range_fac).astype(float32)

        if opt.gating_file:
            logging.info("Gating glitches")
            gate_params = numpy.loadtxt(opt.gating_file)
            if len(gate_params.shape) == 1:
                gate_params = [gate_params]
            strain = gate_data(strain, gate_params,
                            data_start_time=opt.gps_start_time-opt.pad_data)

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
        strain._epoch = lal.LIGOTimeGPS(opt.gps_start_time)

        if opt.injection_file:
            logging.info("Applying injections")
            injections = InjectionSet(opt.injection_file)
            injections.apply(strain, opt.channel_name[0:2])
        
        logging.info("Converting to float32")
        strain = (dyn_range_fac * strain).astype(float32)

    if opt.injection_file:
        strain.injections = injections
    
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
    
    data_reading_group = parser.add_argument_group("Options for obtaining h(t)",
                  "These options are used for generating h(t) either by "
                  "reading from a file or by generating it. This is only "
                  "needed if the PSD is to be estimated from the data, ie. "
                  " if the --psd-estimation option is given.")

    # Required options
    data_reading_group.add_argument("--gps-start-time", 
                            help="The gps start time of the data", type=int)
    data_reading_group.add_argument("--gps-end-time", 
                            help="The gps end time of the data", type=int)
    data_reading_group.add_argument("--strain-high-pass", type=float, 
                            help="High pass frequency")
    data_reading_group.add_argument("--pad-data", 
              help="Extra padding to remove highpass corruption (s)", type=int)
    data_reading_group.add_argument("--sample-rate", type=int, 
                            help="The sample rate to use for h(t) generation.")
    data_reading_group.add_argument("--channel-name", type=str, 
                   help="The channel containing the gravitational strain data")
                   
    #Read from cache file              
    data_reading_group.add_argument("--frame-cache", type=str, nargs="+",
                            help="Cache file containing the frame locations.")
    
    #Read from cache file              
    data_reading_group.add_argument("--frame-files",
                            type=str, nargs="+",
                            help="list of frame files")   
    
    #Generate gaussian noise with given psd           
    data_reading_group.add_argument("--fake-strain",
                help="Name of model PSD for generating fake gaussian noise."
                     " Choices are " + str(psd.get_list()) , 
                     choices=psd.get_list())
    data_reading_group.add_argument("--fake-strain-seed", type=int, default=0,
                help="Seed value for the generation of fake colored"
                     " gaussian noise")
                                
    #optional       
    data_reading_group.add_argument("--injection-file", type=str, 
                      help="(optional) Injection file used to add "
                           "waveforms into the strain")                 

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
    ensure_one_opt(opts, parser, ['--frame-cache', 
                                  '--fake-strain', 
                                  '--frame-files'])  
    required_opts(opts, parser, 
                  ['--gps-start-time', '--gps-end-time', '--strain-high-pass',
                   '--pad-data', '--sample-rate', '--channel-name',
                   ])               

def gate_data(data, gate_params, data_start_time):
    def inverted_tukey(M, n_pad):
        midlen = M - 2*n_pad
        if midlen < 1:
            raise ValueError("No zeros left after applying padding.")
        padarr = 0.5*(1.+numpy.cos(numpy.pi*numpy.arange(n_pad)/n_pad))
        return numpy.concatenate((padarr,numpy.zeros(midlen),padarr[::-1]))

    sample_rate = 1./data.delta_t
    temp = data.data
    for glitch_time, glitch_width, pad_width in gate_params:
        t_start = glitch_time - glitch_width - pad_width - data_start_time
        t_end = glitch_time + glitch_width + pad_width - data_start_time
        if t_start > data.duration or t_end < 0.: 
            continue # Skip gate segments that don't overlap
        win_samples = int(2*sample_rate*(glitch_width+pad_width))
        pad_samples = int(sample_rate*pad_width)
        window = inverted_tukey(win_samples, pad_samples)
        offset = int(t_start * sample_rate)
        idx1 = max(0, -offset)
        idx2 = min(len(window), len(data)-offset)
        temp[idx1+offset:idx2+offset] *= window[idx1:idx2]

    return data

class StrainSegments(object):
    """ Class for managing manipulation of strain data for the purpose of 
        matched filtering. This includes methods for segmenting and 
        conditioning.
    """
    def __init__(self, strain, segment_length=None, segment_start_pad=0, 
                 segment_end_pad=0, trigger_start=None, trigger_end=None, 
                 filter_inj_only=False):
        """ Determine how to chop up the strain data into smaller segents
            for analysis.
        """
        self._fourier_segments = None
        self.strain = strain
        
        self.delta_t = strain.delta_t
        self.sample_rate = strain.sample_rate
        
        if segment_length:
            seg_len = segment_length
        else:
            seg_len = strain.duration
            
        self.delta_f = 1.0 / seg_len
        self.time_len = seg_len * self.sample_rate
        self.freq_len = self.time_len / 2 + 1
            
        seg_end_pad = segment_end_pad
        seg_start_pad = segment_start_pad
        
        if not trigger_start:
            trigger_start = int(strain.start_time)
            
        if not trigger_end:
            trigger_end = int(strain.end_time)
        
        throwaway_size = seg_start_pad + seg_end_pad    
        seg_width = seg_len - throwaway_size   
        
        # The amount of time we can actually analyze given the
        # amount of padding that is needed       
        analyzable = strain.duration - throwaway_size
          
        #number of segments we need to analyze this data                    
        num_segs = int(numpy.ceil(float(analyzable) / float(seg_width)))

        # The offset we will use between segments
        seg_offset = int(numpy.ceil(analyzable / float(num_segs)))
        self.segment_slices = []
        self.analyze_slices = []
        
        # Determine how to chop up the strain into smaller segments
        for nseg in range(num_segs-1):
            # boundaries for time slices into the strain
            seg_start = int((nseg*seg_offset) * strain.sample_rate)
            seg_end = int(seg_start + seg_len * strain.sample_rate)
            seg_slice = slice(seg_start, seg_end)
            self.segment_slices.append(seg_slice)
            
            # boundaries for the analyzable portion of the segment
            ana_start = int(seg_start_pad * strain.sample_rate)
            ana_end = int(ana_start + seg_offset * strain.sample_rate)
            ana_slice = slice(ana_start, ana_end)
            self.analyze_slices.append(ana_slice)
        
        # The last segment takes up any integer boundary slop 
        seg_end = len(strain)
        seg_start = int(seg_end - seg_len * strain.sample_rate)
        seg_slice = slice(seg_start, seg_end)
        self.segment_slices.append(seg_slice)

        remaining = (strain.duration - ((num_segs - 1) * seg_offset + seg_start_pad))
        ana_start = int((seg_len - remaining) * strain.sample_rate)
        ana_end = int((seg_len - seg_end_pad) * strain.sample_rate)
        ana_slice = slice(ana_start, ana_end)
        self.analyze_slices.append(ana_slice)
        
        #Remove segments that are outside trig start and end
        segment_slices_red = []
        analyze_slices_red = []
        trig_start_idx = (trigger_start - int(strain.start_time)) * strain.sample_rate
        trig_end_idx = (trigger_end - int(strain.start_time)) * strain.sample_rate
        
        if filter_inj_only and hasattr(strain, 'injections'):
            inj = strain.injections
            end_times= inj.end_times()
            end_times = [time for time in end_times if float(time) < trigger_end and float(time) > trigger_start]
            
            inj_idx = [(float(time) - float(strain.start_time)) * strain.sample_rate for time in end_times]

        for seg, ana in zip(self.segment_slices, self.analyze_slices):
            start = ana.start
            stop = ana.stop
            cum_start = start + seg.start
            cum_end = stop + seg.start  

            if filter_inj_only and hasattr(strain, 'injections'):
                analyze_this = False
                for inj_id in inj_idx:
                    if inj_id < cum_end and inj_id > cum_start:
                        analyze_this = True
                if not analyze_this:
                    continue     

            # adjust first segment
            if trig_start_idx > cum_start:
                start += (trig_start_idx - cum_start)   
            
            # adjust last segment
            if trig_end_idx < cum_end:
                stop -= (cum_end - trig_end_idx)

            if start < stop:
                segment_slices_red.append(seg)  
                analyze_slices_red.append(slice(start, stop))       
                
        self.segment_slices = segment_slices_red
        self.analyze_slices = analyze_slices_red                   
    
    def fourier_segments(self):
        """ Return a list of the FFT'd segments.
        
        Return the list of FrequencySeries. Additional properties are
        added that describe the strain segment. The property 'analyze' 
        is a slice corresponding to the portion of the time domain equivelant
        of the segment to analyze for triggers. The value 'cumulative_index'
        indexes from the beginning of the original strain series. 
        """
        if not self._fourier_segments:
            self._fourier_segments = []
            for seg_slice, ana in zip(self.segment_slices, self.analyze_slices):
                freq_seg = make_frequency_series(self.strain[seg_slice])
                freq_seg.analyze = ana
                freq_seg.cumulative_index = seg_slice.start + ana.start
                self._fourier_segments.append(freq_seg)  
   
        return self._fourier_segments
        
    @classmethod
    def from_cli(cls, opt, strain):
        """Calculate the segmentation of the strain data for analysis from
        the command line options.
        """
        return cls(strain, segment_length=opt.segment_length, 
                   segment_start_pad=opt.segment_start_pad,
                   segment_end_pad=opt.segment_end_pad,
                   trigger_start=opt.trig_start_time,
                   trigger_end=opt.trig_end_time,
                   filter_inj_only=opt.filter_inj_only)
        
    @classmethod
    def insert_segment_option_group(cls, parser):
        segment_group = parser.add_argument_group(
                                   "Options for segmenting the strain",
                                  "These options are used to determine how to "
                                  "segment the strain into smaller chunks, "
                                  "and for determining the portion of each to "
                                  "analyze for triggers. ")
        segment_group.add_argument("--trig-start-time", type=int, default=0, 
                    help="(optional) The gps time to start recording triggers")
        segment_group.add_argument("--trig-end-time", type=int, default=0,
                    help="(optional) The gps time to stop recording triggers")
        segment_group.add_argument("--segment-length", type=int, 
                          help="The length of each strain segment in seconds.")
        segment_group.add_argument("--segment-start-pad", type=int, 
                          help="The time in seconds to ignore of the "
                               "beginning of each segment in seconds. ")
        segment_group.add_argument("--segment-end-pad", type=int, 
                          help="The time in seconds to ignore at the "
                               "end of each segment in seconds.")
        segment_group.add_argument("--filter-inj-only", action='store_true',
                          help="Analaze only segements that contain an injection.")
        
        
    @classmethod
    def verify_segment_options(cls, opt, parser):
        required_opts(opt, parser, 
                  ['--segment-length',
                   '--segment-start-pad',
                   '--segment-end-pad',
                   ])  
