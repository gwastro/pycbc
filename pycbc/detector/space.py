# -*- coding: UTF-8 -*-
#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
"""
This module provides utilities for simulating the GW response of space-based
observatories.
"""
from abc import ABC, abstractmethod
from pycbc.coordinates.space import TIME_OFFSET_20_DEGREES
from pycbc.types import TimeSeries
import numpy
from numpy import cos, sin
from astropy import constants
import logging

def get_available_space_detectors():
    """List the available space detectors"""
    dets = list(_space_detectors.keys())
    aliases = []
    for i in dets:
        aliases.extend(_space_detectors[i]['aliases'])
    return dets + aliases

def parse_det_name(detector_name):
    """Parse a string into a detector name and TDI channel.
       The input is assumed to look like '{detector name}_{channel name}.'"""
    out = detector_name.split('_', 1)
    det = out[0]
    try:
        chan = out[1]
    except IndexError:
        # detector_name is just the detector, so save channel name as None
        chan = None
    return det, chan

def apply_polarization(hp, hc, polarization):
    """
    Apply polarization rotation matrix.

    Parameters
    ----------
    hp : array
        The plus polarization of the GW.

    hc : array
        The cross polarization of the GW.

    polarization : float
        The SSB polarization angle of the GW in radians.

    Returns
    -------
    (array, array)
        The plus and cross polarizations of the GW rotated by the
        polarization angle.
    """
    cphi = cos(2*polarization)
    sphi = sin(2*polarization)

    hp_ssb = hp*cphi - hc*sphi
    hc_ssb = hp*sphi + hc*cphi

    return hp_ssb, hc_ssb

def check_signal_times(hp, hc, orbit_start_time, orbit_end_time,
                       offset=TIME_OFFSET_20_DEGREES, pad_data=False, t0=1e4):
    """
    Ensure that input signal lies within the provided orbital window. This
    assumes that the start times of hp and hc are relative to the detector
    mission start time.

    Parameters
    ----------
    hp : pycbc.types.TimeSeries
        The plus polarization of the GW.

    hc : pycbc.types.TimeSeries
        The cross polarization of the GW.

    orbit_start_time : float
        SSB start time in seconds of the orbital data. By convention,
        t = 0 corresponds to the mission start time of the detector.

    orbit_end_time : float
        SSB end time in seconds of the orbital data.

    polarization : float (optional)
        The polarization in radians of the GW. Default 0.

    offset : float (optional)
        Time offset in seconds to apply to SSB times to ensure proper
        orientation of the constellation at t=0. Default 7365189.431698299
        for a 20 degree offset from Earth.

    pad_data : bool (optional)
        Flag whether to pad the input GW data with time length t0
        worth of zeros. Default False.

    t0 : float (optional)
        Time duration in seconds by which to pad the data if pad_data
        is True. Default 1e4.

    Returns
    -------
    (pycbc.types.TimeSeries, pycbc.types.TimeSeries)
        The plus and cross polarizations of the GW in the SSB frame,
        padded as requested and/or truncated to fit in the orbital window.
    """
    dt = hp.delta_t

    # apply offsets to wfs
    hp.start_time += offset
    hc.start_time += offset

    # pad the data with zeros
    if pad_data:
        pad_idx = int(t0/dt)
        hp.prepend_zeros(pad_idx)
        hp.append_zeros(pad_idx)
        hc.prepend_zeros(pad_idx)
        hc.append_zeros(pad_idx)

    # make sure signal lies within orbit length
    if hp.duration + hp.start_time > orbit_end_time:
        logging.warning('Time of signal end is greater than end of orbital ' +
                        f'data. Cutting signal at {orbit_end_time}.')
        # cut off data succeeding orbit end time
        end_idx = numpy.argwhere(hp.sample_times.numpy() <= orbit_end_time)[-1][0]
        hp = hp[:end_idx]
        hc = hc[:end_idx]

    if hp.start_time < orbit_start_time:
        logging.warning('Time of signal start is less than start of orbital ' +
                        f'data. Cutting signal at {orbit_start_time}.')
        # cut off data preceding orbit start time
        start_idx = numpy.argwhere(hp.sample_times.numpy() >= orbit_start_time)[0][0]
        hp = hp[start_idx:]
        hc = hc[start_idx:]

    return hp, hc

def cut_channels(tdi_dict, remove_garbage=False, t0=1e4):
    """
    Cut TDI channels if needed.

    Parameters
    ----------
    tdi_dict : dict
        The TDI channels, formatted as a dictionary of TimeSeries arrays
        keyed by the channel label.

    remove_garbage : bool, str (optional)
        Flag whether to remove data from the edges of the channels. If True,
        time length t0 is cut from the start and end. If 'zero', time length
        t0 is zeroed at the start and end. If False, channels are unmodified.
        Default False.

    t0 : float (optional)
        Time in seconds to cut/zero from data if remove_garbage is True/'zero'.
        Default 1e4.
    """
    for chan in tdi_dict.keys():
        if remove_garbage:
            dt = tdi_dict[chan].delta_t
            pad_idx = int(t0/dt)
            if remove_garbage == 'zero':
                # zero the edge data
                tdi_dict[chan][:pad_idx] = 0
                tdi_dict[chan][-pad_idx:] = 0
            elif type(remove_garbage) == bool:
                # cut the edge data
                slc = slice(pad_idx, -pad_idx)
                tdi_dict[chan] = tdi_dict[chan][slc]
            else:
                raise ValueError('remove_garbage arg must be a bool or "zero"')

    return tdi_dict

_space_detectors = {'LISA': {'armlength': 2.5e9,
                             'aliases': ['LISA_A', 'LISA_E', 'LISA_T',
                                         'LISA_X', 'LISA_Y', 'LISA_Z'],
                            },
                   }

class AbsSpaceDet(ABC):
    """
    Abstract base class to set structure for space detector classes.

    Parameters
    ----------
    detector_name : str
        The name of the detector. Accepts any output from
        `get_available_space_detectors`.
    
    reference_time : float (optional)
        The reference time in seconds of the signal in the SSB frame. This is
        defined such that the detector mission start time corresponds to 0.
        Default None.
    """
    def __init__(self, detector_name, reference_time=None, **kwargs):
        self.det, self.chan = parse_det_name(detector_name)
        if detector_name not in get_available_space_detectors():
            raise NotImplementedError('Unrecognized detector. ',
                                      'Currently accepts: ',
                                      f'{get_available_space_detectors()}')
        self.reference_time = reference_time

    @property
    @abstractmethod
    def sky_coords(self):
        """
        List the sky coordinate names for the detector class.
        """
        return

    @abstractmethod
    def project_wave(self, hp, hc, lamb, beta):
        """
        Placeholder for evaluating the TDI channels from the GW projections.
        """
        return


class _LDC_detector(AbsSpaceDet):
    """
    LISA detector modeled using LDC software. Constellation orbits are
    generated using LISA Orbits (https://pypi.org/project/lisaorbits/).
    Link projections are generated using LISA GW Response
    (10.5281/zenodo.6423435). TDI channels are generated using pyTDI
    (10.5281/zenodo.6351736).

    Parameters
    ----------
    detector_name : str
        The name of the detector. Accepts any output from
        `get_available_space_detectors`.
    
    reference_time : float (optional)
        The reference time in seconds of the signal in the SSB frame. This is
        defined such that the detector mission start time corresponds to 0.
        Default None.

    apply_offset : bool (optional)
        Flag whether to shift the times of the input waveforms by
        a given value. Some backends require this such that the
        detector is oriented correctly at t = 0. Default False.

    offset : float (optional)
        The time in seconds by which to offset the input waveform if
        apply_offset is True. Default 7365189.431698299.
    
    orbits : str (optional)
        The constellation orbital data used for generating projections
        and TDI. See self.orbits_init for accepted inputs. Default
        'EqualArmlength'.
    """
    def __init__(self, detector_name, reference_time=None, apply_offset=False,
                 offset=TIME_OFFSET_20_DEGREES,
                 orbits='EqualArmlength', **kwargs):
        super().__init__(detector_name, reference_time, **kwargs)
        assert self.det == 'LISA', 'LDC backend only works with LISA detector'

        # specify whether to apply offsets to GPS times
        if apply_offset:
            self.offset = offset
        else:
            self.offset = 0.

        # orbits properties
        self.orbits = orbits
        self.orbits_start_time = None
        self.orbits_end_time = None

        # waveform properties
        self.dt = None
        self.sample_times = None
        self.start_time = None

        # pre- and post-processing
        self.pad_data = False
        self.remove_garbage = False
        self.t0 = 1e4

        # class initialization
        self.proj_init = None
        self.tdi_init = None
        self.tdi_chan = 'AET'
        if self.chan is not None and self.chan in 'XYZ':
            self.tdi_chan = 'XYZ'

    @property
    def sky_coords(self):
        return 'eclipticlongitude', 'eclipticlatitude'

    def orbits_init(self, orbits, size=316, dt=100000.0, t_init=0.0):
        """
        Initialize the orbital information for the constellation. Defualt args
        generate roughly 1 year worth of data starting at LISA mission
        start time.

        Parameters
        ----------
        orbits : str
            The type of orbit to read in. If "EqualArmlength" or "Keplerian",
            a file is generating using the corresponding method from LISA
            Orbits. Else, the input is treated as a file path following LISA
            Orbits format. Default "EqualArmlength".

        length : int (optional)
            The number of samples to generate if creating a new orbit file.
            Default 316.

        dt : float (optional)
            The time step in seconds to use if generating a new orbit file.
            Default 100000.

        t_init : float (optional)
            The start time in seconds to use if generating a new orbit file.
            Default 0.
        """
        defaults = ['EqualArmlength', 'Keplerian']
        assert type(orbits) == str, ('Must input either a file path as ',
                                     'str, "EqualArmlength", or "Keplerian"')

        # generate a new file
        if orbits in defaults:
            try:
                import lisaorbits
            except ImportError:
                raise ImportError('lisaorbits not found')
            if orbits == 'EqualArmlength':
                o = lisaorbits.EqualArmlengthOrbits()
            if orbits == 'Keplerian':
                o = lisaorbits.KeplerianOrbits()
            o.write('orbits.h5', dt=dt, size=size, t0=t_init, mode='w')
            ofile = 'orbits.h5'
            self.orbits_start_time = t_init
            self.orbits_end_time = t_init + size*dt
            self.orbits = ofile

        # read in from an existing file path
        else:
            import h5py
            ofile = orbits
            with h5py.File(ofile, 'r') as f:
                self.orbits_start_time = f.attrs['t0']
                self.orbits_end_time = self.orbit_start_time + \
                                       f.attrs['dt']*f.attrs['size']

        # add light travel buffer times
        lisa_arm = _space_detectors['LISA']['armlength']
        ltt_au = constants.au.value / constants.c.value
        ltt_arm = lisa_arm / constants.c.value
        self.orbits_start_time += ltt_arm + ltt_au
        self.orbits_end_time += ltt_au

    def strain_container(self, response, orbits=None):
        """
        Read in the necessary link and orbit information for generating TDI
        channels. Replicates the functionality of `pyTDI.Data.from_gws()`.

        Parameters
        ----------
        response : array
            The laser link projections of the GW. Uses get_links output format.

        orbits : str, optional
            The path to the file containing orbital information for the LISA
            constellation. Default to orbits class attribute.

        Returns
        -------
        dict, array
            The arguments and measurements associated with the link and orbital
            data.
        """
        try:
            from pytdi import Data
        except ImportError:
            raise ImportError('pyTDI required for TDI combinations')

        links = ['12', '23', '31', '13', '32', '21']

        # format the measurements from link data
        measurements = {}
        for i, link in enumerate(links):
            measurements[f'isi_{link}'] = response[:, i]
            measurements[f'isi_sb_{link}'] = response[:, i]
            measurements[f'tmi_{link}'] = 0.
            measurements[f'rfi_{link}'] = 0.
            measurements[f'rfi_sb_{link}'] = 0.

        df = 1/self.dt
        t_init = self.orbits_start_time

        # call in the orbital data using pyTDI
        if orbits is None:
            orbits = self.orbits
        return Data.from_orbits(orbits, df, t_init, 'tcb/ltt', **measurements)

    def get_links(self, hp, hc, lamb, beta, polarization):
        """
        Project a radiation frame waveform to the LISA constellation.

        Parameters
        ----------
        hp : pycbc.types.TimeSeries
            The plus polarization of the GW in the radiation frame.

        hc : pycbc.types.TimeSeries
            The cross polarization of the GW in the radiation frame.

        lamb : float
            The ecliptic longitude of the source in the SSB frame.

        beta : float
            The ecliptic latitude of the source in the SSB frame.

        polarization : float (optional)
            The polarization angle of the GW in radians. Default 0.

        Returns
        -------
        ndarray
            The waveform projected to the LISA laser links. Shape is (6, N)
            for input waveforms with N total samples.
        """
        try:
            from lisagwresponse import ReadStrain
        except ImportError:
            raise ImportError('LISA GW Response not found')

        if self.dt is None:
            self.dt = hp.delta_t

        # configure orbits and signal
        self.orbits_init(orbits=self.orbits)
        hp, hc = check_signal_times(hp, hc, self.orbits_start_time,
                                    self.orbits_end_time, offset=self.offset,
                                    pad_data=self.pad_data, t0=self.t0)
        self.start_time = hp.start_time - self.offset
        self.sample_times = hp.sample_times.numpy()

        # apply polarization
        hp, hc = apply_polarization(hp, hc, polarization)

        if self.proj_init is None:
            # initialize the class
            self.proj_init = ReadStrain(self.sample_times, hp, hc,
                                        gw_beta=beta, gw_lambda=lamb,
                                        orbits=self.orbits)
        else:
            # update params in the initialized class
            self.proj_init.gw_beta = beta
            self.proj_init.gw_lambda = lamb
            self.proj_init.set_strain(self.sample_times, hp, hc)

        # project the signal
        wf_proj = self.proj_init.compute_gw_response(self.sample_times,
                                                     self.proj_init.LINKS)

        return wf_proj

    def project_wave(self, hp, hc, lamb, beta, polarization=0,
                     tdi=1.5, tdi_chan=None, pad_data=False,
                     remove_garbage=False, t0=1e4, **kwargs):
        """
        Evaluate the TDI observables.

        The TDI generation requires some startup time at the start and end of
        the waveform, creating erroneous ringing or "garbage" at the edges of
        the signal. By default, this method will cut off a time length t0 from
        the start and end to remove this garbage, which may delete sensitive
        data at the edges of the input strains (e.g., the late inspiral and
        ringdown of a binary merger). Thus, the default output will be shorter
        than the input by (2*t0) seconds. See pad_data and remove_garbage to
        modify this behavior.

        Parameters
        ----------
        hp : pycbc.types.TimeSeries
            The plus polarization of the GW in the radiation frame.

        hc : pycbc.types.TimeSeries
            The cross polarization of the GW in the radiation frame.

        lamb : float
            The ecliptic longitude in the SSB frame.

        beta : float
            The ecliptic latitude in the SSB frame.

        polarization : float
            The polarization angle of the GW in radians.

        tdi : float (optional)
            TDI channel configuration. Accepts 1.5 for 1st generation TDI or
            2 for 2nd generation TDI. Default 1.5.

        tdi_chan : str (optional)
            The TDI observables to calculate. Accepts 'XYZ', 'AET', or 'AE'.
            Default 'AET'.

        pad_data : bool (optional)
            Flag whether to pad the data with time length t0 of zeros at the
            start and end. Default False.

        remove_garbage : bool, str (optional)
            Flag whether to remove gaps in TDI from start and end. If True,
            time length t0 worth of data at the start and end of the waveform
            will be cut from TDI channels. If 'zero', time length t0 worth of
            edge data will be zeroed. If False, TDI channels will not be
            modified. Default False.

        t0 : float (optional)
            Time length in seconds to pad/cut from the start and end of
            the data if pad_data/remove_garbage is True. Default 1e4.

        Returns
        -------
        dict ({str: pycbc.types.TimeSeries})
            The TDI observables as TimeSeries objects keyed by their
            corresponding TDI channel name.
        """
        try:
            from pytdi import michelson
        except ImportError:
            raise ImportError('pyTDI not found')

        # set TDI generation
        if tdi == 1.5:
            X, Y, Z = michelson.X1, michelson.Y1, michelson.Z1
        elif tdi == 2:
            X, Y, Z = michelson.X2, michelson.Y2, michelson.Z2
        else:
            raise ValueError('Unrecognized TDI generation input. ' +
                             'Please input either 1 or 2.')

        # set TDI channels
        if tdi_chan is None:
            tdi_chan = self.tdi_chan

        # generate the Doppler time series
        self.pad_data = pad_data
        self.remove_garbage = remove_garbage
        self.t0 = t0
        response = self.get_links(hp, hc, lamb, beta,
                                  polarization=polarization)

        # load in data using response measurements
        self.tdi_init = self.strain_container(response, self.orbits)

        # generate the XYZ TDI channels
        chanx = X.build(**self.tdi_init.args)(self.tdi_init.measurements)
        chany = Y.build(**self.tdi_init.args)(self.tdi_init.measurements)
        chanz = Z.build(**self.tdi_init.args)(self.tdi_init.measurements)

        # convert to AET if specified
        if tdi_chan == 'XYZ':
            tdi_dict = {'LISA_X': TimeSeries(chanx, delta_t=self.dt,
                                             epoch=self.start_time),
                        'LISA_Y': TimeSeries(chany, delta_t=self.dt,
                                             epoch=self.start_time),
                        'LISA_Z': TimeSeries(chanz, delta_t=self.dt,
                                             epoch=self.start_time)}
        elif tdi_chan == 'AET':
            chana = (chanz - chanx)/numpy.sqrt(2)
            chane = (chanx - 2*chany + chanz)/numpy.sqrt(6)
            chant = (chanx + chany + chanz)/numpy.sqrt(3)
            tdi_dict = {'LISA_A': TimeSeries(chana, delta_t=self.dt,
                                             epoch=self.start_time),
                        'LISA_E': TimeSeries(chane, delta_t=self.dt,
                                             epoch=self.start_time),
                        'LISA_T': TimeSeries(chant, delta_t=self.dt,
                                             epoch=self.start_time)}
        else:
            raise ValueError('Unrecognized TDI channel input. ' +
                             'Please input either "XYZ" or "AET".')

        # processing
        tdi_dict = cut_channels(tdi_dict, remove_garbage=self.remove_garbage,
                                t0=self.t0)
        return tdi_dict


class _FLR_detector(AbsSpaceDet):
    """
    LISA detector modeled using FastLISAResponse. Constellation orbits are
    generated using LISA Analysis Tools (10.5281/zenodo.10930979). Link
    projections and TDI channels are generated using FastLISAResponse
    (https://arxiv.org/abs/2204.06633).

    Parameters
    ----------
    detector_name : str
        The name of the detector. Accepts any output from
        `get_available_space_detectors`.
    
    reference_time : float (optional)
        The reference time in seconds of the signal in the SSB frame. This is
        defined such that the detector mission start time corresponds to 0.
        Default None.

    apply_offset : bool (optional)
        Flag whether to shift the times of the input waveforms by
        a given value. Some backends require this such that the
        detector is oriented correctly at t = 0. Default False.

    offset : float (optional)
        The time in seconds by which to offset the input waveform if
        apply_offset is True. Default 7365189.431698299.

    use_gpu : bool (optional)
        Specify whether to run class on GPU support via CuPy. Default False.

    orbits : str (optional)
        The constellation orbital data used for generating projections
        and TDI. See self.orbits_init for accepted inputs. Default
        'EqualArmlength'.
    """
    def __init__(self, detector_name, reference_time=None, apply_offset=False,
                 offset=TIME_OFFSET_20_DEGREES,
                 orbits='EqualArmlength', use_gpu=False, **kwargs):
        logging.warning('WARNING: FastLISAResponse TDI implementation is a ',
                        'work in progress. Currently unable to reproduce LDC ',
                        'or BBHx waveforms.')
        self.use_gpu = use_gpu
        super().__init__(detector_name, reference_time, **kwargs)
        assert self.det == 'LISA', 'FLR backend only works with LISA detector'

        # specify whether to apply offsets to GPS times
        if apply_offset:
            self.offset = offset
        else:
            self.offset = 0.

        # orbits properties
        self.orbits = orbits
        self.orbits_start_time = None
        self.orbits_end_time = None

        # waveform properties
        self.dt = None
        self.sample_times = None
        self.start_time = None

        # pre- and post-processing
        self.pad_data = False
        self.remove_garbage = False
        self.t0 = 1e4

        # class initialization
        self.tdi_init = None
        self.tdi_chan = 'AET'
        if 'tdi_chan' in kwargs.keys():
            if kwargs['tdi_chan'] is not None and kwargs['tdi_chan'] in 'XYZ':
                self.tdi_chan = 'XYZ'

    @property
    def sky_coords(self):
        return 'eclipticlongitude', 'eclipticlatitude'

    def orbits_init(self, orbits):
        """
        Initialize the orbital information for the constellation.

        Parameters
        ----------
        orbits : str
            The type of orbit to read in. If "EqualArmlength" or "ESA", the
            corresponding Orbits class from LISA Analysis Tools is called.
            Else, the input is treated as a file path following LISA
            Orbits format.
        """
        # if orbits are already a class instance, skip this
        if type(self.orbits) is not (str or None):
            return

        try:
            from lisatools import detector
        except ImportError:
            raise ImportError("LISA Analysis Tools required for FLR orbits")

        # load an orbit from lisatools
        defaults = ['EqualArmlength', 'ESA']
        if orbits in defaults:
            if orbits == 'EqualArmlength':
                o = detector.EqualArmlengthOrbits()
            if orbits == 'ESA':
                o = detector.ESAOrbits()

        # create a new orbits instance for file input
        else:
            class CustomOrbits(detector.Orbits):
                def __init__(self):
                    super().__init__(orbits)
            o = CustomOrbits()

        self.orbits = o
        self.orbits_start_time = self.orbits.t_base[0]
        self.orbits_end_time = self.orbits.t_base[-1]

    def get_links(self, hp, hc, lamb, beta, polarization=0, use_gpu=None):
        """
        Project a radiation frame waveform to the LISA constellation.

        Parameters
        ----------
        hp : pycbc.types.TimeSeries
            The plus polarization of the GW in the radiation frame.

        hc : pycbc.types.TimeSeries
            The cross polarization of the GW in the radiation frame.

        lamb : float
            The ecliptic longitude of the source in the SSB frame.

        beta : float
            The ecliptic latitude of the source in the SSB frame.

        polarization : float (optional)
            The polarization angle of the GW in radians. Default 0.

        use_gpu : bool (optional)
            Flag whether to use GPU support. Default to class input.
            CuPy is required if use_gpu is True.

        Returns
        -------
        ndarray
            The waveform projected to the LISA laser links. Shape is (6, N)
            for input waveforms with N total samples.
        """
        try:
            from fastlisaresponse import pyResponseTDI
        except ImportError:
            raise ImportError('FastLISAResponse not found')

        if self.dt is None:
            self.dt = hp.delta_t

        # configure the orbit and signal
        self.orbits_init(orbits=self.orbits)
        hp, hc = check_signal_times(hp, hc, self.orbits_start_time,
                                    self.orbits_end_time, offset=self.offset,
                                    pad_data=self.pad_data, t0=self.t0)
        self.start_time = hp.start_time - self.offset
        self.sample_times = hp.sample_times.numpy()
        
        # apply polarization
        hp, hc = apply_polarization(hp, hc, polarization)

        # interpolate orbital data to signal sample times
        self.orbits.configure(t_arr=self.sample_times)

        # format wf to hp + i*hc
        hp = hp.numpy()
        hc = hc.numpy()
        wf = hp + 1j*hc

        if use_gpu is None:
            use_gpu = self.use_gpu

        # convert to cupy if needed
        if use_gpu:
            import cupy
            wf = cupy.asarray(wf)

        if self.tdi_init is None:
            # initialize the class
            self.tdi_init = pyResponseTDI(1/self.dt, len(wf),
                                          orbits=self.orbits,
                                          use_gpu=use_gpu)
        else:
            # update params in the initialized class
            self.tdi_init.sampling_frequency = 1/self.dt
            self.tdi_init.num_pts = len(wf)
            self.tdi_init.orbits = self.orbits
            self.tdi_init.use_gpu = use_gpu

        # project the signal
        self.tdi_init.get_projections(wf, lamb, beta, t0=self.t0)
        wf_proj = self.tdi_init.y_gw

        return wf_proj

    def project_wave(self, hp, hc, lamb, beta, polarization=0,
                     tdi=1.5, tdi_chan=None, use_gpu=None, pad_data=False,
                     remove_garbage=False, t0=1e4, **kwargs):
        """
        Evaluate the TDI observables.

        The TDI generation requires some startup time at the start and end of
        the waveform, creating erroneous ringing or "garbage" at the edges of
        the signal. By default, this method will cut off a time length t0 from
        the start and end to remove this garbage, which may delete sensitive
        data at the edges of the input strains (e.g., the late inspiral and
        ringdown of a binary merger). Thus, the default output will be shorter
        than the input by (2*t0) seconds. See pad_data and remove_garbage to
        modify this behavior.

        Parameters
        ----------
        hp : pycbc.types.TimeSeries
            The plus polarization of the GW in the radiation frame.

        hc : pycbc.types.TimeSeries
            The cross polarization of the GW in the radiation frame.

        lamb : float
            The ecliptic longitude in the SSB frame.

        beta : float
            The ecliptic latitude in the SSB frame.

        polarization : float (optional)
            The polarization angle of the GW in radians.

        tdi : float(optional)
            TDI channel configuration. Accepts 1.5 for 1st generation TDI or
            2 for 2nd generation TDI. Default 1.5.

        tdi_chan : str (optional)
            The TDI observables to calculate. Accepts 'XYZ', 'AET', or 'AE'.
            Default 'AET'.

        use_gpu : bool (optional)
            Flag whether to use GPU support. Default False.

        pad_data : bool (optional)
            Flag whether to pad the data with time length t0 of zeros at the
            start and end. Default False.

        remove_garbage : bool, str (optional)
            Flag whether to remove gaps in TDI from start and end. If True,
            time length t0 worth of data at the start and end of the waveform
            will be cut from TDI channels. If 'zero', time length t0 worth of
            edge data will be zeroed. If False, TDI channels will not be
            modified. Default False.

        t0 : float (optional)
            Time length in seconds to pad/cut from the start and end of
            the data if pad_data/remove_garbage is True. Default 1e4.

        Returns
        -------
        dict ({str: pycbc.types.TimeSeries})
            The TDI observables as TimeSeries objects keyed by their
            corresponding TDI channel name.
        """
        # set use_gpu
        if use_gpu is None:
            use_gpu = self.use_gpu

        # generate the Doppler time series
        self.pad_data = pad_data
        self.remove_garbage = remove_garbage
        self.t0 = t0
        self.get_links(hp, hc, lamb, beta, polarization=polarization,
                       use_gpu=use_gpu)

        # set TDI configuration (let FLR handle if not 1 or 2)
        if tdi == 1.5:
            tdi_opt = '1st generation'
        elif tdi == 2:
            tdi_opt = '2nd generation'
        else:
            tdi_opt = tdi

        if tdi_opt != self.tdi_init.tdi:
            # update TDI in existing tdi_init class
            self.tdi_init.tdi = tdi_opt
            self.tdi_init._init_TDI_delays()

        # set TDI channels
        if tdi_chan is None:
            tdi_chan = self.tdi_chan

        if tdi_chan in ['XYZ', 'AET', 'AE']:
            self.tdi_init.tdi_chan = tdi_chan
        else:
            raise ValueError('TDI channels must be one of: XYZ, AET, AE')

        # generate the TDI channels
        tdi_obs = self.tdi_init.get_tdi_delays()

        # processing
        tdi_dict = {}
        for i, chan in enumerate(tdi_chan):
            # save as TimeSeries
            tdi_dict[f'LISA_{chan}'] = TimeSeries(tdi_obs[i], delta_t=self.dt,
                                           epoch=self.start_time)

        tdi_dict = cut_channels(tdi_dict, remove_garbage=self.remove_garbage, 
                                t0=self.t0)
        return tdi_dict


_backends = {'LISA': {'LDC': _LDC_detector,
                      'FLR': _FLR_detector,
                     },
            }

class SpaceDetector(AbsSpaceDet):
    """
    Space-based detector.

    Parameters
    ----------
    detector_name : str
        The name of the detector. Accepts any output from
        `get_available_space_detectors`.
    
    reference_time : float (optional)
        The reference time in seconds of the signal in the SSB frame. This is
        defined such that the detector mission start time corresponds to 0.
        Default None.

    backend : str (optional)
        The backend architecture to use for generating TDI. Accepts 'LDC'
        or 'FLR'. Default 'LDC'.
    """
    def __init__(self, detector_name, reference_time=None, backend='LDC',
                 **kwargs):
        super().__init__(detector_name, reference_time, **kwargs)
        if backend in _backends[self.det].keys():
            c = _backends[self.det][backend]
            self.backend = c(detector_name, reference_time, **kwargs)
        else:
            raise ValueError(f'Detector {self.det} does not support backend ',
                             f'{backend}.This detector accepts: '
                             f'{_backends[self.det].keys()}')

    @property
    def sky_coords(self):
        return self.backend.sky_coords

    def get_links(self, hp, hc, lamb, beta, *args, **kwargs):
        return self.backend.get_links(hp, hc, lamb, beta, *args, **kwargs)

    def project_wave(self, hp, hc, lamb, beta, *args, **kwargs):
        return self.backend.project_wave(hp, hc, lamb, beta, *args, **kwargs)


__all__ = ['get_available_space_detectors', 'SpaceDetector',
           '_space_detectors',]