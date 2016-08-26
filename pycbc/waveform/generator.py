# Copyright (C) 2016  Collin Capano, Alex Nitz
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
This modules provides classes for generating waveforms.
"""

import functools
import waveform
import ringdown
from pycbc.waveform.utils import apply_fd_time_shift
from pycbc.detector import Detector
import lal as _lal

#
#   Generator for CBC waveforms
#

# utility functions/class

class BaseGenerator(object):
    """A wrapper class to call a waveform generator with a set of frozen
    parameters and a set of variable parameters. The frozen parameters and
    values, along with a list of variable parameter names, are set at
    initialization. This way, repeated calls can be made to the underlying
    generator by simply passing a list of values for the variable parameters
    to this class's generate function.

    Parameters
    ----------
    generator : function
        The function that is called for waveform generation.
    variable_args : {(), list}
        A tuple or list of strings giving the names and order of variable
        parameters that will be passed to the waveform generator when the
        generate function is called.
    \**frozen_params :
        These keyword arguments are the ones that will be frozen in the
        waveform generator. For a list of possible parameters, see
        pycbc.waveform.cbc_parameters.

    Attributes
    ----------
    generator : function
        The function that is called for waveform generation.
    variable_args : tuple
        The list of names of variable arguments. Values passed to the generate
        function must be in the same order as the arguments in this list.
    frozen_params : dict
        A dictionary of the frozen keyword arguments that are always passed
        to the waveform generator function.
    current_params : dict
        A dictionary of the frozen keyword arguments and variable arguments
        that were last passed to the waveform generator.

    Methods
    -------
    generate(variable_values)
        Generates a waveform using the variable arguments and the frozen
        arguments.
    """
    def __init__(self, generator, variable_args=(), **frozen_params):
        self.generator = generator
        self.variable_args = tuple(variable_args)
        self.frozen_params = frozen_params
        # we'll keep a dictionary of the current parameters for fast
        # generation
        self.current_params = frozen_params.copy()

    @property
    def static_args(self):
        """Returns a dictionary of the static arguments."""
        return self.frozen_params

    def generate(self, *args):
        """Generates a waveform. The list of arguments must be in the same
        order as self's variable_args attribute.
        """
        if len(args) != len(self.variable_args):
            raise ValueError("variable argument length mismatch")
        self.current_params.update(dict(zip(self.variable_args, args)))
        return self.generator(**self.current_params)

    def generate_from_kwargs(self, **kwargs):
        """Generates a waveform from the keyword args. The current params
        are updated with the given kwargs, then the generator is called.
        """
        self.current_params.update(kwargs)
        return self.generator(**self.current_params)


class FDomainCBCGenerator(BaseGenerator):
    """Uses waveform.get_fd_waveform as a generator function to create frequency-
    domain CBC waveforms in the radiation frame; i.e., with no detector response
    function applied. For more details, see BaseGenerator.

    Examples
    --------
    Initialize a generator:
    >>> generator = waveform.FDomainCBCGenerator(variable_args=['mass1', 'mass2'], delta_f=1./32, f_lower=30., approximant='TaylorF2')

    Create a waveform with the variable arguments (in this case, mass1, mass2):
    >>> generator.generate(1.4, 1.4)
    (<pycbc.types.frequencyseries.FrequencySeries at 0x1110c1450>,
     <pycbc.types.frequencyseries.FrequencySeries at 0x1110c1510>)
    """
    def __init__(self, variable_args=(), **frozen_params):
        super(FDomainCBCGenerator, self).__init__(waveform.get_fd_waveform,
            variable_args=variable_args, **frozen_params)

class TDomainCBCGenerator(BaseGenerator):
    """Uses waveform.get_td_waveform as a generator function to create time-
    domain CBC waveforms in the radiation frame; i.e., with no detector response
    function applied. For more details, see BaseGenerator.

    Examples
    --------
    Initialize a generator:
    >>> generator = waveform.TDomainCBCGenerator(variable_args=['mass1', 'mass2'], delta_t=1./4096, f_lower=30., approximant='TaylorT4')

    Create a waveform with the variable arguments (in this casee, mass1, mass2):
    >>> generator.generate(2., 1.3)
    (<pycbc.types.timeseries.TimeSeries at 0x10e546710>,
     <pycbc.types.timeseries.TimeSeries at 0x115f37690>)
    """
    def __init__(self, variable_args=(), **frozen_params):
        super(TDomainCBCGenerator, self).__init__(waveform.get_td_waveform,
            variable_args=variable_args, **frozen_params)

class FDomainRingdownGenerator(BaseGenerator):
    """Uses ringdown.get_fd_qnm as a generator function to create frequency-
    domain ringdown waveforms in the radiation frame; i.e., with no detector response
    function applied. For more details, see BaseGenerator.

    Examples
    --------
    Initialize a generator:
    >>> generator = waveform.FDomainRingdownGenerator(variable_args=['tau', 'f_0'], delta_f=1./32, f_lower=30., f_final=500)

    Create a ringdown with the variable arguments (in this case, tau, f_0):
    >>> generator.generate(5, 100)
    (<pycbc.types.frequencyseries.FrequencySeries at 0x1110c1450>,
     <pycbc.types.frequencyseries.FrequencySeries at 0x1110c1510>)
    """
    def __init__(self, variable_args=(), **frozen_params):
        super(FDomainRingdownGenerator, self).__init__(ringdown.get_fd_qnm,
            variable_args=variable_args, **frozen_params)

class FDomainMultiModeRingdownGenerator(BaseGenerator):
    """Uses ringdown.get_fd_lm_allmodes as a generator function to create 
    frequency-domain ringdown waveforms with higher modes in the radiation 
    frame; i.e., with no detector response function applied. 
    For more details, see BaseGenerator.

    Examples
    --------
    Initialize a generator:
    >>> generator = waveform.FDomainMultiModeRingdownGenerator(
            variable_args=['final_mass', 'final_spin', 'lmns','amp220','amp210','phi220','phi210'],
            delta_f=1./32, f_lower=30., f_final=500)

    Create a ringdown with the variable arguments:
    >>> generator.generate(65., 0.7, ['221','211'], 1e-21, 0.1e-21, 0., 0.)
    (<pycbc.types.frequencyseries.FrequencySeries at 0x51614d0>,
    <pycbc.types.frequencyseries.FrequencySeries at 0x5161550>)
    """
    def __init__(self, variable_args=(), **frozen_params):
        super(FDomainMultiModeRingdownGenerator, self).__init__(ringdown.get_fd_lm_allmodes,
            variable_args=variable_args, **frozen_params)

class FDomainDetFrameGenerator(object):
    """Generates a waveform using the given radiation frame generator class,
    and applies the detector response function and appropriate time offset.

    Parameters
    ----------
    rFrameGeneratorClass : class
        The class to use for generating the waveform in the radiation frame,
        e.g., FDomainCBCGenerator. This should be the class, not an
        instance of the class (the class will be initialized with the
        appropriate arguments internally).
    detectors : {None, list of strings}
        The names of the detectors to use. If provided, all location parameters
        must be included in either the variable args or the frozen params. If
        None, the generate function will just return the plus polarization
        returned by the rFrameGeneratorClass shifted by any desired time shift.
    epoch : {float, lal.LIGOTimeGPS
        The epoch start time to set the waveform to. A time shift = tc - epoch is
        applied to waveforms before returning.
    variable_args : {(), list or tuple}
        A list or tuple of strings giving the names and order of parameters
        that will be passed to the generate function.
    \**frozen_params
        Keyword arguments setting the parameters that will not be changed from
        call-to-call of the generate function.

    Class Attributes
    ----------------
    location_args : set(['tc', 'ra', 'dec', 'polarization'])
        The set of location parameters. These are not passed to the rFrame
        generator class; instead, they are used to apply the detector response
        function and/or shift the waveform in time. The parameters are:
          * tc: The GPS time of coalescence (should be geocentric time).
          * ra: Right ascension.
          * dec: declination
          * polarization: polarization.
        All of these must be provided in either the variable args or the
        frozen params if detectors is not None. If detectors
        is None, tc may optionally be provided.

    Attributes
    ----------
    detectors : dict
        The dictionary of detectors that antenna patterns are calculated for
        on each call of generate. If no detectors were provided, will be
        ``{'RF': None}``, where "RF" means "radiation frame".
    detector_names : list
        The list of detector names. If no detectors were provided, then this
        will be ['RF'] for "radiation frame".
    epoch : lal.LIGOTimeGPS
        The GPS start time of the frequency series returned by the generate function.
        A time shift is applied to the waveform equal to tc-epoch. Update by using
        ``set_epoch``.
    current_params : dict
        A dictionary of name, value pairs of the arguments that were last
        used by the generate function.
    rframe_generator : instance of rFrameGeneratorClass
        The instance of the radiation-frame generator that is used for waveform
        generation. All parameters in current_params except for the
        location params are passed to this class's generate function.
    frozen_location_args : dict
        Any location parameters that were included in the frozen_params.
    variable_args : tuple
        The list of names of arguments that are passed to the generate
        function.

    Examples
    --------
    Initialize a generator:
    >>> generator = waveform.FDomainDetFrameGenerator(waveform.FDomainCBCGenerator, 0., variable_args=['mass1', 'mass2', 'spin1z', 'spin2z', 'tc', 'ra', 'dec', 'polarization'], detectors=['H1', 'L1'], delta_f=1./64, f_lower=20., approximant='SEOBNRv2_ROM_DoubleSpin')

    Generate a waveform:
    >>> generator.generate(38.6, 29.3, 0.33, -0.94, 2.43, 1.37, -1.26, 2.76)
    {'H1': <pycbc.types.frequencyseries.FrequencySeries at 0x116637350>,
     'L1': <pycbc.types.frequencyseries.FrequencySeries at 0x116637a50>}
    """

    location_args = set(['tc', 'ra', 'dec', 'polarization'])

    def __init__(self, rFrameGeneratorClass, epoch, detectors=None,
            variable_args=(), **frozen_params):
        # initialize frozen & current parameters:
        self.current_params = frozen_params.copy()
        self._static_args = frozen_params.copy()
        # we'll separate out frozen location parameters from the frozen
        # parameters that are sent to the rframe generator
        self.frozen_location_args = {}
        loc_params = set(frozen_params.keys()) & self.location_args
        for param in loc_params:
            self.frozen_location_args[param] = frozen_params.pop(param)
        # set the order of the variable parameters
        self.variable_args = tuple(variable_args)
        # variables that are sent to the rFrame generator
        rframe_variables = list(set(self.variable_args) - self.location_args)
        # initialize the radiation frame generator
        self.rframe_generator = rFrameGeneratorClass(
            variable_args=rframe_variables, **frozen_params)
        self.set_epoch(epoch)
        # if detectors are provided, convert to detector type; also ensure that
        # location variables are specified
        if detectors is not None:
            # FIXME: use the following when we switch to 2.7
            #self.detectors = {det: Detector(det) for det in detectors}
            self.detectors = dict([(det, Detector(det)) for det in detectors])
            missing_args = [arg for arg in self.location_args if not
                (arg in self.current_params or arg in self.variable_args)]
            if any(missing_args):
                raise ValueError("detectors provided, but missing location "
                    "parameters %s. " %(', '.join(missing_args)) +
                    "These must be either in the frozen params or the "
                    "variable args.")
        else:
            self.detectors = {'RF': None}
        self.detector_names = sorted(self.detectors.keys())

    def set_epoch(self, epoch):
        """Sets the epoch; epoch should be a float or a LIGOTimeGPS."""
        self._epoch = float(epoch)

    @property
    def static_args(self):
        """Returns a dictionary of the static arguments."""
        return self._static_args

    @property
    def epoch(self):
        return _lal.LIGOTimeGPS(self._epoch)

    def generate(self, *args):
        """Generates a waveform, applies a time shift and the detector response
        function."""
        self.current_params.update(dict(zip(self.variable_args, args)))
        # FIXME: use the following when we switch to 2.7
        #rfparams = {param: self.current_params[param]
        #    for param in self.rframe_generator.variable_args}
        rfparams = dict([(param, self.current_params[param])
            for param in self.rframe_generator.variable_args])
        hp, hc = self.rframe_generator.generate_from_kwargs(**rfparams)
        hp._epoch = hc._epoch = self._epoch
        h = {}
        if 'tc' in self.current_params:
            try:
                kmin = int(self.current_params['f_lower']/hp.delta_f)
            except KeyError:
                kmin = 0
        if self.detector_names != ['RF']:
            for detname, det in self.detectors.items():
                # apply detector response function
                fp, fc = det.antenna_pattern(self.current_params['ra'],
                            self.current_params['dec'],
                            self.current_params['polarization'],
                            self.current_params['tc'])
                thish = fp*hp + fc*hc
                # apply the time shift
                tc = self.current_params['tc'] + \
                    det.time_delay_from_earth_center(self.current_params['ra'],
                                                     self.current_params['dec'],
                                                     self.current_params['tc'])
                h[detname] = apply_fd_time_shift(thish, tc, kmin=kmin, copy=False)
        else:
            # no detector response, just use the + polarization
            if 'tc' in self.current_params:
                hp = apply_fd_time_shift(hp, self.current_params['tc'],
                            kmin=kmin, copy=False)
            h['RF'] = hp
        return h
