# Copyright (C) 2016  Collin Capano, Alex Nitz, Christopher Biwer
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
import os
import logging

from abc import (ABCMeta, abstractmethod)

from . import waveform
from .waveform import (FailedWaveformError)
from . import ringdown
from . import supernovae
from . import waveform_modes
from pycbc import transforms
from pycbc.types import TimeSeries
from pycbc.waveform import parameters
from pycbc.waveform.utils import apply_fd_time_shift, taper_timeseries, \
                                 ceilpow2
from pycbc.detector import Detector
from pycbc.pool import use_mpi
import lal as _lal
from pycbc import strain


# utility functions/class
failed_counter = 0

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
    record_failures : boolean
        Store output files containing the parameters of failed waveform
        generation. Default is False.
    \**frozen_params :
        These keyword arguments are the ones that will be frozen in the
        waveform generator. For a list of possible parameters, see
        pycbc.waveform.cbc_parameters.

    Attributes
    ----------
    generator : function
        The function that is called for waveform generation.
    variable_args : tuple
        The list of names of variable arguments.
    frozen_params : dict
        A dictionary of the frozen keyword arguments that are always passed
        to the waveform generator function.
    current_params : dict
        A dictionary of the frozen keyword arguments and variable arguments
        that were last passed to the waveform generator.
    """
    def __init__(self, generator, variable_args=(), record_failures=False,
                 **frozen_params):
        self.generator = generator
        self.variable_args = tuple(variable_args)
        self.frozen_params = frozen_params
        # we'll keep a dictionary of the current parameters for fast
        # generation
        self.current_params = frozen_params.copy()
        # keep a list of functions to call before waveform generation
        self._pregenerate_functions = []

        # If we are under mpi, then failed waveform will be stored by
        # mpi rank to avoid file writing conflicts. We'll check for this
        # upfront
        self.record_failures = (record_failures or
                                ('PYCBC_RECORD_FAILED_WAVEFORMS' in os.environ))
        self.mpi_enabled, _, self.mpi_rank = use_mpi()

    @property
    def static_args(self):
        """Returns a dictionary of the static arguments."""
        return self.frozen_params

    def generate(self, **kwargs):
        """Generates a waveform from the keyword args. The current params
        are updated with the given kwargs, then the generator is called.
        """
        self.current_params.update(kwargs)
        return self._generate_from_current()

    def _add_pregenerate(self, func):
        """ Adds a function that will be called by the generator function
        before waveform generation.
        """
        self._pregenerate_functions.append(func)

    def _postgenerate(self, res):
        """Allows the waveform returned by the generator function to be
        manipulated before returning.
        """
        return res

    def _gdecorator(generate_func):
        """A decorator that allows for seemless pre/post manipulation of
        the waveform generator function.
        """
        def dostuff(self):
            for func in self._pregenerate_functions:
                self.current_params = func(self.current_params)
            res = generate_func(self) # pylint:disable=not-callable
            return self._postgenerate(res)
        return dostuff

    @_gdecorator
    def _generate_from_current(self):
        """Generates a waveform from the current parameters.
        """
        try:
            new_waveform = self.generator(**self.current_params)
            return new_waveform
        except RuntimeError as e:
            if self.record_failures:
                import h5py
                from pycbc.io.hdf import dump_state

                global failed_counter

                if self.mpi_enabled:
                    outname = 'failed/params_%s.hdf' % self.mpi_rank
                else:
                    outname = 'failed/params.hdf'

                if not os.path.exists('failed'):
                    os.makedirs('failed')

                with h5py.File(outname) as f:
                    dump_state(self.current_params, f,
                               dsetname=str(failed_counter))
                    failed_counter += 1

            # we'll get a RuntimeError if lalsimulation failed to generate
            # the waveform for whatever reason
            strparams = ' | '.join(['{}: {}'.format(
                p, str(val)) for p, val in self.current_params.items()])
            raise FailedWaveformError("Failed to generate waveform with "
                                      "parameters:\n{}\nError was: {}"
                                      .format(strparams, e))


class BaseCBCGenerator(BaseGenerator):
    """Adds ability to convert from various derived parameters to parameters
    needed by the waveform generators.
    """

    possible_args = set(parameters.td_waveform_params +
                        parameters.fd_waveform_params +
                        ['taper'])
    """set: The set of names of arguments that may be used in the
        `variable_args` or `frozen_params`.
    """

    def __init__(self, generator, variable_args=(), **frozen_params):
        super(BaseCBCGenerator, self).__init__(generator,
            variable_args=variable_args, **frozen_params)
        # decorate the generator function with a list of functions that convert
        # parameters to those used by the waveform generation interface
        all_args = set(list(self.frozen_params.keys()) +
                       list(self.variable_args))
        # check that there are no unused (non-calibration) parameters
        calib_args = set([a for a in self.variable_args if
                          a.startswith('calib_')])
        all_args = all_args - calib_args
        unused_args = all_args - self.possible_args
        if len(unused_args):
            logging.warning("WARNING: The following parameters are generally "
                            "not used by CBC waveform generators: %s. If you "
                            "have provided a transform that converted these "
                            "into known parameters (e.g., mchirp, q to "
                            "mass1, mass2) or you are using a custom model "
                            "that uses these parameters, you can safely "
                            "ignore this message.", ', '.join(unused_args))


class FDomainCBCGenerator(BaseCBCGenerator):
    """Generates frequency-domain CBC waveforms in the radiation frame.

    Uses `waveform.get_fd_waveform` as a generator function to create
    frequency- domain CBC waveforms in the radiation frame; i.e., with no
    detector response function applied. For more details, see `BaseGenerator`.

    Examples
    --------
    Initialize a generator:

    >>> from pycbc.waveform.generator import FDomainCBCGenerator
    >>> generator = FDomainCBCGenerator(variable_args=['mass1', 'mass2'], delta_f=1./32, f_lower=30., approximant='TaylorF2')

    Create a waveform with the variable arguments (in this case, mass1, mass2):

    >>> generator.generate(mass1=1.4, mass2=1.4)
        (<pycbc.types.frequencyseries.FrequencySeries at 0x1110c1450>,
         <pycbc.types.frequencyseries.FrequencySeries at 0x1110c1510>)

    """
    def __init__(self, variable_args=(), **frozen_params):
        super(FDomainCBCGenerator, self).__init__(waveform.get_fd_waveform,
            variable_args=variable_args, **frozen_params)


class FDomainCBCModesGenerator(BaseCBCGenerator):
    """Generates frequency-domain CBC waveform modes.

    Uses :py:func:`waveform_modes.get_fd_waveform_modes` as a generator
    function to create frequency-domain CBC waveforms mode-by-mode, without
    applying spherical harmonics.

    For details, on methods and arguments, see :py:class:`BaseGenerator`.
    """
    def __init__(self, variable_args=(), **frozen_params):
        super(FDomainCBCModesGenerator, self).__init__(
            waveform_modes.get_fd_waveform_modes,
            variable_args=variable_args, **frozen_params)


class TDomainCBCGenerator(BaseCBCGenerator):
    """Create time domain CBC waveforms in the radiation frame.

    Uses waveform.get_td_waveform as a generator function to create time-
    domain CBC waveforms in the radiation frame; i.e., with no detector
    response function applied. For more details, see `BaseGenerator`.

    Examples
    --------
    Initialize a generator:

    >>> from pycbc.waveform.generator import TDomainCBCGenerator
    >>> generator = TDomainCBCGenerator(variable_args=['mass1', 'mass2'], delta_t=1./4096, f_lower=30., approximant='TaylorT4')

    Create a waveform with the variable arguments (in this case, mass1, mass2):

    >>> generator.generate(mass1=2., mass2=1.3)
        (<pycbc.types.timeseries.TimeSeries at 0x10e546710>,
         <pycbc.types.timeseries.TimeSeries at 0x115f37690>)

    """
    def __init__(self, variable_args=(), **frozen_params):
        super(TDomainCBCGenerator, self).__init__(waveform.get_td_waveform,
            variable_args=variable_args, **frozen_params)

    def _postgenerate(self, res):
        """Applies a taper if it is in current params.
        """
        hp, hc = res
        try:
            hp = taper_timeseries(hp, tapermethod=self.current_params['taper'])
            hc = taper_timeseries(hc, tapermethod=self.current_params['taper'])
        except KeyError:
            pass
        return hp, hc


class TDomainCBCModesGenerator(BaseCBCGenerator):
    """Generates time domain CBC waveform modes.

    Uses :py:func:`waveform_modes.get_td_waveform_modes` as a generator
    function to create time-domain CBC waveforms mode-by-mode, without applying
    spherical harmonics. The ``generate`` function returns a dictionary of
    modes -> (real, imag) part of the complex time series.

    For details, on methods and arguments, see :py:class:`BaseGenerator`.
    """
    def __init__(self, variable_args=(), **frozen_params):
        super(TDomainCBCModesGenerator, self).__init__(
            waveform_modes.get_td_waveform_modes,
            variable_args=variable_args, **frozen_params)

    def _postgenerate(self, res):
        """Applies a taper if it is in current params.
        """
        if 'taper' in self.current_params:
            tapermethod = self.current_params['taper']
            for mode in res:
                ulm, vlm = res[mode]
                ulm = taper_timeseries(ulm, tapermethod=tapermethod)
                vlm = taper_timeseries(vlm, tapermethod=tapermethod)
                res[mode] = (ulm, vlm)
        return res


class FDomainMassSpinRingdownGenerator(BaseGenerator):
    """Uses ringdown.get_fd_from_final_mass_spin as a generator function to
    create frequency-domain ringdown waveforms with higher modes in the
    radiation frame; i.e., with no detector response function applied.
    For more details, see BaseGenerator.

    Examples
    --------
    Initialize a generator:

    >>> from pycbc.waveform.generator import FDomainMassSpinRingdownGenerator
    >>> generator = FDomainMassSpinRingdownGenerator(variable_args=['final_mass',
                    'final_spin','amp220','amp210','phi220','phi210'], lmns=['221','211'],
                    delta_f=1./32, f_lower=30., f_final=500)

    Create a ringdown with the variable arguments:

    >>> generator.generate(final_mass=65., final_spin=0.7,
                           amp220=1e-21, amp210=1./10, phi220=0., phi210=0.)
        (<pycbc.types.frequencyseries.FrequencySeries at 0x51614d0>,
         <pycbc.types.frequencyseries.FrequencySeries at 0x5161550>)

    """
    def __init__(self, variable_args=(), **frozen_params):
        super(FDomainMassSpinRingdownGenerator, self).__init__(ringdown.get_fd_from_final_mass_spin,
            variable_args=variable_args, **frozen_params)


class FDomainFreqTauRingdownGenerator(BaseGenerator):
    """Uses ringdown.get_fd_from_freqtau as a generator function to
    create frequency-domain ringdown waveforms with higher modes in the
    radiation frame; i.e., with no detector response function applied.
    For more details, see BaseGenerator.

    Examples
    --------
    Initialize a generator:

    >>> from pycbc.waveform.generator import FDomainFreqTauRingdownGenerator
    >>> generator = FDomainFreqTauRingdownGenerator(variable_args=['f_220',
                    'tau_220','f_210','tau_210','amp220','amp210','phi220','phi210'],
                    lmns=['221','211'], delta_f=1./32, f_lower=30., f_final=500)

    Create a ringdown with the variable arguments:

    >>> generator.generate(f_220=317., tau_220=0.003, f_210=274., tau_210=0.003,
                           amp220=1e-21, amp210=1./10, phi220=0., phi210=0.)
        (<pycbc.types.frequencyseries.FrequencySeries at 0x51614d0>,
         <pycbc.types.frequencyseries.FrequencySeries at 0x5161550>)

    """
    def __init__(self, variable_args=(), **frozen_params):
        super(FDomainFreqTauRingdownGenerator, self).__init__(ringdown.get_fd_from_freqtau,
            variable_args=variable_args, **frozen_params)


class TDomainMassSpinRingdownGenerator(BaseGenerator):
    """Uses ringdown.get_td_from_final_mass_spin as a generator function to
    create time-domain ringdown waveforms with higher modes in the
    radiation frame; i.e., with no detector response function applied.
    For more details, see BaseGenerator.

    Examples
    --------
    Initialize a generator:

    >>> from pycbc.waveform.generator import TDomainMassSpinRingdownGenerator
    >>> generator = TDomainMassSpinRingdownGenerator(variable_args=['final_mass',
                    'final_spin','amp220','amp210','phi220','phi210'], lmns=['221','211'],
                    delta_t=1./2048)

    Create a ringdown with the variable arguments:

    >>> generator.generate(final_mass=65., final_spin=0.7,
                           amp220=1e-21, amp210=1./10, phi220=0., phi210=0.)
        (<pycbc.types.frequencyseries.FrequencySeries at 0x51614d0>,
         <pycbc.types.frequencyseries.FrequencySeries at 0x5161550>)

    """
    def __init__(self, variable_args=(), **frozen_params):
        super(TDomainMassSpinRingdownGenerator, self).__init__(ringdown.get_td_from_final_mass_spin,
            variable_args=variable_args, **frozen_params)


class TDomainFreqTauRingdownGenerator(BaseGenerator):
    """Uses ringdown.get_td_from_freqtau as a generator function to
    create time-domain ringdown waveforms with higher modes in the
    radiation frame; i.e., with no detector response function applied.
    For more details, see BaseGenerator.

    Examples
    --------
    Initialize a generator:

    >>> from pycbc.waveform.generator import FDomainFreqTauRingdownGenerator
    >>> generator = TDomainFreqTauRingdownGenerator(variable_args=['f_220',
                    'tau_220','f_210','tau_210','amp220','amp210','phi220','phi210'],
                    lmns=['221','211'], delta_t=1./2048)

    Create a ringdown with the variable arguments:

    >>> generator.generate(f_220=317., tau_220=0.003, f_210=274., tau_210=0.003,
                           amp220=1e-21, amp210=1./10, phi220=0., phi210=0.)
        (<pycbc.types.frequencyseries.FrequencySeries at 0x51614d0>,
         <pycbc.types.frequencyseries.FrequencySeries at 0x5161550>)

    """
    def __init__(self, variable_args=(), **frozen_params):
        super(TDomainFreqTauRingdownGenerator, self).__init__(ringdown.get_td_from_freqtau,
            variable_args=variable_args, **frozen_params)


class TDomainSupernovaeGenerator(BaseGenerator):
    """Uses supernovae.py to create time domain core-collapse supernovae waveforms
    using a set of Principal Components provided in a .hdf file.
    """
    def __init__(self, variable_args=(), **frozen_params):
        super(TDomainSupernovaeGenerator,
              self).__init__(supernovae.get_corecollapse_bounce,
           variable_args=variable_args, **frozen_params)


#
# =============================================================================
#
#                        Detector-frame generators
#
# =============================================================================
#


class BaseFDomainDetFrameGenerator(metaclass=ABCMeta):
    """Base generator for frquency-domain waveforms in a detector frame.

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

    Attributes
    ----------
    detectors : dict
        The dictionary of detectors that antenna patterns are calculated for
        on each call of generate. If no detectors were provided, will be
        ``{'RF': None}``, where "RF" means "radiation frame".
    detector_names : list
        The list of detector names. If no detectors were provided, then this
        will be ['RF'] for "radiation frame".
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

    """

    location_args = set([])
    """Set: Should be overriden by children classes with a set of parameters
        that set the binary's location.
    """

    def __init__(self, rFrameGeneratorClass, epoch, detectors=None,
                 variable_args=(), recalib=None, gates=None, **frozen_params):
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
        # set calibration model
        self.recalib = recalib
        # if detectors are provided, convert to detector type; also ensure that
        # location variables are specified
        if detectors is not None:
            self.detectors = {det: Detector(det) for det in detectors}
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
        self.gates = gates

    def set_epoch(self, epoch):
        """Sets the epoch; epoch should be a float or a LIGOTimeGPS."""
        self._epoch = float(epoch)

    @property
    def static_args(self):
        """Returns a dictionary of the static arguments."""
        return self._static_args

    @property
    def epoch(self):
        """The GPS start time of the frequency series returned by the generate
        function. A time shift is applied to the waveform equal to tc-epoch.
        Update by using ``set_epoch``
        """
        return _lal.LIGOTimeGPS(self._epoch)

    @abstractmethod
    def generate(self, **kwargs):
        """The function that generates the waveforms.
        """
        pass

    @abstractmethod
    def select_rframe_generator(self, approximant):
        """Method to select waveform generator based on an approximant."""
        pass



class FDomainDetFrameGenerator(BaseFDomainDetFrameGenerator):
    """Generates frequency-domain waveform in a specific frame.

    Generates a waveform using the given radiation frame generator class,
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

    >>> from pycbc.waveform.generator import FDomainDetFrameGenerator
    >>> generator = FDomainDetFrameGenerator(waveform.FDomainCBCGenerator, 0., variable_args=['mass1', 'mass2', 'spin1z', 'spin2z', 'tc', 'ra', 'dec', 'polarization'], detectors=['H1', 'L1'], delta_f=1./64, f_lower=20., approximant='SEOBNRv2_ROM_DoubleSpin')

    Generate a waveform:

    >>> generator.generate(mass1=38.6, mass2=29.3, spin1z=0.33, spin2z=-0.94, tc=2.43, ra=1.37, dec=-1.26, polarization=2.76)
    {'H1': <pycbc.types.frequencyseries.FrequencySeries at 0x116637350>,
     'L1': <pycbc.types.frequencyseries.FrequencySeries at 0x116637a50>}

    """

    location_args = set(['tc', 'ra', 'dec', 'polarization'])
    """set(['tc', 'ra', 'dec', 'polarization']):
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
    """

    def generate(self, **kwargs):
        """Generates a waveform, applies a time shift and the detector response
        function from the given kwargs.
        """
        self.current_params.update(kwargs)
        rfparams = {param: self.current_params[param]
            for param in kwargs if param not in self.location_args}
        hp, hc = self.rframe_generator.generate(**rfparams)
        if isinstance(hp, TimeSeries):
            df = self.current_params['delta_f']
            hp = hp.to_frequencyseries(delta_f=df)
            hc = hc.to_frequencyseries(delta_f=df)
            # time-domain waveforms will not be shifted so that the peak amp
            # happens at the end of the time series (as they are for f-domain),
            # so we add an additional shift to account for it
            tshift = 1./df - abs(hp._epoch)
        else:
            tshift = 0.
        hp._epoch = hc._epoch = self._epoch
        h = {}
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
                         self.current_params['dec'], self.current_params['tc'])
                h[detname] = apply_fd_time_shift(thish, tc+tshift, copy=False)
                if self.recalib:
                    # recalibrate with given calibration model
                    h[detname] = \
                        self.recalib[detname].map_to_adjust(h[detname],
                            **self.current_params)
        else:
            # no detector response, just use the + polarization
            if 'tc' in self.current_params:
                hp = apply_fd_time_shift(hp, self.current_params['tc']+tshift,
                                         copy=False)
            h['RF'] = hp
        if self.gates is not None:
            # resize all to nearest power of 2
            for d in h.values():
                d.resize(ceilpow2(len(d)-1) + 1)
            h = strain.apply_gates_to_fd(h, self.gates)
        return h

    @staticmethod
    def select_rframe_generator(approximant):
        """Returns a radiation frame generator class based on the approximant
        string.
        """
        return select_waveform_generator(approximant)


class FDomainDetFrameTwoPolGenerator(BaseFDomainDetFrameGenerator):
    """Generates frequency-domain waveform in a specific frame.

    Generates both polarizations of a waveform using the given radiation frame
    generator class, and applies the time shift. Detector response functions
    are not applied.

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

    """
    location_args = set(['tc', 'ra', 'dec'])
    """ set(['tc', 'ra', 'dec']):
        The set of location parameters. These are not passed to the rFrame
        generator class; instead, they are used to apply the detector response
        function and/or shift the waveform in time. The parameters are:

          * tc: The GPS time of coalescence (should be geocentric time).
          * ra: Right ascension.
          * dec: declination

        All of these must be provided in either the variable args or the
        frozen params if detectors is not None. If detectors
        is None, tc may optionally be provided.
    """

    def generate(self, **kwargs):
        """Generates a waveform polarizations and applies a time shift.

        Returns
        -------
        dict :
            Dictionary of ``detector names -> (hp, hc)``, where ``hp, hc`` are
            the plus and cross polarization, respectively.
        """
        self.current_params.update(kwargs)
        rfparams = {param: self.current_params[param]
            for param in kwargs if param not in self.location_args}
        hp, hc = self.rframe_generator.generate(**rfparams)
        if isinstance(hp, TimeSeries):
            df = self.current_params['delta_f']
            hp = hp.to_frequencyseries(delta_f=df)
            hc = hc.to_frequencyseries(delta_f=df)
            # time-domain waveforms will not be shifted so that the peak amp
            # happens at the end of the time series (as they are for f-domain),
            # so we add an additional shift to account for it
            tshift = 1./df - abs(hp._epoch)
        else:
            tshift = 0.
        hp._epoch = hc._epoch = self._epoch
        h = {}
        if self.detector_names != ['RF']:
            for detname, det in self.detectors.items():
                # apply the time shift
                tc = self.current_params['tc'] + \
                    det.time_delay_from_earth_center(self.current_params['ra'],
                         self.current_params['dec'], self.current_params['tc'])
                dethp = apply_fd_time_shift(hp, tc+tshift, copy=True)
                dethc = apply_fd_time_shift(hc, tc+tshift, copy=True)
                if self.recalib:
                    # recalibrate with given calibration model
                    dethp = self.recalib[detname].map_to_adjust(
                        dethp, **self.current_params)
                    dethc = self.recalib[detname].map_to_adjust(
                        dethc, **self.current_params)
                h[detname] = (dethp, dethc)
        else:
            # no detector response, just use the + polarization
            if 'tc' in self.current_params:
                hp = apply_fd_time_shift(hp, self.current_params['tc']+tshift,
                                         copy=False)
                hc = apply_fd_time_shift(hc, self.current_params['tc']+tshift,
                                         copy=False)
            h['RF'] = (hp, hc)
        if self.gates is not None:
            # resize all to nearest power of 2
            hps = {}
            hcs = {}
            for det in h:
                hp = h[det]
                hc = h[det]
                hp.resize(ceilpow2(len(hp)-1) + 1)
                hc.resize(ceilpow2(len(hc)-1) + 1)
                hps[det] = hp
                hcs[det] = hc
            hps = strain.apply_gates_to_fd(hps, self.gates)
            hcs = strain.apply_gates_to_fd(hps, self.gates)
            h = {det: (hps[det], hcs[det]) for det in h}
        return h

    @staticmethod
    def select_rframe_generator(approximant):
        """Returns a radiation frame generator class based on the approximant
        string.
        """
        return select_waveform_generator(approximant)

class FDomainDetFrameTwoPolNoRespGenerator(BaseFDomainDetFrameGenerator):
    """Generates frequency-domain waveform in a specific frame.

    Generates both polarizations of a waveform using the given radiation frame
    generator class, and applies the time shift. Detector response functions
    are not applied.

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

    """

    def generate(self, **kwargs):
        """Generates a waveform polarizations

        Returns
        -------
        dict :
            Dictionary of ``detector names -> (hp, hc)``, where ``hp, hc`` are
            the plus and cross polarization, respectively.
        """
        self.current_params.update(kwargs)
        hp, hc = self.rframe_generator.generate(**self.current_params)
        if isinstance(hp, TimeSeries):
            df = self.current_params['delta_f']
            hp = hp.to_frequencyseries(delta_f=df)
            hc = hc.to_frequencyseries(delta_f=df)
            # time-domain waveforms will not be shifted so that the peak amp
            # happens at the end of the time series (as they are for f-domain),
            # so we add an additional shift to account for it
            tshift = 1./df - abs(hp._epoch)
            hp = apply_fd_time_shift(hp, tshift, copy=True)
            hc = apply_fd_time_shift(hc, tshift, copy=True)

        hp._epoch = hc._epoch = self._epoch
        h = {}

        for detname in self.detectors:
            if self.recalib:
                # recalibrate with given calibration model
                hp = self.recalib[detname].map_to_adjust(
                    hp, **self.current_params)
                hc = self.recalib[detname].map_to_adjust(
                    hc, **self.current_params)
            h[detname] = (hp.copy(), hc.copy())
        return h

    @staticmethod
    def select_rframe_generator(approximant):
        """Returns a radiation frame generator class based on the approximant
        string.
        """
        return select_waveform_generator(approximant)

class FDomainDetFrameModesGenerator(BaseFDomainDetFrameGenerator):
    """Generates frequency-domain waveform modes in a specific frame.

    Generates both polarizations of every waveform mode using the given
    radiation frame generator class, and applies the time shift. Detector
    response functions are not applied.

    Parameters
    ----------
    rFrameGeneratorClass : class
        The class to use for generating the waveform modes in the radiation
        frame, e.g., :py:class:`FDomainCBCModesGenerator`. This should be the
        class, not an instance of the class (the class will be initialized with
        the appropriate arguments internally). The class should have a generate
        function that returns a dictionary of waveforms keyed by the modes.
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
        The GPS start time of the frequency series returned by the generate
        function. A time shift is applied to the waveform equal to tc-epoch.
        Update by using ``set_epoch``.
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

    """
    location_args = set(['tc', 'ra', 'dec'])
    """ set(['tc', 'ra', 'dec']):
        The set of location parameters. These are not passed to the rFrame
        generator class; instead, they are used to apply the detector response
        function and/or shift the waveform in time. The parameters are:

          * tc: The GPS time of coalescence (should be geocentric time).
          * ra: Right ascension.
          * dec: declination

        All of these must be provided in either the variable args or the
        frozen params if detectors is not None. If detectors
        is None, tc may optionally be provided.
    """

    def generate(self, **kwargs):
        """Generates and returns a waveform decompsed into separate modes.

        Returns
        -------
        dict :
            Dictionary of ``detector names -> modes -> (ulm, vlm)``, where
            ``ulm, vlm`` are the frequency-domain representations of the real
            and imaginary parts, respectively, of the complex time series
            representation of the ``hlm``.
        """
        self.current_params.update(kwargs)
        rfparams = {param: self.current_params[param]
            for param in kwargs if param not in self.location_args}
        hlms = self.rframe_generator.generate(**rfparams)
        h = {det: {} for det in self.detectors}
        for mode in hlms:
            ulm, vlm = hlms[mode]
            if isinstance(ulm, TimeSeries):
                df = self.current_params['delta_f']
                ulm = ulm.to_frequencyseries(delta_f=df)
                vlm = vlm.to_frequencyseries(delta_f=df)
                # time-domain waveforms will not be shifted so that the peak
                # amplitude happens at the end of the time series (as they are
                # for f-domain), so we add an additional shift to account for
                # it
                tshift = 1./df - abs(ulm._epoch)
            else:
                tshift = 0.
            ulm._epoch = vlm._epoch = self._epoch
            if self.detector_names != ['RF']:
                for detname, det in self.detectors.items():
                    # apply the time shift
                    tc = self.current_params['tc'] + \
                        det.time_delay_from_earth_center(
                            self.current_params['ra'],
                            self.current_params['dec'],
                            self.current_params['tc'])
                    detulm = apply_fd_time_shift(ulm, tc+tshift, copy=True)
                    detvlm = apply_fd_time_shift(vlm, tc+tshift, copy=True)
                    if self.recalib:
                        # recalibrate with given calibration model
                        detulm = self.recalib[detname].map_to_adjust(
                            detulm, **self.current_params)
                        detvlm = self.recalib[detname].map_to_adjust(
                            detvlm, **self.current_params)
                    h[detname][mode] = (detulm, detvlm)
            else:
                # no detector response, just apply time shift
                if 'tc' in self.current_params:
                    ulm = apply_fd_time_shift(ulm,
                                              self.current_params['tc']+tshift,
                                              copy=False)
                    vlm = apply_fd_time_shift(vlm,
                                              self.current_params['tc']+tshift,
                                              copy=False)
                h['RF'][mode] = (ulm, vlm)
            if self.gates is not None:
                # resize all to nearest power of 2
                ulms = {}
                vlms = {}
                for det in h:
                    ulm, vlm = h[det][mode]
                    ulm.resize(ceilpow2(len(ulm)-1) + 1)
                    vlm.resize(ceilpow2(len(vlm)-1) + 1)
                    ulms[det] = ulm
                    vlms[det] = vlm
                ulms = strain.apply_gates_to_fd(ulms, self.gates)
                vlms = strain.apply_gates_to_fd(ulms, self.gates)
                for det in ulms:
                    h[det][mode] = (ulms[det], vlms[det])
        return h

    @staticmethod
    def select_rframe_generator(approximant):
        """Returns a radiation frame generator class based on the approximant
        string.
        """
        return select_waveform_modes_generator(approximant)


#
# =============================================================================
#
#                        Helper functions
#
# =============================================================================
#


def select_waveform_generator(approximant):
    """Returns the single-IFO generator for the approximant.

    Parameters
    ----------
    approximant : str
        Name of waveform approximant. Valid names can be found using
        ``pycbc.waveform`` methods.

    Returns
    -------
    generator : (PyCBC generator instance)
        A waveform generator object.

    Examples
    --------
    Get a list of available approximants:
    >>> from pycbc import waveform
    >>> waveform.fd_approximants()
    >>> waveform.td_approximants()
    >>> from pycbc.waveform import ringdown
    >>> ringdown.ringdown_fd_approximants.keys()

    Get generator object:
    >>> from pycbc.waveform.generator import select_waveform_generator
    >>> select_waveform_generator(waveform.fd_approximants()[0])
    """
    # check if frequency-domain CBC waveform
    if approximant in waveform.fd_approximants():
        return FDomainCBCGenerator
    # check if time-domain CBC waveform
    elif approximant in waveform.td_approximants():
        return TDomainCBCGenerator
    # check if frequency-domain ringdown waveform
    elif approximant in ringdown.ringdown_fd_approximants:
        if approximant == 'FdQNMfromFinalMassSpin':
            return FDomainMassSpinRingdownGenerator
        elif approximant == 'FdQNMfromFreqTau':
            return FDomainFreqTauRingdownGenerator
    elif approximant in ringdown.ringdown_td_approximants:
        if approximant == 'TdQNMfromFinalMassSpin':
            return TDomainMassSpinRingdownGenerator
        elif approximant == 'TdQNMfromFreqTau':
            return TDomainFreqTauRingdownGenerator
    # check if supernovae waveform:
    elif approximant in supernovae.supernovae_td_approximants:
        if approximant == 'CoreCollapseBounce':
            return TDomainSupernovaeGenerator
    # otherwise waveform approximant is not supported
    else:
        raise ValueError("%s is not a valid approximant." % approximant)


def select_waveform_modes_generator(approximant):
    """Returns the single-IFO modes generator for the approximant.

    Parameters
    ----------
    approximant : str
        Name of waveform approximant. Valid names can be found using
        ``pycbc.waveform`` methods.

    Returns
    -------
    generator : (PyCBC generator instance)
        A waveform generator object.
    """
    # check if frequency-domain CBC waveform
    if approximant in waveform.fd_approximants():
        return FDomainCBCModesGenerator
    # check if time-domain CBC waveform
    elif approximant in waveform.td_approximants():
        return TDomainCBCModesGenerator
    # otherwise waveform approximant is not supported
    raise ValueError("%s is not a valid approximant." % approximant)
