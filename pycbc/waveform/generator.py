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

import waveform
import ringdown
from pycbc import filter
from pycbc import transforms
from pycbc.types import TimeSeries
from pycbc.waveform import parameters
from pycbc.waveform.utils import apply_fd_time_shift, taper_timeseries
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
        The list of names of variable arguments. Values passed to the
        `generate_from_args` function must be in the same order as the
        arguments in this list.
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
        # keep a list of functions to call before waveform generation
        self._pregenerate_functions = []

    @property
    def static_args(self):
        """Returns a dictionary of the static arguments."""
        return self.frozen_params

    def generate_from_args(self, *args):
        """Generates a waveform. The list of arguments must be in the same
        order as self's variable_args attribute.
        """
        if len(args) != len(self.variable_args):
            raise ValueError("variable argument length mismatch")
        return self.generate(**dict(zip(self.variable_args, args)))

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
        return self.generator(**self.current_params)


class BaseCBCGenerator(BaseGenerator):
    """Adds ability to convert from various derived parameters to parameters
    needed by the waveform generators.

    Attributes
    ----------
    possible_args : set
        The set of names of arguments that may be used in the `variable_args`
        or `frozen_params`.
    """
    possible_args = set(parameters.td_waveform_params +
                        parameters.fd_waveform_params +
                        ['taper'])
    def __init__(self, generator, variable_args=(), **frozen_params):
        super(BaseCBCGenerator, self).__init__(generator,
            variable_args=variable_args, **frozen_params)
        # decorate the generator function with a list of functions that convert
        # parameters to those used by the waveform generation interface
        all_args = set(self.frozen_params.keys() + list(self.variable_args))
        # compare a set of all args of the generator to the input parameters
        # of the functions that do conversions and adds to list of pregenerate
        # functions if it is needed
        params_used, cs = transforms.get_common_cbc_transforms(
                                       list(self.possible_args), variable_args)
        for c in cs:
            self._add_pregenerate(c)
        # check that there are no unused parameters
        unused_args = all_args.difference(params_used) \
                              .difference(self.possible_args)
        if len(unused_args):
            raise ValueError("The following args are not being used: "
                             "{opts}".format(opts=unused_args))


class FDomainCBCGenerator(BaseCBCGenerator):
    """Generates frequency-domain CBC waveforms in the radiation frame.

    Uses `waveform.get_fd_waveform` as a generator function to create
    frequency- domain CBC waveforms in the radiation frame; i.e., with no
    detector response function applied. For more details, see `BaseGenerator`.

    Derived parameters not understood by `get_fd_waveform` may be used as
    variable args and/or frozen parameters, as long as they can be converted
    into parameters that `get_fd_waveform` can use. For example, `mchirp` and
    `eta` (currently, the only supported derived parameters) may be used as
    variable/frozen params; these are converted to `mass1` and `mass2` prior to
    calling the waveform generator function.

    Examples
    --------
    Initialize a generator:

    >>> generator = waveform.FDomainCBCGenerator(variable_args=['mass1', 'mass2'], delta_f=1./32, f_lower=30., approximant='TaylorF2')

    Create a waveform with the variable arguments (in this case, mass1, mass2):

    >>> generator.generate(mass1=1.4, mass2=1.4)
        (<pycbc.types.frequencyseries.FrequencySeries at 0x1110c1450>,
         <pycbc.types.frequencyseries.FrequencySeries at 0x1110c1510>)

    Initialize a generator using mchirp, eta as the variable args, and generate
    a waveform:

    >>> generator = waveform.FDomainCBCGenerator(variable_args=['mchirp', 'eta'], delta_f=1./32, f_lower=30., approximant='TaylorF2')
    >>> generator.generate(mchirp=1.5, eta=0.25)
        (<pycbc.types.frequencyseries.FrequencySeries at 0x109a104d0>,
         <pycbc.types.frequencyseries.FrequencySeries at 0x109a10b50>)

    Note that the `current_params` contains the mchirp and eta values, along
    with the mass1 and mass2 they were converted to:

    >>> generator.current_params
        {'approximant': 'TaylorF2',
         'delta_f': 0.03125,
         'eta': 0.25,
         'f_lower': 30.0,
         'mass1': 1.7230475324955525,
         'mass2': 1.7230475324955525,
         'mchirp': 1.5}

    """
    def __init__(self, variable_args=(), **frozen_params):
        super(FDomainCBCGenerator, self).__init__(waveform.get_fd_waveform,
            variable_args=variable_args, **frozen_params)


class TDomainCBCGenerator(BaseCBCGenerator):
    """Create time domain CBC waveforms in the radiation frame.

    Uses waveform.get_td_waveform as a generator function to create time-
    domain CBC waveforms in the radiation frame; i.e., with no detector
    response function applied. For more details, see `BaseGenerator`.

    Derived parameters not understood by `get_td_waveform` may be used as
    variable args and/or frozen parameters, as long as they can be converted
    into parameters that `get_td_waveform` can use. For example, `mchirp` and
    `eta` (currently, the only supported derived parameters) may be used as
    variable/frozen params; these are converted to `mass1` and `mass2` prior to
    calling the waveform generator function.

    Examples
    --------
    Initialize a generator:

    >>> generator = waveform.TDomainCBCGenerator(variable_args=['mass1', 'mass2'], delta_t=1./4096, f_lower=30., approximant='TaylorT4')

    Create a waveform with the variable arguments (in this case, mass1, mass2):

    >>> generator.generate(mass1=2., mass2=1.3)
        (<pycbc.types.timeseries.TimeSeries at 0x10e546710>,
         <pycbc.types.timeseries.TimeSeries at 0x115f37690>)

    Initialize a generator using mchirp, eta as the variable args, and generate
    a waveform:

    >>> generator = waveform.TDomainCBCGenerator(variable_args=['mchirp', 'eta'], delta_t=1./4096, f_lower=30., approximant='TaylorT4')
    >>> generator.generate(mchirp=1.75, eta=0.2)
        (<pycbc.types.timeseries.TimeSeries at 0x116ac6050>,
         <pycbc.types.timeseries.TimeSeries at 0x116ac6950>)

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


class FDomainRingdownGenerator(BaseGenerator):
    """Uses ringdown.get_fd_qnm as a generator function to create frequency-
    domain ringdown waveforms in the radiation frame; i.e., with no detector response
    function applied. For more details, see BaseGenerator.

    Examples
    --------
    Initialize a generator:

    >>> generator = waveform.FDomainRingdownGenerator(variable_args=['tau', 'f_0', 'amp', 'phi'], delta_f=1./32, f_lower=30., f_final=500)

    Create a ringdown with the variable arguments (in this case, tau, f_0):

    >>> generator.generate(tau=5., f_0=100., amp=1., phi=0.)
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

    >>> generator = waveform.FDomainMultiModeRingdownGenerator(variable_args=['final_mass', 'final_spin', 'lmns','amp220','amp210','phi220','phi210'], delta_f=1./32, f_lower=30., f_final=500)

    Create a ringdown with the variable arguments:

    >>> generator.generate(final_mass=65., final_spin=0.7, lmns=['221','211'], amp220=1e-21, amp210=1./10, phi220=0., phi210=0.)
        (<pycbc.types.frequencyseries.FrequencySeries at 0x51614d0>,
         <pycbc.types.frequencyseries.FrequencySeries at 0x5161550>)

    """
    def __init__(self, variable_args=(), **frozen_params):
        super(FDomainMultiModeRingdownGenerator, self).__init__(ringdown.get_fd_lm_allmodes,
            variable_args=variable_args, **frozen_params)

class FDomainDetFrameGenerator(object):
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

    >>> generator.generate(mass1=38.6, mass2=29.3, spin1z=0.33, spin2z=-0.94, tc=2.43, ra=1.37, dec=-1.26, polarization=2.76)
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

    def generate_from_args(self, *args):
        """Generates a waveform, applies a time shift and the detector response
        function from the given args.

        The args are assumed to be in the same order as the variable args.
        """
        return self.generate(**dict(zip(self.variable_args, args)))

    def generate(self, **kwargs):
        """Generates a waveform, applies a time shift and the detector response
        function from the given kwargs.
        """
        self.current_params.update(kwargs)
        rfparams = {param: self.current_params[param]
            for param in self.rframe_generator.variable_args}
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
        else:
            # no detector response, just use the + polarization
            if 'tc' in self.current_params:
                hp = apply_fd_time_shift(hp, self.current_params['tc']+tshift,
                                         copy=False)
            h['RF'] = hp
        return h


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
    >>> waveform.select_waveform_generator(waveform.fd_approximants()[0])
    """

    # check if frequency-domain CBC waveform
    if approximant in waveform.fd_approximants():
        return FDomainCBCGenerator

    # check if time-domain CBC waveform
    elif approximant in waveform.td_approximants():
        return TDomainCBCGenerator

    # check if frequency-domain ringdown waveform
    elif approximant in ringdown.ringdown_fd_approximants:
        if approximant == 'FdQNM':
            return FDomainRingdownGenerator
        elif approximant == 'FdQNMmultiModes':
            return FDomainMultiModeRingdownGenerator

    # otherwise waveform approximant is not supported
    elif approximant in ringdown.ringdown_td_approximants:
        raise ValueError("Time domain ringdowns not supported")
    else:
        raise ValueError("%s is not a valid approximant." % approximant)
