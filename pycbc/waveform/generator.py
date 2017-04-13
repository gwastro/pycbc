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

import functools
import waveform
import ringdown
from pycbc import coordinates
from pycbc import filter
from pycbc.types import TimeSeries
from pycbc.waveform import parameters
from pycbc.waveform.utils import apply_fd_time_shift, taper_timeseries, \
    frequency_from_polarizations
from pycbc.detector import Detector
from pycbc import pnutils
import lal as _lal
import numpy
from scipy import signal

#
#   Pregenerator functions for generator
#


def add_attrs(**func_attrs):
    """ Decorator for adding attributes to a function.
    """
    def add_attr_decorator(fn):
        # wraps is a decorator for updating the attributes of the
        # wrapping function to those of the original function
        # ie. __name__, __doc__, and __module__ will point to the
        # original function instead of the wrapped function
        @functools.wraps(fn)
        def wrapper(*args, **kwargs):
            return fn(*args, **kwargs)
        # update with attributes
        for attr, value in func_attrs.iteritems():
            setattr(wrapper, attr, value)
        return wrapper
    return add_attr_decorator


@add_attrs(input_params=[parameters.mchirp, parameters.eta],
           output_params=[parameters.mass1, parameters.mass2])
def generator_mchirp_eta_to_mass1_mass2(generator):
    """Converts mchirp and eta in `current_params`, to mass1 and mass2.
    """
    mchirp = generator.current_params['mchirp']
    eta = generator.current_params['eta']
    m1, m2 = pnutils.mchirp_eta_to_mass1_mass2(mchirp, eta)
    generator.current_params['mass1'] = m1
    generator.current_params['mass2'] = m2


@add_attrs(input_params=[parameters.mtotal, parameters.eta],
           output_params=[parameters.mass1, parameters.mass2])
def generator_mtotal_eta_to_mass1_mass2(generator):
    """Converts mtotal and eta in `current_params`, to mass1 and mass2.
    """
    mtotal = generator.current_params['mtotal']
    eta = generator.current_params['eta']
    m1, m2 = pnutils.mtotal_eta_to_mass1_mass2(mtotal, eta)
    generator.current_params['mass1'] = m1
    generator.current_params['mass2'] = m2


@add_attrs(input_params=[parameters.mchirp, parameters.q],
           output_params=[parameters.mass1, parameters.mass2])
def generator_mchirp_q_to_mass1_mass2(generator):
    """Converts mchirp and q in `current_params`, to mass1 and mass2.
    """
    mchirp = generator.current_params['mchirp']
    q = generator.current_params['q']
    m1, m2 = pnutils.mchirp_q_to_mass1_mass2(mchirp, q)
    generator.current_params['mass1'] = m1
    generator.current_params['mass2'] = m2


@add_attrs(input_params=[parameters.spin1_a, parameters.spin2_a,
                         parameters.spin1_azimuthal,
                         parameters.spin2_azimuthal,
                         parameters.spin1_polar, parameters.spin2_polar],
           output_params=[parameters.spin1x, parameters.spin2x,
                          parameters.spin1y, parameters.spin2y,
                          parameters.spin1z, parameters.spin2z])
def generator_spin_spherical_to_spin_cartesian(generator):
    """Converts spherical spin magnitude and angles in `current_params`,
    to cartesian component spins.
    """
    x, y, z = coordinates.spherical_to_cartesian(
                                   generator.current_params["spin1_a"],
                                   generator.current_params["spin1_azimuthal"],
                                   generator.current_params["spin1_polar"])
    generator.current_params["spin1x"] = x
    generator.current_params["spin1y"] = y
    generator.current_params["spin1z"] = z
    x, y, z = coordinates.spherical_to_cartesian(
                                   generator.current_params["spin2_a"],
                                   generator.current_params["spin2_azimuthal"],
                                   generator.current_params["spin2_polar"])
    generator.current_params["spin2x"] = x
    generator.current_params["spin2y"] = y
    generator.current_params["spin2z"] = z


# a list of all generator functions
generator_functions = [
    generator_mchirp_eta_to_mass1_mass2,
    generator_mtotal_eta_to_mass1_mass2,
    generator_mchirp_q_to_mass1_mass2,
    generator_spin_spherical_to_spin_cartesian,
]

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
        # keep a list of functions to call before waveform generation
        self._pregenerate_functions = []

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
        return self.generate_from_kwargs(**dict(zip(self.variable_args, args)))

    def generate_from_kwargs(self, **kwargs):
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
                func(self)
            res = generate_func(self)
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
                        ['t_final'])
    def __init__(self, generator, variable_args=(), **frozen_params):
        super(BaseCBCGenerator, self).__init__(generator,
            variable_args=variable_args, **frozen_params)
        # decorate the generator function with a list of functions that convert
        # parameters to those used by the waveform generation interface
        all_args = set(self.frozen_params.keys() + list(self.variable_args))
        # compare a set of all args of the generator to the input parameters
        # of the functions that do conversions and adds to list of pregenerate
        # functions if it is needed
        params_used = set([])
        for func in generator_functions:
            if set(func.input_params).issubset(all_args):
                self._add_pregenerate(func)
                params_used.update(func.input_params)
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

    >>> generator.generate(1.4, 1.4)
        (<pycbc.types.frequencyseries.FrequencySeries at 0x1110c1450>,
         <pycbc.types.frequencyseries.FrequencySeries at 0x1110c1510>)

    Initialize a generator using mchirp, eta as the variable args, and generate
    a waveform:

    >>> generator = waveform.FDomainCBCGenerator(variable_args=['mchirp', 'eta'], delta_f=1./32, f_lower=30., approximant='TaylorF2')
    >>> generator.generate(1.5, 0.25)
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

    >>> generator.generate(2., 1.3)
        (<pycbc.types.timeseries.TimeSeries at 0x10e546710>,
         <pycbc.types.timeseries.TimeSeries at 0x115f37690>)

    Initialize a generator using mchirp, eta as the variable args, and generate
    a waveform:

    >>> generator = waveform.TDomainCBCGenerator(variable_args=['mchirp', 'eta'], delta_t=1./4096, f_lower=30., approximant='TaylorT4')
    >>> generator.generate(1.75, 0.2)
        (<pycbc.types.timeseries.TimeSeries at 0x116ac6050>,
         <pycbc.types.timeseries.TimeSeries at 0x116ac6950>)

    """
    def __init__(self, taper=None,
                 taper_function=None, taper_function_param=None,
                 taper_duration=None, taper_whitened=False, psd=None,
                 variable_args=(), **frozen_params):
        super(TDomainCBCGenerator, self).__init__(waveform.get_td_waveform,
            variable_args=variable_args, **frozen_params)
        self.use_end_taper = ('t_final' in variable_args or
                              't_final' in self.frozen_params or
                              'f_final' in variable_args or
                              'f_final' in self.frozen_params or
                              'f_final_func' in variable_args or
                              'f_final_func' in self.frozen_params)
        self.taper = taper
        self.taper_size = self.window = None
        if taper_whitened and not (taper_whitened == 1 or taper_whitened == 2):
            raise ValueError("taper_whitened must be either False (taper "
                             "before whitening), 1 (taper after whitening) "
                             "or 2 (taper after overwhitening)")
        self.taper_whitened = taper_whitened
        self.psd = psd
        self.asd = None
        if self.use_end_taper:
            # construct the windowing function
            if taper_duration is None:
                raise ValueError('must provide a taper duration if using '
                                 't_final, f_final, or f_final_func')
            self.taper_size = int(taper_duration /
                                  self.frozen_params['delta_t'])
            if taper_function is None:
                raise ValueError("must provide a taper function if using "
                                 "t_final, f_final, or f_final_func")
            if taper_function_param is not None:
                taper_function = (taper_function, taper_function_param)
            win = signal.get_window(taper_function, 2*self.taper_size)
            self.window = win[self.taper_size:]
            if self.taper_whitened:
                if psd is None:
                    raise ValueError("must provide a psd if tapering "
                                    "(over-)whitened waveform")
                if self.taper_whitened == 1:
                    self.asd = psd**0.5
                    nzidx = numpy.nonzero(self.asd.data)[0]
                else:
                    nzidx = numpy.nonzero(self.psd.data)[0]
                self.whkmin = nzidx[0]
                self.whkmax = nzidx[-1] + 1
            

    def _postgenerate(self, res):
        """Applies a taper if it is in current params.
        """
        hp, hc = res
        if self.use_end_taper:
            startidx, endidx = self.get_end_taper_range(hp, hc)
        if self.taper is not None:
            hp = taper_timeseries(hp, tapermethod=self.taper)
            hc = taper_timeseries(hc, tapermethod=self.taper)
        if self.use_end_taper and self.taper_whitened:
            hp = hp.to_frequencyseries(delta_f=self.psd.delta_f)
            if self.taper_whitened == 1:
                hp[self.whkmin:self.whkmax] /= \
                    self.asd[self.whkmin:self.whkmax]
            else:
                hp[self.whkmin:self.whkmax] /= \
                    self.psd[self.whkmin:self.whkmax]
            hp.data[:self.whkmin] = 0.
            hp.data[self.whkmax:] = 0.
            hp = hp.to_timeseries()
            # hc
            hc = hc.to_frequencyseries(delta_f=self.psd.delta_f)
            if self.taper_whitened == 1:
                hc[self.whkmin:self.whkmax] /= \
                    self.asd[self.whkmin:self.whkmax]
            else:
                hc[self.whkmin:self.whkmax] /= \
                    self.psd[self.whkmin:self.whkmax]
            hc.data[:self.whkmin] = 0.
            hc.data[self.whkmax:] = 0.
            hc = hc.to_timeseries()
        if self.use_end_taper:
            self.apply_end_taper(hp, hc, startidx, endidx)
        return hp, hc

    def get_end_taper_range(self, hp, hc):
        endidx = len(hp)
        if 'f_final_func' in self.current_params:
            ffunc = self.current_params['f_final_func']
            self.current_params['f_final'] = \
                pnutils.named_frequency_cutoffs[ffunc](self.current_params)
        if 'f_final' in self.current_params:
            # estimate frequency as a function of time
            foft = abs(frequency_from_polarizations(hp, hc))
            endidx = numpy.where(foft.data >= self.current_params['f_final']
                                )[0][0]
        # evaluate taper time
        if 't_final' in self.current_params:
            # epoch gives time until coalescence time
            tfinal = self.current_params['t_final'] - hp.start_time
            tendidx = int(numpy.ceil(tfinal/hp.delta_t))
            # pick t or f final, whichever comes first
            endidx = max(min(tendidx, endidx), 0)
        startidx = max(0, endidx - self.taper_size)
        return startidx, endidx

    def apply_end_taper(self, hp, hc, startidx, endidx):
        # cut the waveform off at a specific frequency, if desired
        getlen = endidx - startidx
        hp.data[startidx:endidx] *= self.window[-getlen:]
        hp.data[endidx:] = 0.
        hc.data[startidx:endidx] *= self.window[-getlen:]
        hc.data[endidx:] = 0.


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

    >>> generator.generate(65., 0.7, ['221','211'], 1e-21, 1./10, 0., 0.)
        (<pycbc.types.frequencyseries.FrequencySeries at 0x51614d0>,
         <pycbc.types.frequencyseries.FrequencySeries at 0x5161550>)

    """
    def __init__(self, variable_args=(), **frozen_params):
        super(FDomainMultiModeRingdownGenerator, self).__init__(
            ringdown.get_fd_lm_allmodes,
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
