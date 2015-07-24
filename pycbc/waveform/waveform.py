# Copyright (C) 2012  Alex Nitz
#
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
"""Convenience functions to genenerate gravitational wave templates and
waveforms.
"""
import lal, lalsimulation
from pycbc.types import TimeSeries, FrequencySeries
from pycbc.types import real_same_precision_as
import pycbc.scheme as _scheme
import inspect
from pycbc.fft import fft
from pycbc import pnutils
from pycbc.waveform import utils as wfutils
import pycbc

default_args = {'spin1x':0, 'spin1y':0, 'spin1z':0, 'spin2x':0, 'spin2y':0,
                'spin2z':0, 'lambda1':0, 'lambda2':0,
                'inclination':0, 'distance':1, 'f_final':0, 'f_ref':0,
                'coa_phase':0, 'amplitude_order':-1, 'phase_order':-1,
                'spin_order':-1, 'tidal_order':-1}

default_sgburst_args = {'eccentricity':0, 'polarization':0}

base_required_args = ['mass1','mass2','f_lower']
td_required_args = base_required_args + ['delta_t']
fd_required_args = base_required_args + ['delta_f']
sgburst_required_args = ['q','frequency','hrss']

# td, fd, filter waveforms generated on the CPU
_lalsim_td_approximants = {}
_lalsim_fd_approximants = {}
_lalsim_enum = {}
_lalsim_sgburst_approximants = {}

def _lalsim_td_waveform(**p):
    flags = lalsimulation.SimInspiralCreateWaveformFlags()
    lalsimulation.SimInspiralSetSpinOrder(flags, p['spin_order'])
    lalsimulation.SimInspiralSetTidalOrder(flags, p['tidal_order'])

    hp, hc = lalsimulation.SimInspiralChooseTDWaveform(float(p['coa_phase']),
               float(p['delta_t']),
               float(pnutils.solar_mass_to_kg(p['mass1'])),
               float(pnutils.solar_mass_to_kg(p['mass2'])),
               float(p['spin1x']), float(p['spin1y']), float(p['spin1z']),
               float(p['spin2x']), float(p['spin2y']), float(p['spin2z']),
               float(p['f_lower']), float(p['f_ref']),
               pnutils.megaparsecs_to_meters(float(p['distance'])),
               float(p['inclination']),
               float(p['lambda1']),  float(p['lambda2']), flags, None,
               int(p['amplitude_order']), int(p['phase_order']),
               _lalsim_enum[p['approximant']])

    hp = TimeSeries(hp.data.data[:]*1, delta_t=hp.deltaT, epoch=hp.epoch)
    hc = TimeSeries(hc.data.data[:]*1, delta_t=hc.deltaT, epoch=hc.epoch)

    return hp, hc

def _lalsim_fd_waveform(**p):
    flags = lalsimulation.SimInspiralCreateWaveformFlags()
    lalsimulation.SimInspiralSetSpinOrder(flags, p['spin_order'])
    lalsimulation.SimInspiralSetTidalOrder(flags, p['tidal_order'])

    hp, hc = lalsimulation.SimInspiralChooseFDWaveform(float(p['coa_phase']),
               float(p['delta_f']),
               float(pnutils.solar_mass_to_kg(p['mass1'])),
               float(pnutils.solar_mass_to_kg(p['mass2'])),
               float(p['spin1x']), float(p['spin1y']), float(p['spin1z']),
               float(p['spin2x']), float(p['spin2y']), float(p['spin2z']),
               float(p['f_lower']), float(p['f_final']), float(p['f_ref']),
               pnutils.megaparsecs_to_meters(float(p['distance'])),
               float(p['inclination']),
               float(p['lambda1']), float(p['lambda2']), flags, None,
               int(p['amplitude_order']), int(p['phase_order']),
               _lalsim_enum[p['approximant']])

    hp = FrequencySeries(hp.data.data[:]*1, delta_f=hp.deltaF,
                            epoch=hp.epoch)

    hc = FrequencySeries(hc.data.data[:]*1, delta_f=hc.deltaF,
                            epoch=hc.epoch)                        
    return hp, hc

def _lalsim_sgburst_waveform(**p):
    hp, hc = lalsimulation.SimBurstSineGaussian(float(p['q']),
               float(p['frequency']),
               float(p['hrss']),
               float(p['eccentricity']),
               float(p['polarization']),
               float(p['delta_t']))

    hp = TimeSeries(hp.data.data[:], delta_t=hp.deltaT, epoch=hp.epoch)
    hc = TimeSeries(hc.data.data[:], delta_t=hc.deltaT, epoch=hc.epoch)

    return hp, hc

for approx_enum in xrange(0, lalsimulation.NumApproximants):
    if lalsimulation.SimInspiralImplementedTDApproximants(approx_enum):
        approx_name = lalsimulation.GetStringFromApproximant(approx_enum)
        _lalsim_enum[approx_name] = approx_enum
        _lalsim_td_approximants[approx_name] = _lalsim_td_waveform

for approx_enum in xrange(0, lalsimulation.NumApproximants):
    if lalsimulation.SimInspiralImplementedFDApproximants(approx_enum):
        approx_name = lalsimulation.GetStringFromApproximant(approx_enum)
        _lalsim_enum[approx_name] = approx_enum
        _lalsim_fd_approximants[approx_name] = _lalsim_fd_waveform

# sine-Gaussian burst
for approx_enum in xrange(0, lalsimulation.NumApproximants):
    if lalsimulation.SimInspiralImplementedFDApproximants(approx_enum):
        approx_name = lalsimulation.GetStringFromApproximant(approx_enum)
        _lalsim_enum[approx_name] = approx_enum
        _lalsim_sgburst_approximants[approx_name] = _lalsim_sgburst_waveform

cpu_sgburst = _lalsim_sgburst_approximants

cpu_td = dict(_lalsim_td_approximants.items())
cpu_fd = _lalsim_fd_approximants

# Waveforms written in CUDA
_cuda_td_approximants = {}
_cuda_fd_approximants = {}

if pycbc.HAVE_CUDA:
    from pycbc.waveform.TaylorF2 import taylorf2 as cuda_taylorf2
    from pycbc.waveform.pycbc_phenomC_tmplt import imrphenomc_tmplt
    from pycbc.waveform.SpinTaylorF2 import spintaylorf2 as cuda_spintaylorf2
    _cuda_fd_approximants["IMRPhenomC"] = imrphenomc_tmplt
    _cuda_fd_approximants["SpinTaylorF2"] = cuda_spintaylorf2

cuda_td = dict(_lalsim_td_approximants.items() + _cuda_td_approximants.items())
cuda_fd = dict(_lalsim_fd_approximants.items() + _cuda_fd_approximants.items())

td_wav = _scheme.ChooseBySchemeDict()
fd_wav = _scheme.ChooseBySchemeDict()
td_wav.update({_scheme.CPUScheme:cpu_td,_scheme.CUDAScheme:cuda_td})
fd_wav.update({_scheme.CPUScheme:cpu_fd,_scheme.CUDAScheme:cuda_fd})
sgburst_wav = {_scheme.CPUScheme:cpu_sgburst}

# List the various available approximants ####################################

def print_td_approximants():
    print("Lalsimulation Approximants")
    for approx in _lalsim_td_approximants.keys():
        print "  " + approx
    print("CUDA Approximants")
    for approx in _cuda_td_approximants.keys():
        print "  " + approx

def print_fd_approximants():
    print("Lalsimulation Approximants")
    for approx in _lalsim_fd_approximants.keys():
        print "  " + approx
    print("CUDA Approximants")
    for approx in _cuda_fd_approximants.keys():
        print "  " + approx

def print_sgburst_approximants():
    print("Lalsimulation Approximants")
    for approx in _lalsim_sgburst_approximants.keys():
        print "  " + approx

def td_approximants(scheme=_scheme.mgr.state):
    """Return a list containing the available time domain approximants for 
       the given processing scheme.
    """
    return td_wav[type(scheme)].keys()

def fd_approximants(scheme=_scheme.mgr.state):
    """Return a list containing the available fourier domain approximants for 
       the given processing scheme.
    """
    return fd_wav[type(scheme)].keys()

def sgburst_approximants(scheme=_scheme.mgr.state):
    """Return a list containing the available time domain sgbursts for
       the given processing scheme.
    """
    return sgburst_wav[type(scheme)].keys()

def filter_approximants(scheme=_scheme.mgr.state):
    """Return a list of fourier domain approximants including those
       written specifically as templates.
    """
    return filter_wav[type(scheme)].keys()

# Input parameter handling ###################################################

def props(obj, **kwargs):
    """ DOCUMENT ME !!
    """
    pr = {}
    if obj is not None:
        if hasattr(obj, '__dict__'):
            pr = obj.__dict__
        elif hasattr(obj, '__slots__'):
            for slot in obj.__slots__:
                if hasattr(obj, slot):
                    pr[slot] = getattr(obj, slot)
        else:
            for name in dir(obj):
                try:
                    value = getattr(obj, name)
                    if not name.startswith('__') and not inspect.ismethod(value):
                        pr[name] = value
                except:
                    continue

    # Get the parameters to generate the waveform
    # Note that keyword arguments override values in the template object
    input_params = default_args.copy()
    input_params.update(pr)
    input_params.update(kwargs)

    return input_params

# Input parameter handling for bursts ########################################

def props_sgburst(obj, **kwargs):
    pr = {}
    if obj is not None:
        for name in dir(obj):
            try:
                value = getattr(obj, name)
                if not name.startswith('__') and not inspect.ismethod(value):
                    pr[name] = value
            except:
                continue

    # Get the parameters to generate the waveform
    # Note that keyword arguments override values in the template object
    input_params = default_sgburst_args.copy()
    input_params.update(pr)
    input_params.update(kwargs)

    return input_params

# Waveform generation ########################################################

def get_td_waveform(template=None, **kwargs):
    """Return the plus and cross polarizations of a time domain waveform. 

    Parameters
    ----------
    template: object
        An object that has attached properties. This can be used to subsitute
        for keyword arguments. A common example would be a row in an xml table. 
    approximant : string
        A string that indicates the chosen approximant. See `td_approximants` 
        for available options. 
    mass1 : float
        The mass of the first component object in the binary in solar masses.
    mass2 : 
        The mass of the second component object in the binary in solar masses.
    delta_t :
        The time step used to generate the waveform. 
    f_lower :
        The starting frequency of the waveform.
    f_ref : {float}, optional
        The reference frequency
    distance : {1, float}, optional
        The distance from the observer to the source in megaparsecs.
    inclination : {0, float}, optional
        The inclination angle of the source. 
    coa_phase : {0, float}, optional
        The final phase or phase at the peak of the wavform. See documentation
        on specific approximants for exact usage. 
    spin1x : {0, float}, optional
        The x component of the first binary component's spin vector.
    spin1y : {0, float}, optional
        y component of the first binary component's spin.
    spin1z : {0, float}, optional
        z component of the first binary component's spin.
    spin2x : {0, float}, optional
        The x component of the second binary component's spin vector.
    spin2y : {0, float}, optional
        y component of the second binary component's spin.
    spin2z : {0, float}, optional
        z component of the second binary component's spin.
    lambda1: {0, float}, optional
        The tidal deformability parameter of object 1.
    lambda2: {0, float}, optional
        The tidal deformability parameter of object 2.
    phase_order: {-1, int}, optional
        The pN order of the orbital phase. The default of -1 indicates that 
        all implemented orders are used.
    spin_order: {-1, int}, optional
        The pN order of the spin corrections. The default of -1 indicates that 
        all implemented orders are used.
    tidal_order: {-1, int}, optional
        The pN order of the tidal corrections. The default of -1 indicates that 
        all implemented orders are used.
    amplitude_order: {-1, int}, optional
        The pN order of the amplitude. The default of -1 indicates that 
        all implemented orders are used.

    Returns
    -------
    hplus: TimeSeries
        The plus polarization of the waveform.
    hcross: TimeSeries
        The cross polarization of the waveform.
    """
    input_params = props(template,**kwargs)
    wav_gen = td_wav[type(_scheme.mgr.state)]

    if 'approximant' not in input_params or input_params['approximant'] is None:
        raise ValueError("Please provide an approximant name")
    elif input_params['approximant'] not in wav_gen:
        raise ValueError("Approximant %s not available" %
                            (input_params['approximant']))

    for arg in td_required_args:
        if arg in input_params:
            pass
        else:
            raise ValueError("Please provide " + str(arg) )

    return wav_gen[input_params['approximant']](**input_params)

def get_fd_waveform(template=None, **kwargs):
    """Return a frequency domain gravitational waveform.

    Parameters
    ----------
    template: object
        An object that has attached properties. This can be used to substitute
        for keyword arguments. A common example would be a row in an xml table. 
    approximant : string
        A string that indicates the chosen approximant. See `fd_approximants` 
        for available options. 
    mass1 : float
        The mass of the first component object in the binary in solar masses.
    mass2 : float
        The mass of the second component object in the binary in solar masses.
    delta_f : float
        The frequency step used to generate the waveform. 
    f_lower : float
        The starting frequency of the waveform.
    f_final : {-1, float}, optional
        The ending frequency of the waveform. The default indicates that the
        choice is made by the respective approximant. 
    f_ref : {float}, optional
        The reference frequency.
    distance : {1, float}, optional
        The distance from the observer to the source in megaparsecs.
    inclination : {0, float}, optional
        The inclination angle of the source. 
    coa_phase : {0, float}, optional
        The final phase or phase at the peak of the waveform. See documentation
        on specific approximants for exact usage. 
    spin1x : {0, float}, optional
        The x component of the first binary component's spin vector.
    spin1y : {0, float}, optional
        y component of the first binary component's spin.
    spin1z : {0, float}, optional
        z component of the first binary component's spin.
    spin2x : {0, float}, optional
        The x component of the second binary component's spin vector.
    spin2y : {0, float}, optional
        y component of the second binary component's spin.
    spin2z : {0, float}, optional
        z component of the second binary component's spin.
    lambda1: {0, float}, optional
        The tidal deformability parameter of object 1.
    lambda2: {0, float}, optional
        The tidal deformability parameter of object 2.
    phase_order: {-1, int}, optional
        The pN order of the orbital phase. The default of -1 indicates that 
        all implemented orders are used.
    spin_order: {-1, int}, optional
        The pN order of the spin corrections. The default of -1 indicates that 
        all implemented orders are used.
    tidal_order: {-1, int}, optional
        The pN order of the tidal corrections. The default of -1 indicates that 
        all implemented orders are used.
    amplitude_order: {-1, int}, optional
        The pN order of the amplitude. The default of -1 indicates that 
        all implemented orders are used.

    Returns
    -------
    hplustilde: FrequencySeries
        The plus phase of the waveform in frequency domain.
    hcrosstilde: FrequencySeries
        The cross phase of the waveform in frequency domain.
    """

    input_params = props(template,**kwargs)
    wav_gen = fd_wav[type(_scheme.mgr.state)]

    if 'approximant' not in input_params:
        raise ValueError("Please provide an approximant name")
    elif input_params['approximant'] not in wav_gen:
        raise ValueError("Approximant %s not available" %
                            (input_params['approximant']))

    for arg in fd_required_args:
        if arg in input_params:
            pass
        else:
            raise ValueError("Please provide " + str(arg) )

    return wav_gen[input_params['approximant']](**input_params)


def get_sgburst_waveform(template=None, **kwargs):
    """Return the plus and cross polarizations of a time domain
    sine-Gaussian burst waveform.

    Parameters
    ----------
    template: object
        An object that has attached properties. This can be used to subsitute
        for keyword arguments. A common example would be a row in an xml table.
    approximant : string
        A string that indicates the chosen approximant. See `td_approximants`
        for available options.
    q : float
        The quality factor of a sine-Gaussian burst
    frequency : float
        The centre-frequency of a sine-Gaussian burst
    delta_t : float
        The time step used to generate the waveform
    hrss : float
        The strain rss
    amplitude: float
        The strain amplitude

    Returns
    -------
    hplus: TimeSeries
        The plus polarization of the waveform.
    hcross: TimeSeries
        The cross polarization of the waveform.
    """
    input_params = props_sgburst(template,**kwargs)

    for arg in sgburst_required_args:
        if arg in input_params:
            pass
        else:
            raise ValueError("Please provide " + str(arg))

    return _lalsim_sgburst_waveform(**input_params)

# Waveform filter routines ###################################################

# Organize Filter Generators
_inspiral_fd_filters = {}
_cuda_fd_filters = {}

from spa_tmplt import spa_tmplt
_cuda_fd_filters['SPAtmplt'] = spa_tmplt
_inspiral_fd_filters['SPAtmplt'] = spa_tmplt

filter_wav = _scheme.ChooseBySchemeDict()
filter_wav.update( {_scheme.CPUScheme:_inspiral_fd_filters,
                    _scheme.CUDAScheme:_cuda_fd_filters,
                   } )

# Organize functions for function conditioning/precalculated values
_filter_norms = {}
_filter_ends = {}
_filter_preconditions = {}
_template_amplitude_norms = {}
_filter_time_lengths = {}

from spa_tmplt import spa_tmplt_norm, spa_tmplt_end, spa_tmplt_precondition, spa_amplitude_factor, spa_length_in_time
_filter_norms["SPAtmplt"] = spa_tmplt_norm
_filter_preconditions["SPAtmplt"] = spa_tmplt_precondition
_filter_ends["SPAtmplt"] = spa_tmplt_end
_template_amplitude_norms["SPAtmplt"] = spa_amplitude_factor
_filter_time_lengths["SPAtmplt"] = spa_length_in_time


def seobnrrom_length_in_time(**kwds):
    """
    This is a stub for holding the calculation for getting length of the ROM
    waveforms.
    """
    mass1 = kwds['mass1']
    mass2 = kwds['mass2']
    spin1z = kwds['spin1z']
    spin2z = kwds['spin2z']
    fmin = kwds['f_lower']
    chi = lalsimulation.SimIMRPhenomBComputeChi(mass1, mass2, spin1z, spin2z)
    time_length = lalsimulation.SimIMRSEOBNRv2ChirpTimeSingleSpin(
                               mass1*lal.MSUN_SI, mass2*lal.MSUN_SI, chi, fmin)
    # FIXME: This is still approximate so add a 10% error margin
    time_length = time_length * 1.1
    return time_length

_filter_time_lengths["SEOBNRv1_ROM_SingleSpin"] = seobnrrom_length_in_time
_filter_time_lengths["SEOBNRv1_ROM_DoubleSpin"] = seobnrrom_length_in_time
_filter_time_lengths["SEOBNRv2_ROM_SingleSpin"] = seobnrrom_length_in_time
_filter_time_lengths["SEOBNRv2_ROM_DoubleSpin"] = seobnrrom_length_in_time


def get_waveform_filter(out, template=None, **kwargs):
    """Return a frequency domain waveform filter for the specified approximant
    """
    n = len(out)

    input_params = props(template, **kwargs)

    if input_params['approximant'] in filter_approximants(_scheme.mgr.state):
        wav_gen = filter_wav[type(_scheme.mgr.state)]
        htilde = wav_gen[input_params['approximant']](out=out, **input_params)
        htilde.resize(n)
        htilde.chirp_length = get_waveform_filter_length_in_time(**input_params)
        htilde.length_in_time = htilde.chirp_length
        return htilde

    if input_params['approximant'] in fd_approximants(_scheme.mgr.state):
        wav_gen = fd_wav[type(_scheme.mgr.state)]
        hp, hc = wav_gen[input_params['approximant']](**input_params)
        hp.resize(n)
        out[0:len(hp)] = hp[:]
        hp = FrequencySeries(out, delta_f=hp.delta_f, copy=False)
        hp.chirp_length = get_waveform_filter_length_in_time(**input_params)
        hp.length_in_time = hp.chirp_length
        return hp

    elif input_params['approximant'] in td_approximants(_scheme.mgr.state):
        # N: number of time samples required
        N = (n-1)*2
        delta_f = 1.0 / (N * input_params['delta_t'])
        wav_gen = td_wav[type(_scheme.mgr.state)]
        hp, hc = wav_gen[input_params['approximant']](**input_params)
        # taper the time series hp if required
        if ('taper' in input_params.keys() and \
            input_params['taper'] is not None):
            hp = wfutils.taper_timeseries(hp, input_params['taper'],
                                          return_lal=False)
        # total duration of the waveform
        tmplt_length = len(hp) * hp.delta_t
        # for IMR templates the zero of time is at max amplitude (merger)
        # thus the start time is minus the duration of the template from
        # lower frequency cutoff to merger, i.e. minus the 'chirp time'
        tChirp = - float( hp.start_time )  # conversion from LIGOTimeGPS
        hp.resize(N)
        k_zero = int(hp.start_time / hp.delta_t)
        hp.roll(k_zero)
        htilde = FrequencySeries(out, delta_f=delta_f, copy=False)
        fft(hp.astype(real_same_precision_as(htilde)), htilde)
        htilde.length_in_time = tmplt_length
        htilde.chirp_length = tChirp
        return htilde

    else:
        raise ValueError("Approximant %s not available" %
                            (input_params['approximant']))

def waveform_norm_exists(approximant):
    if approximant in _filter_norms:
        return True
    else:
        return False

def get_template_amplitude_norm(template=None, **kwargs):
    """ Return additional constant template normalization. This only affects
        the effective distance calculation. Returns None for all templates with a
        physically meaningful amplitude. 
    """
    input_params = props(template,**kwargs)
    approximant = kwargs['approximant']

    if approximant in _template_amplitude_norms:
        return _template_amplitude_norms[approximant](**input_params)
    else:
        return None

def get_waveform_filter_precondition(approximant, length, delta_f):
    """Return the data preconditioning factor for this approximant.
    """
    if approximant in _filter_preconditions:
        return _filter_preconditions[approximant](length, delta_f)
    else:
        return None

def get_waveform_filter_norm(approximant, psd, length, delta_f, f_lower):
    """ Return the normalization vector for the approximant 
    """
    if approximant in _filter_norms:
        return _filter_norms[approximant](psd, length, delta_f, f_lower)
    else:
        return None

def get_waveform_end_frequency(template=None, **kwargs):
   """Return the stop frequency of a template
   """
   input_params = props(template,**kwargs)
   approximant = kwargs['approximant']

   if approximant in _filter_ends:
        return _filter_ends[approximant](**input_params)
   else:
        return None

def get_waveform_filter_length_in_time(approximant,**kwargs):
    """For filter templates, return the length in time of the template.
    """
    if approximant in _filter_time_lengths:
        return _filter_time_lengths[approximant](**kwargs)
    else:
        return None

__all__ = ["get_td_waveform", "get_fd_waveform",
           "print_td_approximants", "print_fd_approximants",
           "td_approximants", "fd_approximants",
           "get_waveform_filter", "filter_approximants",
           "get_waveform_filter_norm", "get_waveform_end_frequency",
           "waveform_norm_exists", "get_template_amplitude_norm",
           "get_waveform_filter_length_in_time", "get_sgburst_waveform",
           "print_sgburst_approximants", "sgburst_approximants"]
