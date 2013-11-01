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
import sys
import lal
import lalsimulation
from pycbc.types import TimeSeries,FrequencySeries,zeros,Array,complex_same_precision_as
from pycbc.types import complex64, float32, complex128
from pycbc import HAVE_CUDA,HAVE_OPENCL
from pycbc.scheme import mgr
from pycbc.types import real_same_precision_as
import pycbc.scheme as _scheme
import inspect
from pycbc.fft import fft
import pycbc


def solar_mass_to_kg(solar_masses):
    return solar_masses * lal.LAL_MSUN_SI
    
def megaparsecs_to_meters(distance):
    return distance * lal.LAL_PC_SI * 1e6

default_args = {'spin1x':0,'spin1y':0,'spin1z':0,
                'spin2x':0,'spin2y':0,'spin2z':0,'lambda1':0, 'lambda2':0,
                'inclination':0,'distance':1,'f_final':0, 'f_ref':0, 'coa_phase':0,
                'amplitude_order':-1,'phase_order':-1,'spin_order':-1,
                'tidal_order':-1}

base_required_args = ['mass1','mass2','f_lower']
td_required_args = base_required_args + ['delta_t']
fd_required_args = base_required_args + ['delta_f']


# td, fd, filter waveforms generated on the CPU
_lalsim_td_approximants = {}
_lalsim_fd_approximants = {}
_lalsim_enum = {}

def _imrphenombfreq(**p):
    import lal, lalinspiral, lalsimulation
    from pycbc import pnutils
    params = lalinspiral.InspiralTemplate()
    m1 = p['mass1']
    m2 = p['mass2']
    
    mc, et = pnutils.mass1_mass2_to_mchirp_eta(m1, m2)
    params.approximant = lalsimulation.IMRPhenomB
    params.fLower = p['f_lower']
    params.eta = et
    params.distance = p['distance'] * lal.LAL_PC_SI * 1e6
    params.mass1 = m1
    params.mass2 = m2
    params.spin1[2] = p['spin1z']
    params.spin2[2] = p['spin2z']
    params.startPhase = p['coa_phase']*2 - 3.141592653
    params.startTime = 0
    
    params.tSampling = 8192
    N = int(params.tSampling / p['delta_f'])
    n = N / 2
    
    # Create temporary memory to hold the results and call the generator
    hpt = zeros(N, dtype=float32)
    hct = zeros(N, dtype=float32)    
    hpt=hpt.lal()
    hct=hct.lal()    
    lalinspiral.BBHPhenWaveBFreqDomTemplates(hpt, hct, params)
    
    # Copy the results to a complex frequencyseries format 
    hctc = FrequencySeries(zeros(n, dtype=complex64), delta_f=p['delta_f'])
    hptc = FrequencySeries(zeros(n, dtype=complex64), delta_f=p['delta_f'])
       
    hptc.data += hpt.data[0:n]
    hptc.data[1:n] += hpt.data[N:N-n:-1] * 1j
    
    hctc.data += hct.data[0:n]
    hctc.data[1:n] += hct.data[N:N-n:-1] * 1j
    
    return hptc.astype(complex128),  hctc.astype(complex128)

def _get_waveform_from_inspiral(**p):
    def get_string_from_order(order):
        if order == 0:
            return 'newtonian'
        if order == 1:
            return 'oneHalfPN'
        if order == 2:
            return 'onePN'
        if order == 3:
            return 'onePointFivePN'
        if order == 4:
            return 'twoPN'
        if order == 5:
            return 'twoPointFivePN'
        if order == 6:
            return 'threePN'
        if (order == 7) or (order == -1):
            return 'threePointFivePN'
        if order == -8:
            return 'pseudoFourPN'

    import lalmetaio
    import lal
    import lalinspiral
    # prefix with 'Inspiral-'
    name = p['approximant'][9:]
    
    if name.startswith('EOB'):
        p['phase_order'] = -8

    params = lalmetaio.SimInspiralTable()
    params.waveform = name + get_string_from_order(p['phase_order'])
    params.mass1= p['mass1']
    params.mass2= p['mass2']
    params.f_lower = p['f_lower']
    params.spin1x = p['spin1x']
    params.spin1y = p['spin1y']
    params.spin1z = p['spin1z']
    params.spin2x = p['spin2x']
    params.spin2y = p['spin2y']
    params.spin2z = p['spin2z']
    params.inclination = p['inclination']
    params.distance = p['distance']
    params.coa_phase = p['coa_phase']
    
    guess_length =  lalinspiral.FindChirpChirpTime(params.mass1, params.mass2, 
                                                        params.f_lower, 7)
    guess_length = max(guess_length, 3)
   
    params.geocent_end_time = guess_length * 1.5
    params.taper = 'TAPER_NONE'
    bufferl = guess_length * 2
    dt = p['delta_t']
    df = 1.0 / bufferl
    sample_rate = int(1.0 / dt)
    epoch = lal.LIGOTimeGPS(0, 0)
    N = bufferl * sample_rate
    n = N / 2 + 1
    
    resp = pycbc.types.FrequencySeries(zeros(n), delta_f = df, 
                                       epoch=epoch, dtype=complex64) + 1
    out = pycbc.types.TimeSeries(zeros(N), delta_t = dt, 
                                 epoch=epoch, dtype=float32)
    outl = out.lal()
    outl.sampleUnits = lal.lalADCCountUnit
    
    out2 = pycbc.types.TimeSeries(zeros(N), delta_t = dt, 
                                 epoch=epoch, dtype=float32)
    outl2 = out.lal()
    outl2.sampleUnits = lal.lalADCCountUnit
    
    respl = resp.lal()
    respl.sampleUnites = lal.lalDimensionlessUnit

    lalinspiral.FindChirpInjectSignals(outl, params, respl)  
    
    params.coa_phase -= lal.LAL_PI / 4
    lalinspiral.FindChirpInjectSignals(outl2, params, respl)
    seriesp = TimeSeries(outl.data.data, delta_t=dt, 
                         epoch=epoch - params.geocent_end_time)
                         
    seriesc = TimeSeries(outl2.data.data, delta_t=dt,
                         epoch=epoch - params.geocent_end_time)
    
    return seriesp, seriesc

def _lalsim_td_waveform(**p):
    flags = lalsimulation.SimInspiralCreateWaveformFlags()
    lalsimulation.SimInspiralSetSpinOrder(flags, p['spin_order'])
    lalsimulation.SimInspiralSetTidalOrder(flags, p['tidal_order'])

    hp, hc = lalsimulation.SimInspiralChooseTDWaveform(float(p['coa_phase']),
               float(p['delta_t']),
               float(solar_mass_to_kg(p['mass1'])),
               float(solar_mass_to_kg(p['mass2'])),
               float(p['spin1x']), float(p['spin1y']), float(p['spin1z']),
               float(p['spin2x']), float(p['spin2y']), float(p['spin2z']),
               float(p['f_lower']), float(p['f_ref']),
               megaparsecs_to_meters(float(p['distance'])),
               float(p['inclination']),
               float(p['lambda1']),  float(p['lambda2']), flags, None,
               int(p['amplitude_order']), int(p['phase_order']),
               _lalsim_enum[p['approximant']])

    hp = TimeSeries(hp.data.data[:], delta_t=hp.deltaT, epoch=hp.epoch)
    hc = TimeSeries(hc.data.data[:], delta_t=hc.deltaT, epoch=hc.epoch)

    return hp, hc

def _lalsim_fd_waveform(**p):
    flags = lalsimulation.SimInspiralCreateWaveformFlags()
    lalsimulation.SimInspiralSetSpinOrder(flags, p['spin_order'])
    lalsimulation.SimInspiralSetTidalOrder(flags, p['tidal_order'])

    hp, hc = lalsimulation.SimInspiralChooseFDWaveform(float(p['coa_phase']),
               float(p['delta_f']),
               float(solar_mass_to_kg(p['mass1'])),
               float(solar_mass_to_kg(p['mass2'])),
               float(p['spin1x']), float(p['spin1y']), float(p['spin1z']),
               float(p['spin2x']), float(p['spin2y']), float(p['spin2z']),
               float(p['f_lower']), float(p['f_final']), float(p['f_ref']),
               megaparsecs_to_meters(float(p['distance'])),
               float(p['inclination']),
               float(p['lambda1']), float(p['lambda2']), flags, None,
               int(p['amplitude_order']), int(p['phase_order']),
               _lalsim_enum[p['approximant']])

    hp = FrequencySeries(hp.data.data[:], delta_f=hp.deltaF,
                            epoch=hp.epoch)
    hc = FrequencySeries(hc.data.data[:], delta_f=hc.deltaF,
                            epoch=hc.epoch)                        
    
    return hp, hc

for approx_enum in xrange(0,lalsimulation.NumApproximants):
    if lalsimulation.SimInspiralImplementedTDApproximants(approx_enum):
        approx_name =  lalsimulation.GetStringFromApproximant(approx_enum)
        _lalsim_enum[approx_name] = approx_enum
        _lalsim_td_approximants[approx_name] = _lalsim_td_waveform

for approx_enum in xrange(0,lalsimulation.NumApproximants):
    if lalsimulation.SimInspiralImplementedFDApproximants(approx_enum):
        approx_name =  lalsimulation.GetStringFromApproximant(approx_enum)
        _lalsim_enum[approx_name] = approx_enum
        _lalsim_fd_approximants[approx_name] = _lalsim_fd_waveform

#Add lalinspiral approximants
insp_td = {}
for apx in ['EOB']:
    name = 'Inspiral-' + apx
    insp_td[name] = _get_waveform_from_inspiral


cpu_td = dict(_lalsim_td_approximants.items() + insp_td.items())
cpu_fd = _lalsim_fd_approximants
cpu_fd['Inspiral-IMRPhenomB'] = _imrphenombfreq

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

# Waveforms written in OpenCL
_opencl_td_approximants = {}
_opencl_fd_approximants = {}

opencl_td = dict(_lalsim_td_approximants.items() + _opencl_td_approximants.items())
opencl_fd = dict(_lalsim_fd_approximants.items() + _opencl_fd_approximants.items())

td_wav = _scheme.ChooseBySchemeDict()
fd_wav = _scheme.ChooseBySchemeDict()
td_wav.update({_scheme.CPUScheme:cpu_td,_scheme.CUDAScheme:cuda_td,_scheme.OpenCLScheme:opencl_td})
fd_wav.update({_scheme.CPUScheme:cpu_fd,_scheme.CUDAScheme:cuda_fd,_scheme.OpenCLScheme:opencl_fd})

# List the various available approximants ####################################

def print_td_approximants():
    print("Lalsimulation Approximants")
    for approx in _lalsim_td_approximants.keys():
        print "  " + approx
    print("CUDA Approximants")
    for approx in _cuda_td_approximants.keys():
        print "  " + approx
    print("OpenCL Approximants")
    for approx in _opencl_td_approximants.keys():
        print "  " + approx
    
def print_fd_approximants():
    print("Lalsimulation Approximants")
    for approx in _lalsim_fd_approximants.keys():
        print "  " + approx
    print("CUDA Approximants")
    for approx in _cuda_fd_approximants.keys():
        print "  " + approx
    print("OpenCL Approximants")
    for approx in _opencl_fd_approximants.keys():
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

def filter_approximants(scheme=_scheme.mgr.state):
    """Return a list of fourier domain approximants including those
       written specifically as templates.
    """
    return filter_wav[type(scheme)].keys()

# Input parameter handling ###################################################

def props(obj, **kwargs):
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
        The mass of teh second component object in the binary in solar masses.
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
        The x component of the first component objects spin vector. 
    spin1y : {0, float}, optional
        The y component of the first component objects spin vector. 
    spin1z : {0, float}, optional
        The z component of the first component objects spin vector. 
    spin2x : {0, float}, optional
        The x component of the second component objects spin vector. 
    spin2y : {0, float}, optional
        The y component of the second component objects spin vector. 
    spin2z : {0, float}, optional
        The z component of the second component objects spin vector. 
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
    wav_gen = td_wav[type(mgr.state)] 

    if 'approximant' not in input_params or input_params['approximant'] is None:
        raise ValueError("Please provide an approximant name")
    elif input_params['approximant'].startswith('Inspiral-'):
        pass
    elif input_params['approximant'] not in wav_gen:
        raise ValueError("Approximant %s not available" % \
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
        An object that has attached properties. This can be used to subsitute
        for keyword arguments. A common example would be a row in an xml table. 
    approximant : string
        A string that indicates the chosen approximant. See `fd_approximants` 
        for available options. 
    mass1 : float
        The mass of the first component object in the binary in solar masses.
    mass2 : 
        The mass of teh second component object in the binary in solar masses.
    delta_f :
        The frequency step used to generate the waveform. 
    f_lower :
        The starting frequency of the waveform.
    f_final : {-1, float}, optional
        The ending frequency of the waveform. The default indicates that the
        choice is made by the respective approximant. 
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
        The x component of the first component objects spin vector. 
    spin1y : {0, float}, optional
        The y component of the first component objects spin vector. 
    spin1z : {0, float}, optional
        The z component of the first component objects spin vector. 
    spin2x : {0, float}, optional
        The x component of the second component objects spin vector. 
    spin2y : {0, float}, optional
        The y component of the second component objects spin vector. 
    spin2z : {0, float}, optional
        The z component of the second component objects spin vector. 
    lambda1: {0, float}, optional
        The tidal deformability parameter of object 1. 
    lambda2: {0, flaot}, optional
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
    wav_gen = fd_wav[type(mgr.state)] 

    if 'approximant' not in input_params:
        raise ValueError("Please provide an approximant name")
    elif input_params['approximant'] not in wav_gen:
        raise ValueError("Approximant %s not available" % \
                         (input_params['approximant']))

    for arg in fd_required_args:
        if arg in input_params:
            pass
        else:
            raise ValueError("Please provide " + str(arg) )

    return wav_gen[input_params['approximant']](**input_params)

    
# Waveform filter routines ###################################################

# Organize Filter Generators
_inspiral_fd_filters = {}
_cuda_fd_filters = {}
opencl_fd_filter = {}

from fctmplt import findchirp_template
from spa_tmplt import spa_tmplt
_inspiral_fd_filters['SPAtmplt'] = spa_tmplt
#_inspiral_fd_filters['FindChirpSP'] = findchirp_template

_cuda_fd_filters['SPAtmplt'] = spa_tmplt
opencl_fd_filter['SPAtmplt'] = spa_tmplt


filter_wav = _scheme.ChooseBySchemeDict()
filter_wav.update( {_scheme.CPUScheme:_inspiral_fd_filters, 
              _scheme.CUDAScheme:_cuda_fd_filters, 
              _scheme.OpenCLScheme:opencl_fd_filter} )

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
#_filter_time_lengths['FindChirpSP'] = spa_length_in_time

def get_waveform_filter(out, template=None, **kwargs):
    """Return a frequency domain waveform filter for the specified approximant
    """
    n = len(out)
    N = (n-1)*2
    
    input_params = props(template,**kwargs)

    if input_params['approximant'] in filter_approximants(mgr.state):
        wav_gen = filter_wav[type(mgr.state)] 
        htilde = wav_gen[input_params['approximant']](out=out, **input_params)
        htilde.resize(n)
        htilde.length_in_time = get_waveform_filter_length_in_time(**input_params)
        return htilde
    if input_params['approximant'] in fd_approximants(mgr.state):
        wav_gen = fd_wav[type(mgr.state)] 
        hp, hc = wav_gen[input_params['approximant']](**input_params)
        hp.resize(n)
        hp.length_in_time = get_waveform_filter_length_in_time(**input_params)
        return hp
    elif input_params['approximant'] in td_approximants(mgr.state):
        delta_f = 1.0 / (N * input_params['delta_t'])
        wav_gen = td_wav[type(mgr.state)] 
        hp, hc = wav_gen[input_params['approximant']](**input_params)
        tmplt_length = len(hp)
        hp.resize(N)
        k_zero = int(hp.start_time / hp.delta_t)
        hp.roll(k_zero)
        htilde = FrequencySeries(out, delta_f=delta_f, copy=False)
        fft(hp.astype(real_same_precision_as(htilde)), htilde)
        htilde.length_in_time = tmplt_length
        return htilde
    else:
        raise ValueError("Approximant %s not available" % \
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

__all__ = ["get_td_waveform","get_fd_waveform","print_td_approximants",
           "print_fd_approximants","td_approximants","fd_approximants", 
           "get_waveform_filter", 
           "filter_approximants", "get_waveform_filter_norm", "get_waveform_end_frequency",
            "waveform_norm_exists", "get_template_amplitude_norm",
           "get_waveform_filter_length_in_time"]
