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
from pycbc import HAVE_CUDA,HAVE_OPENCL
from pycbc.scheme import mgr
import pycbc.scheme as _scheme
import inspect
from pycbc.fft import fft
import pycbc


def solar_mass_to_kg(solar_masses):
    return solar_masses * lal.LAL_MSUN_SI
    
def kiloparsecs_to_meters(distance):
    return distance * lal.LAL_PC_SI * 1000



default_args = {'spin1x':0,'spin1y':0,'spin1z':0,
                'spin2x':0,'spin2y':0,'spin2z':0,'lambda1':0, 'lambda2':0,
                'inclination':0,'distance':1e3,'f_final':0,'phi0':0,
                'amplitude_order':-1,'phase_order':-1,'spin_order':-1,
                'tidal_order':-1}

base_required_args = ['mass1','mass2','f_lower']
td_required_args = base_required_args + ['delta_t']
fd_required_args = base_required_args + ['delta_f']


# td, fd, filter waveforms generated on the CPU
_lalsim_td_approximants = {}
_lalsim_enum = {}

def _lalsim_td_waveform(**p):
    flags = lalsimulation.SimInspiralCreateWaveformFlags()
    lalsimulation.SimInspiralSetSpinOrder(flags, p['spin_order'])
    lalsimulation.SimInspiralSetTidalOrder(flags, p['tidal_order'])

    hp,hc = lalsimulation.SimInspiralChooseTDWaveform(float(p['phi0']),
               float(p['delta_t']),
               float(solar_mass_to_kg(p['mass1'])),
               float(solar_mass_to_kg(p['mass2'])),
               float(p['spin1x']), float(p['spin1y']), float(p['spin1z']),
               float(p['spin2x']), float(p['spin2y']), float(p['spin2z']),
               float(p['f_lower']), 0,
               kiloparsecs_to_meters(float(p['distance'])),
               float(p['inclination']),
               float(p['lambda1']),  float(p['lambda2']), flags, None,
               int(p['amplitude_order']), int(p['phase_order']),
               _lalsim_enum[p['approximant']])

    hp = TimeSeries(hp.data.data,delta_t=hp.deltaT,epoch=hp.epoch)
    hc = TimeSeries(hc.data.data,delta_t=hc.deltaT,epoch=hc.epoch)

    return hp,hc

def _lalsim_fd_waveform(**p):
    flags = lalsimulation.SimInspiralCreateWaveformFlags()
    lalsimulation.SimInspiralSetSpinOrder(flags, p['spin_order'])
    lalsimulation.SimInspiralSetTidalOrder(flags, p['tidal_order'])

    htilde = lalsimulation.SimInspiralChooseFDWaveform(float(p['phi0']),
               float(p['delta_f']),
               float(solar_mass_to_kg(p['mass1'])),
               float(solar_mass_to_kg(p['mass2'])),
               float(p['spin1x']), float(p['spin1y']), float(p['spin1z']),
               float(p['spin2x']), float(p['spin2y']), float(p['spin2z']),
               float(p['f_lower']), float(p['f_final']),
               kiloparsecs_to_meters(float(p['distance'])),
               float(p['inclination']),
               float(p['lambda1']), float(p['lambda2']), flags, None,
               int(p['amplitude_order']), int(p['phase_order']),
               _lalsim_enum[p['approximant']])

    htilde = FrequencySeries(htilde.data.data,delta_f=htilde.deltaF,
                            epoch=htilde.epoch)
    return htilde

for approx_enum in xrange(0,lalsimulation.NumApproximants):
    if lalsimulation.SimInspiralImplementedTDApproximants(approx_enum):
        approx_name =  lalsimulation.GetStringFromApproximant(approx_enum)
        _lalsim_enum[approx_name] = approx_enum
        _lalsim_td_approximants[approx_name] = _lalsim_td_waveform

_lalsim_fd_approximants = {}
for approx_enum in xrange(0,lalsimulation.NumApproximants):
    if lalsimulation.SimInspiralImplementedFDApproximants(approx_enum):
        approx_name =  lalsimulation.GetStringFromApproximant(approx_enum)
        _lalsim_enum[approx_name] = approx_enum
        _lalsim_fd_approximants[approx_name] = _lalsim_fd_waveform

def _filter_td_waveform(p):
    raise NotImplementedError

def _filter_fd_waveform(p):
    raise NotImplementedError


# Waveforms that are optimized to work as filters
_filter_fd_approximants = {}
_filter_td_approximants = {}

cpu_td = _lalsim_td_approximants
cpu_fd = _lalsim_fd_approximants
cpu_td_filter =  dict(_filter_td_approximants.items() + \
              _filter_td_approximants.items() + _lalsim_td_approximants.items())
cpu_fd_filter =  dict(_lalsim_fd_approximants.items() + \
              _filter_fd_approximants.items())

# Waveforms written in CUDA
_cuda_td_approximants = {}
_cuda_fd_approximants = {}

if pycbc.HAVE_CUDA:
    from pycbc.waveform.TaylorF2 import taylorf2 as cuda_taylorf2
    from pycbc.waveform.pycbc_spa_tmplt import spa_tmplt
    from pycbc.waveform.pycbc_phenomC_tmplt import imrphenomc_tmplt
    _cuda_fd_approximants["IMRPhenomC"] = imrphenomc_tmplt
    _cuda_fd_approximants["SPAtmplt"] = spa_tmplt
    _cuda_fd_approximants['TaylorF2'] = cuda_taylorf2

cuda_td = dict(_lalsim_td_approximants.items() + _cuda_td_approximants.items())
cuda_fd = dict(_lalsim_fd_approximants.items() + _cuda_fd_approximants.items())
cuda_td_filter = dict(cpu_td_filter.items() + cuda_td.items())
cuda_fd_filter = dict(cpu_fd_filter.items() + cuda_fd.items())

# Waveforms written in OpenCL
_opencl_td_approximants = {}
_opencl_fd_approximants = {}

opencl_td = dict(_lalsim_td_approximants.items() + _opencl_td_approximants.items())
opencl_fd = dict(_lalsim_fd_approximants.items() + _opencl_fd_approximants.items())
opencl_td_filter = dict(cpu_td_filter.items() + opencl_td.items())
opencl_fd_filter = dict(cpu_fd_filter.items() + opencl_fd.items())

td_wav = {type(None):cpu_td,_scheme.CUDAScheme:cuda_td,_scheme.OpenCLScheme:opencl_td}

fd_wav = {type(None):cpu_fd,_scheme.CUDAScheme:cuda_fd,_scheme.OpenCLScheme:opencl_fd}

td_filter = {type(None):cpu_td_filter,
            _scheme.CUDAScheme:cuda_td_filter,
            _scheme.OpenCLScheme:opencl_td_filter}
fd_filter = {type(None):cpu_fd_filter,
            _scheme.CUDAScheme:cuda_fd_filter,
            _scheme.OpenCLScheme:opencl_fd_filter}

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

def td_approximants(scheme=None):
    """Return a list containing the available time domain approximants for 
       the given processing scheme.
    """
    return td_wav[type(scheme)].keys()

def fd_approximants(scheme=None):
    """Return a list containing the available fourier domain approximants for 
       the given processing scheme.
    """
    return fd_wav[type(scheme)].keys()    

def list_filter_approximants():
    pass

def props(obj, **kwargs):
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
    input_params = default_args.copy()
    input_params.update(pr)
    input_params.update(kwargs)

    return input_params

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
        The mass of the first component object in the binary. 
    mass2 : 
        The mass of teh second component object in the binary.
    delta_t :
        The time step used to generate the waveform. 
    f_lower :
        The starting frequency of the waveform.
    distance : {1e3, float}, optional
        The distance from the observer to the source in kiloparsecs.
    inclination : {0, float}, optional
        The inclination angle of the source. 
    phi0 : {0, float}, optional
        The final phase or phase at the peak of the wavform. See documentation on
    specific approximants for exact usage. 
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
    hplus: TimeSeries
        The plus polarization of the waveform.
    hcross: TimeSeries
        The cross polarization of the waveform.    
    """
    input_params = props(template,**kwargs)
    wav_gen = td_wav[type(mgr.state)] 

    if 'approximant' not in input_params:
        raise ValueError("Please provide an approximant name")
    elif input_params['approximant'] not in wav_gen:
        raise ValueError("Approximant not available")

    for arg in td_required_args:
        if arg in input_params:
            pass
        else:
            raise ValueError("Please provide " + str(arg) )

    hp, hc = wav_gen[input_params['approximant']](**input_params)
    return (hp, hc)

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
        The mass of the first component object in the binary. 
    mass2 : 
        The mass of teh second component object in the binary.
    delta_f :
        The frequency step used to generate the waveform. 
    f_lower :
        The starting frequency of the waveform.
    f_final : {-1, float}, optional
        The ending frequency of the waveform. The default indicates that the
    choice is made by the respective approximant. 
    distance : {1e3, float}, optional
        The distance from the observer to the source in kiloparsecs.
    inclination : {0, float}, optional
        The inclination angle of the source. 
    phi0 : {0, float}, optional
        The final phase or phase at the peak of the wavform. See documentation on
    specific approximants for exact usage. 
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
    htilde: FrequencySeries
        The cosine phase of the waveform in frequency domain.
    """

    input_params = props(template,**kwargs)
    wav_gen = fd_wav[type(mgr.state)] 

    if 'approximant' not in input_params:
        raise ValueError("Please provide an approximant name")
    elif input_params['approximant'] not in wav_gen:
        raise ValueError("Approximant not available")

    for arg in fd_required_args:
        if arg in input_params:
            pass
        else:
            raise ValueError("Please provide " + str(arg) )

    htilde = wav_gen[input_params['approximant']](**input_params)
    return htilde

def get_waveform_filter(length, template=None, **kwargs):
    """Return a frequency domain waveform filter for the specified approximant
    """
    n = length
    N = (n-1)*2

    input_params = props(template,**kwargs)

    if input_params['approximant'] in fd_approximants(mgr.state):
        wav_gen = fd_wav[type(mgr.state)] 
        htilde = wav_gen[input_params['approximant']](**input_params)
        htilde.resize(n)
        return htilde
    elif input_params['approximant'] in td_approximants(mgr.state):
        delta_f = 1.0 / (N * input_params['delta_t'])
        wav_gen = td_wav[type(mgr.state)] 
        hp, hc = wav_gen[input_params['approximant']](**input_params)
        hp.resize(N)
        k_zero = int(hp.start_time / hp.delta_t)
        hp.roll(k_zero)
        htilde = FrequencySeries(zeros(n), delta_f=delta_f,
                                            dtype=complex_same_precision_as(hp))
        fft(hp, htilde)
        return htilde
    else:
        raise ValueError("Approximant Not Available")

    
def get_waveform_filter_precondition(approximant):
    raise NotImplementedError


__all__ = ["get_td_waveform","get_fd_waveform","print_td_approximants",
           "print_fd_approximants","td_approximants","fd_approximants", "get_waveform_filter"]
