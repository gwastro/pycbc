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


#FIXME######################################### RELOCATE THIS TO ANOTHER MODULE

def solar_mass_to_kg(solar_masses):
    return solar_masses * lal.LAL_MSUN_SI
    
def parsecs_to_meters(distance):
    return distance *lal.LAL_PC_SI



default_args = {'spin1x':0,'spin1y':0,'spin1z':0,
                'spin2x':0,'spin2y':0,'spin2z':0,'lambda1':0, 'lambda2':0,
                'inclination':0,'distance':1e6,'f_final':0,'phi0':0,
                'amplitude_order':-1,'phase_order':-1}

base_required_args = ['mass1','mass2','f_lower']
td_required_args = base_required_args + ['delta_t']
fd_required_args = base_required_args + ['delta_f']


# td, fd, filter waveforms generated on the CPU
_lalsim_td_approximants = {}
_lalsim_enum = {}

def _lalsim_td_waveform(**p):
    hp,hc = lalsimulation.SimInspiralChooseTDWaveform(float(p['phi0']),
               float(p['delta_t']),
               float(solar_mass_to_kg(p['mass1'])),
               float(solar_mass_to_kg(p['mass2'])),
               float(p['spin1x']), float(p['spin1y']), float(p['spin1z']),
               float(p['spin2x']), float(p['spin2y']), float(p['spin2z']),
               float(p['f_lower']), 0,
               parsecs_to_meters(float(p['distance'])),
               float(p['inclination']),
               float(p['lambda1']),  float(p['lambda2']), None, None,
               int(p['amplitude_order']), int(p['phase_order']),
               _lalsim_enum[p['approximant']])

    hp = TimeSeries(hp.data.data,delta_t=hp.deltaT,epoch=hp.epoch)
    hc = TimeSeries(hc.data.data,delta_t=hc.deltaT,epoch=hc.epoch)

    return hp,hc

def _lalsim_fd_waveform(**p):
    htilde = lalsimulation.SimInspiralChooseFDWaveform(float(p['phi0']),
               float(p['delta_f']),
               float(solar_mass_to_kg(p['mass1'])),
               float(solar_mass_to_kg(p['mass2'])),
               float(p['spin1x']),float(p['spin1y']),float(p['spin1z']),
               float(p['spin2x']),float(p['spin2y']),float(p['spin2z']),
               float(p['f_lower']),float(p['f_final']),
               parsecs_to_meters(float(p['distance'])),
               float(p['inclination']),
               float(p['lambda1']),float(p['lambda2']),None,None,
               int(p['amplitude_order']),int(p['phase_order']),
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
    """Return the waveform specified by the attributes of the template with 
       overrides given by keyword argument
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

    hp,hc = wav_gen[input_params['approximant']](**input_params)
    return (hp,hc)

def get_fd_waveform(template=None, **kwargs):
    """Return the frequency domain waveform specified by the attributes
       of the template with overrides given by keyword argument
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
    raise NotImplementedError

    
def get_waveform_filter_precondition(approximant):
    raise NotImplementedError


__all__ = ["get_td_waveform","get_fd_waveform","print_td_approximants",
           "print_fd_approximants","td_approximants","fd_approximants"]
