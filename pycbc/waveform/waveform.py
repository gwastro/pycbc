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
"""
"""
import swiglal as lal
import swiglalsimulation as lalsim
from pycbc.types import TimeSeries
from math import sin, cos

def mass_to_massSI(mass):
    return mass * lal.LAL_MSUN_SI

def parsecs_to_meters(distance):
    return distance *lal.LAL_PC_SI

def default_SpinTaylorT4_coordinates(inclination):
    LNhatx = sin(inclination)
    LNhaty = 0
    LNhatz = cos(inclination)
    E1x = cos(inclination)
    E1y = 0
    E1z = - sin(inclination)
    return (LNhatx,LNhaty,LNhatz,E1x,E1y,E1z)

WAVEFORM_DEFS = {
'TaylorEt' : 
    {'FUNC':lalsim.XLALSimInspiralTaylorEtPNGenerator,
     'ARGS':[('phi0',0),('v0',1),('delta_t',None),('m1SI',None),('m2SI',None),('f_lower',None),('distance',None),('inclination',0),('amp_order',0),('phase_order',7)]} ,
'TaylorT4' : 
    {'FUNC':lalsim.XLALSimInspiralTaylorT4PNGenerator,
     'ARGS':[('phi0',0),('v0',1),('delta_t',None),('m1SI',None),('m2SI',None),('f_lower',None),('distance',None),('inclination',0),('amp_order',0),('phase_order',7)]} ,
'TaylorT3' : 
    {'FUNC':lalsim.XLALSimInspiralTaylorT3PNGenerator,
     'ARGS':[('phi0',0),('v0',1),('delta_t',None),('m1SI',None),('m2SI',None),('f_lower',None),('distance',None),('inclination',0),('amp_order',0),('phase_order',7)]} ,
'TaylorT2' : 
    {'FUNC':lalsim.XLALSimInspiralTaylorT2PNGenerator,
     'ARGS':[('phi0',0),('v0',1),('delta_t',None),('m1SI',None),('m2SI',None),('f_lower',None),('distance',None),('inclination',0),('amp_order',0),('phase_order',7)]} ,
'TaylorT1' : 
    {'FUNC':lalsim.XLALSimInspiralTaylorT1PNGenerator,
     'ARGS':[('phi0',0),('v0',1),('delta_t',None),('m1SI',None),('m2SI',None),('f_lower',None),('distance',None),('inclination',0),('amp_order',0),('phase_order',7)]} ,
'SpinTaylorT4':
      {'FUNC':lalsim.XLALSimInspiralSpinTaylorT4,
     'ARGS':[('phi0',0),('v0',1),('delta_t',None),('m1SI',None),('m2SI',None),('f_lower',None),('distance',None),('spin1x',0),
             ('spin1y',0),('spin1z',0),('spin2x',0),('spin2y',0),('spin2z',0),('lambda1',0),('lambda2',0),
             ('SpinTaylorT4_coordinates',(0,0,1,1,0,0)),
             ('interactionFlags',lalsim.LAL_SIM_INSPIRAL_INTERACTION_SPIN_ORBIT_15PN | lalsim.LAL_SIM_INSPIRAL_INTERACTION_SPIN_SPIN_2PN),
             ('phase_order',7),('amp_order',0)]} ,
'EOBNRv2HM'  : 
    {'FUNC':lalsim.XLALSimIMREOBNRv2DominantMode,
     'ARGS':['phi0', 'delta_t', 'm1SI', 'm2SI', 'f_lower', 'distance', 'inclination'] ,
     'ARGS':[('phi0',0),('delta_t',None),('m1SI',None),('m2SI',None),('f_lower',None),('distance',None),('inclination',0)]},
'EOBNRv2'  : 
    {'FUNC':lalsim.XLALSimIMREOBNRv2DominantMode,
     'ARGS':['phi0', 'delta_t', 'm1SI', 'm2SI', 'f_lower', 'distance', 'inclination'] ,
     'ARGS':[('phi0',0),('delta_t',None),('m1SI',None),('m2SI',None),('f_lower',None),('distance',None),('inclination',0)]},
'SEOBNRv1' : 
    {'FUNC':lalsim.XLALSimIMRSpinAlignedEOBWaveform, 
     'ARGS':[('phi0',0),('delta_t',None),('m1SI',None),('m2SI',None),('f_lower',None),('distance',None),('inclination',0),('spin1z',0),('spin2z',0)]},
}

CONVERSIONS = {'m1SI':{('mass1',):mass_to_massSI},
               'm2SI':{('mass2',):mass_to_massSI},
               'distance':  {('distance_in_parsecs',):parsecs_to_meters},
               'SpinTaylorT4_coordinates': {('inclination',):default_SpinTaylorT4_coordinates},
              }

def list_available_approximants():
    for apx in WAVEFORM_DEFS:
        print(apx) 

def list_approximant_args(approx):
    for arg in WAVEFORM_DEFS[approx]['ARGS']:
        print(arg)
        if arg in CONVERSIONS:
            for alt in CONVERSIONS[arg]:
                print( str("     or " +  str(alt) ) )


def make_tuple(value):
    if isinstance(value,tuple):
        return value
    else:
        return value,

def get_waveform(template=None,**kwargs):
    """Return the waveform specified by the attributes of the template with 
       overrides given by keyword argument
    """
    # Get the parameters to generate the waveform
    # Note that keyword arguments override values in the template object
    input_params = {}
    if template is not None:
        input_params.update(template.__dict__)
    input_params.update(kwargs)

    if 'approximant' not in input_params:
        raise RuntimeError("The waveform approximant must be chosen.")

    definition = WAVEFORM_DEFS[input_params['approximant']]

    waveform_args = ()

    for arg,default_value in definition['ARGS']:
        try:
            waveform_args += make_tuple(input_params[arg])
        except KeyError:
            try:
                value = find_conversion(input_params,arg,default_value)
                waveform_args += make_tuple(value)
            except KeyError:
                err_str = "Please provide a value for " + str(arg)
                if arg in CONVERSIONS:
                    for alt in CONVERSIONS[arg]:
                        err_str += " or " +  str(alt)
                raise KeyError( err_str )

    print waveform_args
    lal_hp,lal_hc = definition['FUNC'].__call__(*waveform_args)
    hp = TimeSeries(lal_hp.data.data,delta_t=lal_hp.deltaT)
    hc = TimeSeries(lal_hc.data.data,delta_t=lal_hc.deltaT)

    return (hp,hc)



def find_conversion(input_params, arg,default):
    value = None
    try:
        for alt in CONVERSIONS[arg]:
            try:
                conv_args = ()
                for field in alt:  
                    try:                             
                        conv_arg = (input_params[field],)
                    except KeyError:                                  
                        conv_arg = find_conversion(input_params,field,None)                            
                    conv_args += make_tuple(conv_arg)            
                value = CONVERSIONS[arg][alt].__call__(*conv_args)
            except KeyError:
                continue
            return value   
    except KeyError:
        pass

    if  default is None:
        raise KeyError(arg)
    else:
        return default               




































