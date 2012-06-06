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

def mass_to_massSI(mass):
    return mass * lal.LAL_MSUN_SI

def parsecs_to_meters(distance):
    return distance *lal.LAL_PC_SI

WAVEFORM_DEFS = {'TaylorT4' : {'FUNC':lalsim.XLALSimInspiralTaylorT4PNGenerator,'ARGS':('phi0','v0','delta_t','m1SI','m2SI','f_lower','distance','inclination','amp_order','phase_order')} ,
                'TaylorT3' : {'FUNC':lalsim.XLALSimInspiralTaylorT3PNGenerator,'ARGS':('phi0','v0','delta_t','m1SI','m2SI','f_lower','distance','inclination','amp_order','phase_order')} ,
                'TaylorT2' : {'FUNC':lalsim.XLALSimInspiralTaylorT2PNGenerator,'ARGS':('phi0','v0','delta_t','m1SI','m2SI','f_lower','distance','inclination','amp_order','phase_order')} ,
                'TaylorT1' : {'FUNC':lalsim.XLALSimInspiralTaylorT1PNGenerator,'ARGS':('phi0','v0','delta_t','m1SI','m2SI','f_lower','distance','inclination','amp_order','phase_order')} ,
                'EOBNRv2'  : {'FUNC':lalsim.XLALSimIMREOBNRv2DominantMode,'ARGS':('phi0', 'delta_t', 'm1SI', 'm2SI', 'f_lower', 'distance', 'inclination')}
                }

CONVERSIONS = {'m1SI':{('mass1',):mass_to_massSI},
               'm2SI':{('mass2',):mass_to_massSI},
               'distance':  {('distance_in_parsecs',):parsecs_to_meters},
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


def get_waveform(template=None,**kwargs):
    """
    """
    # Get the parameters to generate the waveform
    # Note that keyword arguments override values in the template object
    input_params = {}
    if template is not None:
        input_params.update(template.__dict__)
    input_params.update(kwargs)

    if 'waveform' not in input_params:
        raise RuntimeError("A waveform must be chosen.")

    definition = WAVEFORM_DEFS[input_params['waveform']]

    waveform_args = ()
    for arg in definition['ARGS']:
       
        try:
            waveform_args += (input_params[arg],)
        except KeyError:
            try:
                for alt in CONVERSIONS[arg]:
                    found_conversion = False
                    try:
                        conv_args = ()
                        for field in alt:                     
                            conv_args += (input_params[field],)
                        waveform_args += (CONVERSIONS[arg][alt].__call__(*conv_args),)
                        found_conversion = True
                    except KeyError:                  
                        continue
                if not found_conversion:
                    raise KeyError(arg)
            except:
                err_str = "Please provide the argument " + str(arg)
                if arg in CONVERSIONS:
                    for alt in CONVERSIONS[arg]:
                        err_str += " or " +  str(alt)
                raise KeyError( err_str )
                

    lal_hp,lal_hc = definition['FUNC'].__call__(*waveform_args)
    hp = TimeSeries(lal_hp.data.data,delta_t=lal_hp.deltaT)
    hc = TimeSeries(lal_hc.data.data,delta_t=lal_hc.deltaT)

    return hp,hc






































