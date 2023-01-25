import astropy.units as u

# # Define the parameters of waveform
#
# # masses are given in solar mass units
# m1 = 20.*u.solMass
# m2 = 20.*u.solMass
#
# # the spins are dimensionless
# s1x = 0.*u.dimensionless_unscaled
# s1y = 0.*u.dimensionless_unscaled
# s1z = 0.*u.dimensionless_unscaled
# s2x = 0.*u.dimensionless_unscaled
# s2y = 0.*u.dimensionless_unscaled
# s2z = 0.*u.dimensionless_unscaled
#
# # we specify the distance in Mpc
# distance = 1000.*u.Mpc
#
# # other parameters are given in SI
# deltaT = 1./1024.*u.s
# deltaF = 1./32.*u.Hz
# f_min = 20.*u.Hz
# f_ref = 20.*u.Hz
# inclination = 0.*u.rad
#
# phiRef = 0.*u.rad
# eccentricity = 0.*u.dimensionless_unscaled
# longAscNodes = 0.*u.rad
# meanPerAno = 0.*u.rad
#
# python_dict = {'mass1' : m1,
#               'mass2' : m2,
#               'spin1x' : s1x,
#               'spin1y' : s1y,
#               'spin1z' : s1z,
#               'spin2x' : s2x,
#               'spin2y' : s2y,
#               'spin2z' : s2z,
#               'deltaT' : deltaT,
#               'deltaF' : deltaF,
#               'f22_start' : f_min,
#               'f22_ref': f_ref,
#               'phi_ref' : phiRef,
#               'distance' : distance,
#               'inclination' : inclination,
#               'eccentricity' : eccentricity,
#               'longAscNodes' : longAscNodes,
#               'meanPerAno' : meanPerAno}


pycbc_to_gws = {
        'delta_t': 'deltaT',
        'delta_f': 'deltaF',
        'f_lower': 'f22_start',
        'f_ref': 'f22_ref',
        'coa_phase': 'phi_ref',
        'long_asc_nodes': 'longAscNodes',
        'mean_per_ano': 'meanPerAno'
        }

gws_units = {'mass1': u.solMass,
             'mass2': u.solMass,
             'spin1x': u.dimensionless_unscaled,
             'spin1y': u.dimensionless_unscaled,
             'spin1z': u.dimensionless_unscaled,
             'spin2x': u.dimensionless_unscaled,
             'spin2y': u.dimensionless_unscaled,
             'spin2z': u.dimensionless_unscaled,
             'deltaT': u.s,
             'deltaF': u.Hz,
             'f22_stat': u.Hz,
             'f22_ref': u.Hz,
             'phi_ref': u.rad,
             'distance': u.Mpc,
             'inclination': u.rad,
             'eccentricity': u.dimensionless_unscaled,
             'longAscNodes': u.rad,
             'meanPerAno': u.rad}

def to_gwsignal_dict(p):
    params = {}
    for k in p.keys() & gws_units.keys():
        if k in pycbc_to_gws:
            knew = pycbc_to_gws.get(k)
        else:
            knew = k


    print('debugging')
    print('input params', p)
    print('modified params', params)
    return params
