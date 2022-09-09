import numpy as np
import pickle

import pycbc.waveform

signals = []

def s_ldc2pycbc(mag, pol):
    return mag*np.cos(pol)

with open('../param_files/MBHB_params_v2.pkl', 'rb') as f:
    pmbhb = pickle.load(f)

    pMBHB = pmbhb[0]

    modes = [(2,2)]

    params = {'approximant': 'connor_bbhx',
        'mass1': pMBHB['Mass1'],
        'mass2': pMBHB['Mass2'],
        'delta_f':1/31536000,
        'inclination': np.pi/3,
        'tc': pMBHB['CoalescenceTime'],
        'polarization': np.pi/2,
        'spin1z': s_ldc2pycbc(pMBHB['Spin1'], pMBHB['PolarAngleOfSpin1']),
        'spin2z': s_ldc2pycbc(pMBHB['Spin2'], pMBHB['PolarAngleOfSpin2']),
        'coa_phase' : pMBHB['PhaseAtCoalescence'],
        'distance': pMBHB['Distance'],
        'eclipticlatitude': pMBHB['EclipticLatitude'],
        'eclipticlongitude': pMBHB['EclipticLongitude'],
        'mode_array':modes}

    signals.append(params)

params = signals[0]



A = pycbc.waveform.get_fd_waveform(f_lower=0.0001, **params)

import matplotlib.pylab as plt


td_A = A[0].to_timeseries()

plt.figure(figsize=(12,6))
plt.plot(td_A.sample_times,td_A)
plt.xlim(params['tc']-1000, params['tc']+1000)
plt.savefig('test_waveform_plugin.jpg')
