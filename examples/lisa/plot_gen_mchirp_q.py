import pickle
import numpy as np

def s_ldc2pycbc(mag, pol):
    return mag*np.cos(pol)

def mr(m1,m2):
    return m1/m2

def mchi(m1, m2):
    return ((m1*m2)**(3/5))/(m1+m2)**(1/5)

def plt(index):

    with open('param_files/MBHB_params_v2.pkl', 'rb') as f:
        pmbhb = pickle.load(f)

    p_index = index
    pMBHB = pmbhb[p_index]

    print(pMBHB)

    from ldc.common import tools

    modes = [(2,2)]#, (2,1), (3,3), (3,2), (4,4), (4,3)]

    psi, incl = tools.aziPolAngleL2PsiIncl(pMBHB["EclipticLatitude"],
                                           pMBHB["EclipticLongitude"],
                                           pMBHB['InitialPolarAngleL'],
                                           pMBHB['InitialAzimuthalAngleL'])

    if psi < 0.:
        psi = psi + 2*np.pi

    q = mr(pMBHB['Mass1'], pMBHB['Mass2'])
    mchirp = mchi(pMBHB['Mass1'],pMBHB['Mass2'])

    params = {'approximant': 'BBHX_PhenomD',
         'mass1': pMBHB['Mass1'],
         'mass2': pMBHB['Mass2'],
         'inclination': incl,
         'tc': pMBHB['CoalescenceTime'],
         'polarization': psi,
         'spin1z': s_ldc2pycbc(pMBHB['Spin1'], pMBHB['PolarAngleOfSpin1']),
         'spin2z': s_ldc2pycbc(pMBHB['Spin2'], pMBHB['PolarAngleOfSpin2']),
         'coa_phase' : pMBHB['PhaseAtCoalescence'],
         'distance': pMBHB['Distance'],
         'eclipticlatitude': pMBHB['EclipticLatitude'],
         'eclipticlongitude': pMBHB['EclipticLongitude'],
         'mchirp':mchirp,
         'q':q,
         'mode_array':modes}

    plot_code = f"""
    pycbc_inference_plot_posterior \
        --input-file simple_test.hdf \
        --output-file simple_run_post.png \
        --z-arg snr --plot-scatter --plot-marginal
        --parameters \
            mass1_from_mchirp_q(mchirp,q):mass1\
            mass2_from_mchirp_q(mchirp,q):mass2\
            tc\
        --expected-parameters  mass1_from_mchirp_q(mchirp,q):{params['mass1']}\
        mass2_from_mchirp_q(mchirp,q):{params['mass2']}\
        tc:{params['tc']}
            """
    return plot_code

import subprocess

p = [0]

for i in p:#range(15):
    process = subprocess.Popen(plt(i).split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    print('rel{} image created'.format(i))


