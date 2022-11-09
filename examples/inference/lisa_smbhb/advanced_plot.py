import subprocess
import pickle
import numpy as np
from pycbc.conversions import q_from_mass1_mass2, mchirp_from_mass1_mass2


def spin_ldc2pycbc(mag, pol):
    return mag*np.cos(pol)

def plt(index):

    with open('./MBHB_params_v2_LISA_frame.pkl', 'rb') as f:
        params_true_all = pickle.load(f)

    p_index = index
    params_true = params_true_all[p_index]
    print(params_true)

    modes = [(2,2)]

    q = q_from_mass1_mass2(params_true['Mass1'], params_true['Mass2'])
    mchirp = mchirp_from_mass1_mass2(params_true['Mass1'],params_true['Mass2'])

    params = {'approximant': 'BBHX_PhenomD',
              'mass1': params_true['Mass1'],
              'mass2': params_true['Mass2'],
              'inclination': params_true['Inclination'],
              'tc_lisa': params_true['CoalescenceTime_LISA'],
              'polarization_lisa': params_true['Polarization_LISA'],
              'spin1z': spin_ldc2pycbc(params_true['Spin1'], params_true['PolarAngleOfSpin1']),
              'spin2z': spin_ldc2pycbc(params_true['Spin2'], params_true['PolarAngleOfSpin2']),
              'coa_phase': params_true['PhaseAtCoalescence'],
              'distance': params_true['Distance'],
              'eclipticlatitude_lisa': params_true['EclipticLatitude_LISA'],
              'eclipticlongitude_lisa': params_true['EclipticLongitude_LISA'],
              'mchirp': mchirp,
              'q': q,
              'mode_array': modes
            }

    plot_code = f"""
                pycbc_inference_plot_posterior \
                --input-file lisa_smbhb.hdf \
                --output-file lisa_smbhb_mass_tc_{p_index}.png \
                --z-arg snr --plot-scatter --plot-marginal \
                --plot-contours --contour-color black \
                --parameters \
                    mass1_from_mchirp_q(mchirp,q):mass1 \
                    mass2_from_mchirp_q(mchirp,q):mass2 \
                    tc \
                --expected-parameters \
                    mass1_from_mchirp_q(mchirp,q):{params['mass1']} \
                    mass2_from_mchirp_q(mchirp,q):{params['mass2']} \
                    tc:{params['tc_lisa']} \
                """
    return plot_code

# The index of first SMBHB in LDC Sangria (0-14) is 0.
p = [0]

for i in p:
    process = subprocess.Popen(plt(i).split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    print('rel{} image created'.format(i))
