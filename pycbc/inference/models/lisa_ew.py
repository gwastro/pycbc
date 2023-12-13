# Copyright (C) 2018  Collin Capano
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

"""
This modules provides models that have analytic solutions for the
log likelihood.
"""

import copy
import logging
import os
import numpy
import numpy.random

from .base import BaseModel

import pycbc.psd

from pycbc.waveform.early_warning_wform import PSDFirKernel, generate_early_warning_psds, generate_data_lisa_ew, generate_waveform_lisa_ew
from pycbc.waveform.waveform import parse_mode_array, props
from .tools import marginalize_likelihood

waveform_params1 = {'approximant': 'BBHX_PhenomD',
    'mass1': 1E6,
    'mass2': 1E6,
    't_obs_start':2592000, # This is setting the data length.
    'f_lower':-1.,
    'inclination': 0.,
    'tc':2580000 , # This is the coalescence time.
    'polarization': 0.,
    'spin1z': 0.,
    'spin2z': 0.,
    'coa_phase' : numpy.pi/2.,
    'distance': 15000,
    'eclipticlatitude': 0.628,
    'eclipticlongitude': 0.785,
    'low-frequency-cutoff': 0.000001, # Why this *and* f_lower??? Don't think this one is used at all.
    'f_final': 0.1,
    'mode_array':[(2,2)],
    'delta_f': 1./2592000,
    'run_phenomd':False}

waveform_params2 = {'approximant': 'BBHX_PhenomD',
    'mass1': 2E6,
    'mass2': 5E5,
    't_obs_start':2592000, # This is setting the data length.
    'f_lower':-1.,
    'inclination': 0.,
    'tc':2560000 , # This is the coalescence time.
    'polarization': 1.3,
    'spin1z': 0.,
    'spin2z': 0.,
    'coa_phase' : 2.67,
    'distance': 36000,
    'eclipticlatitude': 1.322,
    'eclipticlongitude': 2.318,
    'low-frequency-cutoff': 0.000001, # Why this *and* f_lower??? Don't think this one is used at all.
    'f_final': 0.1,
    'mode_array':[(2,2)],
    'run_phenomd':False}

waveform_params3 = {'approximant': 'BBHX_PhenomD',
    'mass1': 1E6,
    'mass2': 7E5,
    't_obs_start':2592000, # This is setting the data length.
    'f_lower':-1.,
    'inclination': 0.2,
    'tc':2570000 , # This is the coalescence time.
    'polarization': 0.5,
    'spin1z': 0.4,
    'spin2z': -0.3,
    'coa_phase' : 1.93,
    'distance': 24000,
    'eclipticlatitude': 0.162,
    'eclipticlongitude': 1.272,
    'low-frequency-cutoff': 0.000001, # Why this *and* f_lower??? Don't think this one is used at all.
    'f_final': 0.1,
    'mode_array':[(2,2)],
    'run_phenomd':False}

waveform_params4 = {'approximant': 'BBHX_PhenomD',
    'mass1': 2.5E6,
    'mass2': 2.5E6,
    't_obs_start':2592000, # This is setting the data length.
    'f_lower':-1.,
    'inclination': 0.,
    'tc':2550000 , # This is the coalescence time.
    'polarization': 2.3,
    'spin1z': -0.9,
    'spin2z': 0.9,
    'coa_phase' : 0.01,
    'distance': 62000,
    'eclipticlatitude': 0.312,
    'eclipticlongitude': numpy.pi/2.,
    'low-frequency-cutoff': 0.000001, # Why this *and* f_lower??? Don't think this one is used at all.
    'f_final': 0.1,
    'mode_array':[(2,2)],
    'run_phenomd':False}

waveform_params5 = {'approximant': 'BBHX_PhenomD',
    'mass1': 6E5,
    'mass2': 5.5E5,
    't_obs_start':2592000, # This is setting the data length.
    'f_lower':-1.,
    'inclination': 0.245,
    'tc':2540000 , # This is the coalescence time.
    'polarization': 0.34,
    'spin1z': 0.1,
    'spin2z': 0.5,
    'coa_phase' : 0.11,
    'distance': 80000,
    'eclipticlatitude': 0.790,
    'eclipticlongitude': 1.865,
    'low-frequency-cutoff': 0.000001, # Why this *and* f_lower??? Don't think this one is used at all.
    'f_final': 0.1,
    'mode_array':[(2,2)],
    'run_phenomd':False}

# As per https://stackoverflow.com/questions/715417 this is
# distutils.util.strtobool, but it's being removed in python3.12 so recommend
# to copy source code as we've done here.
def strtobool(val):
    """Convert a string representation of truth to true (1) or false (0).
    True values are 'y', 'yes', 't', 'true', 'on', and '1'; false values
    are 'n', 'no', 'f', 'false', 'off', and '0'.  Raises ValueError if
    'val' is anything else.
    """
    val = val.lower()
    if val in ('y', 'yes', 't', 'true', 'on', '1'):
        return 1
    elif val in ('n', 'no', 'f', 'false', 'off', '0'):
        return 0
    else:
        raise ValueError("invalid truth value %r" % (val,))



class LISAEarlyWarningModel(BaseModel):
    r"""Ian is messing around

    Parameters
    ----------
    variable_params : (tuple of) string(s)
        A tuple of parameter names that will be varied.
    **kwargs :
        All other keyword arguments are passed to ``BaseModel``.


    """
    name = "lisa_ew"

    def __init__(
            self,
            variable_params,
            static_params=None,
            phase_marginalization=True,
            psd_path=None,
            **kwargs
        ):
        # Pop relevant values from kwargs
        cutoff_time = int(kwargs.pop('cutoff_time'))
        seed = int(kwargs.pop('seed'))
        kernel_length = int(kwargs.pop('kernel_length'))
        window_length = int(kwargs.pop('window_length'))
        tlen = kwargs.pop('tlen')
        inj_keys = [item for item in kwargs.keys() if item.startswith('injparam')]
        inj_params = {}
        for key in inj_keys:
            value = kwargs.pop(key)
            # Type conversion needed ... Ugly!!
            if key in ['injparam_run_phenomd']:
                value = strtobool(value)
            elif key in ['injparam_approximant']:
                pass  # Convert to string, so do nothing
            elif key in ['injparam_t_obs_start']:
                value = int(value)
            elif key in ['injparam_mode_array']:
                # If value is
                if "(" in value:
                    raise NotImplementedError
                elif "[" in value:
                    raise NotImplementedError
                else:
                    # e.g '22 33 44'
                    pass
            else:
                value = float(value)
            inj_params[key.replace('injparam_', '')] = value

        inj_params = parse_mode_array(inj_params)
        
        if isinstance(phase_marginalization, str):
            phase_marginalization = strtobool(phase_marginalization)
        self.phase_marginalization = phase_marginalization

        # set up base likelihood parameters
        super().__init__(variable_params, **kwargs)
        if "mode_array" not in static_params:
            static_params["mode_array"] = [(2,2)]

        run_phenomd = static_params.get("run_phenomd", "False")
        if isinstance(run_phenomd, str):
            run_phenomd = strtobool(run_phenomd)
        static_params["run_phenomd"] = run_phenomd
        self.static_params = parse_mode_array(static_params)

        tlen = 3140000
        sample_rate = 0.2
        length = int(tlen * sample_rate)
        flen = length // 2 + 1

        if psd_path is None:
            psd_path = os.getcwd()

        LISA_A_PSD = pycbc.psd.from_txt(os.path.join(psd_path, 'A_psd_4_smooth.txt'), flen, 1./tlen, 1./tlen, is_asd_file=False)
        LISA_E_PSD = pycbc.psd.from_txt(os.path.join(psd_path, 'E_psd_4_smooth.txt'), flen, 1./tlen, 1./tlen, is_asd_file=False)
        self.psds_for_datagen = {}
        self.psds_for_datagen['LISA_A'] = LISA_A_PSD
        self.psds_for_datagen['LISA_E'] = LISA_E_PSD
        psds_outs = generate_early_warning_psds()
        self.whitening_psds = {}
        self.whitening_psds['LISA_A'] = psds_outs[0][0]
        self.whitening_psds['LISA_E'] = psds_outs[1][0]
        #self.whitening_psds['LISA_A'].save('LISA_A_psd_ew.txt')
        #self.whitening_psds['LISA_E'].save('LISA_E_psd_ew.txt')

        # Store data for doing likelihoods.
        curr_params = inj_params
        self.kernel_length = kernel_length
        self.window_length = window_length
        self.cutoff_time = cutoff_time

        # Want to remove this!
        cutoff_time = self.cutoff_time + (curr_params['t_obs_start'] - curr_params['tc'])
        self.lisa_a_strain, self.lisa_e_strain = generate_data_lisa_ew(curr_params, self.psds_for_datagen, self.whitening_psds, seed, self.window_length, cutoff_time)
        self.lisa_a_strain_fd = pycbc.strain.strain.execute_cached_fft(
            self.lisa_a_strain,
            copy_output=True,
            uid=3223965
        )
        self.lisa_e_strain_fd = pycbc.strain.strain.execute_cached_fft(
            self.lisa_e_strain,
            copy_output=True,
            uid=3223967
        )

        #self.lisa_a_strain.save('LISA_A_strain.txt')
        #self.lisa_e_strain.save('LISA_E_strain.txt')


    def _loglikelihood(self):
        """Returns the log pdf of the multivariate normal.
        """
        cparams = copy.deepcopy(self.static_params)
        cparams.update(self.current_params)
        cutoff_time = self.cutoff_time + (cparams['t_obs_start'] - cparams['tc'])

        ws = generate_waveform_lisa_ew(cparams, self.whitening_psds, self.window_length, cutoff_time, self.kernel_length)
        wform_lisa_a = ws['LISA_A']
        wform_lisa_e = ws['LISA_E']
        #wform_td_A = pycbc.strain.strain.execute_cached_ifft(
        #    wform_lisa_a,
        #    copy_output=False,
        #    uid=12239
        #)
        #wform_td_E = pycbc.strain.strain.execute_cached_ifft(
        #    wform_lisa_e,
        #    copy_output=False,
        #    uid=122391
        #)

        #wform_td_A.save('LISA_A_waveform.txt')
        #wform_td_E.save('LISA_E_waveform.txt')
        snr_A = pycbc.filter.overlap_cplx(wform_lisa_a, self.lisa_a_strain_fd,
                                     normalized=False)
        snr_E = pycbc.filter.overlap_cplx(wform_lisa_e, self.lisa_e_strain_fd,
                                     normalized=False)
        a_norm = pycbc.filter.sigmasq(wform_lisa_a)
        e_norm = pycbc.filter.sigmasq(wform_lisa_e)

        hs = snr_A + snr_E
        hh = (a_norm + e_norm)

        #print ((snr_A[0].real)**2 + (snr_E[0].real)**2)
        #print(snr_A + snr_E, (a_norm + e_norm)/2., (a_norm_test + e_norm_test)/2., cparams['distance'])
        #quit_here

        return marginalize_likelihood(hs, hh, phase=self.phase_marginalization)
