# Copyright (C) 2023 Shichao Wu, Alex Nitz
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

import numpy
from pycbc.coordinates import space
import unittest
from utils import simple_exit


seed = 8202
numpy.random.seed(seed)

# ranges to draw random numbers for each parameter
RANGES = {
    "tc_geo" : (1126259462.43, 1426259462.43),
    "ra" : (0.0, 2 * numpy.pi),
    "dec" : (-numpy.pi / 2, numpy.pi / 2),
    "polarization_geo" : (0.0, 2 * numpy.pi),
    "tc_ssb" : (1126259462.43, 1426259462.43),
    "eclipticlongitude_ssb" : (0.0, 2 * numpy.pi),
    "eclipticlatitude_ssb" : (-numpy.pi / 2, numpy.pi / 2),
    "polarization_ssb" : (0.0, 2 * numpy.pi),
    "tc_lisa" : (1126259462.43, 1426259462.43),
    "eclipticlongitude_lisa" : (0.0, 2 * numpy.pi),
    "eclipticlatitude_lisa" : (-numpy.pi / 2, numpy.pi / 2),
    "polarization_lisa" : (0.0, 2 * numpy.pi),
}

def almost_equal(derived_val, check_val, precision=1e-2):
    """Checks whether the difference in the derived and check values are less
    than the given precision.
    """
    allpass = numpy.allclose(derived_val, check_val, atol=precision)
    if not allpass:
        absdiff = abs(derived_val - check_val)
        maxidx = absdiff.argmax()
        maxdiff = absdiff[maxidx]
    else:
        maxdiff = maxidx = None
    return allpass, maxdiff, maxidx

def curve_similarity(curve_1, curve_2):
    """Using the Euclidean distance to check
    the similarity between two curves.
    """
    return numpy.linalg.norm(curve_1 - curve_2)


class TestParams(unittest.TestCase):
    def setUp(self, *args):
        self.numtests = 1000
        self.precision = 1e-8

        # generate some random points
        random_params = {
            p : numpy.random.uniform(*RANGES[p], size=self.numtests)
            for p in RANGES.keys()}
        self.tc_ssb = random_params['tc_ssb']
        self.eclipticlongitude_ssb = random_params['eclipticlongitude_ssb']
        self.eclipticlatitude_ssb = random_params['eclipticlatitude_ssb']
        self.polarization_ssb = random_params['polarization_ssb']
        self.tc_lisa = random_params['tc_lisa']
        self.eclipticlongitude_lisa = random_params['eclipticlongitude_lisa']
        self.eclipticlatitude_lisa = random_params['eclipticlatitude_lisa']
        self.polarization_lisa = random_params['polarization_lisa']
        self.tc_geo = random_params['tc_geo']
        self.ra = random_params['ra']
        self.dec = random_params['dec']
        self.polarization_geo = random_params['polarization_geo']

        # calculate derived parameters from each

        self.tc_lisa_derived, self.eclipticlongitude_lisa_derived, \
        self.eclipticlatitude_lisa_derived, self.polarization_lisa_derived = \
            space.ssb_to_lisa(
                self.tc_ssb, self.eclipticlongitude_ssb,
                self.eclipticlatitude_ssb, self.polarization_ssb)

        self.tc_geo_derived, self.ra_derived, \
        self.dec_derived, self.polarization_geo_derived = \
            space.lisa_to_geo(
                self.tc_lisa, self.eclipticlongitude_lisa,
                self.eclipticlatitude_lisa, self.polarization_lisa)

        self.tc_ssb_derived, self.eclipticlongitude_ssb_derived, \
        self.eclipticlatitude_ssb_derived, self.polarization_ssb_derived = \
            space.geo_to_ssb(
                self.tc_geo, self.ra,
                self.dec, self.polarization_geo)

    def test_round_robin(self):
        """Computes inverse transformations to get original parameters from
        derived, then compares them to the original.
        """
        msg = '{} does not recover same {}; max difference: {}; inputs: {}'
        # following lists (function to check,
        #                  arguments to pass to the function,
        #                  name of self's attribute to compare to)
        fchecks = [
            (space.ssb_to_geo,
                (self.tc_ssb_derived, self.eclipticlongitude_ssb_derived,
                 self.eclipticlatitude_ssb_derived,
                 self.polarization_ssb_derived),
                ('tc_geo', 'ra',
                'dec', 'polarization_geo')),
            (space.geo_to_lisa,
                (self.tc_geo_derived, self.ra_derived,
                 self.dec_derived, self.polarization_geo_derived),
                ('tc_lisa', 'eclipticlongitude_lisa',
                'eclipticlatitude_lisa', 'polarization_lisa')),
            (space.lisa_to_ssb,
                (self.tc_lisa_derived, self.eclipticlongitude_lisa_derived,
                 self.eclipticlatitude_lisa_derived,
                 self.polarization_lisa_derived),
                ('tc_ssb', 'eclipticlongitude_ssb',
                'eclipticlatitude_ssb', 'polarization_ssb')),
            ]

        for func, args, compval in fchecks:
            func_tuple = func(*args)
            attr_tuple = [getattr(self, i) for i in compval]
            for i in range(len(func_tuple)):
                derived_val = func_tuple[i]
                check_val = attr_tuple[i]
                passed, maxdiff, maxidx = almost_equal(derived_val, check_val,
                                            self.precision)
                if not passed:
                    failinputs = [p[maxidx] for p in args]
                else:
                    failinputs = None
                self.assertTrue(passed, msg.format(
                    func, compval, maxdiff, failinputs))

    def test_polarization(self):
        """Set up a custom V-shape ground-based detector which is almost 
        co-aligned and co-located to LISA-Z channel. When the long-wavelength
        approximation is valid, the detector responses of those two should be
        very similar. The transform functions in `pycbc.coordinates.space` are
        used in the calculation.
        """
        from pycbc.waveform import get_fd_det_waveform

        from pycbc.detector import Detector, add_detector_on_earth
        from pycbc.waveform import get_fd_waveform
        import importlib
        
        YRSID_SI = 31558149.763545603

        def is_module_installed(module_name):
            try:
                importlib.import_module(module_name)
                return True
            except ImportError:
                return False

        def strain_generator(det='D1', model='IMRPhenomD', fs=4096, df=None,
                             flow=2, fref=2, tc=0, params=None,
                             fd_taper=False):

            # Generate a waveform at the detector-frame. Use the same
            # `delta_f` as the LISA-side waveform (passed in via `df`), so
            # both sides land on the same frequency grid; a mismatched grid
            # means the "nearest bin to f_ref" on each side is not really
            # the same frequency, which becomes significant for a fast-
            # evolving phase near f_ref for these masses.
            hp, hc = get_fd_waveform(approximant=model,
                        mass1=params['mass1'], mass2=params['mass2'],
                        spin1x=0, spin1y=0,
                        spin1z=params['spin1z'], spin2x=0,
                        spin2y=0, spin2z=params['spin2z'],
                        distance=params['distance'],
                        coa_phase=params['coa_phase'],
                        inclination=params['inclination'], f_lower=flow,
                        f_ref=fref, delta_t=1.0/fs, delta_f=df)

            # Set merger time to 'tc'.
            hp.start_time += tc
            hc.start_time += tc

            # Project GW waveform onto GW detectors.
            ra = params['ra']
            dec = params['dec']
            psi = params['polarization']
            time = hp.start_time

            det_1 = Detector("D1")
            fp_1, fc_1 = det_1.antenna_pattern(
                            right_ascension=ra, declination=dec,
                            polarization=psi, t_gps=tc)

            # Not take the rotation of the earth into account.
            ht_1 = fp_1*hp + fc_1*hc
            ht_list = [ht_1]

            return ht_list


        # Checking if bbhx is installed, if not, then ignore this test
        # and raise a warning.
        # TODO: we need to implement a better way to install bbhx package.
        bbhx_installed = is_module_installed('bbhx')
        if not bbhx_installed:
            passed = True
            print("Ignore polarization test, because bbhx is not installed.")
        else:
            # cunstom E3
            # All of those hard-coded numbers are carefully chosen,
            # in order to almost co-align and co-locate both 2 detectors.
            fine_tunning = 7992204.094271309
            OMEGA_0 = 1.99098659277e-7
            yangle = numpy.pi / 2 + fine_tunning * OMEGA_0
            add_detector_on_earth(name='D1', longitude=1.8895427761465164,
                                latitude=0.11450614784814996, yangle=yangle,
                                xangle=yangle+numpy.pi/3, height=0)

            # set parameters
            params = {}
            params['tdi'] = '1.5'
            params['ref_frame'] = 'SSB'
            params['approximant'] = 'BBHX_PhenomD'
            params['base_approximant'] = 'BBHX_PhenomD'
            params['coa_phase'] = 0.0
            params['mass1'] = 1e6
            params['mass2'] = 8e5
            params['spin1z'] = 0.0
            params['spin2z'] = 0.0
            params['distance'] = 410
            params['inclination'] = numpy.pi/2  # edge-on
            params['t_obs_start'] = 1*YRSID_SI
            params['delta_f'] = 1./params['t_obs_start']
            params['f_lower'] = 1e-4
            params['f_ref'] = 8e-4
            params['f_final'] = 0.1
            params['delta_t'] = 1/0.2
            params['t_offset'] = 9206958.120016199 # 0 degrees
            t_lisa = YRSID_SI - params['t_offset'] + fine_tunning

            nx = 100
            longitude_array_high_res = numpy.linspace(
                            0, 2*numpy.pi, num=nx, endpoint=False)

            amp_E3_psi = []
            amp_LISA3_psi = []
            phase_diff_psi = []

            for polarization_lisa in numpy.linspace(
                    0, 2*numpy.pi, endpoint=False, num=nx):
                params['tc'], params['eclipticlongitude'], \
                params['eclipticlatitude'], params['polarization'] = \
                    space.lisa_to_ssb(t_lisa, 0, numpy.pi/4,
                                polarization_lisa, params['t_offset'])

                nparams = {'mass1':params['mass1'], 'mass2':params['mass2'],
                            'spin1z':params['spin1z'],
                            'spin2z':params['spin2z'],
                            'f_lower':params['f_lower']}

                bbhx_fd = get_fd_det_waveform(
                            ifos=['LISA_A','LISA_E','LISA_T'], **params)
                lisa_X_fd = -numpy.sqrt(2)/2 * bbhx_fd['LISA_A'] + \
                            numpy.sqrt(6)/6 * bbhx_fd['LISA_E'] + \
                            numpy.sqrt(3)/3 * bbhx_fd['LISA_T']
                lisa_Y_fd = -numpy.sqrt(6)/3 * bbhx_fd['LISA_E'] + \
                            numpy.sqrt(3)/3 * bbhx_fd['LISA_T']
                lisa_Z_fd = numpy.sqrt(2)/2 * bbhx_fd['LISA_A'] + \
                            numpy.sqrt(6)/6 * bbhx_fd['LISA_E'] + \
                            numpy.sqrt(3)/3 * bbhx_fd['LISA_T']
                channel_xyz_fd = {'LISA_X':lisa_X_fd,
                                'LISA_Y':lisa_Y_fd,
                                'LISA_Z':lisa_Z_fd}

                t_geo, ra, dec, polarization_geo = \
                    space.ssb_to_geo(
                        params['tc'],
                        params['eclipticlongitude'], 
                        params['eclipticlatitude'],
                        params['polarization']
                    )

                params_3g = params.copy()
                params_3g['tc'] = t_geo
                params_3g['ra'] = ra
                params_3g['dec'] = dec
                params_3g['polarization'] = polarization_geo
                params_3g['coa_phase'] = numpy.pi/2 - params['coa_phase']

                E3_signal = strain_generator(det='D1', model='IMRPhenomD',
                                            fs=1./params_3g['delta_t'],
                                            df=params_3g['delta_f'],
                                            flow=params_3g['f_lower'],
                                            fref=params_3g['f_ref'],
                                            tc=params_3g['tc'],
                                            params=params_3g,
                                            fd_taper=False)
                E3_signal_fd = {
                    'DIY_E3':E3_signal[0].to_frequencyseries()
                }

                # Use the nearest frequency bin to f_ref on each side, rather
                # than the first bin within an absolute tolerance: with
                # `delta_f` now matched between E3 and LISA that tolerance
                # window can still span several bins, and phase evolves fast
                # enough at f_ref for these masses that which of those bins
                # gets used matters.
                index_E3 = numpy.argmin(numpy.abs(
                            E3_signal_fd['DIY_E3'].sample_frequencies -
                            params_3g['f_ref']))
                index_LISA3 = numpy.argmin(numpy.abs(
                            bbhx_fd['LISA_A'].sample_frequencies -
                            params['f_ref']))
                val_E3 = E3_signal_fd['DIY_E3'][index_E3]
                val_LISA3 = channel_xyz_fd['LISA_Z'][index_LISA3]

                amp_E3_psi.append(numpy.abs(val_E3))
                amp_LISA3_psi.append(numpy.abs(val_LISA3))
                # Compare phase through the angle of E3 relative to LISA-Z,
                # instead of differencing two independently-computed
                # `numpy.angle` values in degrees: two phases that are a
                # fraction of a degree apart, but fall on either side of the
                # +-180 deg branch cut, would otherwise show up as a
                # spurious ~360 deg difference.
                phase_diff_psi.append(numpy.angle(val_E3 * numpy.conj(val_LISA3)))

            dist_amp = curve_similarity(amp_E3_psi/numpy.mean(amp_E3_psi),
                                        amp_LISA3_psi/numpy.mean(
                                            amp_LISA3_psi))
            # E3 and LISA-Z are generated from different waveform models
            # (IMRPhenomD vs BBHX_PhenomD) with different reference-phase
            # conventions, so a constant offset between their phases is
            # expected and not itself meaningful. What the frame/
            # polarization transforms in `pycbc.coordinates.space` are
            # responsible for is that this offset stays constant as the
            # polarization angle changes -- a bug there (e.g. a sign error
            # or a missing factor of 2 in the polarization angle) would show
            # up as the relative phase oscillating with the polarization
            # angle instead of sitting at a fixed value. So we unwrap the
            # relative phase across the polarization-angle sweep and check
            # how much it scatters around its mean, rather than comparing
            # its absolute value to a hard-coded number.
            phase_diff_psi = numpy.unwrap(phase_diff_psi)
            phase_resid_std = numpy.rad2deg(
                numpy.std(phase_diff_psi - numpy.mean(phase_diff_psi)))
            # Here we compare the amplitude as a function of the
            # polarization angle in the LISA frame, because E3 and LISA-Z are
            # almost co-aligned and co-located, their detector response
            # should be very similar (when long-wavelength approximation
            # holds); `dist_amp` should stay close to the hard-coded
            # reference value below (obtained on the author's laptop) if the
            # frame transform functions work correctly. Users can modify
            # this script to plot amp_E3_psi, amp_LISA3_psi, and
            # phase_diff_psi to inspect the curves directly.
            if (numpy.abs(dist_amp - 0.2838670151034317) < 1e-2) and \
                (phase_resid_std < 2.0):
                passed = True
            else:
                passed = False

        self.assertTrue(passed) 


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestParams))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)