import copy
import unittest
import random
import os
import numpy
from numpy import complex128, real, sqrt, sin, cos, angle, ceil, log
from numpy import zeros, argmax, array
from astropy.utils.data import download_file
from pycbc import DYN_RANGE_FAC
from pycbc.waveform import get_td_waveform, get_fd_waveform, td_approximants, fd_approximants
from pycbc.pnutils import nearest_larger_binary_number
from pycbc.types import FrequencySeries, TimeSeries, complex_same_precision_as
from pycbc.types import load_frequencyseries
from pycbc.filter import sigmasq, overlap_cplx, matched_filter_core
from pycbc.filter import compute_max_snr_over_sky_loc_stat
from pycbc.filter import compute_max_snr_over_sky_loc_stat_no_phase
from pycbc.filter import compute_u_val_for_sky_loc_stat_no_phase
from pycbc.filter import compute_u_val_for_sky_loc_stat
from pycbc import psd
from pycbc import vetoes
from utils import parse_args_all_schemes, simple_exit

_scheme, _context = parse_args_all_schemes("correlate")

expected_results = {}
for idx in range(4):
    expected_results[idx] = {}
    for jdx in range(4):
        expected_results[idx][jdx] = {}

expected_results[0][0]['Ip_snr'] = 100.0
expected_results[0][0]['Ip_angle'] = -1.88619488652e-18
expected_results[0][0]['Ip_argmax'] = 0
expected_results[0][0]['Ic_snr'] = 98.7349100759
expected_results[0][0]['Ic_angle'] = 1.60960753393
expected_results[0][0]['Ic_argmax'] = 3

expected_results[0][1]['Ip_snr'] = 96.3390579783
expected_results[0][1]['Ip_angle'] = 0.511744420131
expected_results[0][1]['Ip_argmax'] = 1387
expected_results[0][1]['Ic_snr'] = 96.3390579783
expected_results[0][1]['Ic_angle'] = 2.08254074693
expected_results[0][1]['Ic_argmax'] = 1387

expected_results[0][2]['Ip_snr'] = 98.8434423546
expected_results[0][2]['Ip_angle'] = 0.566451407787
expected_results[0][2]['Ip_argmax'] = 1100
expected_results[0][2]['Ic_snr'] = 98.8485523538
expected_results[0][2]['Ic_angle'] = 2.10418718318
expected_results[0][2]['Ic_argmax'] = 1099

expected_results[0][3]['Ip_snr'] = 96.4239530554
expected_results[0][3]['Ip_angle'] = -0.946889162447
expected_results[0][3]['Ip_argmax'] = 1447
expected_results[0][3]['Ic_snr'] = 95.6528566731
expected_results[0][3]['Ic_angle'] = 1.10466400896
expected_results[0][3]['Ic_argmax'] = 1484

expected_results[1][0]['Ip_snr'] = 96.3390579783
expected_results[1][0]['Ip_angle'] = -0.511744420131
expected_results[1][0]['Ip_argmax'] = 260757
expected_results[1][0]['Ic_snr'] = 94.4282712604
expected_results[1][0]['Ic_angle'] = 0.957548192904
expected_results[1][0]['Ic_argmax'] = 260753

expected_results[1][1]['Ip_snr'] = 100.0
expected_results[1][1]['Ip_angle'] = -1.73241699772e-18
expected_results[1][1]['Ip_argmax'] = 0
expected_results[1][1]['Ic_snr'] = 100.0
expected_results[1][1]['Ic_angle'] = 1.57079632679
expected_results[1][1]['Ic_argmax'] = 0

expected_results[1][2]['Ip_snr'] = 97.6397283701
expected_results[1][2]['Ip_angle'] = 0.156578884423
expected_results[1][2]['Ip_argmax'] = 261862
expected_results[1][2]['Ic_snr'] = 97.7584101045
expected_results[1][2]['Ic_angle'] = 1.69481552225
expected_results[1][2]['Ic_argmax'] = 261861

expected_results[1][3]['Ip_snr'] = 99.4434573331
expected_results[1][3]['Ip_angle'] = -1.50330148916
expected_results[1][3]['Ip_argmax'] = 58
expected_results[1][3]['Ic_snr'] = 98.0771994342
expected_results[1][3]['Ic_angle'] = 0.521399663782
expected_results[1][3]['Ic_argmax'] = 94

expected_results[2][0]['Ip_snr'] = 98.8434423546
expected_results[2][0]['Ip_angle'] = -0.566451407787
expected_results[2][0]['Ip_argmax'] = 261044
expected_results[2][0]['Ic_snr'] = 96.8727372348
expected_results[2][0]['Ic_angle'] = 0.921151545297
expected_results[2][0]['Ic_argmax'] = 261041

expected_results[2][1]['Ip_snr'] = 97.6397283701
expected_results[2][1]['Ip_angle'] = -0.156578884423
expected_results[2][1]['Ip_argmax'] = 282
expected_results[2][1]['Ic_snr'] = 97.6397283701
expected_results[2][1]['Ic_angle'] = 1.41421744237
expected_results[2][1]['Ic_argmax'] = 282

expected_results[2][2]['Ip_snr'] = 100.0
expected_results[2][2]['Ip_angle'] = 2.41820532326e-18
expected_results[2][2]['Ip_argmax'] = 0
expected_results[2][2]['Ic_snr'] = 99.988543227
expected_results[2][2]['Ic_angle'] = 1.5377922043
expected_results[2][2]['Ic_argmax'] = 262143

expected_results[2][3]['Ip_snr'] = 97.3725007917
expected_results[2][3]['Ip_angle'] = -1.61086911817
expected_results[2][3]['Ip_argmax'] = 342
expected_results[2][3]['Ic_snr'] = 96.2035744912
expected_results[2][3]['Ic_angle'] = 0.442670138337
expected_results[2][3]['Ic_argmax'] = 379

expected_results[3][0]['Ip_snr'] = 96.4239530554
expected_results[3][0]['Ip_angle'] = 0.946889162447
expected_results[3][0]['Ip_argmax'] = 260697
expected_results[3][0]['Ic_snr'] = 94.4958639934
expected_results[3][0]['Ic_angle'] = 2.41579357775
expected_results[3][0]['Ic_argmax'] = 260693

expected_results[3][1]['Ip_snr'] = 99.4434573331
expected_results[3][1]['Ip_angle'] = 1.50330148916
expected_results[3][1]['Ip_argmax'] = 262086
expected_results[3][1]['Ic_snr'] = 99.4434573331
expected_results[3][1]['Ic_angle'] = 3.07409781595
expected_results[3][1]['Ic_argmax'] = 262086

expected_results[3][2]['Ip_snr'] = 97.3725007917
expected_results[3][2]['Ip_angle'] = 1.61086911817
expected_results[3][2]['Ip_argmax'] = 261802
expected_results[3][2]['Ic_snr'] = 97.3906866656
expected_results[3][2]['Ic_angle'] = -3.11216854627
expected_results[3][2]['Ic_argmax'] = 261802

expected_results[3][3]['Ip_snr'] = 100.0
expected_results[3][3]['Ip_angle'] = 9.03608368726e-19
expected_results[3][3]['Ip_argmax'] = 0
expected_results[3][3]['Ic_snr'] = 99.4335056063
expected_results[3][3]['Ic_angle'] = 2.02876392072
expected_results[3][3]['Ic_argmax'] = 36

def generate_detector_strain(template_params, h_plus, h_cross):
    polarization = 0

    if hasattr(template_params, 'polarization'):
        polarization = template_params.polarization

    f_plus = cos(polarization)
    f_cross = sin(polarization)

    return h_plus * f_plus + h_cross * f_cross

def make_padded_frequency_series(vec, filter_N=None, delta_f=None):
    """Convert vec (TimeSeries or FrequencySeries) to a FrequencySeries. If
    filter_N and/or delta_f are given, the output will take those values. If
    not told otherwise the code will attempt to pad a timeseries first such that
    the waveform will not wraparound. However, if delta_f is specified to be
    shorter than the waveform length then wraparound *will* be allowed.
    """
    if filter_N is None:
        power = ceil(log(len(vec), 2)) + 1
        N = 2 ** power
    else:
        N = filter_N
    n = N / 2 + 1

    if isinstance(vec, FrequencySeries):
        vectilde = FrequencySeries(zeros(n, dtype=complex_same_precision_as(vec)),
                                   delta_f=1.0, copy=False)
        if len(vectilde) < len(vec):
            cplen = len(vectilde)
        else:
            cplen = len(vec)
        vectilde[0:cplen] = vec[0:cplen]
        delta_f = vec.delta_f

    elif isinstance(vec, TimeSeries):
        # First determine if the timeseries is too short for the specified df
        # and increase if necessary
        curr_length = len(vec)
        new_length = int(nearest_larger_binary_number(curr_length))
        while new_length * vec.delta_t < 1./delta_f:
            new_length = new_length * 2
        vec.resize(new_length)
        # Then convert to frequencyseries
        v_tilde = vec.to_frequencyseries()
        # Then convert frequencyseries to required length and spacing by keeping
        # only every nth sample if delta_f needs increasing, and cutting at
        # Nyquist if the max frequency is too high.
        # NOTE: This assumes that the input and output data is using binary
        #       lengths.
        i_delta_f = v_tilde.get_delta_f()
        v_tilde = v_tilde.numpy()
        df_ratio = int(delta_f / i_delta_f)
        n_freq_len = int((n-1) * df_ratio +1)
        assert(n <= len(v_tilde))
        df_ratio = int(delta_f / i_delta_f)
        v_tilde = v_tilde[:n_freq_len:df_ratio]
        vectilde = FrequencySeries(v_tilde, delta_f=delta_f, dtype=complex128)

    return FrequencySeries(vectilde * DYN_RANGE_FAC, delta_f=delta_f,
                           dtype=complex128)

def get_waveform(wf_params, start_frequency, sample_rate, length,
                 filter_rate, sky_max_template=False):
    delta_f = filter_rate / float(length)
    if wf_params.approximant in fd_approximants():
        hp, hc = get_fd_waveform(wf_params, delta_f=delta_f,
                                 f_lower=start_frequency)

    elif wf_params.approximant in td_approximants():
        hp, hc = get_td_waveform(wf_params, delta_t=1./sample_rate,
                                 f_lower=start_frequency)

    if not sky_max_template:
        hvec = generate_detector_strain(wf_params, hp, hc)
        return make_padded_frequency_series(hvec, length, delta_f=delta_f)
    else:
        return make_padded_frequency_series(hp, length, delta_f=delta_f), \
            make_padded_frequency_series(hc, length, delta_f=delta_f)

class DummyClass(object):
    pass

class TestChisq(unittest.TestCase):
    __test__ = False
    def setUp(self, *args):
        # Where are my data files?
        self.context = _context
        self.scheme = _scheme
        self.tolerance = 1e-6
        self.filter_t_length = 16
        self.low_freq_filter = 30.
        self.sample_rate = 16384
        self.filter_N = int(self.filter_t_length * self.sample_rate)
        self.filter_n = int(self.filter_N / 2 + 1)
        self.filter_delta_f = 1.0 / self.filter_t_length
        self.psd = psd.from_string('aLIGOZeroDetHighPowerGWINC',
                                   self.filter_n, self.filter_delta_f,
                                   self.low_freq_filter)
        self.psd *= DYN_RANGE_FAC*DYN_RANGE_FAC
        wps1 = DummyClass()
        wps1.mass1 = 123.7627
        wps1.mass2 = 72.55471
        wps1.inclination = 1.125029
        wps1.coa_phase = 2.906049
        wps1.approximant='EOBNRv2HM'
        self.wps1 = wps1

        wps2 = DummyClass()
        wps2.mass1 = 131.460647583
        wps2.mass2 = 69.0030059814
        wps2.inclination = 0.8432287
        wps2.coa_phase = 0.2
        wps2.approximant = 'SEOBNRv4_ROM'
        self.wps2 = wps2

        wps3 = copy.deepcopy(wps2)
        wps3.approximant = 'EOBNRv2HM_ROM'
        self.wps3 = wps3

        wps4 = copy.deepcopy(wps2)
        wps4.spin1x = 0.8
        wps4.spin2y = -0.9
        wps4.approximant = 'IMRPhenomPv2'
        self.wps4 = wps4

        self.wps_list = [wps1, wps2, wps3, wps4]

        self.sm_power_chisq = vetoes.SingleDetSkyMaxPowerChisq(num_bins='1')
        self.sm_power_chisq2 = vetoes.SingleDetSkyMaxPowerChisq(num_bins='10')
        self.power_chisq = vetoes.SingleDetPowerChisq(num_bins='1')
        self.power_chisq2 = vetoes.SingleDetPowerChisq(num_bins='10')


    def test_filtering(self):
        idx = self.idx
        jdx = self.jdx
        # Uncomment these lines if needing to regenerate data files
        #w1 = self.wps_list[idx]
        #w2 = self.wps_list[jdx]
        #stilde = get_waveform(w1, self.low_freq_filter-1,
        #                      self.sample_rate, self.filter_N,
        #                      self.sample_rate)
        #try:
        #    stilde.save('data/skymaxtest_stilde_%d.hdf' % idx)
        #except:
        #    pass
        url = ('https://github.com/gwastro/pycbc-config/raw/master/'
               'test_data_files/{}')
        fname = f'skymaxtest_stilde_{idx}.hdf'
        apy_fname = download_file(url.format(fname), cache=False)
        # Astropy will not download with the .hdf extension, which we need,
        # so symlink
        os.symlink(apy_fname, fname)

        stilde = load_frequencyseries(fname)
        os.unlink(fname)
        s_norm = sigmasq(stilde, psd=self.psd,
                         low_frequency_cutoff=self.low_freq_filter)
        stilde /= sqrt(float(s_norm))
        stilde *= 100
        # Uncomment these lines if needing to regenerate data files
        #hplus, hcross = get_waveform(w2, self.low_freq_filter-1,
        #                             self.sample_rate, self.filter_N,
        #                             self.sample_rate, sky_max_template=True)
        #try:
        #    hplus.save('data/skymaxtest_hplus_%d.hdf' % jdx)
        #    hcross.save('data/skymaxtest_hcross_%d.hdf' % jdx)
        #except:
        #    pass
        fname = f'skymaxtest_hplus_{jdx}.hdf'
        apy_fname = download_file(url.format(fname), cache=False)
        # Astropy will not download with the .hdf extension, which we need,
        # so symlink
        os.symlink(apy_fname, fname)
        hplus = load_frequencyseries(fname)
        os.unlink(fname)

        fname = f'skymaxtest_hcross_{jdx}.hdf'
        apy_fname = download_file(url.format(fname), cache=False)
        # Astropy will not download with the .hdf extension, which we need,
        # so symlink
        os.symlink(apy_fname, fname)
        hcross = load_frequencyseries(fname)
        os.unlink(fname)

        hplus.f_lower = self.low_freq_filter
        hplus.params = random.randint(0,100000000000)
        hcross.f_lower = self.low_freq_filter
        hcross.params = random.randint(0,100000000000)
        hp_norm = sigmasq(hplus, psd=self.psd,
                          low_frequency_cutoff=self.low_freq_filter)
        hc_norm = sigmasq(hcross, psd=self.psd,
                          low_frequency_cutoff=self.low_freq_filter)
        hplus /= sqrt(float(hp_norm))
        hcross /= sqrt(float(hc_norm))
        hpc_corr = overlap_cplx(hplus, hcross, psd=self.psd,
                                low_frequency_cutoff=self.low_freq_filter,
                                normalized=False)
        hpc_corr_R = real(hpc_corr)
        I_plus, corr_plus, n_plus = matched_filter_core\
            (hplus, stilde, psd=self.psd,
             low_frequency_cutoff=self.low_freq_filter, h_norm=1.)
        # FIXME: Remove the deepcopies before merging with master
        I_plus = copy.deepcopy(I_plus)
        corr_plus = copy.deepcopy(corr_plus)
        I_cross, corr_cross, n_cross = matched_filter_core\
            (hcross, stilde, psd=self.psd,
             low_frequency_cutoff=self.low_freq_filter, h_norm=1.)
        I_cross = copy.deepcopy(I_cross)
        corr_cross = copy.deepcopy(corr_cross)
        I_plus = I_plus * n_plus
        I_cross = I_cross * n_cross
        IPM = abs(I_plus.data).argmax()
        ICM = abs(I_cross.data).argmax()
        self.assertAlmostEqual(abs(I_plus[IPM]),
                               expected_results[idx][jdx]['Ip_snr'])
        self.assertAlmostEqual(angle(I_plus[IPM]),
                               expected_results[idx][jdx]['Ip_angle'])
        self.assertEqual(IPM, expected_results[idx][jdx]['Ip_argmax'])
        self.assertAlmostEqual(abs(I_cross[ICM]),
                               expected_results[idx][jdx]['Ic_snr'])
        self.assertAlmostEqual(angle(I_cross[ICM]),
                               expected_results[idx][jdx]['Ic_angle'])
        self.assertEqual(ICM, expected_results[idx][jdx]['Ic_argmax'])

        #print "expected_results[{}][{}]['Ip_snr'] = {}" .format(idx,jdx,abs(I_plus[IPM]))
        #print "expected_results[{}][{}]['Ip_angle'] = {}".format(idx,jdx,angle(I_plus[IPM]))
        #print "expected_results[{}][{}]['Ip_argmax'] = {}".format(idx,jdx, IPM)
        #print "expected_results[{}][{}]['Ic_snr'] = {}" .format(idx,jdx,abs(I_cross[ICM]))
        #print "expected_results[{}][{}]['Ic_angle'] = {}".format(idx,jdx,angle(I_cross[ICM]))
        #print "expected_results[{}][{}]['Ic_argmax'] = {}".format(idx,jdx, ICM)

        det_stat_prec = compute_max_snr_over_sky_loc_stat\
            (I_plus, I_cross, hpc_corr_R, hpnorm=1., hcnorm=1.,
             thresh=0.1, analyse_slice=slice(0,len(I_plus.data)))
        det_stat_hom = compute_max_snr_over_sky_loc_stat_no_phase\
            (I_plus, I_cross, hpc_corr_R, hpnorm=1., hcnorm=1.,
             thresh=0.1, analyse_slice=slice(0,len(I_plus.data)))
        idx_max_prec = argmax(det_stat_prec.data)
        idx_max_hom = argmax(det_stat_hom.data)
        max_ds_prec = det_stat_prec[idx_max_prec]
        max_ds_hom = det_stat_hom[idx_max_hom]

        uvals_prec, _ = compute_u_val_for_sky_loc_stat\
            (I_plus.data, I_cross.data, hpc_corr_R, indices=[idx_max_prec],
             hpnorm=1., hcnorm=1.)
        with numpy.errstate(divide="ignore"):
            uvals_hom, _ = compute_u_val_for_sky_loc_stat_no_phase\
                (I_plus.data, I_cross.data, hpc_corr_R, indices=[idx_max_hom],
                 hpnorm=1., hcnorm=1.)

        ht = hplus * uvals_hom[0] + hcross
        ht_norm = sigmasq(ht, psd=self.psd,
                          low_frequency_cutoff=self.low_freq_filter)
        ht /= sqrt(float(ht_norm))
        ht.f_lower = self.low_freq_filter
        ht.params = random.randint(0,100000000000)
        I_t, corr_t, n_t = matched_filter_core\
            (ht, stilde, psd=self.psd,
             low_frequency_cutoff=self.low_freq_filter, h_norm=1.)
        I_t = I_t * n_t
        self.assertAlmostEqual(abs(real(I_t.data[idx_max_hom])), max_ds_hom)
        self.assertEqual(abs(real(I_t.data[idx_max_hom])),
                         max(abs(real(I_t.data))))
        with numpy.errstate(invalid='ignore', divide='ignore'):
            chisq, _ = self.power_chisq.values\
                (corr_t, array([max_ds_hom]) / n_plus, n_t,
                 self.psd, array([idx_max_hom]), ht)

        ht = hplus * uvals_prec[0] + hcross
        ht_norm = sigmasq(ht, psd=self.psd,
                          low_frequency_cutoff=self.low_freq_filter)
        ht /= sqrt(float(ht_norm))
        ht.f_lower = self.low_freq_filter
        ht.params = random.randint(0,100000000000)
        I_t, corr_t, n_t = matched_filter_core\
            (ht, stilde, psd=self.psd,
             low_frequency_cutoff=self.low_freq_filter, h_norm=1.)
        I_t = I_t * n_t

        with numpy.errstate(divide="ignore", invalid='ignore'):
            chisq, _ = self.power_chisq.values\
                (corr_t, array([max_ds_prec]) / n_plus, n_t, self.psd,
                 array([idx_max_prec]), ht)

        self.assertAlmostEqual(abs(I_t.data[idx_max_prec]), max_ds_prec)
        self.assertEqual(idx_max_prec, abs(I_t.data).argmax())
        self.assertTrue(chisq < 1E-4)



def skymax_test_maker(class_name, idx, jdx):
    class Test(class_name):
        __test__ = True
        idx = idx
        jdx = jdx

    Test.__name__ = "Test %s" % '_'.join([str(idx),str(jdx)])
    return Test


suite = unittest.TestSuite()
for idx in range(4):
    for jdx in range(4):
        curr_cls = skymax_test_maker(TestChisq, idx, jdx)
        vars()[curr_cls.__name__] = curr_cls
        suite.addTest(unittest.TestLoader().loadTestsFromTestCase(curr_cls))
        del curr_cls

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)

