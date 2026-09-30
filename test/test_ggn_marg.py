import os
import unittest
import shutil
from pycbc.inference import models
from pycbc.workflow import WorkflowConfigParser
import numpy
import tempfile
import subprocess
import sys
from utils import simple_exit

TESTDIR = os.path.dirname(os.path.abspath(__file__))
DATADIR = os.path.join(TESTDIR, 'data', 'ggn_marg')
CREATE_INJECTIONS = os.path.join(TESTDIR, '..', 'bin',
                                 'pycbc_create_injections')

MAXL_PHI220 = 0.66926183
TEMPLATE_PARAMS = dict(ra = 3.5,
                       dec = 0.73,
                       delta_tc = 0.01525955,
                       inclination = 0.52359877,
                       final_mass = 305.33635657,
                       final_spin = 0.49861983,
                       polarization = 1.62913401,
                       amp220 = 7.41039482e-21,
                       amp330 = 0.14645439,
                       phi330 = 0.00250140 - MAXL_PHI220,
                       phi220 = 0.,
                       delta_phi330 = 0.00250140 - MAXL_PHI220, # brute-force
                       amp220_snr = 91.2 # approx. SNR for inj
                       )

class TestMargModels(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # create the injection
        cls.tmpdir = tempfile.mkdtemp()
        cls.injfile = os.path.join(cls.tmpdir, 'injection.hdf')
        subprocess.run(
            [sys.executable, CREATE_INJECTIONS,
             '--config-files', os.path.join(DATADIR, 'injection.ini'),
             '--ninjections', '1', '--seed', '10',
             '--output-file', cls.injfile,
             '--variable-params-section', 'variable_params',
             '--static-params-section', 'static_params',
             '--dist-section', 'prior', '--force'],
            check=True)
        
        # evaluated margphase likelihood
        cls.margphase_margl = -5.326111182154365

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir)
        
    @classmethod
    def load_model(cls, model_name, fid_amp=None, phase_samples=500000):
        # load one of the models from given configs
        filename = os.path.join(DATADIR, f'{model_name}.ini')
        cp = WorkflowConfigParser(configFiles=[filename])
        cp.set('data', 'injection-file', cls.injfile)
        
        # opt to change fiducial amplitude for multimargphase
        if fid_amp is not None:
            cp.set('model', 'fiducial-amp-value', str(fid_amp))
            
        # opt to change # of phase sample points for brute-force marg
        if fid_amp is not None:
            cp.set('model', 'phase-samples', str(phase_samples))
        
        # load in nominal params to model
        model = models.read_from_config(cp)
        model.update(**TEMPLATE_PARAMS)
        maxl = model.loglikelihood
        return model, maxl
    
    def _test_marginalization(self, model, marglogl):
        # generic marginalization tests
        maxl = model.current_stats['maxl_logl']
        maxl_phase = model.current_stats['maxl_phase']
        
        # test that the maxL loglikelihood is ~0
        self.assertAlmostEqual(maxl, 0., places = 1,
                              msg=f"calculated maximum loglikelihood {maxl} "
                                   "is not close to zero")
        
        # test that the marginalization recovers the right phase
        self.assertAlmostEqual(maxl_phase, MAXL_PHI220, delta = 0.001,
                              msg=f"maximum likelihood phase {maxl_phase} "
                                   "is not close to injected value "
                                   f"{MAXL_PHI220}")
        
        # test that the marginalized likelihood is strictly less than the maxL
        self.assertTrue(marglogl < maxl,
                       msg="total marginalized likelihood is greater than "
                           "maximum loglikelihood")

    def test_margphase(self):
        '''Test the marginalized phase model.
        '''
        model, marglogl = self.load_model("gated_gaussian_margphase")
        self._test_marginalization(model, marglogl)
        
        # save likelihood value for comparisons
        self.margphase_margl = marglogl
    
    def test_multimargphase_amps(self):
        '''Test the multimargphase model, sampling in amplitude. This should
        give the exact same results as the margphase model when provided
        with the same parameters.
        '''
        model, marglogl = self.load_model('gated_gaussian_multimargphase_amps')
        self._test_marginalization(model, marglogl)
        
        # test that the likelihood is close to the marginalized phase model
        self.assertAlmostEqual(marglogl, self.margphase_margl,
                               delta = 0.001,
                               msg=f"marginalized loglikelihood {marglogl} "
                                   f"does not match margphase "
                                   f"{self.margphase_margl}")
        
    def _test_multimargphase(self, fid_amp=1e-20):
        # Generic multimargphase tests given fiducial amp
        model, marglogl = self.load_model('gated_gaussian_multimargphase',
                                          fid_amp=fid_amp)
        self._test_marginalization(model, marglogl)

        # test that the likelihood is close to the marginalized phase model
        self.assertAlmostEqual(marglogl, self.margphase_margl,
                               delta = 0.001,
                               msg=f"marginalized loglikelihood {marglogl} "
                                   f"does not match margphase "
                                   f"{self.margphase_margl}")

        # test that the loaded SNR reconstructs the correct amplitude
        scale_factor = model.current_stats['scale_factor_220']
        fid_amp = model.fiducial_amp_value
        self.assertAlmostEqual(scale_factor * fid_amp,
                               TEMPLATE_PARAMS['amp220'],
                               delta = 0.001,
                               msg=f"calculated amp {scale_factor * fid_amp} "
                                   f"does not match injected value "
                                   f"{TEMPLATE_PARAMS['amp220']}")

    def test_multimargphase(self):
        '''Test the multimargphase model under change of fiducial amplitude.'''
        log_fid_amps = numpy.arange(-25, 1).astype(numpy.float64)
        fid_amps = 10.**log_fid_amps
        for fa in fid_amps:
            self._test_multimargphase(fid_amp=fa)
   
    def test_brute_force_marg(self):
        '''Test that the margphase model matches brute force marginalization.
        '''
        from scipy.special import logsumexp
        margphase_model, mp_margl = self.load_model("gated_gaussian_margphase",
                                                   phase_samples=1000)
        bruteforce_model, _ = self.load_model("brute_force")
        brute_samples = numpy.linspace(0, 2*numpy.pi, 1000,
                                       endpoint=False)

        # cycle over sample points and brute-force marginalize
        bf_logls = []
        params = TEMPLATE_PARAMS.copy()
        for phi in brute_samples:
            params['phi220'] = phi
            bruteforce_model.update(**params)
            bf_logls.append(bruteforce_model.loglr)
        bf_margl = logsumexp(bf_logls) + bruteforce_model.lognl - \
                    numpy.log(1000)
        
        self.assertAlmostEqual(bf_margl, mp_margl, delta = 0.01,
                               msg=f"Brute-force likelihood {bf_margl} "
                                   f"does not match marginalized model "
                                   f"{mp_margl}")
        
        # test that maxl occurs at correct point
        maxlidx = numpy.array(bf_logls).argmax()
        bf_phase = brute_samples[maxlidx]
        mp_phase = margphase_model.current_stats['maxl_phase']
        self.assertAlmostEqual(bf_phase, mp_phase, delta=0.01,
                               msg=f"Brute-force maxl phase {bf_phase} "
                                   f"does not match marginalized model "
                                   f"{mp_phase}")

        # test likelihood at injected phase
        params['phi220'] = MAXL_PHI220
        bruteforce_model.update(**params)
        self.assertAlmostEqual(bruteforce_model.loglikelihood, 
                               margphase_model.current_stats['maxl_logl'],
                               delta=0.01,
                               msg=f"Brute-force maxL "
                                   f"{bruteforce_model.loglikelihood} "
                                   f"does not match marginalized model "
                                   f"{margphase_model.current_stats['maxl_logl']}")


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestMargModels))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)