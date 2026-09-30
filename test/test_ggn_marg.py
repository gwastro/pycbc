import os
import unittest
import shutil
from pycbc.inference import io, models
from pycbc.workflow import WorkflowConfigParser
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
                       amp220_snr = 91.2)

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
        
        # placeholders for comparisons
        cls.margphase_margl = 0.

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir)
        
    @classmethod
    def load_model(cls, model_name):
        # load one of the models from given configs
        filename = os.path.join(DATADIR, f'{model_name}.ini')
        cp = WorkflowConfigParser(configFiles=[filename])
        cp.set('data', 'injection-file', filename)
        model = models.read_from_config(cp)
        
        # load in nominal params
        model.update(**TEMPLATE_PARAMS)
        maxl = model.loglikelihood
        return model, maxl
    
    @classmethod
    def test_marginalization(cls, model, marglogl):
        # generic marginalization tests
        maxl = model.current_stats['maxl_logl']
        maxl_phase = model.current_stats['maxl_phase']
        
        # test that the maxL loglikelihood is ~0
        cls.assertAlmostEqual(maxl, 0., places = 1)
        
        # test that the marginalization recovers the right phase
        cls.assertAlmostEqual(maxl_phase, MAXL_PHI220, delta = 0.001)
        
        # test that the marginalized likelihood is strictly less than the maxL
        cls.assertTrue(marglogl < maxl)

    def test_margphase(self):
        # Test the marginalized phase model
        model, marglogl = self.load_model("gated_gaussian_margphase")
        self.test_marginalization(model, marglogl)
        
        # save likelihood value for comparisons
        self.margphase_margl = marglogl
    
    def test_multimargphase_amps(self):
        # Test the multimargphase model, sampling in amplitude
        model, marglogl = self.load_model('gated_gaussian_multimargphase_amps')
        self.test_marginalization(model, marglogl)
        
        # test that the likelihood is close to the marginalized phase model
        self.assertAlmostEqual(marglogl, self.margphase_margl,
                               delta = 0.01)
        
    def test_multimargphase(self):
        # Test the multimargphase model, sampling in SNR
        model, marglogl = self.load_model('gated_gaussian_multimargphase')
        self.test_marginalization(model, marglogl)
        
        # test that the likelihood is close to the marginalized phase model
        self.assertAlmostEqual(marglogl, self.margphase_margl,
                               delta = 0.01)
        
        # test that the loaded SNR reconstructs the correct amplitude
        scale_factor = model.current_stats['scale_factor_220']
        fid_amp = model.fiducial_amp_value
        self.assertAlmostEqual(scale_factor * fid_amp,
                               TEMPLATE_PARAMS['amp220'],
                               delta = 0.01)

suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(
    TestMargModels))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)