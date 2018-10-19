import unittest

from pycbc import strain
from pycbc.types import FrequencySeries
import numpy as np
from pycbc.workflow.configuration import WorkflowConfigParser
from pycbc.strain.recalibrate import Recalibrate
from utils import simple_exit

def strain_array():
    frequency_array = np.linspace(0, 2048, 100)
    delta_f = frequency_array[1] - frequency_array[0]
    strain_array = frequency_array**2
    return FrequencySeries(strain_array, delta_f)

class CalibrateTestBase(unittest.TestCase):
    def test_instantiate(self):
        self.assertRaises(TypeError, Recalibrate, 'test')

    def test_instantiation(self):
        params = dict(ifo_name='test',
                          minimum_frequency=10,
                          maximum_frequency=1024,
                          n_points=5)
        model = strain.models['cubic_spline'](**params)
        self.assertTrue(model.name == 'cubic_spline')
        
    def test_instantiation_from_config(self):
        parameters = dict(minimum_frequency='10', maximum_frequency='1024',
                          n_points='5')
    
        cp = WorkflowConfigParser()
        cp.add_section('test')
        ifo_name = 'ifo'
        cp.set('test', '{}_model'.format(ifo_name), 'cubic_spline')
        for key in parameters:
            cp.set('test', '{}_{}'.format(ifo_name, key), parameters[key])
        from_config = strain.read_model_from_config(cp, ifo_name, 'test')
        self.assertTrue(all([from_config.name == 'cubic_spline',
                             from_config.ifo_name == ifo_name]))

    def test_too_few_spline_points_fails(self):
        self.assertRaises(ValueError, strain.CubicSpline,
                         ifo_name='test', minimum_frequency=10,
                         maximum_frequency=1024, n_points=3)
   
    def test_update_parameters(self):
        init_params = dict(ifo_name='test', minimum_frequency=10,
          maximum_frequency=1024, n_points=5)
        dict_params = dict(recalib_amplitude_test_0=0.1,
          recalib_amplitude_test_1=0.1,
          recalib_amplitude_test_2=0.1,
          recalib_amplitude_test_3=0.1,
          recalib_amplitude_test_4=0.1,
          recalib_phase_test_0=0.1,
          recalib_phase_test_1=0.1,
          recalib_phase_test_2=0.1,
          recalib_phase_test_3=0.1,
          recalib_phase_test_4=0.1,
          )
        model = strain.models['cubic_spline'](**init_params)
        dict_params.update(dict(foo='bar'))
        prefix = 'recalib_'
        model.map_to_adjust(strain_array(), prefix=prefix, **dict_params)
        model_keys = [key[len(prefix):] for key in dict_params if prefix in key]
        self.assertTrue(all([model.params[key] == dict_params[prefix+key]
                    for key in model_keys]))
 

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(unittest.TestLoader().loadTestsFromTestCase(CalibrateTestBase))
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)

