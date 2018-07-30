import pytest

from pycbc import strain
from pycbc.types import FrequencySeries
import numpy as np
from pycbc.workflow.configuration import WorkflowConfigParser
from pycbc.strain.recalibrate import Recalibrate


@pytest.fixture
def strain_array():
    frequency_array = np.linspace(0, 2048, 100)
    delta_f = frequency_array[1] - frequency_array[0]
    strain_array = frequency_array**2
    return FrequencySeries(strain_array, delta_f)


def test_cannot_instantiate_base_class():
    with pytest.raises(TypeError):
        Recalibrate('test')


@pytest.mark.parametrize("name, params", [
    ('cubic_spline', dict(ifo_name='test',
                          minimum_frequency=10,
                          maximum_frequency=1024,
                          n_points=5)),
])
def test_instantiation(name, params):
    model = strain.models[name](**params)
    assert model.name == name


@pytest.mark.parametrize("name, init_params, dict_params", [
    ('cubic_spline',
     dict(ifo_name='test', minimum_frequency=10,
          maximum_frequency=1024, n_points=5),
     dict(recalib_amplitude_test_0=0.1,
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
     )
])
def test_update_parameters(name, init_params, dict_params):
    model = strain.models[name](**init_params)
    dict_params.update(dict(foo='bar'))
    prefix = 'recalib_'
    model.map_to_adjust(strain_array(), prefix=prefix, **dict_params)
    model_keys = [key[len(prefix):] for key in dict_params if prefix in key]
    assert all([model.params[key] == dict_params[prefix+key]
                for key in model_keys])


@pytest.mark.parametrize("name, parameters", [
    ('cubic_spline', dict(minimum_frequency='10', maximum_frequency='1024',
                          n_points='5'))
])
def test_instantiation_from_config(name, parameters):
    cp = WorkflowConfigParser()
    cp.add_section('test')
    ifo_name = 'ifo'
    cp.set('test', '{}_model'.format(ifo_name), name)
    for key in parameters:
        cp.set('test', '{}_{}'.format(ifo_name, key), parameters[key])
    from_config = strain.read_model_from_config(cp, ifo_name, 'test')
    assert all([from_config.name == name, from_config.ifo_name == ifo_name])


def test_too_few_spline_points_fails():
    with pytest.raises(ValueError):
        strain.CubicSpline(
            ifo_name='test', minimum_frequency=10,
            maximum_frequency=1024, n_points=3)
