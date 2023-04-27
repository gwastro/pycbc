import numpy
from pycbc import coordinates
from pycbc import distributions
from pycbc import conversions
import unittest
from utils import simple_exit

seed = 8202
numpy.random.seed(seed)

def almost_equal(derived_val, check_val, precision=1e-8):
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

def angle_almost_equal(derived_val, check_val, precision=1e-8):
    """Checks whether the given angles are almost equal. This is done by
    taking the modulus of each value on [0, 2*pi) before comparing.
    """
    derived_val = numpy.mod(derived_val, 2*numpy.pi)
    check_val = numpy.mod(check_val, 2*numpy.pi)
    return almost_equal(derived_val, check_val, precision=precision)

class TestParams(unittest.TestCase):
    def setUp(self, *args):
        self.numtests = 1000
        self.precision = 1e-8
        self.f_lower = 10.
        # create some component masses to work with
        self.m1 = numpy.random.uniform(1., 100., size=self.numtests)
        self.m2 = numpy.random.uniform(1., 100., size=self.numtests)
        # create some spins to work with
        spin_angledist = distributions.UniformSolidAngle()
        rvals = spin_angledist.rvs(size=self.numtests)
        self.spin1_polar = rvals['theta']
        self.spin1_az = rvals['phi']
        self.spin1_amp = numpy.random.uniform(0., 1., size=self.numtests)
        rvals = spin_angledist.rvs(size=self.numtests)
        self.spin2_polar = rvals['theta']
        self.spin2_az = rvals['phi']
        self.spin2_amp = numpy.random.uniform(0., 1., size=self.numtests)

        # calculate derived parameters from each
        self.mp = conversions.primary_mass(self.m1, self.m2)
        self.ms = conversions.secondary_mass(self.m1, self.m2)
        self.mtotal = conversions.mtotal_from_mass1_mass2(self.m1, self.m2)
        self.q = conversions.q_from_mass1_mass2(self.m1, self.m2)
        self.invq = conversions.invq_from_mass1_mass2(self.m1, self.m2)
        self.mchirp = conversions.mchirp_from_mass1_mass2(self.m1, self.m2)
        self.eta = conversions.eta_from_mass1_mass2(self.m1, self.m2)
        self.tau0 = conversions.tau0_from_mtotal_eta(self.mtotal, self.eta,
                                                     self.f_lower)
        self.tau3 = conversions.tau3_from_mtotal_eta(self.mtotal, self.eta,
                                                     self.f_lower)
        self.spin1x, self.spin1y, self.spin1z = \
            coordinates.spherical_to_cartesian(self.spin1_amp, self.spin1_az,
                                self.spin1_polar)
        self.spin2x, self.spin2y, self.spin2z = \
            coordinates.spherical_to_cartesian(self.spin2_amp, self.spin2_az,
                                self.spin2_polar)
        self.effective_spin = conversions.chi_eff(self.m1, self.m2,
                                self.spin1z, self.spin2z)
        self.chi_p = conversions.chi_p(self.m1, self.m2, self.spin1x,
            self.spin1y, self.spin2x, self.spin2y)
        self.primary_spinx = conversions.primary_spin(self.m1, self.m2,
                                self.spin1x, self.spin2x)
        self.primary_spiny = conversions.primary_spin(self.m1, self.m2,
                                self.spin1y, self.spin2y)
        self.primary_spinz = conversions.primary_spin(self.m1, self.m2,
                                self.spin1z, self.spin2z)
        self.secondary_spinx = conversions.secondary_spin(self.m1, self.m2,
                                    self.spin1x, self.spin2x)
        self.secondary_spiny = conversions.secondary_spin(self.m1, self.m2,
                                    self.spin1y, self.spin2y)
        self.secondary_spinz = conversions.secondary_spin(self.m1, self.m2,
                                    self.spin1z, self.spin2z)


    def test_physical_consistency(self):
        """Tests whether derived parameters pass physical checks; e.g., eta <=
        0.25.
        """
        self.assertTrue((self.mp >= self.ms).all(),
                        'primary mass not >= secondary mass')
        self.assertTrue((self.q >= 1.).all(),
                        'mass ratio not >= 1')
        self.assertTrue((self.invq <= 1.).all(),
                        'inverse mass ratio not <= 1')
        self.assertTrue((self.eta <= 0.25).all(),
                        'eta not <= 0.25')
        self.assertTrue((abs(self.effective_spin) <= 1.).all(),
                        'abs(effective spin) not <= 1')
        for which_comp in ['primary', 'secondary']:
            for coord in ['x', 'y', 'z']:
                spinparam = '{}_spin{}'.format(which_comp, coord)
                self.assertTrue((abs(getattr(self, spinparam)) <= 1.).all(),
                                '{} not <= 1.'.format(spinparam))

    def test_round_robin(self):
        """Computes inverse transformations to get original parameters from
        derived, then compares them to the original.
        """
        msg = '{} does not recover same {}; max difference: {}; inputs: {}'
        # following lists (function to check,
        #                  arguments to pass to the function,
        #                  name of self's attribute to compare to)
        fchecks = [
            (conversions.mass1_from_mtotal_q, (self.mtotal, self.q), 'mp'),
            (conversions.mass2_from_mtotal_q, (self.mtotal, self.q), 'ms'),
            (conversions.mass1_from_mtotal_eta, (self.mtotal, self.eta), 'mp'),
            (conversions.mass2_from_mtotal_eta, (self.mtotal, self.eta), 'ms'),
            (conversions.mtotal_from_mchirp_eta, (self.mchirp, self.eta), 'mtotal'),
            (conversions.mass1_from_mchirp_eta, (self.mchirp, self.eta), 'mp'),
            (conversions.mass2_from_mchirp_eta, (self.mchirp, self.eta), 'ms'),
            (conversions.mass2_from_mchirp_mass1, (self.mchirp, self.mp), 'ms'),
            (conversions.mass2_from_mass1_eta, (self.mp, self.eta), 'ms'),
            (conversions.mass1_from_mass2_eta, (self.ms, self.eta), 'mp'),
            (conversions.eta_from_q, (self.q,), 'eta'),
            (conversions.mass1_from_mchirp_q, (self.mchirp, self.q), 'mp'),
            (conversions.mass2_from_mchirp_q, (self.mchirp, self.q), 'ms'),
            (conversions.tau0_from_mchirp, (self.mchirp, self.f_lower), 'tau0'),
            (conversions.tau0_from_mass1_mass2, (self.m1, self.m2, self.f_lower), 'tau0'),
            (conversions.tau3_from_mass1_mass2, (self.m1, self.m2, self.f_lower), 'tau3'),
            (conversions.mchirp_from_tau0, (self.tau0, self.f_lower), 'mchirp'),
            (conversions.mtotal_from_tau0_tau3, (self.tau0, self.tau3, self.f_lower), 'mtotal'),
            (conversions.eta_from_tau0_tau3, (self.tau0, self.tau3, self.f_lower), 'eta'),
            (conversions.mass1_from_tau0_tau3, (self.tau0, self.tau3, self.f_lower), 'mp'),
            (conversions.mass2_from_tau0_tau3, (self.tau0, self.tau3, self.f_lower), 'ms'),
            (conversions.chi_eff_from_spherical,
                (self.m1, self.m2, self.spin1_amp, self.spin1_polar,
                 self.spin2_amp, self.spin2_polar), 'effective_spin'),
            (conversions.chi_p_from_spherical,
                (self.m1, self.m2, self.spin1_amp, self.spin1_az,
                 self.spin1_polar, self.spin2_amp, self.spin2_az,
                 self.spin2_polar), 'chi_p'),
            ]

        for func, args, compval in fchecks:
            passed, maxdiff, maxidx = almost_equal(func(*args), getattr(self, compval),
                                         self.precision)
            if not passed:
                failinputs = [p[maxidx] for p in args]
            else:
                failinputs = None
            self.assertTrue(passed, msg.format(func, compval, maxdiff, failinputs))

    def test_chip_compare_lalsuite(self):
        """Compares effective precession parameter bewteen
        the pycbc implementation and the lalsuite implementation.
        """
        import lal
        import lalsimulation as lalsim

        msg = '{} does not recover same {}; max difference: {}; inputs: {}'

        f_ref = self.f_lower
        chip_lal = []
        for i in range(len(self.m1)):
            _,_,tmp,_,_,_,_ = lalsim.SimIMRPhenomPCalculateModelParametersFromSourceFrame(
                self.m1[i]*lal.MSUN_SI, self.m2[i]*lal.MSUN_SI,
                f_ref, 0., 0.,
                self.spin1x[i], self.spin1y[i], self.spin1z[i],
                self.spin2x[i], self.spin2y[i], self.spin2z[i], lalsim.IMRPhenomPv2_V)
            chip_lal.append(tmp)

        chip_pycbc = conversions.chi_p(
            self.m1,self.m2,self.spin1x,self.spin1y,self.spin2x,self.spin2y)

        passed, maxdiff, maxidx = almost_equal(chip_lal, chip_pycbc, self.precision)
        failinputs = (
            self.m1[maxidx], self.m2[maxidx],
            self.spin1x[maxidx], self.spin1y[maxidx],
            self.spin2x[maxidx],self.spin2y[maxidx]
            )
        self.assertTrue(passed, msg.format("conversions.chi_p", "chi_p", maxdiff, failinputs))

suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestParams))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
