# Copyright (C) 2013  Alex Nitz
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

#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
"""
These are simple unit tests for lalsimulation
"""
import unittest
import copy

import numpy

import lal, lalsimulation
from pycbc.filter import match, overlap, sigma, make_frequency_series
from pycbc.waveform import td_approximants, fd_approximants, \
        get_td_waveform, get_fd_waveform, TimeSeries

import optparse
from utils import simple_exit, _check_scheme_cpu

parser = optparse.OptionParser()
parser.add_option('--scheme','-s', action='callback', type = 'choice',
                   choices = ('cpu','cuda'),
                   default = 'cpu', dest = 'scheme', callback = _check_scheme_cpu,
                   help = optparse.SUPPRESS_HELP)
parser.add_option('--device-num','-d', action='store', type = 'int',
                   dest = 'devicenum', default=0,
                   help = optparse.SUPPRESS_HELP)

parser.add_option('--show-plots', action='store_true',
                   help = 'show the plots generated in this test suite')
parser.add_option('--save-plots', action='store_true',
                   help = 'save the plots generated in this test suite')

parser.add_option('--approximant', type = 'choice', choices = td_approximants() + fd_approximants(),
                  help = "Choices are %s" % str(td_approximants() + fd_approximants()))

parser.add_option('--mass1', type = float, default=10, help = "[default: %default]")
parser.add_option('--mass2', type = float, default=9, help = "[default: %default]")
parser.add_option('--spin1x', type = float, default=0, help = "[default: %default]")
parser.add_option('--spin1y', type = float, default=0, help = "[default: %default]")
parser.add_option('--spin1z', type = float, default=0, help = "[default: %default]")
parser.add_option('--spin2x', type = float, default=0, help = "[default: %default]")
parser.add_option('--spin2y', type = float, default=0, help = "[default: %default]")
parser.add_option('--spin2z', type = float, default=0, help = "[default: %default]")
parser.add_option('--lambda1', type = float, default=0, help = "[default: %default]")
parser.add_option('--lambda2', type = float, default=0, help = "[default: %default]")
parser.add_option('--coa-phase', type = float, default=0, help = "[default: %default]")
parser.add_option('--inclination', type = float, default=0, help = "[default: %default]")

parser.add_option('--delta-t', type = float, default=1.0/8192,  help = "[default: %default]")
parser.add_option('--delta-f', type = float, default=1.0/256,  help = "[default: %default]")
parser.add_option('--f-lower', type = float, default=30, help = "[default: %default]")

parser.add_option('--phase-order', type = int, default=-1, help = "[default: %default]")
parser.add_option('--amplitude-order', type = int, default=-1, help = "[default: %default]")
parser.add_option('--spin-order', type = int, default=-1, help = "[default: %default]")
parser.add_option('--tidal-order', type = int, default=-1, help = "[default: %default]")

(opt, args) = parser.parse_args()

print(72*'=')
print("Running {0} unit tests for {1}:".format('CPU', "Lalsimulation Waveforms"))

import matplotlib
if not opt.show_plots:
    matplotlib.use('Agg')
import pylab

def get_waveform(p, **kwds):
    """ Given the input parameters get me the waveform, whether it is TD or
        FD
    """
    params = copy.copy(p.__dict__)
    params.update(kwds)

    if params['approximant'] in td_approximants():
        return get_td_waveform(**params)
    else:
        return get_fd_waveform(**params)

class TestLALSimulation(unittest.TestCase):
    def setUp(self,*args):
        self.save_plots = opt.save_plots
        self.show_plots = opt.show_plots
        self.plot_dir = "."

        class params(object):
            pass

        self.p = params()

        # Overide my parameters with the program input arguments
        self.p.__dict__.update(vars(opt))

        if 'approximant' in self.kwds:
            self.p.approximant = self.kwds['approximant']

        from pycbc import version
        self.version_txt = "pycbc: %s  %s\n" % (version.git_hash, version.date) + \
                           "lalsimulation: %s  %s" % (lalsimulation.SimulationVCSIdentInfo.vcsId, lalsimulation.SimulationVCSIdentInfo.vcsDate)


    def test_varying_orbital_phase(self):
        #"""Check that the waveform is consistent under phase changes
        #"""

        if self.p.approximant in td_approximants():
            sample_attr = 'sample_times'
        else:
            sample_attr = 'sample_frequencies'

        f = pylab.figure()
        pylab.axes([.1, .2, 0.8, 0.70])
        hp_ref, hc_ref = get_waveform(self.p, coa_phase=0)
        pylab.plot(getattr(hp_ref, sample_attr), hp_ref.real(), label="phiref")

        hp, hc = get_waveform(self.p, coa_phase=lal.PI/4)
        m, i = match(hp_ref, hp)
        self.assertAlmostEqual(1, m, places=2)
        o = overlap(hp_ref, hp)
        pylab.plot(getattr(hp, sample_attr), hp.real(), label="$phiref \pi/4$")

        hp, hc = get_waveform(self.p, coa_phase=lal.PI/2)
        m, i = match(hp_ref, hp)
        o = overlap(hp_ref, hp)
        self.assertAlmostEqual(1, m, places=7)
        self.assertAlmostEqual(-1, o, places=7)
        pylab.plot(getattr(hp, sample_attr), hp.real(), label="$phiref \pi/2$")

        hp, hc = get_waveform(self.p, coa_phase=lal.PI)
        m, i = match(hp_ref, hp)
        o = overlap(hp_ref, hp)
        self.assertAlmostEqual(1, m, places=7)
        self.assertAlmostEqual(1, o, places=7)
        pylab.plot(getattr(hp, sample_attr), hp.real(), label="$phiref \pi$")

        pylab.xlim(min(getattr(hp, sample_attr)), max(getattr(hp, sample_attr)))
        pylab.title("Vary %s oribital phiref, h+" % self.p.approximant)

        if self.p.approximant in td_approximants():
            pylab.xlabel("Time to coalescence (s)")
        else:
            pylab.xlabel("GW Frequency (Hz)")

        pylab.ylabel("GW Strain (real part)")
        pylab.legend(loc="upper left")

        info = self.version_txt
        pylab.figtext(0.05, 0.05, info)

        if self.save_plots:
            pname = self.plot_dir + "/%s-vary-phase.png" % self.p.approximant
            pylab.savefig(pname)
        if self.show_plots:
            pylab.show()
        else:
            pylab.close(f)


    def test_distance_scaling(self):
        #""" Check that the waveform is consistent under distance changes
        #"""
        distance = 1e6
        tolerance = 1e-5
        fac = 10

        hpc, hcc = get_waveform(self.p, distance=distance)
        hpm, hcm = get_waveform(self.p, distance=distance*fac)
        hpf, hcf = get_waveform(self.p, distance=distance*fac*fac)
        hpn, hcn = get_waveform(self.p, distance=distance/fac)

        f = pylab.figure()
        pylab.axes([.1, .2, 0.8, 0.70])
        htilde = make_frequency_series(hpc)
        pylab.loglog(htilde.sample_frequencies, abs(htilde), label="D")

        htilde = make_frequency_series(hpm)
        pylab.loglog(htilde.sample_frequencies, abs(htilde), label="D * %s" %fac)

        htilde = make_frequency_series(hpf)
        pylab.loglog(htilde.sample_frequencies, abs(htilde), label="D * %s" %(fac*fac))

        htilde = make_frequency_series(hpn)
        pylab.loglog(htilde.sample_frequencies, abs(htilde), label="D / %s" %fac)

        pylab.title("Vary %s distance, $\\tilde{h}$+" % self.p.approximant)
        pylab.xlabel("GW Frequency (Hz)")
        pylab.ylabel("GW Strain")
        pylab.legend()
        pylab.xlim(xmin=self.p.f_lower)

        info = self.version_txt
        pylab.figtext(0.05, .05, info)

        if self.save_plots:
            pname = self.plot_dir + "/%s-distance-scaling.png" % self.p.approximant
            pylab.savefig(pname)

        if self.show_plots:
            pylab.show()
        else:
            pylab.close(f)

        self.assertTrue(hpc.almost_equal_elem(hpm * fac, tolerance, relative=True))
        self.assertTrue(hpc.almost_equal_elem(hpf * fac * fac, tolerance, relative=True))
        self.assertTrue(hpc.almost_equal_elem(hpn / fac, tolerance, relative=True))

    def test_nearby_waveform_agreement(self):
        #""" Check that the overlaps are consistent for nearby waveforms
        #"""
        def nearby(params):
            tol = 1e-7

            from numpy.random import uniform
            nearby_params = copy.copy(params)
            nearby_params.mass1 *= uniform(low=1-tol, high=1+tol)
            nearby_params.mass2 *= uniform(low=1-tol, high=1+tol)
            nearby_params.spin1x *= uniform(low=1-tol, high=1+tol)
            nearby_params.spin1y *= uniform(low=1-tol, high=1+tol)
            nearby_params.spin1z *= uniform(low=1-tol, high=1+tol)
            nearby_params.spin2x *= uniform(low=1-tol, high=1+tol)
            nearby_params.spin2y *= uniform(low=1-tol, high=1+tol)
            nearby_params.spin2z *= uniform(low=1-tol, high=1+tol)
            nearby_params.inclination *= uniform(low=1-tol, high=1+tol)
            nearby_params.coa_phase *= uniform(low=1-tol, high=1+tol)
            return nearby_params

        hp, hc = get_waveform(self.p)

        for i in range(10):
            p_near = nearby(self.p)
            hpn, hcn = get_waveform(p_near)

            maxlen = max(len(hpn), len(hp))
            hp.resize(maxlen)
            hpn.resize(maxlen)
            o = overlap(hp, hpn)
            self.assertAlmostEqual(1, o, places=5)

    def test_almost_equal_mass_waveform(self):
        #""" Check that the overlaps are consistent for nearby waveforms
        #"""
        def nearby(params):
            tol = 1e-7

            from numpy.random import uniform
            nearby_params = copy.copy(params)
            nearby_params.mass2 = nearby_params.mass1 * \
                uniform(low=1-tol, high=1+tol)
            nearby_params.mass1 *= uniform(low=1-tol, high=1+tol)
            nearby_params.spin1x *= uniform(low=1-tol, high=1+tol)
            nearby_params.spin1y *= uniform(low=1-tol, high=1+tol)
            nearby_params.spin1z *= uniform(low=1-tol, high=1+tol)
            nearby_params.spin2x *= uniform(low=1-tol, high=1+tol)
            nearby_params.spin2y *= uniform(low=1-tol, high=1+tol)
            nearby_params.spin2z *= uniform(low=1-tol, high=1+tol)
            nearby_params.inclination *= uniform(low=1-tol, high=1+tol)
            nearby_params.coa_phase *= uniform(low=1-tol, high=1+tol)
            return nearby_params

        for i in range(10):
            p_near = nearby(self.p)
            hpn, hcn = get_waveform(p_near)


    def test_varying_inclination(self):
        #""" Test that the waveform is consistent for changes in inclination
        #"""
        sigmas = []
        incs = numpy.arange(0, 21, 1.0) * lal.PI / 10.0

        for inc in incs:
            # WARNING: This does not properly handle the case of SpinTaylor*
            # where the spin orientation is not relative to the inclination
            hp, hc = get_waveform(self.p, inclination=inc)
            s = sigma(hp, low_frequency_cutoff=self.p.f_lower)
            sigmas.append(s)

        f = pylab.figure()
        pylab.axes([.1, .2, 0.8, 0.70])
        pylab.plot(incs, sigmas)
        pylab.title("Vary %s inclination, $\\tilde{h}$+" % self.p.approximant)
        pylab.xlabel("Inclination (radians)")
        pylab.ylabel("sigma (flat PSD)")

        info = self.version_txt
        pylab.figtext(0.05, 0.05, info)

        if self.save_plots:
            pname = self.plot_dir + "/%s-vary-inclination.png" % self.p.approximant
            pylab.savefig(pname)

        if self.show_plots:
            pylab.show()
        else:
            pylab.close(f)

        self.assertAlmostEqual(sigmas[-1], sigmas[0], places=7)
        self.assertAlmostEqual(max(sigmas), sigmas[0], places=7)
        self.assertTrue(sigmas[0] > sigmas[5])

    def test_swapping_constituents(self):
        #""" Test that waveform remains unchanged under swapping both objects
        #"""

        hp, hc = get_waveform(self.p)
        hpswap, hcswap = get_waveform(self.p, mass1=self.p.mass2, mass2=self.p.mass1,
                spin1x=self.p.spin2x, spin1y=self.p.spin2y, spin1z=self.p.spin2z,
                spin2x=self.p.spin1x, spin2y=self.p.spin1y, spin2z=self.p.spin1z,
                lambda1=self.p.lambda2, lambda2=self.p.lambda1)
        op = overlap(hp, hpswap)
        self.assertAlmostEqual(1, op, places=7)
        oc = overlap(hc, hcswap)
        self.assertAlmostEqual(1, oc, places=7)

    def test_change_rate(self):
        #""" Test that waveform remains unchanged under changing rate
        #"""
        hp, hc = get_waveform(self.p)
        hp2dec, hc2dec = get_waveform(self.p, delta_t=self.p.delta_t*2.)

        hpdec=numpy.zeros(len(hp2dec.data))
        hcdec=numpy.zeros(len(hp2dec.data))

        for idx in range(min(len(hp2dec.data),int(len(hp.data)/2))):
            hpdec[idx]=hp.data[2*idx]
            hcdec[idx]=hc.data[2*idx]

        hpTS=TimeSeries(hpdec, delta_t=self.p.delta_t*2.,epoch=hp.start_time)
        hcTS=TimeSeries(hcdec, delta_t=self.p.delta_t*2.,epoch=hc.start_time)

        f = pylab.figure()
        pylab.plot(hp.sample_times, hp.data,label="rate %s Hz" %"{:.0f}".format(1./self.p.delta_t))
        pylab.plot(hp2dec.sample_times, hp2dec.data, label="rate %s Hz" %"{:.0f}".format(1./(self.p.delta_t*2.)))

        pylab.title("Halving %s rate, $\\tilde{h}$+" % self.p.approximant)
        pylab.xlabel("time (sec)")
        pylab.ylabel("amplitude")
        pylab.legend()

        info = self.version_txt
        pylab.figtext(0.05, 0.05, info)

        if self.save_plots:
            pname = self.plot_dir + "/%s-vary-rate.png" % self.p.approximant
            pylab.savefig(pname)

        if self.show_plots:
            pylab.show()
        else:
            pylab.close(f)

        op=overlap(hpTS,hp2dec)
        self.assertAlmostEqual(1., op, places=2)
        oc=overlap(hcTS,hc2dec)
        self.assertAlmostEqual(1., oc, places=2)

def test_maker(class_name, name, **kwds):
    class Test(class_name):
        def __init__(self, *args):
            self.kwds = kwds
            class_name.__init__(self, *args)

    Test.__name__ = "Test %s" % name
    return Test

suite = unittest.TestSuite()

if opt.approximant:
    apxs = [opt.approximant]
else:
    apxs = td_approximants() + fd_approximants()

# These waveforms fail the current sanity checks, and are not used in current
# analyses. Tracking down reasons for each of these failures is a lot of work,
# so for now I just exclude these from tests.
fail_list = ['EOBNRv2', 'HGimri', 'SEOBNRv1', 'SpinDominatedWf',
             'PhenSpinTaylor', 'PhenSpinTaylorRD', 'EccentricTD',
             'EccentricFD', 'Lackey_Tidal_2013_SEOBNRv2_ROM']

for apx in apxs:
    # The inspiral wrapper is only single precision we won't bother checking
    # it here. It may need different tolerances and some special care.
    if apx.startswith("Inspiral-"):
        continue

    # The INTERP waveforms are designed only for filters
    if apx.endswith('_INTERP') and not opt.approximant:
        continue

    if apx in fail_list and not opt.approximant:
        # These waveforms segfault and prints debugging to screen
        # Only test this is specifically told to do so
        continue
    if apx in ['NR_hdf5']:
        # We'll need an example file for this. Also it will need a special
        # set of tests.
        continue

    vars()[apx] = test_maker(TestLALSimulation, apx, approximant=apx)
    suite.addTest( unittest.TestLoader().loadTestsFromTestCase(vars()[apx]) )

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
