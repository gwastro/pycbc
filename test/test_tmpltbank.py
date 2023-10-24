# Copyright (C) 2013 Ian Harry
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
These are the unittests for the pycbc.tmpltbank module
"""

import math
import gzip
import numpy
from astropy.utils.data import download_file
from astropy.utils.data import conf
import pycbc.tmpltbank
# Old LigoLW output functions are not imported at tmpltbank level
import pycbc.tmpltbank.bank_output_utils as llw_output
import pycbc.psd
import pycbc.pnutils
from pycbc import pnutils
from pycbc.types import Array
from pycbc.filter import match
from pycbc.waveform import get_fd_waveform
import matplotlib
matplotlib.use('Agg')

import unittest
from utils import parse_args_cpu_only, simple_exit

# This will return whatever is appropriate, depending on whether this
# particular instance of the unittest was called for CPU, CUDA, or OpenCL
parse_args_cpu_only("Template bank module")

import argparse
parser = argparse.ArgumentParser()

DATA_FILE_URL = 'https://github.com/gwastro/pycbc-config/raw/master/test_data_files/{}'
# Allow astropy more time before downloads timeout
conf.remote_timeout = 100

def update_mass_parameters(tmpltbank_class):
    """
    Choose various sets of mass parameters for testing.
    """
    num_comp_masses = 3
    min_mass1 = [1,2,6]
    max_mass1 = [5,8,12]
    min_mass2 = [1,1,1]
    max_mass2 = [5,5,5]
    num_tot_masses = 3
    # These *must* be provided
    min_tot_mass = [None, 2.5, 3.5]
    max_tot_mass = [None, 11, 7.5]
    num_chirp_masses = 3
    max_chirp_mass = [None, 2.43, 3.5]
    min_chirp_mass = [None, 1.218, 2.43]
    num_etas = 3
    max_eta = [0.25, 0.24, 0.23]
    min_eta = [None, 0.16, 0.17]

    max_iter_idx = num_comp_masses * num_tot_masses *\
                   num_chirp_masses * num_etas

    for idx in range(max_iter_idx):
        comp_masses_idx = idx % num_comp_masses
        tmpltbank_class.min_mass1 = min_mass1[comp_masses_idx]
        tmpltbank_class.max_mass1 = max_mass1[comp_masses_idx]
        tmpltbank_class.min_mass2 = min_mass2[comp_masses_idx]
        tmpltbank_class.max_mass2 = max_mass2[comp_masses_idx]
        reduced_idx = idx // num_comp_masses
        tot_mass_idx = reduced_idx % num_tot_masses
        tmpltbank_class.min_total_mass = min_tot_mass[tot_mass_idx]
        tmpltbank_class.max_total_mass = max_tot_mass[tot_mass_idx]
        reduced_idx = reduced_idx // num_tot_masses
        chirp_mass_idx = reduced_idx % num_chirp_masses
        tmpltbank_class.min_chirp_mass = min_chirp_mass[chirp_mass_idx]
        tmpltbank_class.max_chirp_mass = max_chirp_mass[chirp_mass_idx]
        reduced_idx = reduced_idx // num_chirp_masses
        eta_idx = reduced_idx
        tmpltbank_class.max_eta = max_eta[eta_idx]
        tmpltbank_class.min_eta = min_eta[eta_idx]
        yield idx

    return

class TmpltbankTestClass(unittest.TestCase):
    def setUp(self):

        self.deltaF = 0.1
        self.f_low = 15
        self.f_upper = 2000
        self.f0 = 70
        self.sampleRate = 4096
        self.pnOrder = 'threePointFivePN'
        self.min_mass1 = 1
        self.min_mass2 = 1
        self.max_mass1 = 5
        self.max_mass2 = 5
        self.max_ns_spin_mag = 0.5
        self.max_bh_spin_mag = 0.9
        self.ns_bh_boundary_mass = 2.0
        self.min_total_mass = 2.5
        self.max_total_mass = 6.0
        self.max_chirp_mass = 2.4375415772291475
        self.min_chirp_mass = 1.2187707886145738
        self.max_eta = 0.24
        self.min_eta = 0.16

        # Sanity check these
        pycbc.tmpltbank.verify_mass_range_options(self, parser=parser)

        # Need to use F2 metric for ethinca
        self.ethincaOrder = 'threePointFivePN'
        self.ethincaCutoff = 'SchwarzISCO'
        self.ethincaFreqStep = 10.

        self.segLen = 1./self.deltaF
        self.psdSize = int(self.segLen * self.sampleRate / 2.) + 1

        apy_fname = download_file(
            DATA_FILE_URL.format('ZERO_DET_high_P.txt.gz'),
            cache=True
        )
        match_psd_size = int(256 * self.sampleRate / 2.) + 1
        with gzip.open(apy_fname) as apy_fp:
            self.psd = pycbc.psd.from_txt(apy_fp, self.psdSize, self.deltaF,
                                          self.f_low, is_asd_file=True)
            apy_fp.seek(0)  # Reset file-pointer to start of file
            self.psd_for_match = pycbc.psd.from_txt(apy_fp, match_psd_size,
                                                    1./256., self.f_low,
                                                    is_asd_file=True)

        metricParams = pycbc.tmpltbank.metricParameters(self.pnOrder,\
                         self.f_low, self.f_upper, self.deltaF, self.f0)
        metricParams.psd = self.psd

        massRangeParams = pycbc.tmpltbank.massRangeParameters(self.min_mass1,\
                            self.max_mass1, self.min_mass2, self.max_mass2,\
                            maxNSSpinMag=self.max_ns_spin_mag,\
                            maxBHSpinMag=self.max_bh_spin_mag,\
                            maxTotMass=self.max_total_mass,\
                            minTotMass=self.min_total_mass,\
                            max_chirp_mass=self.max_chirp_mass,\
                            min_chirp_mass=self.min_chirp_mass,\
                            maxEta=self.max_eta,\
                            minEta=self.min_eta,\
                            ns_bh_boundary_mass=self.ns_bh_boundary_mass)

        # And again with the nsbh flag
        massRangeParams2 = pycbc.tmpltbank.massRangeParameters(self.min_mass1,\
                            self.max_mass1, self.min_mass2, self.max_mass2,\
                            maxNSSpinMag=self.max_ns_spin_mag,\
                            maxBHSpinMag=self.max_bh_spin_mag,\
                            maxTotMass=self.max_total_mass,\
                            minTotMass=self.min_total_mass,\
                            max_chirp_mass=self.max_chirp_mass,\
                            min_chirp_mass=self.min_chirp_mass,\
                            maxEta=self.max_eta,\
                            minEta=self.min_eta,\
                            nsbhFlag=True)

        metricParams = pycbc.tmpltbank.determine_eigen_directions(metricParams)

        vals=pycbc.tmpltbank.estimate_mass_range(100000, massRangeParams,\
               metricParams, self.f_upper, covary=False)

        cov = numpy.cov(vals)
        _,self.evecsCV = numpy.linalg.eig(cov)
        metricParams.evecsCV = {}
        metricParams.evecsCV[self.f_upper] = self.evecsCV

        vals=pycbc.tmpltbank.estimate_mass_range(100000, massRangeParams,\
               metricParams, self.f_upper, covary=False)

        self.metricParams = metricParams
        self.massRangeParams = massRangeParams
        self.massRangeParams2 = massRangeParams2
        self.ethincaParams = pycbc.tmpltbank.ethincaParameters(
            self.ethincaOrder, self.ethincaCutoff, self.ethincaFreqStep,
            full_ethinca=False, time_ethinca=False)

        self.xis = vals

    def test_eigen_directions(self):
        fname='stockEvals.dat.gz'
        apy_fname = download_file(DATA_FILE_URL.format(fname), cache=False)
        with gzip.open(apy_fname) as apy_fp:
            evalsStock = Array(numpy.loadtxt(apy_fp))

        fname='stockEvecs.dat.gz'
        apy_fname = download_file(DATA_FILE_URL.format(fname), cache=False)
        with gzip.open(apy_fname) as apy_fp:
            evecsStock = Array(numpy.loadtxt(apy_fp))

        maxEval = max(evalsStock)
        evalsCurr = Array(self.metricParams.evals[self.f_upper])
        evecsCurr = Array(self.metricParams.evecs[self.f_upper])
        # Uncomment these lines to regenerate the data files
        #numpy.savetxt('newEvals.dat', evalsCurr)
        #numpy.savetxt('newEvecs.dat', evecsCurr)
        errMsg = "pycbc.tmpltbank.determine_eigen_directions has failed "
        errMsg += "sanity check."
        evalsDiff = abs(evalsCurr - evalsStock)/maxEval
        self.assertTrue(not (evalsDiff > 1E-5).any(), msg=errMsg)
        for stock,test in zip(evecsStock.data,evecsCurr.data):
            stockScaled = stock * evalsCurr.data**0.5
            testScaled = test * evalsCurr.data**0.5
            diff = stockScaled - testScaled
            self.assertTrue(not (diff > 1E-4).any(), msg=errMsg)

    def test_get_random_mass(self):
        # Want to do this for a variety of mass combinations
        for i in update_mass_parameters(self):
            curr_min_mass = self.min_total_mass
            curr_max_mass = self.max_total_mass
            try:
                pycbc.tmpltbank.verify_mass_range_options(self, parser=parser)
            except ValueError:
                # Some of the inputs are unphysical and will fail.
                # These cases are known to fail, the inputs are unphysical
                # 35 has inconsistent total mass and eta restrictions
                # 38 Component mass, [upper] chirp mass and [lower] eta limits
                #    rule out the entire space.
                # 41 Same as 38
                # 44 Same as 38
                # 62 From component mass and total mass limits only total masses
                #    between 7 and 7.5 are possible. This range all has eta
                #    lower than the limit of 0.17.
                # 65 Same as 38
                # 68 Same as 38
                # 71 Same as 38
                # 80 Same as 62
                if i in [35,38,41,44,62,65,68,71,80]:
                    continue
                raise

            # Check that if the mass limits have changed, it was right to do so
            # This is not exhaustive, but gets most things
            if (curr_min_mass is None) or \
                    not math.isclose(self.min_total_mass, curr_min_mass,
                                     rel_tol=1e-06):
                min_comp_mass = self.min_mass1 + self.min_mass2
                min_eta = self.min_mass1 * self.min_mass2 /\
                           (min_comp_mass * min_comp_mass)
                min_chirp_mass = min_comp_mass * min_eta**(3./5.)
                if (min_comp_mass is not None) and \
                        math.isclose(self.min_total_mass, min_comp_mass,
                                     rel_tol=1e-06):
                    # Okay, the total mass is changed by the components
                    pass
                elif (self.min_eta and min_eta < self.min_eta) or \
                        (self.max_eta and min_eta > self.max_eta):
                    # Okay, not possible from eta
                    pass
                elif self.min_chirp_mass and \
                        min_chirp_mass < self.min_chirp_mass:
                    # Okay, not possible from chirp mass
                    pass
                else:
                    err_msg = "Minimum total mass changed unexpectedly."
                    self.fail(err_msg)
            if (curr_max_mass is None) or \
                    not math.isclose(self.max_total_mass, curr_max_mass,
                                 rel_tol=1e-06):
                max_comp_mass = self.max_mass1 + self.max_mass2
                max_eta = self.max_mass1 * self.max_mass2 /\
                           (max_comp_mass * max_comp_mass)
                max_chirp_mass = max_comp_mass * max_eta**(3./5.)
                if (max_comp_mass is not None) and \
                        math.isclose(self.max_total_mass, max_comp_mass,
                                     rel_tol=1e-06):
                    # Okay, the total mass is changed by the components
                    pass
                elif (self.min_eta and max_eta < self.min_eta) or\
                        (self.max_eta and max_eta > self.max_eta):
                    # Okay, not possible from eta
                    pass
                elif self.max_chirp_mass and \
                        max_chirp_mass > self.max_chirp_mass:
                    # Okay, not possible from chirp mass
                    pass
                else:
                    err_msg = "Maximum total mass changed unexpectedly."
                    self.fail(err_msg)

            massRangeParams = pycbc.tmpltbank.massRangeParameters(\
                                self.min_mass1,\
                                self.max_mass1, self.min_mass2, self.max_mass2,\
                                maxNSSpinMag=self.max_ns_spin_mag,\
                                maxBHSpinMag=self.max_bh_spin_mag,\
                                maxTotMass=self.max_total_mass,\
                                minTotMass=self.min_total_mass,\
                                max_chirp_mass=self.max_chirp_mass,\
                                min_chirp_mass=self.min_chirp_mass,\
                                maxEta=self.max_eta,\
                                minEta=self.min_eta,\
                                ns_bh_boundary_mass=self.ns_bh_boundary_mass)

            # And again with the nsbh flag
            massRangeParams2 = pycbc.tmpltbank.massRangeParameters(\
                                self.min_mass1,\
                                self.max_mass1, self.min_mass2, self.max_mass2,\
                                maxNSSpinMag=self.max_ns_spin_mag,\
                                maxBHSpinMag=self.max_bh_spin_mag,\
                                maxTotMass=self.max_total_mass,\
                                minTotMass=self.min_total_mass,\
                                max_chirp_mass=self.max_chirp_mass,\
                                min_chirp_mass=self.min_chirp_mass,\
                                maxEta=self.max_eta,\
                                minEta=self.min_eta,\
                                nsbhFlag=True)

            mass1, mass2, spin1z, spin2z = \
                 pycbc.tmpltbank.get_random_mass(100000, massRangeParams)
            mass = mass1 + mass2
            errMsg = "pycbc.tmpltbank.get_random_mass returns invalid ranges."
            self.assertTrue(not (mass < self.min_total_mass).any(),msg=errMsg)
            self.assertTrue(not (mass > self.max_total_mass).any(),msg=errMsg)
            self.assertTrue(not (mass1 > self.max_mass1 * 1.001).any(),
                            msg=errMsg)
            self.assertTrue(not (mass1 < self.min_mass1 * 0.999).any(),
                            msg=errMsg)
            self.assertTrue(not (mass2 > self.max_mass2 * 1.001).any(),
                            msg=errMsg)
            self.assertTrue(not (mass2 < self.min_mass2 * 0.999).any(),
                            msg=errMsg)
            self.assertTrue(not (mass1 < mass2).any(),msg=errMsg)
            # Chirp mass and eta
            mchirp, eta = pnutils.mass1_mass2_to_mchirp_eta(mass1,mass2)
            if self.max_chirp_mass:
                self.assertTrue(not (mchirp > self.max_chirp_mass*1.0001).any(),
                                msg=errMsg)
            if self.min_chirp_mass:
                self.assertTrue(not (mchirp < self.min_chirp_mass*0.9999).any(),
                                msg=errMsg)
            if self.min_eta:
                self.assertTrue(not (eta < self.min_eta*0.9999).any(),
                                msg=errMsg)
                self.assertTrue(not (eta > self.max_eta*1.0001).any(),
                                msg=errMsg)
            nsSpin1 = spin1z[mass1 < self.ns_bh_boundary_mass]
            nsSpin2 = spin2z[mass2 < self.ns_bh_boundary_mass]
            bhSpin1 = spin1z[mass1 > self.ns_bh_boundary_mass]
            bhSpin2 = spin2z[mass2 > self.ns_bh_boundary_mass]
            self.assertTrue(not (abs(nsSpin1) > 0.5).any(), msg=errMsg)
            self.assertTrue(not (abs(nsSpin2) > 0.5).any(), msg=errMsg)
            self.assertTrue(not (abs(bhSpin1) > 0.9).any(), msg=errMsg)
            self.assertTrue(not (abs(bhSpin2) > 0.9).any(), msg=errMsg)
            # Check that *some* spins are bigger than 0.5
            if len(bhSpin1):
                self.assertTrue((abs(bhSpin1) > 0.5).any(), msg=errMsg)
            if len(bhSpin2):
                self.assertTrue((abs(bhSpin2) > 0.5).any(), msg=errMsg)

            # Check nsbh flag
            mass1, mass2, spin1z, spin2z = \
                 pycbc.tmpltbank.get_random_mass(100000, massRangeParams2)
            self.assertTrue(not (abs(spin1z) > 0.9).any(), msg=errMsg)
            self.assertTrue(not (abs(spin2z) > 0.5).any(), msg=errMsg)
            self.assertTrue((abs(spin1z) > 0.5).any(), msg=errMsg)

    def test_metric_match_prediction(self):
        mass1a, mass2a, spin1za, spin2za = \
                 pycbc.tmpltbank.get_random_mass(10, self.massRangeParams)
        mass1b, mass2b, spin1zb, spin2zb = \
                 pycbc.tmpltbank.get_random_mass(10, self.massRangeParams)
        for idx in range(10):
            masses1 = [mass1a[idx], mass2a[idx], spin1za[idx], spin2za[idx]]
            masses2 = [mass1b[idx], mass2b[idx], spin1zb[idx], spin2zb[idx]]
            dist, _, _ = pycbc.tmpltbank.get_point_distance \
                (masses1,  masses2, self.metricParams, self.f_upper)
            opt_dist = 0.02
            while dist > opt_dist * 1.01  or dist < opt_dist * 0.99:
                dist_fac = opt_dist / dist
                dist_fac = dist_fac**0.5
                if dist_fac < 0.01:
                    dist_fac = 0.01
                if dist_fac > 2:
                    dist_fac = 2
                for idx, curr_mass2 in enumerate(masses2):
                    masses2[idx] = masses1[idx] + \
                        (curr_mass2 - masses1[idx]) * dist_fac
                dist, _, _ = pycbc.tmpltbank.get_point_distance \
                    (masses1,  masses2, self.metricParams, self.f_upper)
            self.assertFalse(numpy.isnan(dist))

            htilde1, _ = get_fd_waveform\
                (approximant='TaylorF2', mass1=masses1[0], mass2=masses1[1],
                 spin1z=masses1[2], spin2z=masses1[3], delta_f=1.0/256,
                 f_lower=15, f_final=2000)
            htilde2, _ = get_fd_waveform\
                (approximant='TaylorF2', mass1=masses2[0], mass2=masses2[1],
                 spin1z=masses2[2], spin2z=masses2[3], delta_f=1.0/256,
                 f_lower=15, f_final=2000)
            overlap, _ = match(htilde1, htilde2, psd=self.psd_for_match,
                            low_frequency_cutoff=15)
            self.assertTrue(overlap > 0.97 and overlap < 0.985)


    def test_chirp_params(self):
        chirps=pycbc.tmpltbank.get_chirp_params(2.2, 1.8, 0.2, 0.3,
                              self.metricParams.f0, self.metricParams.pnOrder)
        fname = 'stockChirps.dat.gz'
        apy_fname = download_file(DATA_FILE_URL.format(fname), cache=False)
        with gzip.open(apy_fname) as apy_fp:
            stockChirps = numpy.loadtxt(apy_fp)
        diff = (chirps - stockChirps) / stockChirps
        errMsg = "Calculated chirp params differ from that expected."
        self.assertTrue( not (abs(diff) > 1E-4).any(), msg=errMsg)

    def test_hexagonal_placement(self):
        arrz = pycbc.tmpltbank.generate_hexagonal_lattice(10, 0, 10, 0, 0.03)
        arrz = numpy.array(arrz)
        fname = 'stockHexagonal.dat.gz'
        apy_fname = download_file(DATA_FILE_URL.format(fname), cache=False)
        with gzip.open(apy_fname) as apy_fp:
            stockGrid = numpy.loadtxt(apy_fp)
        diff = arrz - stockGrid
        errMsg = "Calculated lattice differs from that expected."
        self.assertTrue( not (diff > 1E-4).any(), msg=errMsg)

    def test_anstar_placement(self):
        arrz = pycbc.tmpltbank.generate_anstar_3d_lattice(0, 10, 0, 10, 0, \
                                                          10, 0.03)
        arrz = numpy.array(arrz)
        fname = 'stockAnstar3D.dat.gz'
        apy_fname = download_file(DATA_FILE_URL.format(fname), cache=False)
        with gzip.open(apy_fname) as apy_fp:
            stockGrid = numpy.loadtxt(apy_fp)
        # Uncomment this line to regenerate the data file
        #numpy.savetxt("new_example.dat", arrz)
        errMsg = "Calculated lattice differs from that expected."
        self.assertTrue(len(arrz) == len(stockGrid), msg=errMsg)
        diff = arrz - stockGrid
        self.assertTrue( not (diff > 1E-4).any(), msg=errMsg)

    def test_get_mass_distribution(self):
        # Just run the function, no checking output
        pycbc.tmpltbank.get_mass_distribution([1.35,0.239,0.4,-0.2], 2, \
                          self.massRangeParams, self.metricParams, \
                          self.f_upper, \
                          numJumpPoints=123, chirpMassJumpFac=0.0002, \
                          etaJumpFac=0.009, spin1zJumpFac=0.1, \
                          spin2zJumpFac=0.2)

    def test_get_phys_cov_masses(self):
        evecs = self.metricParams.evecs[self.f_upper]
        evals = self.metricParams.evals[self.f_upper]
        masses1 = [2.2,1.8,0.4,0.3]
        masses2 = [2.21,1.79,0.41,0.29]
        xis1 = pycbc.tmpltbank.get_cov_params(masses1[0], masses1[1],
                 masses1[2], masses1[3], self.metricParams, self.f_upper)
        xis2 = pycbc.tmpltbank.get_cov_params(masses2[0], masses2[1],
                 masses2[2], masses2[3], self.metricParams, self.f_upper)

        testXis = [xis1[0],xis1[1]]
        b_mtot, b_eta = pnutils.mass1_mass2_to_mtotal_eta(masses2[0],
                                                          masses2[1])
        bestMasses = [b_mtot, b_eta, masses2[2], masses2[3]]
        bestXis = xis2
        output = pycbc.tmpltbank.get_physical_covaried_masses(testXis, \
                   bestMasses, bestXis, 0.0001, self.massRangeParams, \
                   self.metricParams, self.f_upper)
        # Test that returned xis are close enough
        diff = (output[6][0] - testXis[0])**2
        diff += (output[6][1] - testXis[1])**2
        errMsg = 'pycbc.tmpltbank.get_physical_covaried_masses '
        errMsg += 'failed to find a point within the desired limits.'
        self.assertTrue( diff < 1E-4,msg=errMsg)
        # Test that returned masses and xis agree
        massT = output[0] + output[1]
        etaT = output[0]*output[1] / (massT*massT)
        spinSetT = pycbc.pnutils.get_beta_sigma_from_aligned_spins(\
                     etaT, output[2], output[3])
        xisT = pycbc.tmpltbank.get_cov_params(output[0], output[1],
                 output[2], output[3], self.metricParams, self.f_upper)
        errMsg = "Recovered xis do not agree with those expected."
        self.assertTrue( abs(xisT[0] - output[6][0]) < 1E-5, msg=errMsg)
        self.assertTrue( abs(xisT[1] - output[6][1]) < 1E-5, msg=errMsg)
        self.assertTrue( abs(xisT[2] - output[6][2]) < 1E-5, msg=errMsg)
        self.assertTrue( abs(xisT[3] - output[6][3]) < 1E-5, msg=errMsg)

        # Test again with nsbh flag on
        output = pycbc.tmpltbank.get_physical_covaried_masses(testXis, \
                   bestMasses, bestXis, 0.0001, self.massRangeParams2, \
                   self.metricParams, self.f_upper)
        # Test that returned xis are close enough
        diff = (output[6][0] - testXis[0])**2
        diff += (output[6][1] - testXis[1])**2
        errMsg = 'pycbc.tmpltbank.get_physical_covaried_masses '
        errMsg += 'failed to find a point within the desired limits.'
        self.assertTrue( diff < 1E-4,msg=errMsg)
        # Test that returned masses and xis agree
        xisT = pycbc.tmpltbank.get_cov_params(output[0], output[1],
                 output[2], output[3], self.metricParams, self.f_upper)
        errMsg = "Recovered xis do not agree with those expected."
        self.assertTrue( abs(xisT[0] - output[6][0]) < 1E-5, msg=errMsg)
        self.assertTrue( abs(xisT[1] - output[6][1]) < 1E-5, msg=errMsg)
        self.assertTrue( abs(xisT[2] - output[6][2]) < 1E-5, msg=errMsg)
        self.assertTrue( abs(xisT[3] - output[6][3]) < 1E-5, msg=errMsg)

    def test_stack_xi_direction(self):
        # Just run the function, no checking output
        evecs = self.metricParams.evecs[self.f_upper]
        evals = self.metricParams.evals[self.f_upper]
        masses1 = [2.2,1.8,0.4,0.3]
        masses2 = [2.21,1.79,0.41,0.29]
        xis1 = pycbc.tmpltbank.get_cov_params(masses1[0], masses1[1], \
                 masses1[2], masses1[3], self.metricParams, self.f_upper)
        xis2 = pycbc.tmpltbank.get_cov_params(masses2[0], masses2[1], \
                 masses2[2], masses2[3], self.metricParams, self.f_upper)
        testXis = [xis1[0],xis1[1]]
        b_mtot, b_eta = pnutils.mass1_mass2_to_mtotal_eta(masses2[0],
                                                          masses2[1])
        bestMasses = [b_mtot, b_eta, masses2[2], masses2[3]]
        bestXis = xis2

        depths = pycbc.tmpltbank.stack_xi_direction_brute(testXis, \
              bestMasses, bestXis, 3, 0.03, self.massRangeParams, \
              self.metricParams, self.f_upper, numIterations=50)

    def test_point_distance(self):
        masses1 = [2,2,0.4,0.6]
        masses2 = [2.02,1.97,0.41,0.59]
        dist, xis1, xis2 = pycbc.tmpltbank.get_point_distance(masses1, \
                             masses2, self.metricParams, self.f_upper)
        diff = abs((dist - 23.3681922039) / dist)

        errMsg = "Obtained distance does not agree with expected value."
        self.assertTrue( diff < 1E-5, msg=errMsg)

    def test_conv_to_sngl(self):
        # Just run the function, no checking output
        masses1 = [(2,2,0.4,0.3),(4.01,0.249,0.41,0.29)]
        llw_output.convert_to_sngl_inspiral_table(masses1, "a")

    def test_ethinca_calc(self):
        # Just run the function, no checking output
        m1 = 2.
        m2 = 2.
        s1z = 0.
        s2z = 0.
        # ethinca calc breaks unless f0 = fLow
        self.metricParams.f0 = self.metricParams.fLow
        output = llw_output.calculate_ethinca_metric_comps(
            self.metricParams, self.ethincaParams, m1, m2, s1z, s2z)
        # restore initial f0 value
        self.metricParams.f0 = self.f0

    def tearDown(self):
        pass

suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TmpltbankTestClass))

if  __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)

