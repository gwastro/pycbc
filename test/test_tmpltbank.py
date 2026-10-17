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

import gzip
import math

import matplotlib
import numpy
from astropy.utils.data import conf

import pycbc.pnutils
import pycbc.psd
import pycbc.tmpltbank

# Old LigoLW output functions are not imported at tmpltbank level
import pycbc.tmpltbank.bank_output_utils as llw_output
from pycbc import pnutils
from pycbc.filter import match
from pycbc.io import get_file
from pycbc.types import Array
from pycbc.waveform import get_fd_waveform

matplotlib.use("Agg")

import unittest

from utils import parse_args_cpu_only, simple_exit

# This will return whatever is appropriate, depending on whether this
# particular instance of the unittest was called for CPU, CUDA, or OpenCL
parse_args_cpu_only("Template bank module")

import argparse

parser = argparse.ArgumentParser()

DATA_FILE_URL = "https://github.com/gwastro/pycbc-config/raw/master/test_data_files/{}"
# Allow astropy more time before downloads timeout
conf.remote_timeout = 100


def update_mass_parameters(tmpltbank_class):
    """
    Choose various sets of mass parameters for testing.
    """
    num_comp_masses = 3
    min_mass1 = [1, 2, 6]
    max_mass1 = [5, 8, 12]
    min_mass2 = [1, 1, 1]
    max_mass2 = [5, 5, 5]
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

    max_iter_idx = num_comp_masses * num_tot_masses * num_chirp_masses * num_etas

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


class TmpltbankTestClass(unittest.TestCase):
    def setUp(self):

        self.delta_f = 0.1
        self.f_low = 15
        self.f_upper = 2000
        self.f0 = 70
        self.sample_rate = 4096
        self.pn_order = "threePointFivePN"
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
        self.ethinca_order = "threePointFivePN"
        self.ethinca_cutoff = "SchwarzISCO"
        self.ethinca_freq_step = 10.0

        self.seg_len = 1.0 / self.delta_f
        self.psd_size = int(self.seg_len * self.sample_rate / 2.0) + 1

        apy_fname = get_file(DATA_FILE_URL.format("ZERO_DET_high_P.txt.gz"), cache=True)
        match_psd_size = int(256 * self.sample_rate / 2.0) + 1
        with gzip.open(apy_fname) as apy_fp:
            self.psd = pycbc.psd.from_txt(
                apy_fp, self.psd_size, self.delta_f, self.f_low, is_asd_file=True
            )
            apy_fp.seek(0)  # Reset file-pointer to start of file
            self.psd_for_match = pycbc.psd.from_txt(
                apy_fp, match_psd_size, 1.0 / 256.0, self.f_low, is_asd_file=True
            )

        metric_params = pycbc.tmpltbank.MetricParameters(
            self.pn_order, self.f_low, self.f_upper, self.delta_f, self.f0
        )
        metric_params.psd = self.psd

        mass_range_params = pycbc.tmpltbank.MassRangeParameters(
            self.min_mass1,
            self.max_mass1,
            self.min_mass2,
            self.max_mass2,
            max_ns_spin_mag=self.max_ns_spin_mag,
            max_bh_spin_mag=self.max_bh_spin_mag,
            max_tot_mass=self.max_total_mass,
            min_tot_mass=self.min_total_mass,
            max_chirp_mass=self.max_chirp_mass,
            min_chirp_mass=self.min_chirp_mass,
            max_eta=self.max_eta,
            min_eta=self.min_eta,
            ns_bh_boundary_mass=self.ns_bh_boundary_mass,
        )

        # And again with the nsbh flag
        mass_range_params2 = pycbc.tmpltbank.MassRangeParameters(
            self.min_mass1,
            self.max_mass1,
            self.min_mass2,
            self.max_mass2,
            max_ns_spin_mag=self.max_ns_spin_mag,
            max_bh_spin_mag=self.max_bh_spin_mag,
            max_tot_mass=self.max_total_mass,
            min_tot_mass=self.min_total_mass,
            max_chirp_mass=self.max_chirp_mass,
            min_chirp_mass=self.min_chirp_mass,
            max_eta=self.max_eta,
            min_eta=self.min_eta,
            nsbh_flag=True,
        )

        metric_params = pycbc.tmpltbank.determine_eigen_directions(metric_params)

        vals = pycbc.tmpltbank.estimate_mass_range(
            100000, mass_range_params, metric_params, self.f_upper, covary=False
        )

        cov = numpy.cov(vals)
        _, self.evecs_cv = numpy.linalg.eig(cov)
        metric_params.evecs_cv = {}
        metric_params.evecs_cv[self.f_upper] = self.evecs_cv

        vals = pycbc.tmpltbank.estimate_mass_range(
            100000, mass_range_params, metric_params, self.f_upper, covary=False
        )

        self.metric_params = metric_params
        self.mass_range_params = mass_range_params
        self.mass_range_params2 = mass_range_params2
        self.ethinca_params = pycbc.tmpltbank.EthincaParameters(
            self.ethinca_order,
            self.ethinca_cutoff,
            self.ethinca_freq_step,
            full_ethinca=False,
            time_ethinca=False,
        )

        self.xis = vals

    def test_eigen_directions(self):
        fname = "stockEvals.dat.gz"
        apy_fname = get_file(DATA_FILE_URL.format(fname), cache=True)
        with gzip.open(apy_fname) as apy_fp:
            evals_stock = Array(numpy.loadtxt(apy_fp))

        fname = "stockEvecs.dat.gz"
        apy_fname = get_file(DATA_FILE_URL.format(fname), cache=True)
        with gzip.open(apy_fname) as apy_fp:
            evecs_stock = Array(numpy.loadtxt(apy_fp))

        max_eval = max(evals_stock)
        evals_curr = Array(self.metric_params.evals[self.f_upper])
        evecs_curr = Array(self.metric_params.evecs[self.f_upper])
        # Uncomment these lines to regenerate the data files
        # numpy.savetxt('newEvals.dat', evals_curr)
        # numpy.savetxt('newEvecs.dat', evecs_curr)
        err_msg = "pycbc.tmpltbank.determine_eigen_directions has failed "
        err_msg += "sanity check."
        evals_diff = abs(evals_curr - evals_stock) / max_eval
        self.assertTrue(not (evals_diff > 1e-5).any(), msg=err_msg)
        for stock, test in zip(evecs_stock.data, evecs_curr.data):
            stock_scaled = stock * evals_curr.data**0.5
            test_scaled = test * evals_curr.data**0.5
            diff = stock_scaled - test_scaled
            self.assertTrue(not (diff > 1e-4).any(), msg=err_msg)

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
                if i in [35, 38, 41, 44, 62, 65, 68, 71, 80]:
                    continue
                raise

            # Check that if the mass limits have changed, it was right to do so
            # This is not exhaustive, but gets most things
            if (curr_min_mass is None) or not math.isclose(
                self.min_total_mass, curr_min_mass, rel_tol=1e-06
            ):
                min_comp_mass = self.min_mass1 + self.min_mass2
                min_eta = (
                    self.min_mass1 * self.min_mass2 / (min_comp_mass * min_comp_mass)
                )
                min_chirp_mass = min_comp_mass * min_eta ** (3.0 / 5.0)
                if (min_comp_mass is not None) and math.isclose(
                    self.min_total_mass, min_comp_mass, rel_tol=1e-06
                ):
                    # Okay, the total mass is changed by the components
                    pass
                elif (self.min_eta and min_eta < self.min_eta) or (
                    self.max_eta and min_eta > self.max_eta
                ):
                    # Okay, not possible from eta
                    pass
                elif self.min_chirp_mass and min_chirp_mass < self.min_chirp_mass:
                    # Okay, not possible from chirp mass
                    pass
                else:
                    err_msg = "Minimum total mass changed unexpectedly."
                    self.fail(err_msg)
            if (curr_max_mass is None) or not math.isclose(
                self.max_total_mass, curr_max_mass, rel_tol=1e-06
            ):
                max_comp_mass = self.max_mass1 + self.max_mass2
                max_eta = (
                    self.max_mass1 * self.max_mass2 / (max_comp_mass * max_comp_mass)
                )
                max_chirp_mass = max_comp_mass * max_eta ** (3.0 / 5.0)
                if (max_comp_mass is not None) and math.isclose(
                    self.max_total_mass, max_comp_mass, rel_tol=1e-06
                ):
                    # Okay, the total mass is changed by the components
                    pass
                elif (self.min_eta and max_eta < self.min_eta) or (
                    self.max_eta and max_eta > self.max_eta
                ):
                    # Okay, not possible from eta
                    pass
                elif self.max_chirp_mass and max_chirp_mass > self.max_chirp_mass:
                    # Okay, not possible from chirp mass
                    pass
                else:
                    err_msg = "Maximum total mass changed unexpectedly."
                    self.fail(err_msg)

            mass_range_params = pycbc.tmpltbank.MassRangeParameters(
                self.min_mass1,
                self.max_mass1,
                self.min_mass2,
                self.max_mass2,
                max_ns_spin_mag=self.max_ns_spin_mag,
                max_bh_spin_mag=self.max_bh_spin_mag,
                max_tot_mass=self.max_total_mass,
                min_tot_mass=self.min_total_mass,
                max_chirp_mass=self.max_chirp_mass,
                min_chirp_mass=self.min_chirp_mass,
                max_eta=self.max_eta,
                min_eta=self.min_eta,
                ns_bh_boundary_mass=self.ns_bh_boundary_mass,
            )

            # And again with the nsbh flag
            mass_range_params2 = pycbc.tmpltbank.MassRangeParameters(
                self.min_mass1,
                self.max_mass1,
                self.min_mass2,
                self.max_mass2,
                max_ns_spin_mag=self.max_ns_spin_mag,
                max_bh_spin_mag=self.max_bh_spin_mag,
                max_tot_mass=self.max_total_mass,
                min_tot_mass=self.min_total_mass,
                max_chirp_mass=self.max_chirp_mass,
                min_chirp_mass=self.min_chirp_mass,
                max_eta=self.max_eta,
                min_eta=self.min_eta,
                nsbh_flag=True,
            )

            mass1, mass2, spin1z, spin2z = pycbc.tmpltbank.get_random_mass(
                100000, mass_range_params
            )
            mass = mass1 + mass2
            err_msg = "pycbc.tmpltbank.get_random_mass returns invalid ranges."
            self.assertTrue(not (mass < self.min_total_mass).any(), msg=err_msg)
            self.assertTrue(not (mass > self.max_total_mass).any(), msg=err_msg)
            self.assertTrue(not (mass1 > self.max_mass1 * 1.001).any(), msg=err_msg)
            self.assertTrue(not (mass1 < self.min_mass1 * 0.999).any(), msg=err_msg)
            self.assertTrue(not (mass2 > self.max_mass2 * 1.001).any(), msg=err_msg)
            self.assertTrue(not (mass2 < self.min_mass2 * 0.999).any(), msg=err_msg)
            self.assertTrue(not (mass1 < mass2).any(), msg=err_msg)
            # Chirp mass and eta
            mchirp, eta = pnutils.mass1_mass2_to_mchirp_eta(mass1, mass2)
            if self.max_chirp_mass:
                self.assertTrue(
                    not (mchirp > self.max_chirp_mass * 1.0001).any(), msg=err_msg
                )
            if self.min_chirp_mass:
                self.assertTrue(
                    not (mchirp < self.min_chirp_mass * 0.9999).any(), msg=err_msg
                )
            if self.min_eta:
                self.assertTrue(not (eta < self.min_eta * 0.9999).any(), msg=err_msg)
                self.assertTrue(not (eta > self.max_eta * 1.0001).any(), msg=err_msg)
            ns_spin1 = spin1z[mass1 < self.ns_bh_boundary_mass]
            ns_spin2 = spin2z[mass2 < self.ns_bh_boundary_mass]
            bh_spin1 = spin1z[mass1 > self.ns_bh_boundary_mass]
            bh_spin2 = spin2z[mass2 > self.ns_bh_boundary_mass]
            self.assertTrue(not (abs(ns_spin1) > 0.5).any(), msg=err_msg)
            self.assertTrue(not (abs(ns_spin2) > 0.5).any(), msg=err_msg)
            self.assertTrue(not (abs(bh_spin1) > 0.9).any(), msg=err_msg)
            self.assertTrue(not (abs(bh_spin2) > 0.9).any(), msg=err_msg)
            # Check that *some* spins are bigger than 0.5
            if len(bh_spin1):
                self.assertTrue((abs(bh_spin1) > 0.5).any(), msg=err_msg)
            if len(bh_spin2):
                self.assertTrue((abs(bh_spin2) > 0.5).any(), msg=err_msg)

            # Check nsbh flag
            mass1, mass2, spin1z, spin2z = pycbc.tmpltbank.get_random_mass(
                100000, mass_range_params2
            )
            self.assertTrue(not (abs(spin1z) > 0.9).any(), msg=err_msg)
            self.assertTrue(not (abs(spin2z) > 0.5).any(), msg=err_msg)
            self.assertTrue((abs(spin1z) > 0.5).any(), msg=err_msg)

    def test_metric_match_prediction(self):
        mass1a, mass2a, spin1za, spin2za = pycbc.tmpltbank.get_random_mass(
            10, self.mass_range_params
        )
        mass1b, mass2b, spin1zb, spin2zb = pycbc.tmpltbank.get_random_mass(
            10, self.mass_range_params
        )
        for idx in range(10):
            masses1 = [mass1a[idx], mass2a[idx], spin1za[idx], spin2za[idx]]
            masses2 = [mass1b[idx], mass2b[idx], spin1zb[idx], spin2zb[idx]]
            dist, _, _ = pycbc.tmpltbank.get_point_distance(
                masses1, masses2, self.metric_params, self.f_upper
            )
            opt_dist = 0.02
            while dist > opt_dist * 1.01 or dist < opt_dist * 0.99:
                dist_fac = opt_dist / dist
                dist_fac = dist_fac**0.5
                dist_fac = max(dist_fac, 0.01)
                dist_fac = min(dist_fac, 2)
                for idx, curr_mass2 in enumerate(masses2):
                    masses2[idx] = masses1[idx] + (curr_mass2 - masses1[idx]) * dist_fac
                dist, _, _ = pycbc.tmpltbank.get_point_distance(
                    masses1, masses2, self.metric_params, self.f_upper
                )
            self.assertFalse(numpy.isnan(dist))

            htilde1, _ = get_fd_waveform(
                approximant="TaylorF2",
                mass1=masses1[0],
                mass2=masses1[1],
                spin1z=masses1[2],
                spin2z=masses1[3],
                delta_f=1.0 / 256,
                f_lower=15,
                f_final=2000,
            )
            htilde2, _ = get_fd_waveform(
                approximant="TaylorF2",
                mass1=masses2[0],
                mass2=masses2[1],
                spin1z=masses2[2],
                spin2z=masses2[3],
                delta_f=1.0 / 256,
                f_lower=15,
                f_final=2000,
            )
            overlap, _ = match(
                htilde1, htilde2, psd=self.psd_for_match, low_frequency_cutoff=15
            )
            self.assertTrue(overlap > 0.97 and overlap < 0.985)

    def test_chirp_params(self):
        chirps = pycbc.tmpltbank.get_chirp_params(
            2.2, 1.8, 0.2, 0.3, self.metric_params.f0, self.metric_params.pn_order
        )
        fname = "stockChirps.dat.gz"
        apy_fname = get_file(DATA_FILE_URL.format(fname), cache=True)
        with gzip.open(apy_fname) as apy_fp:
            stock_chirps = numpy.loadtxt(apy_fp)
        diff = (chirps - stock_chirps) / stock_chirps
        err_msg = "Calculated chirp params differ from that expected."
        self.assertTrue(not (abs(diff) > 1e-4).any(), msg=err_msg)

    def test_hexagonal_placement(self):
        arrz = pycbc.tmpltbank.generate_hexagonal_lattice(10, 0, 10, 0, 0.03)
        arrz = numpy.array(arrz)
        fname = "stockHexagonal.dat.gz"
        apy_fname = get_file(DATA_FILE_URL.format(fname), cache=True)
        with gzip.open(apy_fname) as apy_fp:
            stock_grid = numpy.loadtxt(apy_fp)
        diff = arrz - stock_grid
        err_msg = "Calculated lattice differs from that expected."
        self.assertTrue(not (diff > 1e-4).any(), msg=err_msg)

    def test_anstar_placement(self):
        arrz = pycbc.tmpltbank.generate_anstar_3d_lattice(0, 10, 0, 10, 0, 10, 0.03)
        arrz = numpy.array(arrz)
        fname = "stockAnstar3D.dat.gz"
        apy_fname = get_file(DATA_FILE_URL.format(fname), cache=True)
        with gzip.open(apy_fname) as apy_fp:
            stock_grid = numpy.loadtxt(apy_fp)
        # Uncomment this line to regenerate the data file
        # numpy.savetxt("new_example.dat", arrz)
        err_msg = "Calculated lattice differs from that expected."
        self.assertTrue(len(arrz) == len(stock_grid), msg=err_msg)
        diff = arrz - stock_grid
        self.assertTrue(not (diff > 1e-4).any(), msg=err_msg)

    def test_get_mass_distribution(self):
        # Just run the function, no checking output
        pycbc.tmpltbank.get_mass_distribution(
            [1.35, 0.239, 0.4, -0.2],
            2,
            self.mass_range_params,
            self.metric_params,
            self.f_upper,
            num_jump_points=123,
            chirp_mass_jump_fac=0.0002,
            eta_jump_fac=0.009,
            spin1z_jump_fac=0.1,
            spin2z_jump_fac=0.2,
        )

    def test_get_phys_cov_masses(self):
        evecs = self.metric_params.evecs[self.f_upper]
        evals = self.metric_params.evals[self.f_upper]
        masses1 = [2.2, 1.8, 0.4, 0.3]
        masses2 = [2.21, 1.79, 0.41, 0.29]
        xis1 = pycbc.tmpltbank.get_cov_params(
            masses1[0],
            masses1[1],
            masses1[2],
            masses1[3],
            self.metric_params,
            self.f_upper,
        )
        xis2 = pycbc.tmpltbank.get_cov_params(
            masses2[0],
            masses2[1],
            masses2[2],
            masses2[3],
            self.metric_params,
            self.f_upper,
        )

        test_xis = [xis1[0], xis1[1]]
        b_mtot, b_eta = pnutils.mass1_mass2_to_mtotal_eta(masses2[0], masses2[1])
        best_masses = [b_mtot, b_eta, masses2[2], masses2[3]]
        best_xis = xis2
        output = pycbc.tmpltbank.get_physical_covaried_masses(
            test_xis,
            best_masses,
            best_xis,
            0.0001,
            self.mass_range_params,
            self.metric_params,
            self.f_upper,
        )
        # Test that returned xis are close enough
        diff = (output[6][0] - test_xis[0]) ** 2
        diff += (output[6][1] - test_xis[1]) ** 2
        err_msg = "pycbc.tmpltbank.get_physical_covaried_masses "
        err_msg += "failed to find a point within the desired limits."
        self.assertTrue(diff < 1e-4, msg=err_msg)
        # Test that returned masses and xis agree
        mass = output[0] + output[1]
        eta = output[0] * output[1] / (mass * mass)
        spin_set = pycbc.pnutils.get_beta_sigma_from_aligned_spins(
            eta, output[2], output[3]
        )
        xis = pycbc.tmpltbank.get_cov_params(
            output[0], output[1], output[2], output[3], self.metric_params, self.f_upper
        )
        err_msg = "Recovered xis do not agree with those expected."
        self.assertTrue(abs(xis[0] - output[6][0]) < 1e-5, msg=err_msg)
        self.assertTrue(abs(xis[1] - output[6][1]) < 1e-5, msg=err_msg)
        self.assertTrue(abs(xis[2] - output[6][2]) < 1e-5, msg=err_msg)
        self.assertTrue(abs(xis[3] - output[6][3]) < 1e-5, msg=err_msg)

        # Test again with nsbh flag on
        output = pycbc.tmpltbank.get_physical_covaried_masses(
            test_xis,
            best_masses,
            best_xis,
            0.0001,
            self.mass_range_params2,
            self.metric_params,
            self.f_upper,
        )
        # Test that returned xis are close enough
        diff = (output[6][0] - test_xis[0]) ** 2
        diff += (output[6][1] - test_xis[1]) ** 2
        err_msg = "pycbc.tmpltbank.get_physical_covaried_masses "
        err_msg += "failed to find a point within the desired limits."
        self.assertTrue(diff < 1e-4, msg=err_msg)
        # Test that returned masses and xis agree
        xis = pycbc.tmpltbank.get_cov_params(
            output[0], output[1], output[2], output[3], self.metric_params, self.f_upper
        )
        err_msg = "Recovered xis do not agree with those expected."
        self.assertTrue(abs(xis[0] - output[6][0]) < 1e-5, msg=err_msg)
        self.assertTrue(abs(xis[1] - output[6][1]) < 1e-5, msg=err_msg)
        self.assertTrue(abs(xis[2] - output[6][2]) < 1e-5, msg=err_msg)
        self.assertTrue(abs(xis[3] - output[6][3]) < 1e-5, msg=err_msg)

    def test_stack_xi_direction(self):
        # Just run the function, no checking output
        evecs = self.metric_params.evecs[self.f_upper]
        evals = self.metric_params.evals[self.f_upper]
        masses1 = [2.2, 1.8, 0.4, 0.3]
        masses2 = [2.21, 1.79, 0.41, 0.29]
        xis1 = pycbc.tmpltbank.get_cov_params(
            masses1[0],
            masses1[1],
            masses1[2],
            masses1[3],
            self.metric_params,
            self.f_upper,
        )
        xis2 = pycbc.tmpltbank.get_cov_params(
            masses2[0],
            masses2[1],
            masses2[2],
            masses2[3],
            self.metric_params,
            self.f_upper,
        )
        test_xis = [xis1[0], xis1[1]]
        b_mtot, b_eta = pnutils.mass1_mass2_to_mtotal_eta(masses2[0], masses2[1])
        best_masses = [b_mtot, b_eta, masses2[2], masses2[3]]
        best_xis = xis2

        depths = pycbc.tmpltbank.stack_xi_direction_brute(
            test_xis,
            best_masses,
            best_xis,
            3,
            0.03,
            self.mass_range_params,
            self.metric_params,
            self.f_upper,
            num_iterations=50,
        )

    def test_point_distance(self):
        masses1 = [2, 2, 0.4, 0.6]
        masses2 = [2.02, 1.97, 0.41, 0.59]
        dist, xis1, xis2 = pycbc.tmpltbank.get_point_distance(
            masses1, masses2, self.metric_params, self.f_upper
        )
        diff = abs((dist - 23.3681922039) / dist)

        err_msg = "Obtained distance does not agree with expected value."
        self.assertTrue(diff < 1e-5, msg=err_msg)

    def test_conv_to_sngl(self):
        # Just run the function, no checking output
        masses1 = [(2, 2, 0.4, 0.3), (4.01, 0.249, 0.41, 0.29)]
        llw_output.convert_to_sngl_inspiral_table(masses1, "a")

    def test_ethinca_calc(self):
        # Just run the function, no checking output
        m1 = 2.0
        m2 = 2.0
        s1z = 0.0
        s2z = 0.0
        # ethinca calc breaks unless f0 = fLow
        self.metric_params.f0 = self.metric_params.f_low
        output = llw_output.calculate_ethinca_metric_comps(
            self.metric_params, self.ethinca_params, m1, m2, s1z, s2z
        )
        # restore initial f0 value
        self.metric_params.f0 = self.f0

    def tearDown(self):
        pass


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TmpltbankTestClass))

if __name__ == "__main__":
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
