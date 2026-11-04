import argparse
import unittest

import numpy as np
from utils import simple_exit

from pycbc.filter import matched_filter_core, sigma
from pycbc.noise import frequency_noise_from_psd
from pycbc.psd import aLIGOZeroDetHighPower
from pycbc.vetoes import SingleDetBankVeto
from pycbc.waveform import get_fd_waveform


class TestBankChisq(unittest.TestCase):
    def setUp(self):
        self.f_low = 30.0
        self.delta_f = 1.0 / 16.0
        self.flen = int(2048.0 / self.delta_f) // 2 + 1
        self.psd = aLIGOZeroDetHighPower(self.flen, self.delta_f, self.f_low)
        kmin = int(self.f_low / self.delta_f)
        self.psd.data[:kmin] = self.psd[kmin]

        self.template = self._make_template(30.0)
        self.filters = [self._make_template(15.0), self._make_template(20.0)]

        self.veto = SingleDetBankVeto(None)
        self.veto.do = True
        self.veto.f_low = self.f_low
        self.veto.filters = self.filters
        self.veto.dof = 2 * len(self.filters)
        self.veto._overlaps_cache = {}
        self.veto._segment_snrs_cache = {}

        self.noise = frequency_noise_from_psd(self.psd, seed=10)
        self.peak = np.array([int(8.0 * 2048)])
        self.signal = self._inject(self.template, snr=15.0)

    def _make_template(self, mass):
        hp, _ = get_fd_waveform(
            approximant="IMRPhenomD",
            mass1=mass,
            mass2=mass,
            f_lower=self.f_low,
            delta_f=self.delta_f,
        )
        hp.resize(self.flen)
        hp.f_lower = self.f_low
        hp.params = argparse.Namespace(mass1=mass, mass2=mass)
        return hp

    def _inject(self, waveform, snr):
        norm = snr / sigma(waveform, psd=self.psd, low_frequency_cutoff=self.f_low)
        return self.noise + waveform.cyclic_time_shift(-8.0) * norm

    def _evaluate_at_peak(self, stilde, overwhitened=False):
        raw = stilde * self.psd if overwhitened else stilde
        snrv, _, norm = matched_filter_core(
            self.template, raw, psd=self.psd, low_frequency_cutoff=self.f_low
        )
        return self.veto.values(
            self.template,
            self.psd,
            stilde,
            snrv.numpy()[self.peak],
            norm,
            self.peak,
            overwhitened=overwhitened,
        )

    def test_matching_signal_in_gaussian_noise(self):
        chisq, dof = self._evaluate_at_peak(self.signal)
        self.assertEqual(dof[0], 4)
        self.assertAlmostEqual(chisq[0] / dof[0], 1.0, delta=1.0)

    def test_mismatched_signal(self):
        mismatched = self._inject(self.filters[0], snr=15.0)
        chisq, dof = self._evaluate_at_peak(mismatched)
        self.assertGreater(chisq[0] / dof[0], 10.0)

    def test_overwhitened_data(self):
        chisq, dof = self._evaluate_at_peak(self.signal, overwhitened=False)
        chisq_ow, dof_ow = self._evaluate_at_peak(
            self.signal / self.psd, overwhitened=True
        )
        np.testing.assert_array_equal(dof_ow, dof)
        np.testing.assert_allclose(chisq_ow, chisq, rtol=1e-5)

    def test_cli_options(self):
        parser = argparse.ArgumentParser()
        SingleDetBankVeto.insert_option_group(parser)
        self.assertFalse(SingleDetBankVeto.from_cli(parser.parse_args([])).do)
        veto = SingleDetBankVeto.from_cli(
            parser.parse_args(["--bank-veto-bank-file", "bank.hdf"])
        )
        self.assertTrue(veto.do)
        self.assertEqual(veto.bank_file, "bank.hdf")


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestBankChisq))

if __name__ == "__main__":
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
