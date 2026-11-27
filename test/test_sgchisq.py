import argparse
import unittest

import numpy as np
from utils import simple_exit

from pycbc.filter import matched_filter_core, sigma
from pycbc.noise import frequency_noise_from_psd
from pycbc.psd import aLIGOZeroDetHighPower
from pycbc.vetoes.sgchisq import SingleDetSGChisq
from pycbc.waveform import get_fd_waveform, sinegauss


class TestSGChisq(unittest.TestCase):
    def setUp(self):
        self.f_low = 30.0
        self.delta_f = 1.0 / 16.0
        self.flen = int(2048.0 / self.delta_f) // 2 + 1
        self.psd = aLIGOZeroDetHighPower(self.flen, self.delta_f, self.f_low)
        kmin = int(self.f_low / self.delta_f)
        self.psd.data[:kmin] = self.psd[kmin]

        hp, _ = get_fd_waveform(
            approximant="IMRPhenomD",
            mass1=30.0,
            mass2=30.0,
            f_lower=self.f_low,
            delta_f=self.delta_f,
        )
        hp.resize(self.flen)
        hp.f_lower = self.f_low
        hp.params = argparse.Namespace(template_hash=1)
        self.template = hp

        self.sg = SingleDetSGChisq(None, num_bins="16")
        self.sg.do = True
        self.sg.params = {1: "20-300,20-350"}

        self.noise = frequency_noise_from_psd(self.psd, seed=42)
        sig_norm = 15.0 / sigma(
            self.template, psd=self.psd, low_frequency_cutoff=self.f_low
        )
        self.signal = self.noise + self.template.cyclic_time_shift(-8.0) * sig_norm
        snrv, _, self.norm = matched_filter_core(
            self.template, self.signal, psd=self.psd, low_frequency_cutoff=self.f_low
        )
        self.peak = np.array([int(np.argmax(abs(snrv.numpy())))])
        self.snrv = snrv.numpy()[self.peak]

    def test_matching_signal_in_gaussian_noise(self):
        chisq = self.sg.values(
            self.signal,
            self.template,
            self.psd,
            self.snrv,
            self.norm,
            indices=self.peak,
            overwhitened=False,
        )
        self.assertAlmostEqual(chisq[0], 1.0, delta=1.0)

    def test_high_frequency_glitch(self):
        fcen = self.sg.peak_frequency(self.template, self.psd) + 300.0
        glitch = sinegauss.fd_sine_gaussian(
            1.0, 20.0, fcen, fcen - 50.0, self.flen * self.delta_f, self.delta_f
        ).cyclic_time_shift(-8.0)
        glitch *= 15.0 / sigma(glitch, psd=self.psd, low_frequency_cutoff=fcen - 50.0)
        chisq = self.sg.values(
            self.signal + glitch,
            self.template,
            self.psd,
            self.snrv,
            self.norm,
            indices=self.peak,
            overwhitened=False,
        )
        self.assertGreater(chisq[0], 10.0)

    def test_overwhitened_data(self):
        chisq = self.sg.values(
            self.signal,
            self.template,
            self.psd,
            self.snrv,
            self.norm,
            indices=self.peak,
            overwhitened=False,
        )
        chisq_ow = self.sg.values(
            self.signal / self.psd,
            self.template,
            self.psd,
            self.snrv,
            self.norm,
            indices=self.peak,
            overwhitened=True,
        )
        np.testing.assert_allclose(chisq_ow, chisq, rtol=1e-5)

    def test_tile_cache_is_bounded_by_memory(self):
        # The same waveform under 24 template hashes: 24 sets of tiles
        first = self.sg.values(
            self.signal,
            self.template,
            self.psd,
            self.snrv,
            self.norm,
            indices=self.peak,
            overwhitened=False,
        )
        tiles = next(iter(self.sg._tile_cache.values()))
        one_set = sum(tile[0].data.nbytes for tile in tiles)
        # Each tile holds only the band it is evaluated on
        for gtem, kmin, kmax, _, _ in tiles:
            self.assertEqual(len(gtem), kmax - kmin)

        # Room for 20 sets: the cache holds 20, and no fewer than 16 with room for none
        for budget, held in ((20 * one_set, 20), (0, 16)):
            sg = SingleDetSGChisq(None, num_bins="16")
            sg.do = True
            sg.params = {}
            sg.tile_cache_bytes = budget
            for template_hash in range(1, 25):
                sg.params[template_hash] = "20-300,20-350"
                self.template.params = argparse.Namespace(template_hash=template_hash)
                chisq = sg.values(
                    self.signal,
                    self.template,
                    self.psd,
                    self.snrv,
                    self.norm,
                    indices=self.peak,
                    overwhitened=False,
                )
                np.testing.assert_array_equal(chisq, first)
            self.assertEqual(len(sg._tile_cache), held)
        self.template.params = argparse.Namespace(template_hash=1)

    def test_cli_options(self):
        parser = argparse.ArgumentParser()
        SingleDetSGChisq.insert_option_group(parser)
        parser.add_argument("--chisq-bins", default="16")
        self.assertFalse(SingleDetSGChisq.from_cli(parser.parse_args([])).do)

    def test_overlapping_tile_subtract_signal(self):
        sg_inband = SingleDetSGChisq(None, num_bins="16")
        sg_inband.do = True
        sg_inband.params = {1: "20--120"}

        unsubtracted = sg_inband.values(
            self.signal,
            self.template,
            self.psd,
            self.snrv,
            self.norm,
            indices=self.peak,
            overwhitened=False,
            subtract_signal=False,
        )
        subtracted = sg_inband.values(
            self.signal,
            self.template,
            self.psd,
            self.snrv,
            self.norm,
            indices=self.peak,
            overwhitened=False,
            subtract_signal=True,
        )
        self.assertGreater(unsubtracted[0], 5.0)
        self.assertAlmostEqual(subtracted[0], 1.0, delta=1.0)


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestSGChisq))

if __name__ == "__main__":
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
