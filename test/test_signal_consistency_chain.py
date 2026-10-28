import argparse
import types
import unittest

import numpy

from pycbc.filter import matched_filter, matched_filter_core, sigma
from pycbc.types import FrequencySeries, TimeSeries
from pycbc.vetoes.autochisq import SingleDetAutoChisq
from pycbc.vetoes.bank_chisq import SingleDetBankVeto
from pycbc.vetoes.chain import SignalConsistencyChain
from pycbc.vetoes.chisq import SingleDetPowerChisq
from pycbc.vetoes.sgchisq import SingleDetSGChisq
from pycbc.waveform import get_fd_waveform


class StubPower(SingleDetPowerChisq):
    def __init__(self):
        super().__init__(num_bins=4)
        self.calls = 0

    def values(self, corr, snrv, snr_norm, psd, indices, template):
        self.calls += 1
        return numpy.array([10 * abs(snrv[0])]), numpy.array([2 * 4 - 2])


PLACEHOLDER = FrequencySeries(numpy.ones(4), delta_f=1.0)


def run(
    chain,
    snr,
    times,
    tids=None,
    template=lambda j: PLACEHOLDER,
    correlation=lambda j, reach: (None, 0),
    **kwargs,
):
    tids = numpy.zeros(len(snr), dtype=int) if tids is None else tids
    return chain.values(
        numpy.array(snr, dtype=complex),
        numpy.array(times),
        numpy.array(tids),
        template,
        correlation,
        **kwargs,
    )


class TestChain(unittest.TestCase):
    def test_window_and_iteration_cap(self):
        power = StubPower()
        chain = SignalConsistencyChain(power, steps=["power"], window=10)
        keep, vals = run(chain, [5, 9, 6, 7], [0, 5, 100, 200])
        self.assertTrue(keep.all())
        numpy.testing.assert_allclose(vals["chisq"], [90, 90, 60, 70])
        self.assertEqual(power.calls, 3)

        # Window sharing is per-template
        _, vals = run(chain, [9, 5], [0, 5], tids=[0, 1])
        numpy.testing.assert_allclose(vals["chisq"], [90, 50])

        # max_iterations caps evaluations, loudest trigger first
        capped = SignalConsistencyChain(
            StubPower(), steps=["power"], window=0, max_iterations=1
        )
        _, vals = run(capped, [5, 9, 7], [0, 50, 100])
        numpy.testing.assert_allclose(vals["chisq"], [0, 90, 0])

    def test_thresholds_and_cuts(self):
        power = StubPower()
        chain = SignalConsistencyChain(
            power, steps=["snr:4:lower:remove", "snr:6:lower:skip", "power"]
        )
        keep, vals = run(chain, [5, 9, 3], [0, 50, 100])
        numpy.testing.assert_array_equal(keep, [True, True, False])
        numpy.testing.assert_allclose(vals["chisq"], [0, 90, 0])
        numpy.testing.assert_array_equal(vals["chisq_dof"], [-100, 6, -100])
        self.assertEqual(power.calls, 1)

        # Cut on a ranking statistic computed from an earlier stage
        chain = SignalConsistencyChain(
            StubPower(), steps=["power", "newsnr:2.4:lower:remove"]
        )
        keep, _ = run(chain, [2.5, 9.0], [0, 50])
        numpy.testing.assert_array_equal(keep, [False, True])

        # Cut on a supplied chisq choice
        chain = SignalConsistencyChain(
            StubPower(),
            steps=["traditional:2:upper:remove", "power"],
            input_columns=["chisq", "chisq_dof"],
        )
        cols = {"chisq": numpy.array([6.0, 18.0]), "chisq_dof": numpy.array([6, 6])}
        keep, _ = run(chain, [5, 9], [0, 50], columns=cols)
        numpy.testing.assert_array_equal(keep, [True, False])

    def test_validation_and_cli(self):
        with self.assertRaises(ValueError):
            SignalConsistencyChain(
                StubPower(), steps=["newsnr:5:lower:remove", "power"]
            )
        with self.assertRaises(ValueError):
            SignalConsistencyChain(StubPower(), steps=["snr:5:lower", "power"])
        with self.assertRaises(ValueError):
            SignalConsistencyChain(StubPower(), steps=["power", "sg"])

        parser = argparse.ArgumentParser()
        SignalConsistencyChain.insert_option_group(parser)
        args = parser.parse_args(
            [
                "--chisq-bins",
                "16",
                "--signal-consistency-chain",
                "snr:6.0:lower:skip",
                "power",
            ]
        )
        chain = SignalConsistencyChain.from_cli(args, None, 1024.0)
        self.assertEqual(chain.description, "snr:6.0:lower:skip power")

    def test_multi_stage_progression(self):
        class StubSG(SingleDetSGChisq):
            do = True
            calls = 0

            def __init__(self):
                pass

            def values(self, *args, **kwargs):
                self.calls += 1
                return numpy.array([2.0])

        stilde = types.SimpleNamespace(delta_f=1.0, psd=None)
        sg = StubSG()
        chain = SignalConsistencyChain(
            StubPower(), sg_chisq=sg, steps=["power", "sg"], window=0, max_iterations=1
        )
        _, vals = run(chain, [5, 9, 7], [0, 50, 100], stilde=stilde)
        self.assertEqual(sg.calls, 1)
        numpy.testing.assert_allclose(vals["sg_chisq"], [1, 2, 1])

    def test_auto_and_bank_stages(self):
        rng = numpy.random.default_rng(1)
        n, sr = 8192, 1024.0
        df = sr / n
        hp, _ = get_fd_waveform(
            approximant="TaylorF2", mass1=10, mass2=10, f_lower=30, delta_f=df
        )
        hp.resize(n // 2 + 1)
        hp.f_lower = 30.0
        hp.params = types.SimpleNamespace(mass1=10.0, mass2=10.0)

        psd = FrequencySeries(numpy.ones(n // 2 + 1), delta_f=df)
        noise = TimeSeries(rng.normal(size=n) * (sr / 2) ** 0.5, delta_t=1 / sr)
        signal = (
            hp.cyclic_time_shift(2.0).to_timeseries()
            * 20
            / sigma(hp, psd=psd, low_frequency_cutoff=30)
        )
        stilde = TimeSeries(
            noise.numpy() + signal.numpy(), delta_t=1 / sr
        ).to_frequencyseries()
        stilde.psd = psd

        snr = matched_filter(hp, stilde, psd=psd, low_frequency_cutoff=30)
        index = int(numpy.argmax(abs(snr)))
        ref_auto = SingleDetAutoChisq(2, 10, twophase=True)
        expect_auto, expect_auto_dof = ref_auto.values(
            snr, numpy.array([index]), hp, psd, 1.0, low_frequency_cutoff=30
        )

        corr = FrequencySeries(numpy.fft.fft(snr.numpy()) / n, delta_f=df)
        hp.psd = psd
        auto = SingleDetAutoChisq(2, 10, twophase=True)
        chain = SignalConsistencyChain(autochisq=auto, steps=["auto"])
        _, vals = run(
            chain,
            [abs(snr[index])],
            [index],
            template=lambda j: hp,
            correlation=lambda j, reach: (corr, index),
        )
        numpy.testing.assert_allclose(vals["cont_chisq"], expect_auto, rtol=1e-6)
        numpy.testing.assert_array_equal(vals["cont_chisq_dof"], [expect_auto_dof])

        bank_tmplt, _ = get_fd_waveform(
            approximant="TaylorF2", mass1=12, mass2=12, f_lower=30, delta_f=df
        )
        bank_tmplt.resize(n // 2 + 1)
        bank_tmplt.f_lower = 30.0
        bank_tmplt.params = types.SimpleNamespace(mass1=12.0, mass2=12.0)

        snrv, _, norm = matched_filter_core(
            hp, stilde, psd=psd, low_frequency_cutoff=30
        )
        bank_veto = SingleDetBankVeto(None)
        bank_veto.do = True
        bank_veto.f_low = 30.0
        bank_veto.filters = [bank_tmplt]
        bank_veto.dof = 2
        bank_veto._overlaps_cache = {}
        bank_veto._segment_snrs_cache = {}

        expect_bank, expect_bank_dof = bank_veto.values(
            hp,
            psd,
            stilde,
            numpy.array([snrv[index]]),
            norm,
            numpy.array([index]),
            overwhitened=False,
        )
        chain = SignalConsistencyChain(
            bank_chisq=bank_veto, steps=["bank", "bank:50:upper:remove"]
        )
        keep, vals = run(
            chain,
            [snrv[index] * norm],
            [index],
            template=lambda j: hp,
            stilde=stilde,
        )
        self.assertTrue(keep.all())
        numpy.testing.assert_allclose(vals["bank_chisq"], expect_bank, rtol=1e-5)
        numpy.testing.assert_array_equal(vals["bank_chisq_dof"], expect_bank_dof)


if __name__ == "__main__":
    unittest.main()
