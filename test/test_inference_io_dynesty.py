import os
import unittest
import tempfile
import numpy as np
from utils import simple_exit, parse_args_cpu_only
from pycbc.inference.io.dynesty import DynestyFile

parse_args_cpu_only("inference.io.dynesty")


def write_dummy_file(path, nsamples=2000, seed=42):
    """Write a minimal file with the datasets DynestyFile needs."""
    rng = np.random.default_rng(seed)
    with DynestyFile(path, "w") as fp:
        group = fp.create_group("samples")
        group["mass1"] = rng.uniform(10., 50., nsamples)
        group["logwt"] = rng.normal(0., 3., nsamples)
        group["loglikelihood"] = rng.normal(100., 5., nsamples)
        fp.attrs["log_evidence"] = 0.0


class TestDynestyIO(unittest.TestCase):

    def test_posterior_is_reproducible(self):
        """Reading the same file twice must give the same posterior.

        The weighted samples are resampled on every read, so the draw has to
        come from the seeded generator rather than numpy's global one.
        """
        with tempfile.TemporaryDirectory() as td:
            path = os.path.join(td, "dynesty.hdf")
            write_dummy_file(path)

            with DynestyFile(path, "r") as fp:
                first = fp.read_raw_samples(["mass1"])
            with DynestyFile(path, "r") as fp:
                second = fp.read_raw_samples(["mass1"])

            for param in ("mass1", "loglikelihood"):
                np.testing.assert_array_equal(first[param], second[param])

    def test_seed_changes_the_posterior(self):
        """A different seed must give a different, but repeatable, draw."""
        with tempfile.TemporaryDirectory() as td:
            path = os.path.join(td, "dynesty.hdf")
            write_dummy_file(path)

            with DynestyFile(path, "r") as fp:
                zero = fp.read_raw_samples(["mass1"])["mass1"]
            with DynestyFile(path, "r") as fp:
                one = fp.read_raw_samples(["mass1"], seed=1)["mass1"]
            with DynestyFile(path, "r") as fp:
                one_again = fp.read_raw_samples(["mass1"], seed=1)["mass1"]

            self.assertFalse(np.array_equal(zero, one))
            np.testing.assert_array_equal(one, one_again)

    def test_raw_samples_are_unweighted(self):
        """``raw_samples`` returns the datasets untouched."""
        with tempfile.TemporaryDirectory() as td:
            path = os.path.join(td, "dynesty.hdf")
            write_dummy_file(path)

            with DynestyFile(path, "r") as fp:
                raw = fp.read_raw_samples(["mass1"], raw_samples=True)["mass1"]
                stored = fp["samples"]["mass1"][:]

            np.testing.assert_array_equal(raw, stored)


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestDynestyIO))


if __name__ == "__main__":
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
