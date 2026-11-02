# Copyright (C) 2026 Alex Nitz
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

"""Unittests for the model stats cache."""

import unittest

import numpy
from utils import parse_args_cpu_only, simple_exit

from pycbc.distributions import JointDistribution, Uniform
from pycbc.inference.models import TestNormal

parse_args_cpu_only("model stats cache")


class TestStatsCache(unittest.TestCase):
    """The stats an evaluation worked out can be asked for afterwards."""

    def model(self, size=None):
        params = ["x", "y"]
        prior = JointDistribution(
            params, Uniform(x=(-5.0, 5.0)), Uniform(y=(-5.0, 5.0))
        )
        return TestNormal(params, prior=prior, stats_cache=size)

    @staticmethod
    def points(n, seed=4):
        rng = numpy.random.RandomState(seed)
        return [
            {"x": float(a), "y": float(b)}
            for a, b in rng.uniform(-3.0, 3.0, size=(n, 2))
        ]

    def evaluate(self, model, points):
        for point in points:
            model.update(**point)
            _ = model.logposterior

    def test_recovers_each_points_own_stats(self):
        """
        Asking for a point returns what that point produced, not what the
        most recent evaluation produced, and the order asked in does not
        matter.
        """
        model = self.model(100)
        points = self.points(20)
        self.evaluate(model, points)
        for point in points[::-1]:
            cached = model.cached_stats(point)
            self.assertIsNotNone(cached)
            model.update(**point)
            _ = model.logposterior
            numpy.testing.assert_array_equal(
                numpy.array(cached, dtype=float),
                numpy.array(model.get_current_stats(), dtype=float),
            )

    def test_a_point_never_evaluated_is_a_miss(self):
        """A miss is None, so a caller can tell it has to evaluate."""
        model = self.model(100)
        self.evaluate(model, self.points(5))
        self.assertIsNone(model.cached_stats({"x": 99.0, "y": 99.0}))

    def test_keeps_everything_by_default(self):
        """
        The default is to keep all of them, which is what the samplers
        that carried their own stats already did. A bound is opt in.
        """
        model = self.model()
        points = self.points(50)
        self.evaluate(model, points)
        # asking moves the evaluation still being held into the store
        self.assertIsNotNone(model.cached_stats(points[0]))
        self.assertEqual(len(model._stats_cache), len(points))
        self.assertIsNotNone(model.cached_stats(points[-1]))

    def test_a_size_of_zero_keeps_nothing(self):
        """Zero is a bound like any other, and it bounds to nothing."""
        model = self.model(0)
        points = self.points(5)
        self.evaluate(model, points)
        self.assertIsNone(model.cached_stats(points[0]))

    def test_oldest_go_when_full(self):
        """
        The size is a bound: a run longer than it keeps the recent points
        and forgets the early ones rather than growing without limit.
        """
        model = self.model(10)
        points = self.points(40)
        self.evaluate(model, points)
        self.assertLessEqual(len(model._stats_cache), 10)
        self.assertIsNone(model.cached_stats(points[0]))
        self.assertIsNotNone(model.cached_stats(points[-1]))

    def test_parameters_need_not_be_numbers(self):
        """
        Static params carry strings -- approximant, for one -- so the key
        must not assume it can convert what it is given.
        """
        params = ["x", "y"]
        prior = JointDistribution(
            params, Uniform(x=(-5.0, 5.0)), Uniform(y=(-5.0, 5.0))
        )
        model = TestNormal(
            params,
            prior=prior,
            static_params={"approximant": "TaylorF2"},
            stats_cache=10,
        )
        point = {"x": 0.5, "y": -0.5}
        model.update(**point)
        _ = model.logposterior
        self.assertIsNotNone(model.cached_stats(point))

    def test_revisit_keeps_stats(self):
        """Going back to a point does not throw away what it worked out."""
        model = self.model()
        points = self.points(3)
        self.evaluate(model, points)
        before = model.cached_stats(points[0])
        # revisit without evaluating, then move on, which stores whatever
        # the revisit left behind
        model.update(**points[0])
        model.update(**points[1])
        self.assertEqual(model.cached_stats(points[0]), before)

    def test_revisit_is_evaluated_again(self):
        """
        The cache only answers cached_stats; the model still evaluates.

        Anything about the model other than the parameters can change
        between visits, ignore_failed_waveforms for one.
        """
        model = self.model()
        calls = []
        evaluate = model._loglikelihood
        model._loglikelihood = lambda: calls.append(1) or evaluate()
        points = self.points(2)
        self.evaluate(model, [points[0], points[1], points[0]])
        self.assertEqual(len(calls), 3)

    def test_array_valued_params_are_not_cached(self):
        """Reconstruction updates with arrays of trial values, not a point."""
        model = self.model()
        point = self.points(1)[0]
        self.evaluate(model, [point])
        trial = {"x": numpy.linspace(-1.0, 1.0, 8), "y": numpy.zeros(8)}
        model.update(**trial)  # must not raise on an unhashable key
        self.assertIsNone(model.cached_stats(trial))
        # and the point it did remember is still there afterwards
        self.assertIsNotNone(model.cached_stats(point))

    def test_key_ignores_parameter_order(self):
        """The same point is the same point whatever order it is given in."""
        model = self.model()
        point = self.points(1)[0]
        self.evaluate(model, [point])
        reversed_point = dict(reversed(list(point.items())))
        self.assertEqual(model.cached_stats(reversed_point), model.cached_stats(point))


class TestStatsCacheInAPool(unittest.TestCase):
    """What a worker worked out has to reach the process that asks."""

    def test_gathers_what_the_workers_hold(self):
        """
        Including the evaluation each worker is still holding: that one
        has not been moved into its store yet, and it is what a sampler's
        most recent points are.
        """
        from pycbc.inference import models
        from pycbc.pool import BroadcastPool

        params = ["x", "y"]
        prior = JointDistribution(
            params, Uniform(x=(-5.0, 5.0)), Uniform(y=(-5.0, 5.0))
        )
        model = TestNormal(params, prior=prior)
        models._global_instance = models.CallModel(
            model, "logposterior", return_all_stats=False
        )
        pool = BroadcastPool(2)
        points = [
            [float(a), float(b)]
            for a, b in numpy.random.RandomState(8).uniform(-3.0, 3.0, size=(24, 2))
        ]
        pool.map(models._call_global_model, points)
        model.gather_stats_cache(pool)
        pool.close_pool()
        found = [model.cached_stats({"x": p[0], "y": p[1]}) for p in points]
        self.assertEqual(sum(f is not None for f in found), len(points))

    def test_gather_without_global_instance(self):
        """Samplers like nessai attach a pool without setting _global_instance."""
        from pycbc.inference import models
        from pycbc.pool import SinglePool

        params = ["x", "y"]
        prior = JointDistribution(
            params, Uniform(x=(-5.0, 5.0)), Uniform(y=(-5.0, 5.0))
        )
        model = TestNormal(params, prior=prior)
        models._global_instance = None
        model.update(x=0.5, y=-0.5)
        _ = model.logposterior
        model.gather_stats_cache(SinglePool())
        self.assertIsNotNone(model.cached_stats({"x": 0.5, "y": -0.5}))

    def test_sampling_transforms_hit_and_miss(self):
        """
        Samplers step in sampling_params but write variable_params; both a
        cache hit and a cache miss must work in variable_params space.
        """
        from pycbc.inference.models.base import SamplingTransforms
        from pycbc.transforms import Mass1Mass2ToMchirpQ

        params = ["mass1", "mass2"]
        prior = JointDistribution(
            params, Uniform(mass1=(10.0, 50.0)), Uniform(mass2=(10.0, 50.0))
        )
        st = SamplingTransforms(
            params, ["mchirp", "q"], params, [Mass1Mass2ToMchirpQ()]
        )
        model = TestNormal(params, prior=prior, sampling_transforms=st)
        sp = st.apply({"mass1": 30.0, "mass2": 20.0})
        model.update(mchirp=sp["mchirp"], q=sp["q"])
        _ = model.logposterior
        vp = {p: model.current_params[p] for p in params}
        self.assertIsNotNone(model.cached_stats(vp))


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestStatsCache))
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestStatsCacheInAPool))

if __name__ == "__main__":
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
