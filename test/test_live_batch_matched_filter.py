# Copyright (C) 2026  Ian Harry
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

"""
Unit tests for pycbc.filter.matchedfilter.LiveBatchMatchedFilter.

LiveBatchMatchedFilter shares a single pair of correlation / SNR workspace
buffers between all template groups. Groups have different slot sizes, so
the negative-frequency part of one group's slots overlaps data written by
other groups (and by the veto correlation). These tests check that every
template's SNR time series matches a direct calculation regardless of that
stale data.
"""
import unittest
import numpy
from pycbc.types import FrequencySeries
from pycbc.filter.matchedfilter import LiveBatchMatchedFilter, correlate
from utils import parse_args_cpu_only, simple_exit

parse_args_cpu_only("LiveBatchMatchedFilter")

SAMPLE_RATE = 256.0
# (duration in s, number of templates). Different durations give groups with
# different slot sizes, so each group's zero-padded region overlaps data
# written by the others.
BINS = [(4, 3), (8, 5), (12, 2), (20, 4)]
PARAMS_DTYPE = numpy.dtype([('mass1', 'f8'), ('mass2', 'f8')])


def _random_fseries(rng, flen, delta_f):
    data = rng.normal(size=flen) + 1j * rng.normal(size=flen)
    data[0] = 0
    return FrequencySeries(data.astype(numpy.complex64), delta_f=delta_f)


def _make_templates(seed=1):
    rng = numpy.random.default_rng(seed)
    templates = []
    for dur, count in BINS:
        flen = int(dur * SAMPLE_RATE) // 2 + 1
        for _ in range(count):
            htilde = _random_fseries(rng, flen, 1.0 / dur)
            htilde.params = numpy.array(
                [(1.4, 1.4)], dtype=PARAMS_DTYPE
            )[0]
            htilde.id = len(templates)
            htilde.sigmasq = lambda psd: 1.0
            templates.append(htilde)
    return templates


class _FakeStrain(object):
    """Minimal stand-in for the StrainBuffer interface used by
    LiveBatchMatchedFilter: whitened data at each template resolution."""
    sample_rate = SAMPLE_RATE
    trim_padding = 16
    blocksize = 1
    start_time = 1000000000

    def __init__(self, seed):
        self.rng = numpy.random.default_rng(seed)
        self.cache = {}

    def overwhitened_data(self, delta_f):
        if delta_f not in self.cache:
            flen = int(round(SAMPLE_RATE / delta_f)) // 2 + 1
            stilde = _random_fseries(self.rng, flen, delta_f)
            stilde.psd = FrequencySeries(
                numpy.ones(flen, dtype=numpy.float32), delta_f=delta_f
            )
            self.cache[delta_f] = stilde
        return self.cache[delta_f]


def _expected_snr(htilde, stilde):
    """Complex SNR time series computed directly: correlate the positive
    frequencies, zero the negative frequencies, unnormalized inverse FFT."""
    flen = len(htilde)
    n = (flen - 1) * 2
    corr = numpy.zeros(n, dtype=numpy.complex128)
    corr[:flen] = numpy.conj(htilde.numpy()) * stilde.numpy()
    return numpy.fft.ifft(corr) * n


class TestLiveBatchMatchedFilter(unittest.TestCase):

    def _check_all_groups(self, mf, strain):
        """Run every template group once and compare each template's SNR
        series with the direct calculation, right after its own batch."""
        mf.set_data(strain)
        ngroups = len(mf.tgroups)
        self.assertGreater(ngroups, 1)
        for _ in range(ngroups):
            tgroup = mf.tgroups[mf.block_id]
            mf._process_batch()
            for htilde in tgroup:
                stilde = strain.overwhitened_data(htilde.delta_f)
                expected = _expected_snr(htilde, stilde)
                got = htilde.out.numpy().astype(numpy.complex128)
                self.assertEqual(len(got), len(expected))
                err = numpy.max(numpy.abs(got - expected))
                self.assertLess(
                    err, 1e-4 * numpy.max(numpy.abs(expected)),
                    msg=f'SNR mismatch for template {htilde.id} '
                        f'(duration {1.0 / htilde.delta_f:g} s)'
                )
        return ngroups

    def _dirty_workspace(self, mf, strain):
        """Write into the shared correlation buffer the way the chisq veto
        does (only the first len(template) samples of each slot)."""
        for tgroup in mf.tgroups:
            for htilde in tgroup:
                stilde = strain.overwhitened_data(htilde.delta_f)
                correlate(htilde, stilde, htilde.cout)

    def _run(self, maxelements):
        mf = LiveBatchMatchedFilter(
            _make_templates(), snr_threshold=0.0, chisq_bins=0,
            sg_chisq=None, maxelements=maxelements
        )
        # First pass over every group, then dirty the workspace, then a
        # second pass with new data: each group then runs on memory that
        # other (differently sized) groups have written into.
        strain = _FakeStrain(seed=2)
        self._check_all_groups(mf, strain)
        self._dirty_workspace(mf, strain)
        self._check_all_groups(mf, _FakeStrain(seed=3))
        return mf

    def test_snr_one_group_per_duration(self):
        """Each duration bin fits in one batch, as on a live worker with
        many MPI ranks."""
        mf = self._run(maxelements=2**30)
        self.assertEqual(len(mf.tgroups), len(BINS))

    def test_snr_split_groups(self):
        """A small batch size splits bins into several groups, so groups of
        the same size reuse one slice of the workspace."""
        mf = self._run(maxelements=4 * 2048)
        self.assertGreater(len(mf.tgroups), len(BINS))

    def test_workspace_is_shared(self):
        """All groups use views of one buffer pair, sized for the largest
        group, instead of one allocation per group."""
        mf = LiveBatchMatchedFilter(
            _make_templates(), snr_threshold=0.0, chisq_bins=0,
            sg_chisq=None, maxelements=2**30
        )
        for mem in (mf.out_mem, mf.cout_mem):
            largest = max(mem.values(), key=len)
            sizes = [len(tgroup) * psize for tgroup, psize
                     in zip(mf.tgroups, mf.chunk_tsamples)]
            self.assertEqual(len(largest), max(sizes))
            for view in mem.values():
                self.assertTrue(
                    numpy.shares_memory(view.numpy(), largest.numpy())
                )


suite = unittest.TestSuite()
suite.addTest(
    unittest.TestLoader().loadTestsFromTestCase(TestLiveBatchMatchedFilter)
)

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
