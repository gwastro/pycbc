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

"""Unittests for the inference sampler file IO helpers."""

import os
import shutil
import tempfile
import unittest

import h5py
import numpy as np
from utils import simple_exit, parse_args_cpu_only

from pycbc.inference.sampler.nautilus import sync_state

parse_args_cpu_only("inference sampler IO")


class TestSyncState(unittest.TestCase):
    """Tests mirroring a sampler's own state into a group of our file."""

    def setUp(self):
        self.dir = tempfile.mkdtemp()
        self.src = h5py.File(os.path.join(self.dir, 'src.h5'), 'w')
        self.dst = h5py.File(os.path.join(self.dir, 'dst.h5'), 'w')

    def tearDown(self):
        self.src.close()
        self.dst.close()
        shutil.rmtree(self.dir)

    def write_src(self, npoints, nbounds):
        """Writes a state shaped the way nautilus writes one: a group of
        points per shell, and a group of weights per bound.
        """
        for name in list(self.src):
            del self.src[name]
        group = self.src.require_group('sampler')
        group.attrs['n_like'] = npoints
        for i in range(nbounds):
            group.create_dataset('points_%d' % i, maxshape=(None, 2),
                                 data=np.arange(2 * npoints).reshape(-1, 2))
            self.src.require_group('bound_%d' % i).create_dataset(
                'weights', data=np.full(npoints, float(i)))

    def assert_mirrors_src(self):
        sync_state(self.src, self.dst)
        paths = []
        self.src.visititems(lambda name, obj: paths.append(name))
        self.assertEqual(sorted(self.dst), sorted(self.src))
        for name in paths:
            self.assertEqual(dict(self.dst[name].attrs),
                             dict(self.src[name].attrs))
            if isinstance(self.src[name], h5py.Dataset):
                np.testing.assert_array_equal(self.dst[name][()],
                                              self.src[name][()])

    def test_mirrors_a_new_state(self):
        self.write_src(4, nbounds=3)
        self.assert_mirrors_src()

    def test_replaces_a_state_that_changed_shape(self):
        self.write_src(4, nbounds=3)
        sync_state(self.src, self.dst)
        self.write_src(9, nbounds=3)
        self.assert_mirrors_src()

    def test_drops_shells_that_are_no_longer_there(self):
        self.write_src(4, nbounds=3)
        sync_state(self.src, self.dst)
        self.write_src(4, nbounds=2)
        self.assert_mirrors_src()
        self.assertNotIn('bound_2', self.dst)
        self.assertNotIn('points_2', self.dst['sampler'])


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestSyncState))

if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
