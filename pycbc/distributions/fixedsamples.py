# Copyright (C) 2020 Alexander Nitz
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
This modules provides classes for evaluating distributions based on a fixed
set of points
"""
import logging, numpy
import numpy.random
from pycbc import VARARGS_DELIM

class FixedSamples(object):

    name = "fixed_samples"

    def __init__(self, params, samples):
        self.params = params
        self.samples = samples

        self.p1 = self.samples[params[0]]
        self.frac = len(self.p1)**0.5 / len(self.p1)
        self.sort = self.p1.argsort()
        self.p1sorted = self.p1[self.sort]

        assert(len(numpy.unique(self.p1)) == len(self.p1))

        if len(params) > 2:
            raise ValueError("Only one or two parameters supported "
                             "for fixed sample distribution")

    def rvs(self, size=1, **kwds):
        "Draw random value"
        i = numpy.random.randint(0, high=len(p1), size=size)
        return {p: self.samples[p][i] for p in self.params}

    def cdfinv(self, **original):
        new = {}

        #First dimension
        u1 = original[self.params[0]]
        i1 = int(round(u1 * len(self.p1)))
        if i1 >= len(self.p1):
            i1 = len(self.p1) - 1
        if i1 < 0:
            i1 = 0
        new[self.params[0]] = p1v = self.p1sorted[i1]
        if len(self.params) == 1:
            return new

        # possible second dimension, probably shouldn't
        # do more dimensions than this
        u2 = original[self.params[1]]
        l = numpy.searchsorted(self.p1sorted, p1v * (1 - self.frac))
        r = numpy.searchsorted(self.p1sorted, p1v * (1 + self.frac))
        if r < l:
            l, r = r, l

        region = numpy.array(self.sort[l:r], ndmin=1)
        p2 = self.samples[self.params[1]]
        p2part = numpy.array(p2[region], ndmin=1)
        l = p2part.argsort()
        p2part = numpy.array(p2part[l], ndmin=1)

        i2 = int(round(u2 * len(p2part)))
        if i2 >= len(p2part):
            i2 = len(p2part) - 1
        if i2 < 0:
            i2 = 0
        new[self.params[1]] = p2part[i2]

        p1part = numpy.array(self.p1[region[l]], ndmin=1)
        new[self.params[0]] = p1part[i2]
        return new

    def apply_boundary_conditions(self, **params):
        return params

    def __call__(self, **kwds):
        raise NotImplementedError("Fixed Sample distributions doesn't support"
                                  "for evaluating a PDF")

    @classmethod
    def from_config(cls, cp, section, tag):
        from pycbc.distributions import read_distributions_from_config
        from pycbc.transforms import (read_transforms_from_config,
                                      apply_transforms, BaseTransform)
        from pycbc.transforms import transforms as global_transforms

        params = tag.split(VARARGS_DELIM)
        subname = cp.get_opt_tag(section, 'subname', tag)
        size = cp.get_opt_tag(section, 'sample-size', tag)

        distsec = '{}_sample'.format(subname)
        dist = read_distributions_from_config(cp, section=distsec)
        if len(dist) > 1:
            raise ValueError("Fixed sample distrubtion only supports a single"
                             " distribution to sample from.")

        logging.info('Drawing samples for fixed sample distribution:{}'.format(params))
        samples = dist[0].rvs(size=int(float(size)))
        samples = {p: samples[p] for p in samples.dtype.names}

        transec = '{}_transform'.format(subname)
        trans = read_transforms_from_config(cp, section=transec)
        if len(trans) > 0:
            trans = trans[0]
            samples = apply_transforms(samples, [trans])
            p1 = samples[params[0]]
            class Thook(BaseTransform):
                name = subname
                _inputs = trans.outputs
                _outputs = trans.inputs
                p1name = params[0]
                sort = p1.argsort()
                p1sorted = p1[sort]
                def transform(self, maps):
                    idx = numpy.searchsorted(self.p1sorted, maps[self.p1name])
                    out = {p: samples[p][self.sort[idx]] for p in self.outputs}
                    return self.format_output(maps, out)
            global_transforms[Thook.name] = Thook
        return cls(params, samples)

__all__ = ['FixedSamples']
