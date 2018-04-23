# Copyright (C) 2018 Miriam Cabero, Collin Capano
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


from pycbc import conversions
from . import uniform

class UniformF0Tau(uniform.Uniform):
    """A distribution uniform in QNM frequency and damping time.

    Constraints may be placed to exclude frequencies and damping times
    corresponding to specific masses and spins.
    """

    name = 'uniform_f0_tau'

    def __init__(self, rdfreq='f_0', damping_time='tau', final_mass=None,
                 final_spin=None, **kwargs):
        self.rdfreq = rdfreq
        self.damping_time = damping_time
        super(UniformF0Tau, self).__init__(**kwargs)
        if final_mass is None:
            final_mass = (1., numpy.inf)
        elif isinstance(final_mass, str) or isinstance(final_mass, unicode):
            final_mass = map(float, final_mass.split(','))
        if final_spin is None:
            final_spin = (-0.996, 0.996)
        elif isinstance(final_spin, str) or isinstance(final_spin, unicode):
            final_spin = map(float, final_spin.split(','))
        self.final_mass_bounds = bounded.boundaries.Bounds(
            min_bound=final_mass[0], max_bound=final_mass[1])
        self.final_spin_bounds = bounded.boundaries.Bounds(
            min_bound=final_spin[0], max_bound=final_spin[1])

    def __contains__(self, params):
        isin = super(UniformF0Tau, self).__contains__(params)
        if isin:
            isin &= self._constraints(params)
        return isin

    def _constraints(self, params):
        f_0 = params[self.rdfreq]
        tau = params[self.damping_time]
        mf = conversions.final_mass_from_f0_tau(f_0, tau)
        sf = conversions.final_spin_from_f0_tau(f_0, tau)
        return (self.final_mass_bounds.__contains__(mf)) & (
                self.final_spin_bounds.__contains__(sf))

    def rvs(self, size=1):
        size = int(size)
        dtype = [(p, float) for p in self.params]
        arr = numpy.zeros(size, dtype=dtype)
        remaining = size
        keepidx = 0
        while remaining:
            draws = super(UniformF0Tau, self).rvs(size=remaining)
            mask = self._constraints(draws)
            addpts = mask.sum()
            arr[keepidx:keepidx+addpts] = draws[mask]
            keepidx += addpts
            remaining = size - keepidx
    return arr
