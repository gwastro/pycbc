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


import numpy
import pycbc
from pycbc import conversions, boundaries
from . import uniform, bounded


class UniformF0Tau(uniform.Uniform):
    """A distribution uniform in QNM frequency and damping time.

    Constraints may be placed to exclude frequencies and damping times
    corresponding to specific masses and spins.

    To ensure a properly normalized pdf that accounts for the constraints
    on final mass and spin, a renormalization factor is calculated upon
    initialization. This is calculated numerically: f0 and tau are drawn
    randomly, then the norm is scaled by the fraction of points that yield
    final masses and spins within the constraints. The `norm_tolerance` keyword
    arguments sets the error on the estimate of the norm from this numerical
    method. If this value is too large, such that no points are found in
    the allowed region, a ValueError is raised.

    Parameters
    ----------
    f0 : tuple or boundaries.Bounds
        The range of QNM frequencies (in Hz).
    tau : tuple or boundaries.Bounds
        The range of QNM damping times (in s).
    final_mass : tuple or boundaries.Bounds, optional
        The range of final masses to allow. Default is [0,inf).
    final_spin : tuple or boundaries.Bounds, optional
        The range final spins to allow. Must be in [-0.996, 0.996], which is
        the default.
    rdfreq : str, optional
        Use the given string as the name for the f0 parameter. Default is 'f0'.
    damping_time : str, optional
        Use the given string as the name for the tau parameter. Default is
        'tau'.
    norm_tolerance : float, optional
        The tolerance on the estimate of the normalization. Default is 1e-3.
    norm_seed : int, optional
        Seed to use for the random number generator when estimating the norm.
        Default is 0. After the norm is estimated, the random number generator
        is set back to the state it was in upon initialization.

    Examples
    --------

    Create a distribution:

    >>> dist = UniformF0Tau(f0=(10., 2048.), tau=(1e-4,1e-2))

    Check that all random samples drawn from the distribution yield final
    masses > 1:

    >>> from pycbc import conversions
    >>> samples = dist.rvs(size=1000)
    >>> (conversions.final_mass_from_f0_tau(samples['f0'],
            samples['tau']) > 1.).all()
    True

    Create a distribution with tighter bounds on final mass and spin:

    >>> dist = UniformF0Tau(f0=(10., 2048.), tau=(1e-4,1e-2),
            final_mass=(20., 200.), final_spin=(0,0.996))

    Check that all random samples drawn from the distribution are in the
    final mass and spin constraints:

    >>> samples = dist.rvs(size=1000)
    >>> (conversions.final_mass_from_f0_tau(samples['f0'],
            samples['tau']) >= 20.).all()
    True
    >>> (conversions.final_mass_from_f0_tau(samples['f0'],
            samples['tau']) < 200.).all()
    True
    >>> (conversions.final_spin_from_f0_tau(samples['f0'],
            samples['tau']) >= 0.).all()
    True
    >>> (conversions.final_spin_from_f0_tau(samples['f0'],
            samples['tau']) < 0.996).all()
    True

    """

    name = 'uniform_f0_tau'

    def __init__(self, f0=None, tau=None, final_mass=None, final_spin=None,
                 rdfreq='f0', damping_time='tau', norm_tolerance=1e-3,
                 norm_seed=0):
        if f0 is None:
            raise ValueError("must provide a range for f0")
        if tau is None:
            raise ValueError("must provide a range for tau")
        self.rdfreq = rdfreq
        self.damping_time = damping_time
        parent_args = {rdfreq: f0, damping_time: tau}
        super(UniformF0Tau, self).__init__(**parent_args)
        if final_mass is None:
            final_mass = (0., numpy.inf)
        if final_spin is None:
            final_spin = (-0.996, 0.996)
        self.final_mass_bounds = boundaries.Bounds(
            min_bound=final_mass[0], max_bound=final_mass[1])
        self.final_spin_bounds = boundaries.Bounds(
            min_bound=final_spin[0], max_bound=final_spin[1])
        # Re-normalize to account for cuts: we'll do this by just sampling
        # a large number of spaces f0 taus, and seeing how many are in the
        # desired range.
        # perseve the current random state
        s = numpy.random.get_state()
        numpy.random.seed(norm_seed)
        nsamples = int(1./norm_tolerance**2)
        draws = super(UniformF0Tau, self).rvs(size=nsamples)
        # reset the random state
        numpy.random.set_state(s)
        num_in = self._constraints(draws).sum()
        # if num_in is 0, than the requested tolerance is too large
        if num_in == 0:
            raise ValueError("the normalization is < then the norm_tolerance; "
                             "try again with a smaller nrom_tolerance")
        self._lognorm += numpy.log(num_in) - numpy.log(nsamples)
        self._norm = numpy.exp(self._lognorm)

    def __contains__(self, params):
        isin = super(UniformF0Tau, self).__contains__(params)
        if isin:
            isin &= self._constraints(params)
        return isin

    def _constraints(self, params):
        f0 = params[self.rdfreq]
        tau = params[self.damping_time]
        # temporarily silence invalid warnings... these will just be ruled out
        # automatically
        orig = numpy.seterr(invalid='ignore')
        mf = conversions.final_mass_from_f0_tau(f0, tau)
        sf = conversions.final_spin_from_f0_tau(f0, tau)
        isin = (self.final_mass_bounds.__contains__(mf)) & (
                self.final_spin_bounds.__contains__(sf))
        numpy.seterr(**orig)
        return isin

    def rvs(self, size=1):
        """Draw random samples from this distribution.

        Parameters
        ----------
        size : int, optional
            The number of draws to do. Default is 1.

        Returns
        -------
        array
            A structured array of the random draws.
        """
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

    @classmethod
    def from_config(cls, cp, section, variable_args):
        """Initialize this class from a config file.

        Bounds on ``f0``, ``tau``, ``final_mass`` and ``final_spin`` should
        be specified by providing ``min-{param}`` and ``max-{param}``. If
        the ``f0`` or ``tau`` param should be renamed, ``rdfreq`` and
        ``damping_time`` should be provided; these must match
        ``variable_args``. If ``rdfreq`` and ``damping_time`` are not
        provided, ``variable_args`` are expected to be ``f0`` and ``tau``.

        Only ``min/max-f0`` and ``min/max-tau`` need to be provided.

        Example:

        .. code-block:: ini

            [{section}-f0+tau]
            name = uniform_f0_tau
            min-f0 = 10
            max-f0 = 2048
            min-tau = 0.0001
            max-tau = 0.010
            min-final_mass = 10

        Parameters
        ----------
        cp : pycbc.workflow.WorkflowConfigParser
            WorkflowConfigParser instance to read.
        section : str
            The name of the section to read.
        variable_args : str
            The name of the variable args. These should be separated by
            ``pycbc.VARARGS_DELIM``.

        Returns
        -------
        UniformF0Tau :
            This class initialized with the parameters provided in the config
            file.
        """
        tag = variable_args
        variable_args = set(variable_args.split(pycbc.VARARGS_DELIM))
        # get f0 and tau
        f0 = bounded.get_param_bounds_from_config(cp, section, tag, 'f0')
        tau = bounded.get_param_bounds_from_config(cp, section, tag, 'tau')
        # see if f0 and tau should be renamed
        if cp.has_option_tag(section, 'rdfreq', tag):
            rdfreq = cp.get_opt_tag(section, 'rdfreq', tag)
        else:
            rdfreq = 'f0'
        if cp.has_option_tag(section, 'damping_time', tag):
            damping_time = cp.get_opt_tag(section, 'damping_time', tag)
        else:
            damping_time = 'tau'
        # check that they match whats in the variable args
        if not variable_args == set([rdfreq, damping_time]):
            raise ValueError("variable args do not match rdfreq and "
                             "damping_time names")
        # get the final mass and spin values, if provided
        final_mass = bounded.get_param_bounds_from_config(
            cp, section, tag, 'final_mass')
        final_spin = bounded.get_param_bounds_from_config(
            cp, section, tag, 'final_spin')
        extra_opts = {}
        if cp.has_option_tag(section, 'norm_tolerance', tag):
            extra_opts['norm_tolerance'] = float(
                cp.get_opt_tag(section, 'norm_tolerance', tag))
        if cp.has_option_tag(section, 'norm_seed', tag):
            extra_opts['norm_seed'] = int(
                cp.get_opt_tag(section, 'norm_seed', tag))
        return cls(f0=f0, tau=tau,
                   final_mass=final_mass, final_spin=final_spin,
                   rdfreq=rdfreq, damping_time=damping_time,
                   **extra_opts)
