# Copyright (C) 2018 Collin Capano
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
This modules provides mass distributions of CBCs.
"""

from pycbc import conversions
from pycbc.distributions.uniform import Uniform


class UniformComponentMasses(Uniform):
    """A distribution uniform in mass1 and mass2.

    When initializing the distribution, the possible range in the component
    masses must be specified. Random samples drawn from this distribution
    always have ``mass1 >= mass2``.

    .. note::

        The ranges of mass1 and mass2 may cover the same range; for instance,
        in the examples below, both mass1 and mass2 have ranges ``[10, 80)``.
        Likewise, the ``(log)pdf`` function does not require that ``mass1 >=
        mass2``. It is only in the ``rvs`` function that this convention is
        enforced.

    Parameters
    ----------
    mass1 : tuple or ``boundaries.Bounds``
        The range for mass1.

    mass2 : tuple of ``boundaries.Bounds``
        The range for mass2.

    Examples
    --------
    Initialize a distribution:

    >>> from pycbc.distributions import UniformComponentMasses
    >>> d = UniformComponentMasses(mass1=(10., 80.), mass2=(10., 80.))

    Draw random variates from the distribution and check that mass1 >= mass2
    for all values:

    >>> vals = d.rvs(size=1000)
    >>> (vals['mass1'] >= vals['mass2']).all()
    True

    Note that the pdf function does not require ``mass1 >= mass2``:

    >>> d.pdf(mass1=60., mass2=20.)
    0.00020408163265306107
    >>> d.pdf(mass1=20., mass2=60.)
    0.00020408163265306107

    """
    name = 'uniform_component_masses'

    def __init__(self, mass1=None, mass2=None):
        # check that parameters for mass1 and mass2 are both provided
        if mass1 is None:
            raise ValueError("must provide limits for mass1")
        if mass2 is None:
            raise ValueError("must provide limits for mass2")
        super(UniformComponentMasses, self).__init__(mass1=mass2, mass2=mass2)

    def rvs(self, size=1, param=None):
        """Gives a set of random values drawn from this distribution.

        In the returned set, mass2 <= mass1.

        Parameters
        ----------
        size : int
            The number of values to generate; default is 1.
        param : str, optional
            If provided, will just return values for the given parameter.

        Returns
        -------
        structured array
            The random values in a numpy structured array.
        """
        arr = super(UniformComponentMasses, self).rvs(size=size)
        # enforce m1 > m2
        m1 = conversions.primary_mass(arr['mass1'], arr['mass2'])
        m2 = conversions.secondary_mass(arr['mass1'], arr['mass2'])
        arr['mass1'][:] = m1
        arr['mass2'][:] = m2
        if param is not None:
            arr = arr[param]
        return arr
