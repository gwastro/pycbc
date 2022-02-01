# Copyright (C) 2017  Christopher M. Biwer
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
This modules provides classes and functions for transforming parameters.
"""

import os
import logging
import numpy
from pycbc import conversions
from pycbc import coordinates
from pycbc import cosmology
from pycbc.io import record
from pycbc.waveform import parameters
from pycbc.boundaries import Bounds
from pycbc import VARARGS_DELIM
from pycbc.pnutils import jframe_to_l0frame


class BaseTransform(object):
    """A base class for transforming between two sets of parameters."""

    name = None
    inverse = None
    _inputs = []
    _outputs = []

    def __init__(self):
        self.inputs = set(self._inputs)
        self.outputs = set(self._outputs)

    def __call__(self, maps):
        return self.transform(maps)

    def transform(self, maps):
        """This function transforms from inputs to outputs."""
        raise NotImplementedError("Not added.")

    def inverse_transform(self, maps):
        """The inverse conversions of transform. This function transforms from
        outputs to inputs.
        """
        raise NotImplementedError("Not added.")

    def jacobian(self, maps):
        """The Jacobian for the inputs to outputs transformation."""
        raise NotImplementedError("Jacobian transform not implemented.")

    def inverse_jacobian(self, maps):
        """The Jacobian for the outputs to inputs transformation."""
        raise NotImplementedError("Jacobian transform not implemented.")

    @staticmethod
    def format_output(old_maps, new_maps):
        """This function takes the returned dict from `transform` and converts
        it to the same datatype as the input.

        Parameters
        ----------
        old_maps : {FieldArray, dict}
            The mapping object to add new maps to.
        new_maps : dict
            A dict with key as parameter name and value is numpy.array.

        Returns
        -------
        {FieldArray, dict}
            The old_maps object with new keys from new_maps.
        """

        # if input is FieldArray then return FieldArray
        if isinstance(old_maps, record.FieldArray):
            keys = new_maps.keys()
            values = [new_maps[key] for key in keys]
            for key, vals in zip(keys, values):
                try:
                    old_maps = old_maps.add_fields([vals], [key])
                except ValueError:
                    old_maps[key] = vals
            return old_maps

        # if input is dict then return dict
        elif isinstance(old_maps, dict):
            out = old_maps.copy()
            out.update(new_maps)
            return out

        # else error
        else:
            raise TypeError("Input type must be FieldArray or dict.")

    @classmethod
    def from_config(cls, cp, section, outputs,
                    skip_opts=None, additional_opts=None):
        """Initializes a transform from the given section.

        Parameters
        ----------
        cp : pycbc.workflow.WorkflowConfigParser
            A parsed configuration file that contains the transform options.
        section : str
            Name of the section in the configuration file.
        outputs : str
            The names of the parameters that are output by this transformation,
            separated by `VARARGS_DELIM`. These must appear in the "tag" part
            of the section header.
        skip_opts : list, optional
            Do not read options in the given list.
        additional_opts : dict, optional
            Any additional arguments to pass to the class. If an option is
            provided that also exists in the config file, the value provided
            will be used instead of being read from the file.

        Returns
        -------
        cls
            An instance of the class.
        """
        tag = outputs
        if skip_opts is None:
            skip_opts = []
        if additional_opts is None:
            additional_opts = {}
        else:
            additional_opts = additional_opts.copy()
        outputs = set(outputs.split(VARARGS_DELIM))
        special_args = ["name"] + skip_opts + list(additional_opts.keys())
        # get any extra arguments to pass to init
        extra_args = {}
        for opt in cp.options("-".join([section, tag])):
            if opt in special_args:
                continue
            # check if option can be cast as a float
            val = cp.get_opt_tag(section, opt, tag)
            try:
                val = float(val)
            except ValueError:
                pass
            # add option
            extra_args.update({opt: val})
        extra_args.update(additional_opts)
        out = cls(**extra_args)
        # check that the outputs matches
        if outputs - out.outputs != set() or out.outputs - outputs != set():
            raise ValueError(
                "outputs of class do not match outputs specified " "in section"
            )
        return out


class CustomTransform(BaseTransform):
    """Allows for any transform to be defined.

    Parameters
    ----------
    input_args : (list of) str
        The names of the input parameters.
    output_args : (list of) str
        The names of the output parameters.
    transform_functions : dict
        Dictionary mapping input args to a string giving a function call;
        e.g., ``{'q': 'q_from_mass1_mass2(mass1, mass2)'}``.
    jacobian : str, optional
        String giving a jacobian function. The function must be in terms of
        the input arguments.

    Examples
    --------
    Create a custom transform that converts mass1, mass2 to mtotal, q:

    >>> t = transforms.CustomTransform(['mass1', 'mass2'], ['mtotal', 'q'], {'mtotal': 'mass1+mass2', 'q': 'mass1/mass2'}, '(mass1 + mass2) / mass2**2')

    Evaluate a pair of masses:

    >>> t.transform({'mass1': 10., 'mass2': 5.})
    {'mass1': 10.0, 'mass2': 5.0, 'mtotal': 15.0, 'q': 2.0}

    The Jacobian for the same pair of masses:

    >>> t.jacobian({'mass1': 10., 'mass2': 5.})
    0.59999999999999998

    """

    name = "custom"

    def __init__(self, input_args, output_args, transform_functions,
                 jacobian=None):
        if isinstance(input_args, str):
            input_args = [input_args]
        if isinstance(output_args, str):
            output_args = [output_args]
        self.inputs = set(input_args)
        self.outputs = set(output_args)
        self.transform_functions = transform_functions
        self._jacobian = jacobian
        # we'll create a scratch FieldArray space to do transforms on
        # we'll default to length 1; this will be changed if a map is passed
        # with more than one value in it
        self._createscratch()

    def _createscratch(self, shape=1):
        """Creates a scratch FieldArray to use for transforms."""
        self._scratch = record.FieldArray(
            shape, dtype=[(p, float) for p in self.inputs]
        )

    def _copytoscratch(self, maps):
        """Copies the data in maps to the scratch space.

        If the maps contain arrays that are not the same shape as the scratch
        space, a new scratch space will be created.
        """
        try:
            for p in self.inputs:
                self._scratch[p][:] = maps[p]
        except ValueError:
            # we'll get a ValueError if the scratch space isn't the same size
            # as the maps; in that case, re-create the scratch space with the
            # appropriate size and try again
            invals = maps[list(self.inputs)[0]]
            if isinstance(invals, numpy.ndarray):
                shape = invals.shape
            else:
                shape = len(invals)
            self._createscratch(shape)
            for p in self.inputs:
                self._scratch[p][:] = maps[p]

    def _getslice(self, maps):
        """Determines how to slice the scratch for returning values."""
        invals = maps[list(self.inputs)[0]]
        if not isinstance(invals, (numpy.ndarray, list)):
            getslice = 0
        else:
            getslice = slice(None, None)
        return getslice

    def transform(self, maps):
        """Applies the transform functions to the given maps object.

        Parameters
        ----------
        maps : dict, or FieldArray

        Returns
        -------
        dict or FieldArray
            A map object containing the transformed variables, along with the
            original variables. The type of the output will be the same as the
            input.
        """
        if self.transform_functions is None:
            raise NotImplementedError("no transform function(s) provided")
        # copy values to scratch
        self._copytoscratch(maps)
        # ensure that we return the same data type in each dict
        getslice = self._getslice(maps)
        # evaluate the functions
        out = {
            p: self._scratch[func][getslice]
            for p, func in self.transform_functions.items()
        }
        return self.format_output(maps, out)

    def jacobian(self, maps):
        if self._jacobian is None:
            raise NotImplementedError("no jacobian provided")
        # copy values to scratch
        self._copytoscratch(maps)
        out = self._scratch[self._jacobian]
        if isinstance(out, numpy.ndarray):
            out = out[self._getslice(maps)]
        return out

    @classmethod
    def from_config(cls, cp, section, outputs):
        """Loads a CustomTransform from the given config file.

        Example section:

        .. code-block:: ini

            [{section}-outvar1+outvar2]
            name = custom
            inputs = inputvar1, inputvar2
            outvar1 = func1(inputs)
            outvar2 = func2(inputs)
            jacobian = func(inputs)
        """
        tag = outputs
        outputs = set(outputs.split(VARARGS_DELIM))
        inputs = map(str.strip,
                     cp.get_opt_tag(section, "inputs", tag).split(","))
        # get the functions for each output
        transform_functions = {}
        for var in outputs:
            # check if option can be cast as a float
            func = cp.get_opt_tag(section, var, tag)
            transform_functions[var] = func
        s = "-".join([section, tag])
        if cp.has_option(s, "jacobian"):
            jacobian = cp.get_opt_tag(section, "jacobian", tag)
        else:
            jacobian = None
        return cls(inputs, outputs, transform_functions, jacobian=jacobian)


#
# =============================================================================
#
#                             Forward Transforms
#
# =============================================================================
#


class MchirpQToMass1Mass2(BaseTransform):
    """Converts chirp mass and mass ratio to component masses."""

    name = "mchirp_q_to_mass1_mass2"

    def __init__(
        self, mass1_param=None, mass2_param=None, mchirp_param=None, q_param=None
    ):
        if mass1_param is None:
            mass1_param = parameters.mass1
        if mass2_param is None:
            mass2_param = parameters.mass2
        if mchirp_param is None:
            mchirp_param = parameters.mchirp
        if q_param is None:
            q_param = parameters.q
        self.mass1_param = mass1_param
        self.mass2_param = mass2_param
        self.mchirp_param = mchirp_param
        self.q_param = q_param
        self._inputs = [self.mchirp_param, self.q_param]
        self._outputs = [self.mass1_param, self.mass2_param]
        super(MchirpQToMass1Mass2, self).__init__()

    def transform(self, maps):
        """This function transforms from chirp mass and mass ratio to component
        masses.

        Parameters
        ----------
        maps : a mapping object

        Examples
        --------
        Convert a dict of numpy.array:

        >>> import numpy
        >>> from pycbc import transforms
        >>> t = transforms.MchirpQToMass1Mass2()
        >>> t.transform({'mchirp': numpy.array([10.]), 'q': numpy.array([2.])})
        {'mass1': array([ 16.4375183]), 'mass2': array([ 8.21875915]),
         'mchirp': array([ 10.]), 'q': array([ 2.])}

        Returns
        -------
        out : dict
            A dict with key as parameter name and value as numpy.array or float
            of transformed values.
        """
        out = {}
        out[self.mass1_param] = conversions.mass1_from_mchirp_q(
            maps[self.mchirp_param], maps[self.q_param]
        )
        out[self.mass2_param] = conversions.mass2_from_mchirp_q(
            maps[self.mchirp_param], maps[self.q_param]
        )
        return self.format_output(maps, out)

    def inverse_transform(self, maps):
        """This function transforms from component masses to chirp mass and
        mass ratio.

        Parameters
        ----------
        maps : a mapping object

        Examples
        --------
        Convert a dict of numpy.array:

        >>> import numpy
        >>> from pycbc import transforms
        >>> t = transforms.MchirpQToMass1Mass2()
        >>> t.inverse_transform({'mass1': numpy.array([16.4]), 'mass2': numpy.array([8.2])})
            {'mass1': array([ 16.4]), 'mass2': array([ 8.2]),
             'mchirp': array([ 9.97717521]), 'q': 2.0}

        Returns
        -------
        out : dict
            A dict with key as parameter name and value as numpy.array or float
            of transformed values.
        """
        out = {}
        m1 = maps[self.mass1_param]
        m2 = maps[self.mass2_param]
        out[self.mchirp_param] = conversions.mchirp_from_mass1_mass2(m1, m2)
        out[self.q_param] = m1 / m2
        return self.format_output(maps, out)

    def jacobian(self, maps):
        """Returns the Jacobian for transforming mchirp and q to mass1 and
        mass2.
        """
        mchirp = maps[self.mchirp_param]
        q = maps[self.q_param]
        return mchirp * ((1.0 + q) / q ** 3.0) ** (2.0 / 5)

    def inverse_jacobian(self, maps):
        """Returns the Jacobian for transforming mass1 and mass2 to
        mchirp and q.
        """
        m1 = maps[self.mass1_param]
        m2 = maps[self.mass2_param]
        return conversions.mchirp_from_mass1_mass2(m1, m2) / m2 ** 2.0


class MchirpEtaToMass1Mass2(BaseTransform):
    """Converts chirp mass and symmetric mass ratio to component masses."""

    name = "mchirp_eta_to_mass1_mass2"
    _inputs = [parameters.mchirp, parameters.eta]
    _outputs = [parameters.mass1, parameters.mass2]

    def transform(self, maps):
        """This function transforms from chirp mass and symmetric mass ratio to
        component masses.

        Parameters
        ----------
        maps : a mapping object

        Examples
        --------
        Convert a dict of numpy.array:

        >>> import numpy
        >>> from pycbc import transforms
        >>> t = transforms.MchirpEtaToMass1Mass2()
        >>> t.transform({'mchirp': numpy.array([10.]), 'eta': numpy.array([0.25])})
        {'mass1': array([ 16.4375183]), 'mass2': array([ 8.21875915]),
         'mchirp': array([ 10.]), 'eta': array([ 0.25])}

        Returns
        -------
        out : dict
            A dict with key as parameter name and value as numpy.array or float
            of transformed values.
        """
        out = {}
        out[parameters.mass1] = conversions.mass1_from_mchirp_eta(
            maps[parameters.mchirp], maps[parameters.eta]
        )
        out[parameters.mass2] = conversions.mass2_from_mchirp_eta(
            maps[parameters.mchirp], maps[parameters.eta]
        )
        return self.format_output(maps, out)

    def inverse_transform(self, maps):
        """This function transforms from component masses to chirp mass and
        symmetric mass ratio.

        Parameters
        ----------
        maps : a mapping object

        Examples
        --------
        Convert a dict of numpy.array:

        >>> import numpy
        >>> from pycbc import transforms
        >>> t = transforms.MchirpQToMass1Mass2()
        >>> t.inverse_transform({'mass1': numpy.array([8.2]), 'mass2': numpy.array([8.2])})
            {'mass1': array([ 8.2]), 'mass2': array([ 8.2]),
             'mchirp': array([ 9.97717521]), 'eta': 0.25}

        Returns
        -------
        out : dict
            A dict with key as parameter name and value as numpy.array or float
            of transformed values.
        """
        out = {}
        m1 = maps[parameters.mass1]
        m2 = maps[parameters.mass2]
        out[parameters.mchirp] = conversions.mchirp_from_mass1_mass2(m1, m2)
        out[parameters.eta] = conversions.eta_from_mass1_mass2(m1, m2)
        return self.format_output(maps, out)

    def jacobian(self, maps):
        """Returns the Jacobian for transforming mchirp and eta to mass1 and
        mass2.
        """
        mchirp = maps[parameters.mchirp]
        eta = maps[parameters.eta]
        m1 = conversions.mass1_from_mchirp_eta(mchirp, eta)
        m2 = conversions.mass2_from_mchirp_eta(mchirp, eta)
        return mchirp * (m1 - m2) / (m1 + m2) ** 3

    def inverse_jacobian(self, maps):
        """Returns the Jacobian for transforming mass1 and mass2 to
        mchirp and eta.
        """
        m1 = maps[parameters.mass1]
        m2 = maps[parameters.mass2]
        mchirp = conversions.mchirp_from_mass1_mass2(m1, m2)
        eta = conversions.eta_from_mass1_mass2(m1, m2)
        return -1.0 * mchirp / eta ** (6.0 / 5)


class ChirpDistanceToDistance(BaseTransform):
    """Converts chirp distance to luminosity distance, given the chirp mass."""

    name = "chirp_distance_to_distance"
    _inputs = [parameters.chirp_distance, parameters.mchirp]
    _outputs = [parameters.distance]

    def __init__(self, ref_mass=1.4):
        self.inputs = set(self._inputs)
        self.outputs = set(self._outputs)
        self.ref_mass = ref_mass

    def transform(self, maps):
        """This function transforms from chirp distance to luminosity distance,
        given the chirp mass.

        Parameters
        ----------
        maps : a mapping object

        Examples
        --------
        Convert a dict of numpy.array:

        >>> import numpy as np
        >>> from pycbc import transforms
        >>> t = transforms.ChirpDistanceToDistance()
        >>> t.transform({'chirp_distance': np.array([40.]), 'mchirp': np.array([1.2])})
        {'mchirp': array([ 1.2]), 'chirp_distance': array([ 40.]), 'distance': array([ 39.48595679])}

        Returns
        -------
        out : dict
            A dict with key as parameter name and value as numpy.array or float
            of transformed values.
        """
        out = {}
        out[parameters.distance] = conversions.distance_from_chirp_distance_mchirp(
            maps[parameters.chirp_distance],
            maps[parameters.mchirp],
            ref_mass=self.ref_mass,
        )
        return self.format_output(maps, out)

    def inverse_transform(self, maps):
        """This function transforms from luminosity distance to chirp distance,
        given the chirp mass.

        Parameters
        ----------
        maps : a mapping object

        Examples
        --------
        Convert a dict of numpy.array:

        >>> import numpy as np
        >>> from pycbc import transforms
        >>> t = transforms.ChirpDistanceToDistance()
        >>> t.inverse_transform({'distance': np.array([40.]), 'mchirp': np.array([1.2])})
        {'distance': array([ 40.]), 'chirp_distance': array([ 40.52073522]), 'mchirp': array([ 1.2])}

        Returns
        -------
        out : dict
            A dict with key as parameter name and value as numpy.array or float
            of transformed values.
        """
        out = {}
        out[parameters.chirp_distance] = conversions.chirp_distance(
            maps[parameters.distance], maps[parameters.mchirp], ref_mass=self.ref_mass
        )
        return self.format_output(maps, out)

    def jacobian(self, maps):
        """Returns the Jacobian for transforming chirp distance to
        luminosity distance, given the chirp mass.
        """
        ref_mass = 1.4
        mchirp = maps["mchirp"]
        return (2.0 ** (-1.0 / 5) * self.ref_mass / mchirp) ** (-5.0 / 6)

    def inverse_jacobian(self, maps):
        """Returns the Jacobian for transforming luminosity distance to
        chirp distance, given the chirp mass.
        """
        ref_mass = 1.4
        mchirp = maps["mchirp"]
        return (2.0 ** (-1.0 / 5) * self.ref_mass / mchirp) ** (5.0 / 6)


class AlignTotalSpin(BaseTransform):
    """Converts angles from total angular momentum J frame to orbital angular
     momentum L (waveform) frame"""

    name = "align_total_spin"
    _inputs = [parameters.thetajn, parameters.spin1x, parameters.spin1y,
               parameters.spin1z, parameters.spin2x, parameters.spin2y,
               parameters.spin2z, parameters.mass1, parameters.mass2,
               parameters.f_ref, "phi_ref"]
    _outputs = [parameters.inclination, parameters.spin1x, parameters.spin1y,
               parameters.spin1z, parameters.spin2x, parameters.spin2y,
               parameters.spin2z]

    def __init__(self):
        self.inputs = set(self._inputs)
        self.outputs = set(self._outputs)
        super(AlignTotalSpin, self).__init__()

    def transform(self, maps):
        """
        Rigidly rotate binary so that the total angular momentum has the given
        inclination (iota) instead of the orbital angular momentum. Return
        the new inclination, s1, and s2. s1 and s2 are dimensionless spin.
        Note: the spins are assumed to be given in the frame defined by the
        orbital angular momentum.
        """

        if isinstance(maps, dict):
            maps = record.FieldArray.from_kwargs(**maps)
        newfields = [n for n in self._outputs if n not in maps.fieldnames]
        newmaps = maps.add_fields([numpy.zeros(len(maps))]*len(newfields),
                                  names=newfields)
        for item in newmaps:
            if not all(s == 0.0 for s in
                       [item[parameters.spin1x], item[parameters.spin1y],
                        item[parameters.spin2x], item[parameters.spin2y]]):

                # Calculate the quantities required by jframe_to_l0frame
                s1_a, s1_az, s1_pol = coordinates.cartesian_to_spherical(
                        item[parameters.spin1x], item[parameters.spin1y],
                        item[parameters.spin1z])
                s2_a, s2_az, s2_pol = coordinates.cartesian_to_spherical(
                        item[parameters.spin2x], item[parameters.spin2y],
                        item[parameters.spin2z])

                out = jframe_to_l0frame(
                    item[parameters.mass1],
                    item[parameters.mass2],
                    item[parameters.f_ref],
                    phiref=item["phi_ref"],
                    thetajn=item[parameters.thetajn],
                    phijl=numpy.pi,
                    spin1_a=s1_a,
                    spin2_a=s2_a,
                    spin1_polar=s1_pol,
                    spin2_polar=s2_pol,
                    spin12_deltaphi=s1_az-s2_az
                )

                for key in out:
                    item[key] = out[key]
            else:
                item[parameters.inclination] = item[parameters.thetajn]

        return newmaps


class SphericalToCartesian(BaseTransform):
    """Converts spherical coordinates to cartesian.

    Parameters
    ----------
    x : str
        The name of the x parameter.
    y : str
        The name of the y parameter.
    z : str
        The name of the z parameter.
    radial : str
        The name of the radial parameter.
    azimuthal : str
        The name of the azimuthal angle parameter.
    polar : str
        The name of the polar angle parameter.
    """

    name = "spherical_to_cartesian"

    def __init__(self, x, y, z, radial, azimuthal, polar):
        self.x = x
        self.y = y
        self.z = z
        self.radial = radial
        self.polar = polar
        self.azimuthal = azimuthal
        self._inputs = [self.radial, self.azimuthal, self.polar]
        self._outputs = [self.x, self.y, self.z]
        super(SphericalToCartesian, self).__init__()

    def transform(self, maps):
        """This function transforms from spherical to cartesian spins.

        Parameters
        ----------
        maps : a mapping object

        Examples
        --------
        Convert a dict of numpy.array:

        >>> import numpy
        >>> from pycbc import transforms
        >>> t = transforms.SphericalToCartesian('x', 'y', 'z',
                                                'a', 'phi', 'theta')
        >>> t.transform({'a': numpy.array([0.1]), 'phi': numpy.array([0.1]),
                        'theta': numpy.array([0.1])})
            {'a': array([ 0.1]), 'phi': array([ 0.1]), 'theta': array([ 0.1]),
             'x': array([ 0.00993347]), 'y': array([ 0.00099667]),
             'z': array([ 0.09950042])}

        Returns
        -------
        out : dict
            A dict with key as parameter name and value as numpy.array or float
            of transformed values.
        """
        a = self.radial
        az = self.azimuthal
        po = self.polar
        x, y, z = coordinates.spherical_to_cartesian(maps[a], maps[az], maps[po])
        out = {self.x: x, self.y: y, self.z: z}
        return self.format_output(maps, out)

    def inverse_transform(self, maps):
        """This function transforms from cartesian to spherical spins.

        Parameters
        ----------
        maps : a mapping object

        Returns
        -------
        out : dict
            A dict with key as parameter name and value as numpy.array or float
            of transformed values.
        """
        x = self.x
        y = self.y
        z = self.z
        a, az, po = coordinates.cartesian_to_spherical(maps[x], maps[y], maps[z])
        out = {self.radial: a, self.azimuthal: az, self.polar: po}
        return self.format_output(maps, out)


class SphericalSpin1ToCartesianSpin1(SphericalToCartesian):
    """Converts spherical spin parameters (radial and two angles) to
    catesian spin parameters. This class only transforms spsins for the first
    component mass.

    **Deprecation Warning:** This will be removed in a future update. Use
    :py:class:`SphericalToCartesian` with spin-parameter names passed in
    instead.
    """

    name = "spherical_spin_1_to_cartesian_spin_1"

    def __init__(self):
        logging.warning(
            "Deprecation warning: the {} transform will be "
            "removed in a future update. Please use {} instead, "
            "passing spin1x, spin1y, spin1z, spin1_a, "
            "spin1_azimuthal, spin1_polar as arguments.".format(
                self.name, SphericalToCartesian.name
            )
        )
        super(SphericalSpin1ToCartesianSpin1, self).__init__(
            "spin1x", "spin1y", "spin1z", "spin1_a",
            "spin1_azimuthal", "spin1_polar"
        )


class SphericalSpin2ToCartesianSpin2(SphericalToCartesian):
    """Converts spherical spin parameters (radial and two angles) to
    catesian spin parameters. This class only transforms spsins for the first
    component mass.

    **Deprecation Warning:** This will be removed in a future update. Use
    :py:class:`SphericalToCartesian` with spin-parameter names passed in
    instead.
    """

    name = "spherical_spin_2_to_cartesian_spin_2"

    def __init__(self):
        logging.warning(
            "Deprecation warning: the {} transform will be "
            "removed in a future update. Please use {} instead, "
            "passing spin2x, spin2y, spin2z, spin2_a, "
            "spin2_azimuthal, spin2_polar as arguments.".format(
                self.name, SphericalToCartesian.name
            )
        )
        super(SphericalSpin2ToCartesianSpin2, self).__init__(
            "spin2x", "spin2y", "spin2z",
            "spin2_a", "spin2_azimuthal", "spin2_polar"
        )


class DistanceToRedshift(BaseTransform):
    """Converts distance to redshift."""

    name = "distance_to_redshift"
    inverse = None
    _inputs = [parameters.distance]
    _outputs = [parameters.redshift]

    def transform(self, maps):
        """This function transforms from distance to redshift.

        Parameters
        ----------
        maps : a mapping object

        Examples
        --------
        Convert a dict of numpy.array:

        >>> import numpy
        >>> from pycbc import transforms
        >>> t = transforms.DistanceToRedshift()
        >>> t.transform({'distance': numpy.array([1000])})
            {'distance': array([1000]), 'redshift': 0.19650987609144363}

        Returns
        -------
        out : dict
            A dict with key as parameter name and value as numpy.array or float
            of transformed values.
        """
        out = {parameters.redshift: cosmology.redshift(maps[parameters.distance])}
        return self.format_output(maps, out)


class AlignedMassSpinToCartesianSpin(BaseTransform):
    """Converts mass-weighted spins to cartesian z-axis spins."""

    name = "aligned_mass_spin_to_cartesian_spin"
    _inputs = [parameters.mass1, parameters.mass2, parameters.chi_eff, "chi_a"]
    _outputs = [
        parameters.mass1,
        parameters.mass2,
        parameters.spin1z,
        parameters.spin2z,
    ]

    def transform(self, maps):
        """This function transforms from aligned mass-weighted spins to
        cartesian spins aligned along the z-axis.

        Parameters
        ----------
        maps : a mapping object

        Returns
        -------
        out : dict
            A dict with key as parameter name and value as numpy.array or float
            of transformed values.
        """
        mass1 = maps[parameters.mass1]
        mass2 = maps[parameters.mass2]
        out = {}
        out[parameters.spin1z] = conversions.spin1z_from_mass1_mass2_chi_eff_chi_a(
            mass1, mass2, maps[parameters.chi_eff], maps["chi_a"]
        )
        out[parameters.spin2z] = conversions.spin2z_from_mass1_mass2_chi_eff_chi_a(
            mass1, mass2, maps[parameters.chi_eff], maps["chi_a"]
        )
        return self.format_output(maps, out)

    def inverse_transform(self, maps):
        """This function transforms from component masses and cartesian spins
        to mass-weighted spin parameters aligned with the angular momentum.

        Parameters
        ----------
        maps : a mapping object

        Returns
        -------
        out : dict
            A dict with key as parameter name and value as numpy.array or float
            of transformed values.
        """
        mass1 = maps[parameters.mass1]
        spin1z = maps[parameters.spin1z]
        mass2 = maps[parameters.mass2]
        spin2z = maps[parameters.spin2z]
        out = {
            parameters.chi_eff:
            conversions.chi_eff(mass1, mass2, spin1z, spin2z),
            "chi_a": conversions.chi_a(mass1, mass2, spin1z, spin2z),
        }
        return self.format_output(maps, out)


class PrecessionMassSpinToCartesianSpin(BaseTransform):
    """Converts mass-weighted spins to cartesian x-y plane spins."""

    name = "precession_mass_spin_to_cartesian_spin"
    _inputs = [parameters.mass1, parameters.mass2,
               "xi1", "xi2", "phi_a", "phi_s"]
    _outputs = [
        parameters.mass1,
        parameters.mass2,
        parameters.spin1x,
        parameters.spin1y,
        parameters.spin2x,
        parameters.spin2y,
    ]

    def transform(self, maps):
        """This function transforms from mass-weighted spins to caretsian spins
        in the x-y plane.

        Parameters
        ----------
        maps : a mapping object

        Returns
        -------
        out : dict
            A dict with key as parameter name and value as numpy.array or float
            of transformed values.
        """

        # find primary and secondary masses
        # since functions in conversions.py map to primary/secondary masses
        m_p = conversions.primary_mass(maps["mass1"], maps["mass2"])
        m_s = conversions.secondary_mass(maps["mass1"], maps["mass2"])

        # find primary and secondary xi
        # can re-purpose spin functions for just a generic variable
        xi_p = conversions.primary_spin(
            maps["mass1"], maps["mass2"], maps["xi1"], maps["xi2"]
        )
        xi_s = conversions.secondary_spin(
            maps["mass1"], maps["mass2"], maps["xi1"], maps["xi2"]
        )

        # convert using convention of conversions.py that is mass1 > mass2
        spinx_p = conversions.spin1x_from_xi1_phi_a_phi_s(
            xi_p, maps["phi_a"], maps["phi_s"]
        )
        spiny_p = conversions.spin1y_from_xi1_phi_a_phi_s(
            xi_p, maps["phi_a"], maps["phi_s"]
        )
        spinx_s = conversions.spin2x_from_mass1_mass2_xi2_phi_a_phi_s(
            m_p, m_s, xi_s, maps["phi_a"], maps["phi_s"]
        )
        spiny_s = conversions.spin2y_from_mass1_mass2_xi2_phi_a_phi_s(
            m_p, m_s, xi_s, maps["phi_a"], maps["phi_s"]
        )

        # map parameters from primary/secondary to indices
        out = {}
        if isinstance(m_p, numpy.ndarray):
            mass1, mass2 = map(numpy.array, [maps["mass1"], maps["mass2"]])
            mask_mass1_gte_mass2 = mass1 >= mass2
            mask_mass1_lt_mass2 = mass1 < mass2
            out[parameters.spin1x] = numpy.concatenate(
                (spinx_p[mask_mass1_gte_mass2], spinx_s[mask_mass1_lt_mass2])
            )
            out[parameters.spin1y] = numpy.concatenate(
                (spiny_p[mask_mass1_gte_mass2], spiny_s[mask_mass1_lt_mass2])
            )
            out[parameters.spin2x] = numpy.concatenate(
                (spinx_p[mask_mass1_lt_mass2], spinx_s[mask_mass1_gte_mass2])
            )
            out[parameters.spin2y] = numpy.concatenate(
                (spinx_p[mask_mass1_lt_mass2], spinx_s[mask_mass1_gte_mass2])
            )
        elif maps["mass1"] > maps["mass2"]:
            out[parameters.spin1x] = spinx_p
            out[parameters.spin1y] = spiny_p
            out[parameters.spin2x] = spinx_s
            out[parameters.spin2y] = spiny_s
        else:
            out[parameters.spin1x] = spinx_s
            out[parameters.spin1y] = spiny_s
            out[parameters.spin2x] = spinx_p
            out[parameters.spin2y] = spiny_p

        return self.format_output(maps, out)

    def inverse_transform(self, maps):
        """This function transforms from component masses and cartesian spins to
        mass-weighted spin parameters perpendicular with the angular momentum.

        Parameters
        ----------
        maps : a mapping object

        Returns
        -------
        out : dict
            A dict with key as parameter name and value as numpy.array or float
            of transformed values.
        """

        # convert
        out = {}
        xi1 = conversions.primary_xi(
            maps[parameters.mass1],
            maps[parameters.mass2],
            maps[parameters.spin1x],
            maps[parameters.spin1y],
            maps[parameters.spin2x],
            maps[parameters.spin2y],
        )
        xi2 = conversions.secondary_xi(
            maps[parameters.mass1],
            maps[parameters.mass2],
            maps[parameters.spin1x],
            maps[parameters.spin1y],
            maps[parameters.spin2x],
            maps[parameters.spin2y],
        )
        out["phi_a"] = conversions.phi_a(
            maps[parameters.mass1],
            maps[parameters.mass2],
            maps[parameters.spin1x],
            maps[parameters.spin1y],
            maps[parameters.spin2x],
            maps[parameters.spin2y],
        )
        out["phi_s"] = conversions.phi_s(
            maps[parameters.spin1x],
            maps[parameters.spin1y],
            maps[parameters.spin2x],
            maps[parameters.spin2y],
        )

        # map parameters from primary/secondary to indices
        if isinstance(xi1, numpy.ndarray):
            mass1, mass2 = map(
                numpy.array, [maps[parameters.mass1], maps[parameters.mass2]]
            )
            mask_mass1_gte_mass2 = mass1 >= mass2
            mask_mass1_lt_mass2 = mass1 < mass2
            out["xi1"] = numpy.concatenate(
                (xi1[mask_mass1_gte_mass2], xi2[mask_mass1_lt_mass2])
            )
            out["xi2"] = numpy.concatenate(
                (xi1[mask_mass1_gte_mass2], xi2[mask_mass1_lt_mass2])
            )
        elif maps["mass1"] > maps["mass2"]:
            out["xi1"] = xi1
            out["xi2"] = xi2
        else:
            out["xi1"] = xi2
            out["xi2"] = xi1

        return self.format_output(maps, out)


class CartesianSpinToChiP(BaseTransform):
    """Converts cartesian spins to `chi_p`."""

    name = "cartesian_spin_to_chi_p"
    _inputs = [
        parameters.mass1,
        parameters.mass2,
        parameters.spin1x,
        parameters.spin1y,
        parameters.spin2x,
        parameters.spin2y,
    ]
    _outputs = ["chi_p"]

    def transform(self, maps):
        """This function transforms from component masses and caretsian spins
        to chi_p.

        Parameters
        ----------
        maps : a mapping object

        Examples
        --------
        Convert a dict of numpy.array:

        Returns
        -------
        out : dict
            A dict with key as parameter name and value as numpy.array or float
            of transformed values.
        """
        out = {}
        out["chi_p"] = conversions.chi_p(
            maps[parameters.mass1],
            maps[parameters.mass2],
            maps[parameters.spin1x],
            maps[parameters.spin1y],
            maps[parameters.spin2x],
            maps[parameters.spin2y],
        )
        return self.format_output(maps, out)


class LambdaFromTOVFile(BaseTransform):
    """Transforms mass values corresponding to Lambda values for a given EOS
    interpolating from the mass-Lambda data for that EOS read in from an
    external ASCII file.

    The interpolation of the mass-Lambda data is a one-dimensional piecewise
    linear interpolation. If the ``redshift_mass`` keyword argument is ``True``
    (the default), the mass values to be transformed are assumed to be detector
    frame masses. In that case, a distance should be provided along with the
    mass for transformation to the source frame mass before the Lambda values
    are extracted from the interpolation. If the transform is read in from a
    config file, an example code block would be:

    .. code-block:: ini

        [{section}-lambda1]
        name = lambda_from_tov_file
        mass_param = mass1
        lambda_param = lambda1
        distance = 40
        mass_lambda_file = {filepath}

    If this transform is used in a parameter estimation analysis where
    distance is a variable parameter, the distance to be used will vary
    with each draw. In that case, the example code block will be:

    .. code-block:: ini

        [{section}-lambda1]
        name = lambda_from_tov_file
        mass_param = mass1
        lambda_param = lambda1
        mass_lambda_file = filepath

    If your prior is in terms of the source-frame masses (``srcmass``), then
    you can shut off the redshifting by adding ``do-not-redshift-mass`` to the
    config file. In this case you do not need to worry about a distance.
    Example:

    .. code-block:: ini

        [{section}-lambda1]
        name = lambda_from_tov_file
        mass_param = srcmass1
        lambda_param = lambda1
        mass_lambda_file = filepath
        do-not-redshift-mass =

    Parameters
    ----------
    mass_param : str
        The name of the mass parameter to transform.
    lambda_param : str
        The name of the tidal deformability parameter that mass_param is to
        be converted to interpolating from the data in the mass-Lambda file.
    mass_lambda_file : str
        Path of the mass-Lambda data file. The first column in the data file
        should contain mass values, and the second column Lambda values.
    distance : float, optional
        The distance (in Mpc) of the source. Used to redshift the mass. Needed
        if ``redshift_mass`` is True and no distance parameter exists If
        None, then a distance must be provided to the transform.
    redshift_mass : bool, optional
        Redshift the mass parameters when computing the lambdas. Default is
        False.
    file_columns : list of str, optional
        The names and order of columns in the ``mass_lambda_file``. Must
        contain at least 'mass' and 'lambda'. If not provided, will assume the
        order is ('mass', 'lambda').
    """

    name = "lambda_from_tov_file"

    def __init__(
        self,
        mass_param,
        lambda_param,
        mass_lambda_file,
        distance=None,
        redshift_mass=True,
        file_columns=None,
    ):
        self._mass_lambda_file = mass_lambda_file
        self._mass_param = mass_param
        self._lambda_param = lambda_param
        self.redshift_mass = redshift_mass
        self._distance = distance
        self._inputs = [mass_param, "distance"]
        self._outputs = [lambda_param]
        if file_columns is None:
            file_columns = ["mass", "lambda"]
        dtype = [(fname, float) for fname in file_columns]
        data = numpy.loadtxt(self._mass_lambda_file, dtype=dtype)
        self._data = data
        super(LambdaFromTOVFile, self).__init__()

    @property
    def mass_param(self):
        """Returns the input mass parameter."""
        return self._mass_param

    @property
    def lambda_param(self):
        """Returns the output lambda parameter."""
        return self._lambda_param

    @property
    def data(self):
        return self._data

    @property
    def mass_data(self):
        """Returns the mass data read from the mass-Lambda data file for
        an EOS.
        """
        return self._data["mass"]

    @property
    def lambda_data(self):
        """Returns the Lambda data read from the mass-Lambda data file for
        an EOS.
        """
        return self._data["lambda"]

    @property
    def distance(self):
        """Returns the fixed distance to transform mass samples from detector
        to source frame if one is specified.
        """
        return self._distance

    @staticmethod
    def lambda_from_tov_data(m_src, mass_data, lambda_data):
        """Returns Lambda corresponding to a given mass interpolating from the
        TOV data.

        Parameters
        ----------
        m : float
            Value of the mass.
        mass_data : array
            Mass array from the Lambda-M curve of an EOS.
        lambda_data : array
            Lambda array from the Lambda-M curve of an EOS.

        Returns
        -------
        lambdav : float
            The Lambda corresponding to the mass `m` for the EOS considered.
        """
        if m_src > mass_data.max():
            # assume black hole
            lambdav = 0.0
        else:
            lambdav = numpy.interp(m_src, mass_data, lambda_data)
        return lambdav

    def transform(self, maps):
        """Computes the transformation of mass to Lambda.

        Parameters
        ----------
        maps : dict or FieldArray
            A dictionary or FieldArray which provides a map between the
            parameter name of the variable to transform and its value(s).

        Returns
        -------
        out : dict or FieldArray
            A map between the transformed variable name and value(s), along
            with the original variable name and value(s).
        """
        m = maps[self._mass_param]
        if self.redshift_mass:
            if self._distance is not None:
                d = self._distance
            else:
                try:
                    d = maps["distance"]
                except KeyError as e:
                    logging.warning(
                        "Either provide distance samples in the "
                        "list of samples to be transformed, or "
                        "provide a fixed distance value as input "
                        "when initializing LambdaFromTOVFile."
                    )
                    raise e
            shift = 1.0 / (1.0 + cosmology.redshift(abs(d)))
        else:
            shift = 1.0
        out = {
            self._lambda_param: self.lambda_from_tov_data(
                m * shift, self._data["mass"], self._data["lambda"]
            )
        }
        return self.format_output(maps, out)

    @classmethod
    def from_config(cls, cp, section, outputs):
        # see if we're redshifting masses
        if cp.has_option("-".join([section, outputs]), "do-not-redshift-mass"):
            additional_opts = {"redshift_mass": False}
            skip_opts = ["do-not-redshift-mass"]
        else:
            additional_opts = None
            skip_opts = None
        return super(LambdaFromTOVFile, cls).from_config(
            cp, section, outputs, skip_opts=skip_opts, additional_opts=additional_opts
        )


class LambdaFromMultipleTOVFiles(BaseTransform):
    """Uses multiple equation of states.

    Parameters
    ----------
    mass_param : str
        The name of the mass parameter to transform.
    lambda_param : str
        The name of the tidal deformability parameter that mass_param is to
        be converted to interpolating from the data in the mass-Lambda file.
    mass_lambda_file : str
        Path of the mass-Lambda data file. The first column in the data file
        should contain mass values, and the second column Lambda values.
    distance : float, optional
        The distance (in Mpc) of the source. Used to redshift the mass. If
        None, then a distance must be provided to the transform.
    file_columns : list of str, optional
        The names and order of columns in the ``mass_lambda_file``. Must
        contain at least 'mass' and 'lambda'. If not provided, will assume the
        order is ('radius', 'mass', 'lambda').
    """

    name = "lambda_from_multiple_tov_files"

    def __init__(
        self,
        mass_param,
        lambda_param,
        map_file,
        distance=None,
        redshift_mass=True,
        file_columns=None,
    ):
        self._map_file = map_file
        self._mass_param = mass_param
        self._lambda_param = lambda_param
        self._distance = distance
        self.redshift_mass = redshift_mass
        self._inputs = [mass_param, "eos", "distance"]
        self._outputs = [lambda_param]
        # create a dictionary of the EOS files from the map_file
        self._eos_files = {}
        with open(self._map_file, "r") as fp:
            for line in fp:
                fname = line.rstrip("\n")
                eosidx = int(os.path.basename(fname).split(".")[0])
                self._eos_files[eosidx] = os.path.abspath(fname)
        # create an eos cache for fast load later
        self._eos_cache = {}
        if file_columns is None:
            file_columns = ("radius", "mass", "lambda")
        self._file_columns = file_columns
        super(LambdaFromMultipleTOVFiles, self).__init__()

    @property
    def mass_param(self):
        """Returns the input mass parameter."""
        return self._mass_param

    @property
    def lambda_param(self):
        """Returns the output lambda parameter."""
        return self._lambda_param

    @property
    def map_file(self):
        """Returns the mass data read from the mass-Lambda data file for
        an EOS.
        """
        return self._map_file

    @property
    def distance(self):
        """Returns the fixed distance to transform mass samples from detector
        to source frame if one is specified.
        """
        return self._distance

    def get_eos(self, eos_index):
        """Gets the EOS for the given index.

        If the index is not in range returns None.
        """
        try:
            eos = self._eos_cache[eos_index]
        except KeyError:
            try:
                fname = self._eos_files[eos_index]
                eos = LambdaFromTOVFile(
                    mass_param=self._mass_param,
                    lambda_param=self._lambda_param,
                    mass_lambda_file=fname,
                    distance=self._distance,
                    redshift_mass=self.redshift_mass,
                    file_columns=self._file_columns,
                )
                self._eos_cache[eos_index] = eos
            except KeyError:
                eos = None
        return eos

    def transform(self, maps):
        """Transforms mass value and eos index into a lambda value"""
        m = maps[self._mass_param]
        # floor
        eos_index = int(maps["eos"])
        eos = self.get_eos(eos_index)
        if eos is not None:
            return eos.transform(maps)
        else:
            # no eos, just return nan
            out = {self._lambda_param: numpy.nan}
            return self.format_output(maps, out)

    @classmethod
    def from_config(cls, cp, section, outputs):
        # see if we're redshifting masses
        if cp.has_option("-".join([section, outputs]), "do-not-redshift-mass"):
            additional_opts = {"redshift_mass": False}
            skip_opts = ["do-not-redshift-mass"]
        else:
            additional_opts = None
            skip_opts = None
        return super(LambdaFromMultipleTOVFiles, cls).from_config(
            cp, section, outputs, skip_opts=skip_opts, additional_opts=additional_opts
        )


class Log(BaseTransform):
    """Applies a log transform from an `inputvar` parameter to an `outputvar`
    parameter. This is the inverse of the exponent transform.

    Parameters
    ----------
    inputvar : str
        The name of the parameter to transform.
    outputvar : str
        The name of the transformed parameter.
    """

    name = "log"

    def __init__(self, inputvar, outputvar):
        self._inputvar = inputvar
        self._outputvar = outputvar
        self._inputs = [inputvar]
        self._outputs = [outputvar]
        super(Log, self).__init__()

    @property
    def inputvar(self):
        """Returns the input parameter."""
        return self._inputvar

    @property
    def outputvar(self):
        """Returns the output parameter."""
        return self._outputvar

    def transform(self, maps):
        r"""Computes :math:`\log(x)`.

        Parameters
        ----------
        maps : dict or FieldArray
            A dictionary or FieldArray which provides a map between the
            parameter name of the variable to transform and its value(s).

        Returns
        -------
        out : dict or FieldArray
            A map between the transformed variable name and value(s), along
            with the original variable name and value(s).
        """
        x = maps[self._inputvar]
        out = {self._outputvar: numpy.log(x)}
        return self.format_output(maps, out)

    def inverse_transform(self, maps):
        r"""Computes :math:`y = e^{x}`.

        Parameters
        ----------
        maps : dict or FieldArray
            A dictionary or FieldArray which provides a map between the
            parameter name of the variable to transform and its value(s).

        Returns
        -------
        out : dict or FieldArray
            A map between the transformed variable name and value(s), along
            with the original variable name and value(s).
        """
        y = maps[self._outputvar]
        out = {self._inputvar: numpy.exp(y)}
        return self.format_output(maps, out)

    def jacobian(self, maps):
        r"""Computes the Jacobian of :math:`y = \log(x)`.

        This is:

        .. math::

            \frac{\mathrm{d}y}{\mathrm{d}x} = \frac{1}{x}.

        Parameters
        ----------
        maps : dict or FieldArray
            A dictionary or FieldArray which provides a map between the
            parameter name of the variable to transform and its value(s).

        Returns
        -------
        float
            The value of the jacobian at the given point(s).
        """
        x = maps[self._inputvar]
        return 1.0 / x

    def inverse_jacobian(self, maps):
        r"""Computes the Jacobian of :math:`y = e^{x}`.

        This is:

        .. math::

            \frac{\mathrm{d}y}{\mathrm{d}x} = e^{x}.

        Parameters
        ----------
        maps : dict or FieldArray
            A dictionary or FieldArray which provides a map between the
            parameter name of the variable to transform and its value(s).

        Returns
        -------
        float
            The value of the jacobian at the given point(s).
        """
        x = maps[self._outputvar]
        return numpy.exp(x)


class Logit(BaseTransform):
    """Applies a logit transform from an `inputvar` parameter to an `outputvar`
    parameter. This is the inverse of the logistic transform.

    Typically, the input of the logit function is assumed to have domain
    :math:`\in (0, 1)`. However, the `domain` argument can be used to expand
    this to any finite real interval.

    Parameters
    ----------
    inputvar : str
        The name of the parameter to transform.
    outputvar : str
        The name of the transformed parameter.
    domain : tuple or distributions.bounds.Bounds, optional
        The domain of the input parameter. Can be any finite
        interval. Default is (0., 1.).
    """

    name = "logit"

    def __init__(self, inputvar, outputvar, domain=(0.0, 1.0)):
        self._inputvar = inputvar
        self._outputvar = outputvar
        self._inputs = [inputvar]
        self._outputs = [outputvar]
        self._bounds = Bounds(domain[0], domain[1],
                              btype_min="open", btype_max="open")
        # shortcuts for quick access later
        self._a = domain[0]
        self._b = domain[1]
        super(Logit, self).__init__()

    @property
    def inputvar(self):
        """Returns the input parameter."""
        return self._inputvar

    @property
    def outputvar(self):
        """Returns the output parameter."""
        return self._outputvar

    @property
    def bounds(self):
        """Returns the domain of the input parameter."""
        return self._bounds

    @staticmethod
    def logit(x, a=0.0, b=1.0):
        r"""Computes the logit function with domain :math:`x \in (a, b)`.

        This is given by:

        .. math::

            \mathrm{logit}(x; a, b) = \log\left(\frac{x-a}{b-x}\right).

        Note that this is also the inverse of the logistic function with range
        :math:`(a, b)`.

        Parameters
        ----------
        x : float
            The value to evaluate.
        a : float, optional
            The minimum bound of the domain of x. Default is 0.
        b : float, optional
            The maximum bound of the domain of x. Default is 1.

        Returns
        -------
        float
            The logit of x.
        """
        return numpy.log(x - a) - numpy.log(b - x)

    @staticmethod
    def logistic(x, a=0.0, b=1.0):
        r"""Computes the logistic function with range :math:`\in (a, b)`.

        This is given by:

        .. math::

            \mathrm{logistic}(x; a, b) = \frac{a + b e^x}{1 + e^x}.

        Note that this is also the inverse of the logit function with domain
        :math:`(a, b)`.

        Parameters
        ----------
        x : float
            The value to evaluate.
        a : float, optional
            The minimum bound of the range of the logistic function. Default
            is 0.
        b : float, optional
            The maximum bound of the range of the logistic function. Default
            is 1.

        Returns
        -------
        float
            The logistic of x.
        """
        expx = numpy.exp(x)
        return (a + b * expx) / (1.0 + expx)

    def transform(self, maps):
        r"""Computes :math:`\mathrm{logit}(x; a, b)`.

        The domain :math:`a, b` of :math:`x` are given by the class's bounds.

        Parameters
        ----------
        maps : dict or FieldArray
            A dictionary or FieldArray which provides a map between the
            parameter name of the variable to transform and its value(s).

        Returns
        -------
        out : dict or FieldArray
            A map between the transformed variable name and value(s), along
            with the original variable name and value(s).
        """
        x = maps[self._inputvar]
        # check that x is in bounds
        isin = self._bounds.__contains__(x)
        if isinstance(isin, numpy.ndarray):
            isin = isin.all()
        if not isin:
            raise ValueError("one or more values are not in bounds")
        out = {self._outputvar: self.logit(x, self._a, self._b)}
        return self.format_output(maps, out)

    def inverse_transform(self, maps):
        r"""Computes :math:`y = \mathrm{logistic}(x; a,b)`.

        The codomain :math:`a, b` of :math:`y` are given by the class's bounds.

        Parameters
        ----------
        maps : dict or FieldArray
            A dictionary or FieldArray which provides a map between the
            parameter name of the variable to transform and its value(s).

        Returns
        -------
        out : dict or FieldArray
            A map between the transformed variable name and value(s), along
            with the original variable name and value(s).
        """
        y = maps[self._outputvar]
        out = {self._inputvar: self.logistic(y, self._a, self._b)}
        return self.format_output(maps, out)

    def jacobian(self, maps):
        r"""Computes the Jacobian of :math:`y = \mathrm{logit}(x; a,b)`.

        This is:

        .. math::

            \frac{\mathrm{d}y}{\mathrm{d}x} = \frac{b -a}{(x-a)(b-x)},

        where :math:`x \in (a, b)`.

        Parameters
        ----------
        maps : dict or FieldArray
            A dictionary or FieldArray which provides a map between the
            parameter name of the variable to transform and its value(s).

        Returns
        -------
        float
            The value of the jacobian at the given point(s).
        """
        x = maps[self._inputvar]
        # check that x is in bounds
        isin = self._bounds.__contains__(x)
        if isinstance(isin, numpy.ndarray) and not isin.all():
            raise ValueError("one or more values are not in bounds")
        elif not isin:
            raise ValueError("{} is not in bounds".format(x))
        return (self._b - self._a) / ((x - self._a) * (self._b - x))

    def inverse_jacobian(self, maps):
        r"""Computes the Jacobian of :math:`y = \mathrm{logistic}(x; a,b)`.

        This is:

        .. math::

            \frac{\mathrm{d}y}{\mathrm{d}x} = \frac{e^x (b-a)}{(1+e^y)^2},

        where :math:`y \in (a, b)`.

        Parameters
        ----------
        maps : dict or FieldArray
            A dictionary or FieldArray which provides a map between the
            parameter name of the variable to transform and its value(s).

        Returns
        -------
        float
            The value of the jacobian at the given point(s).
        """
        x = maps[self._outputvar]
        expx = numpy.exp(x)
        return expx * (self._b - self._a) / (1.0 + expx) ** 2.0

    @classmethod
    def from_config(cls, cp, section, outputs,
                    skip_opts=None, additional_opts=None):
        """Initializes a Logit transform from the given section.

        The section must specify an input and output variable name. The domain
        of the input may be specified using `min-{input}`, `max-{input}`.
        Example:

        .. code-block:: ini

            [{section}-logitq]
            name = logit
            inputvar = q
            outputvar = logitq
            min-q = 1
            max-q = 8

        Parameters
        ----------
        cp : pycbc.workflow.WorkflowConfigParser
            A parsed configuration file that contains the transform options.
        section : str
            Name of the section in the configuration file.
        outputs : str
            The names of the parameters that are output by this transformation,
            separated by `VARARGS_DELIM`. These must appear in the "tag" part
            of the section header.
        skip_opts : list, optional
            Do not read options in the given list.
        additional_opts : dict, optional
            Any additional arguments to pass to the class. If an option is
            provided that also exists in the config file, the value provided
            will be used instead of being read from the file.

        Returns
        -------
        cls
            An instance of the class.
        """
        # pull out the minimum, maximum values of the input variable
        inputvar = cp.get_opt_tag(section, "inputvar", outputs)
        s = "-".join([section, outputs])
        opt = "min-{}".format(inputvar)
        if skip_opts is None:
            skip_opts = []
        if additional_opts is None:
            additional_opts = {}
        else:
            additional_opts = additional_opts.copy()
        if cp.has_option(s, opt):
            a = cp.get_opt_tag(section, opt, outputs)
            skip_opts.append(opt)
        else:
            a = None
        opt = "max-{}".format(inputvar)
        if cp.has_option(s, opt):
            b = cp.get_opt_tag(section, opt, outputs)
            skip_opts.append(opt)
        else:
            b = None
        if a is None and b is not None or b is None and a is not None:
            raise ValueError(
                "if providing a min(max)-{}, must also provide "
                "a max(min)-{}".format(inputvar, inputvar)
            )
        elif a is not None:
            additional_opts.update({"domain": (float(a), float(b))})
        return super(Logit, cls).from_config(
            cp, section, outputs, skip_opts, additional_opts
        )


#
# =============================================================================
#
#                             Inverse Transforms
#
# =============================================================================
#
class Mass1Mass2ToMchirpQ(MchirpQToMass1Mass2):
    """The inverse of MchirpQToMass1Mass2."""

    name = "mass1_mass2_to_mchirp_q"
    inverse = MchirpQToMass1Mass2
    transform = inverse.inverse_transform
    inverse_transform = inverse.transform
    jacobian = inverse.inverse_jacobian
    inverse_jacobian = inverse.jacobian

    def __init__(
        self, mass1_param=None, mass2_param=None, mchirp_param=None, q_param=None
    ):
        if mass1_param is None:
            mass1_param = parameters.mass1
        if mass2_param is None:
            mass2_param = parameters.mass2
        if mchirp_param is None:
            mchirp_param = parameters.mchirp
        if q_param is None:
            q_param = parameters.q
        self.mass1_param = mass1_param
        self.mass2_param = mass2_param
        self.mchirp_param = mchirp_param
        self.q_param = q_param
        self._inputs = [self.mass1_param, self.mass2_param]
        self._outputs = [self.mchirp_param, self.q_param]
        BaseTransform.__init__(self)


class Mass1Mass2ToMchirpEta(MchirpEtaToMass1Mass2):
    """The inverse of MchirpEtaToMass1Mass2."""

    name = "mass1_mass2_to_mchirp_eta"
    inverse = MchirpEtaToMass1Mass2
    _inputs = inverse._outputs
    _outputs = inverse._inputs
    transform = inverse.inverse_transform
    inverse_transform = inverse.transform
    jacobian = inverse.inverse_jacobian
    inverse_jacobian = inverse.jacobian


class DistanceToChirpDistance(ChirpDistanceToDistance):
    """The inverse of ChirpDistanceToDistance."""

    name = "distance_to_chirp_distance"
    inverse = ChirpDistanceToDistance
    _inputs = [parameters.distance, parameters.mchirp]
    _outputs = [parameters.chirp_distance]
    transform = inverse.inverse_transform
    inverse_transform = inverse.transform
    jacobian = inverse.inverse_jacobian
    inverse_jacobian = inverse.jacobian


class CartesianToSpherical(SphericalToCartesian):
    """Converts spherical coordinates to cartesian.

    Parameters
    ----------
    x : str
        The name of the x parameter.
    y : str
        The name of the y parameter.
    z : str
        The name of the z parameter.
    radial : str
        The name of the radial parameter.
    azimuthal : str
        The name of the azimuthal angle parameter.
    polar : str
        The name of the polar angle parameter.
    """

    name = "cartesian_to_spherical"
    inverse = SphericalToCartesian
    transform = inverse.inverse_transform
    inverse_transform = inverse.transform
    jacobian = inverse.inverse_jacobian
    inverse_jacobian = inverse.jacobian

    def __init__(self, *args):
        super(CartesianToSpherical, self).__init__(*args)
        # swap inputs and outputs
        outputs = self._inputs
        inputs = self._outputs
        self._inputs = inputs
        self._outputs = outputs
        self.inputs = set(self._inputs)
        self.outputs = set(self._outputs)


class CartesianSpin1ToSphericalSpin1(CartesianToSpherical):
    """The inverse of SphericalSpin1ToCartesianSpin1.

    **Deprecation Warning:** This will be removed in a future update. Use
    :py:class:`CartesianToSpherical` with spin-parameter names passed in
    instead.
    """

    name = "cartesian_spin_1_to_spherical_spin_1"

    def __init__(self):
        logging.warning(
            "Deprecation warning: the {} transform will be "
            "removed in a future update. Please use {} instead, "
            "passing spin1x, spin1y, spin1z, spin1_a, "
            "spin1_azimuthal, spin1_polar as arguments.".format(
                self.name, CartesianToSpherical.name
            )
        )
        super(CartesianSpin1ToSphericalSpin1, self).__init__(
            "spin1x", "spin1y", "spin1z",
            "spin1_a", "spin1_azimuthal", "spin1_polar"
        )


class CartesianSpin2ToSphericalSpin2(CartesianToSpherical):
    """The inverse of SphericalSpin2ToCartesianSpin2.

    **Deprecation Warning:** This will be removed in a future update. Use
    :py:class:`CartesianToSpherical` with spin-parameter names passed in
    instead.
    """

    name = "cartesian_spin_2_to_spherical_spin_2"

    def __init__(self):
        logging.warning(
            "Deprecation warning: the {} transform will be "
            "removed in a future update. Please use {} instead, "
            "passing spin2x, spin2y, spin2z, spin2_a, "
            "spin2_azimuthal, spin2_polar as arguments.".format(
                self.name, CartesianToSpherical.name
            )
        )
        super(CartesianSpin2ToSphericalSpin2, self).__init__(
            "spin2x", "spin2y", "spin2z",
            "spin2_a", "spin2_azimuthal", "spin2_polar"
        )


class CartesianSpinToAlignedMassSpin(AlignedMassSpinToCartesianSpin):
    """The inverse of AlignedMassSpinToCartesianSpin."""

    name = "cartesian_spin_to_aligned_mass_spin"
    inverse = AlignedMassSpinToCartesianSpin
    _inputs = inverse._outputs
    _outputs = inverse._inputs
    transform = inverse.inverse_transform
    inverse_transform = inverse.transform
    jacobian = inverse.inverse_jacobian
    inverse_jacobian = inverse.jacobian


class CartesianSpinToPrecessionMassSpin(PrecessionMassSpinToCartesianSpin):
    """The inverse of PrecessionMassSpinToCartesianSpin."""

    name = "cartesian_spin_to_precession_mass_spin"
    inverse = PrecessionMassSpinToCartesianSpin
    _inputs = inverse._outputs
    _outputs = inverse._inputs
    transform = inverse.inverse_transform
    inverse_transform = inverse.transform
    jacobian = inverse.inverse_jacobian
    inverse_jacobian = inverse.jacobian


class ChiPToCartesianSpin(CartesianSpinToChiP):
    """The inverse of `CartesianSpinToChiP`."""

    name = "cartesian_spin_to_chi_p"
    inverse = CartesianSpinToChiP
    _inputs = inverse._outputs
    _outputs = inverse._inputs
    transform = inverse.inverse_transform
    inverse_transform = inverse.transform
    jacobian = inverse.inverse_jacobian
    inverse_jacobian = inverse.jacobian


class Exponent(Log):
    """Applies an exponent transform to an `inputvar` parameter.

    This is the inverse of the log transform.

    Parameters
    ----------
    inputvar : str
        The name of the parameter to transform.
    outputvar : str
        The name of the transformed parameter.
    """

    name = "exponent"
    inverse = Log
    transform = inverse.inverse_transform
    inverse_transform = inverse.transform
    jacobian = inverse.inverse_jacobian
    inverse_jacobian = inverse.jacobian

    def __init__(self, inputvar, outputvar):
        super(Exponent, self).__init__(outputvar, inputvar)


class Logistic(Logit):
    """Applies a logistic transform from an `input` parameter to an `output`
    parameter. This is the inverse of the logit transform.

    Typically, the output of the logistic function has range :math:`\in [0,1)`.
    However, the `codomain` argument can be used to expand this to any
    finite real interval.

    Parameters
    ----------
    inputvar : str
        The name of the parameter to transform.
    outputvar : str
        The name of the transformed parameter.
    frange : tuple or distributions.bounds.Bounds, optional
        The range of the output parameter. Can be any finite
        interval. Default is (0., 1.).
    """

    name = "logistic"
    inverse = Logit
    transform = inverse.inverse_transform
    inverse_transform = inverse.transform
    jacobian = inverse.inverse_jacobian
    inverse_jacobian = inverse.jacobian

    def __init__(self, inputvar, outputvar, codomain=(0.0, 1.0)):
        super(Logistic, self).__init__(outputvar, inputvar, domain=codomain)

    @property
    def bounds(self):
        """Returns the range of the output parameter."""
        return self._bounds

    @classmethod
    def from_config(cls, cp, section, outputs,
                    skip_opts=None, additional_opts=None):
        """Initializes a Logistic transform from the given section.

        The section must specify an input and output variable name. The
        codomain of the output may be specified using `min-{output}`,
        `max-{output}`. Example:

        .. code-block:: ini

            [{section}-q]
            name = logistic
            inputvar = logitq
            outputvar = q
            min-q = 1
            max-q = 8

        Parameters
        ----------
        cp : pycbc.workflow.WorkflowConfigParser
            A parsed configuration file that contains the transform options.
        section : str
            Name of the section in the configuration file.
        outputs : str
            The names of the parameters that are output by this transformation,
            separated by `VARARGS_DELIM`. These must appear in the "tag" part
            of the section header.
        skip_opts : list, optional
            Do not read options in the given list.
        additional_opts : dict, optional
            Any additional arguments to pass to the class. If an option is
            provided that also exists in the config file, the value provided
            will be used instead of being read from the file.

        Returns
        -------
        cls
            An instance of the class.
        """
        # pull out the minimum, maximum values of the output variable
        outputvar = cp.get_opt_tag(section, "output", outputs)
        if skip_opts is None:
            skip_opts = []
        if additional_opts is None:
            additional_opts = {}
        else:
            additional_opts = additional_opts.copy()
        s = "-".join([section, outputs])
        opt = "min-{}".format(outputvar)
        if cp.has_option(s, opt):
            a = cp.get_opt_tag(section, opt, outputs)
            skip_opts.append(opt)
        else:
            a = None
        opt = "max-{}".format(outputvar)
        if cp.has_option(s, opt):
            b = cp.get_opt_tag(section, opt, outputs)
            skip_opts.append(opt)
        else:
            b = None
        if a is None and b is not None or b is None and a is not None:
            raise ValueError(
                "if providing a min(max)-{}, must also provide "
                "a max(min)-{}".format(outputvar, outputvar)
            )
        elif a is not None:
            additional_opts.update({"codomain": (float(a), float(b))})
        return super(Logistic, cls).from_config(
            cp, section, outputs, skip_opts, additional_opts
        )


# set the inverse of the forward transforms to the inverse transforms
MchirpQToMass1Mass2.inverse = Mass1Mass2ToMchirpQ
ChirpDistanceToDistance.inverse = DistanceToChirpDistance
SphericalToCartesian.inverse = CartesianToSpherical
SphericalSpin1ToCartesianSpin1.inverse = CartesianSpin1ToSphericalSpin1
SphericalSpin2ToCartesianSpin2.inverse = CartesianSpin2ToSphericalSpin2
AlignedMassSpinToCartesianSpin.inverse = CartesianSpinToAlignedMassSpin
PrecessionMassSpinToCartesianSpin.inverse = CartesianSpinToPrecessionMassSpin
ChiPToCartesianSpin.inverse = CartesianSpinToChiP
Log.inverse = Exponent
Logit.inverse = Logistic


#
# =============================================================================
#
#                      Collections of transforms
#
# =============================================================================
#

# dictionary of all transforms
transforms = {
    CustomTransform.name: CustomTransform,
    MchirpQToMass1Mass2.name: MchirpQToMass1Mass2,
    Mass1Mass2ToMchirpQ.name: Mass1Mass2ToMchirpQ,
    MchirpEtaToMass1Mass2.name: MchirpEtaToMass1Mass2,
    Mass1Mass2ToMchirpEta.name: Mass1Mass2ToMchirpEta,
    ChirpDistanceToDistance.name: ChirpDistanceToDistance,
    DistanceToChirpDistance.name: DistanceToChirpDistance,
    SphericalToCartesian.name: SphericalToCartesian,
    CartesianToSpherical.name: CartesianToSpherical,
    SphericalSpin1ToCartesianSpin1.name: SphericalSpin1ToCartesianSpin1,
    CartesianSpin1ToSphericalSpin1.name: CartesianSpin1ToSphericalSpin1,
    SphericalSpin2ToCartesianSpin2.name: SphericalSpin2ToCartesianSpin2,
    CartesianSpin2ToSphericalSpin2.name: CartesianSpin2ToSphericalSpin2,
    DistanceToRedshift.name: DistanceToRedshift,
    AlignedMassSpinToCartesianSpin.name: AlignedMassSpinToCartesianSpin,
    CartesianSpinToAlignedMassSpin.name: CartesianSpinToAlignedMassSpin,
    PrecessionMassSpinToCartesianSpin.name: PrecessionMassSpinToCartesianSpin,
    CartesianSpinToPrecessionMassSpin.name: CartesianSpinToPrecessionMassSpin,
    ChiPToCartesianSpin.name: ChiPToCartesianSpin,
    CartesianSpinToChiP.name: CartesianSpinToChiP,
    Log.name: Log,
    Exponent.name: Exponent,
    Logit.name: Logit,
    Logistic.name: Logistic,
    LambdaFromTOVFile.name: LambdaFromTOVFile,
    LambdaFromMultipleTOVFiles.name: LambdaFromMultipleTOVFiles,
    AlignTotalSpin.name: AlignTotalSpin,
}

# standard CBC transforms: these are transforms that do not require input
# arguments; they are typically used in CBC parameter estimation to transform
# to coordinates understood by the waveform generator
common_cbc_forward_transforms = [
    MchirpQToMass1Mass2(),
    DistanceToRedshift(),
    SphericalToCartesian(
        parameters.spin1x,
        parameters.spin1y,
        parameters.spin1z,
        parameters.spin1_a,
        parameters.spin1_azimuthal,
        parameters.spin1_polar,
    ),
    SphericalToCartesian(
        parameters.spin2x,
        parameters.spin2y,
        parameters.spin2z,
        parameters.spin2_a,
        parameters.spin2_azimuthal,
        parameters.spin2_polar,
    ),
    AlignedMassSpinToCartesianSpin(),
    PrecessionMassSpinToCartesianSpin(),
    ChiPToCartesianSpin(),
    ChirpDistanceToDistance(),
]
common_cbc_inverse_transforms = [
    _t.inverse()
    for _t in common_cbc_forward_transforms
    if not (_t.inverse is None or _t.name == "spherical_to_cartesian")
]
common_cbc_inverse_transforms.extend(
    [
        CartesianToSpherical(
            parameters.spin1x,
            parameters.spin1y,
            parameters.spin1z,
            parameters.spin1_a,
            parameters.spin1_azimuthal,
            parameters.spin1_polar,
        ),
        CartesianToSpherical(
            parameters.spin2x,
            parameters.spin2y,
            parameters.spin2z,
            parameters.spin2_a,
            parameters.spin2_azimuthal,
            parameters.spin2_polar,
        ),
    ]
)

common_cbc_transforms = common_cbc_forward_transforms \
                        + common_cbc_inverse_transforms


def get_common_cbc_transforms(requested_params, variable_args, valid_params=None):
    """Determines if any additional parameters from the InferenceFile are
    needed to get derived parameters that user has asked for.

    First it will try to add any base parameters that are required to calculate
    the derived parameters. Then it will add any sampling parameters that are
    required to calculate the base parameters needed.

    Parameters
    ----------
    requested_params : list
        List of parameters that user wants.
    variable_args : list
        List of parameters that InferenceFile has.
    valid_params : list
        List of parameters that can be accepted.

    Returns
    -------
    requested_params : list
        Updated list of parameters that user wants.
    all_c : list
        List of BaseTransforms to apply.
    """
    variable_args = (
        set(variable_args) if not isinstance(variable_args, set) else variable_args
    )

    # try to parse any equations by putting all strings together
    # this will get some garbage but ensures all alphanumeric/underscored
    # parameter names are added
    new_params = []
    for opt in requested_params:
        s = ""
        for ch in opt:
            s += ch if ch.isalnum() or ch == "_" else " "
        new_params += s.split(" ")
    requested_params = set(list(requested_params) + list(new_params))

    # can pass a list of valid parameters to remove garbage from parsing above
    if valid_params:
        valid_params = set(valid_params)
        requested_params = requested_params.intersection(valid_params)

    # find all the transforms for the requested derived parameters
    # calculated from base parameters
    from_base_c = []
    for converter in common_cbc_inverse_transforms:
        if converter.outputs.issubset(variable_args) or \
           converter.outputs.isdisjoint(requested_params):
            continue
        intersect = converter.outputs.intersection(requested_params)
        if (
            not intersect
            or intersect.issubset(converter.inputs)
            or intersect.issubset(variable_args)
        ):
            continue
        requested_params.update(converter.inputs)
        from_base_c.append(converter)

    # find all the tranforms for the required base parameters
    # calculated from sampling parameters
    to_base_c = []
    for converter in common_cbc_forward_transforms:
        if (
            converter.inputs.issubset(variable_args)
            and len(converter.outputs.intersection(requested_params)) > 0
        ):
            requested_params.update(converter.inputs)
            to_base_c.append(converter)
            variable_args.update(converter.outputs)

    # get list of transforms that converts sampling parameters to the base
    # parameters and then converts base parameters to the derived parameters
    all_c = to_base_c + from_base_c

    return list(requested_params), all_c


def apply_transforms(samples, transforms, inverse=False):
    """Applies a list of BaseTransform instances on a mapping object.

    Parameters
    ----------
    samples : {FieldArray, dict}
        Mapping object to apply transforms to.
    transforms : list
        List of BaseTransform instances to apply. Nested transforms are assumed
        to be in order for forward transforms.
    inverse : bool, optional
        Apply inverse transforms. In this case transforms will be applied in
        the opposite order. Default is False.

    Returns
    -------
    samples : {FieldArray, dict}
        Mapping object with transforms applied. Same type as input.
    """
    if inverse:
        transforms = transforms[::-1]
    for t in transforms:
        try:
            if inverse:
                samples = t.inverse_transform(samples)
            else:
                samples = t.transform(samples)
        except NotImplementedError:
            continue
    return samples


def compute_jacobian(samples, transforms, inverse=False):
    """Computes the jacobian of the list of transforms at the given sample
    points.

    Parameters
    ----------
    samples : {FieldArray, dict}
        Mapping object specifying points at which to compute jacobians.
    transforms : list
        List of BaseTransform instances to apply. Nested transforms are assumed
        to be in order for forward transforms.
    inverse : bool, optional
        Compute inverse jacobians. Default is False.

    Returns
    -------
    float :
        The product of the jacobians of all fo the transforms.
    """
    j = 1.0
    if inverse:
        for t in transforms:
            j *= t.inverse_jacobian(samples)
    else:
        for t in transforms:
            j *= t.jacobian(samples)
    return j


def order_transforms(transforms):
    """Orders transforms to ensure proper chaining.

    For example, if `transforms = [B, A, C]`, and `A` produces outputs needed
    by `B`, the transforms will be re-rorderd to `[A, B, C]`.

    Parameters
    ----------
    transforms : list
        List of transform instances to order.

    Outputs
    -------
    list :
        List of transformed ordered such that forward transforms can be carried
        out without error.
    """
    # get a set of all inputs and all outputs
    outputs = set().union(*[set(t.outputs)-set(t.inputs) for t in transforms])
    out = []
    remaining = [t for t in transforms]
    while remaining:
        # pull out transforms that have no inputs in the set of outputs
        leftover = []
        for t in remaining:
            if t.inputs.isdisjoint(outputs):
                out.append(t)
                outputs -= t.outputs
            else:
                leftover.append(t)
        remaining = leftover
    return out


def read_transforms_from_config(cp, section="transforms"):
    """Returns a list of PyCBC transform instances for a section in the
    given configuration file.

    If the transforms are nested (i.e., the output of one transform is the
    input of another), the returned list will be sorted by the order of the
    nests.

    Parameters
    ----------
    cp : WorflowConfigParser
        An open config file to read.
    section : {"transforms", string}
        Prefix on section names from which to retrieve the transforms.

    Returns
    -------
    list
        A list of the parsed transforms.
    """
    trans = []
    for subsection in cp.get_subsections(section):
        name = cp.get_opt_tag(section, "name", subsection)
        t = transforms[name].from_config(cp, section, subsection)
        trans.append(t)
    return order_transforms(trans)
