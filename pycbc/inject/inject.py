# Copyright (C) 2012  Alex Nitz, Tito Dal Canton
#
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


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
"""This module provides utilities for injecting signals into data"""

import os
import numpy as np
import lal
import copy
import logging
from abc import ABCMeta, abstractmethod
import h5py
from pycbc import waveform, frame, libutils
from pycbc.opt import LimitedSizeDict
from pycbc.waveform import (get_td_waveform, fd_det,
                            get_td_det_waveform_from_fd_det)
from pycbc.waveform import utils as wfutils
from pycbc.waveform import ringdown_td_approximants
from pycbc.types import float64, float32, TimeSeries, load_timeseries
from pycbc.detector import Detector
from pycbc.conversions import tau0_from_mass1_mass2
from pycbc.filter import resample_to_delta_t
import pycbc.io
from pycbc.io.ligolw import LIGOLWContentHandler
from ligo.lw import utils as ligolw_utils, ligolw, lsctables

sim = libutils.import_optional('lalsimulation')

injection_func_map = {
    np.dtype(float32): lambda *args: sim.SimAddInjectionREAL4TimeSeries(*args),
    np.dtype(float64): lambda *args: sim.SimAddInjectionREAL8TimeSeries(*args),
}

# Map parameter names used in pycbc to names used in the sim_inspiral
# table, if they are different
sim_inspiral_map = {
    'ra': 'longitude',
    'dec': 'latitude',
    'approximant': 'waveform',
    }

def set_sim_data(inj, field, data):
    """Sets data of a SimInspiral instance."""
    try:
        sim_field = sim_inspiral_map[field]
    except KeyError:
        sim_field = field
    # for tc, map to geocentric times
    if sim_field == 'tc':
        inj.geocent_end_time = int(data)
        inj.geocent_end_time_ns = int(1e9*(data % 1))
    # for spin1 and spin2 we need data to be an array
    if sim_field in ['spin1', 'spin2']:
        setattr(inj, sim_field, [0, 0, data])
    else:
        setattr(inj, sim_field, data)


def projector(detector_name, inj, hp, hc, distance_scale=1):
    """ Use the injection row to project the polarizations into the
    detector frame
    """
    detector = Detector(detector_name)

    hp /= distance_scale
    hc /= distance_scale

    try:
        tc = inj.tc
        ra = inj.ra
        dec = inj.dec
    except:
        tc = inj.time_geocent
        ra = inj.longitude
        dec = inj.latitude

    hp.start_time += tc
    hc.start_time += tc

    # taper the polarizations
    try:
        hp_tapered = wfutils.taper_timeseries(hp, inj.taper)
        hc_tapered = wfutils.taper_timeseries(hc, inj.taper)
    except AttributeError:
        hp_tapered = hp
        hc_tapered = hc

    projection_method = 'lal'
    if hasattr(inj, 'detector_projection_method'):
        projection_method = inj.detector_projection_method

    logging.info('Injecting at %s, method is %s', tc, projection_method)

    # compute the detector response and add it to the strain
    signal = detector.project_wave(hp_tapered, hc_tapered,
                                   ra, dec, inj.polarization,
                                   method=projection_method,
                                   reference_time=tc,)
    return signal

def legacy_approximant_name(apx):
    """Convert the old style xml approximant name to a name
    and phase_order. Alex: I hate this function. Please delete this when we
    use Collin's new tables.
    """
    apx = str(apx)
    try:
        order = sim.GetOrderFromString(apx)
    except:
        print("Warning: Could not read phase order from string, using default")
        order = -1
    name = sim.GetStringFromApproximant(sim.GetApproximantFromString(apx))
    return name, order


class _XMLInjectionSet(object):

    """Manages sets of injections: reads injections from LIGOLW XML files
    and injects them into time series.

    Parameters
    ----------
    sim_file : string
        Path to a LIGOLW XML file containing a SimInspiralTable
        with injection definitions.

    Attributes
    ----------
    indoc
    table
    """

    def __init__(self, sim_file, **kwds):
        self.indoc = ligolw_utils.load_filename(
            sim_file, False, contenthandler=LIGOLWContentHandler)
        self.table = lsctables.SimInspiralTable.get_table(self.indoc)
        self.extra_args = kwds

    def apply(self, strain, detector_name, f_lower=None, distance_scale=1,
              simulation_ids=None,
              inj_filter_rejector=None,
              injection_sample_rate=None,):
        """Add injections (as seen by a particular detector) to a time series.

        Parameters
        ----------
        strain : TimeSeries
            Time series to inject signals into, of type float32 or float64.
        detector_name : string
            Name of the detector used for projecting injections.
        f_lower : {None, float}, optional
            Low-frequency cutoff for injected signals. If None, use value
            provided by each injection.
        distance_scale: {1, float}, optional
            Factor to scale the distance of an injection with. The default is
            no scaling.
        simulation_ids: iterable, optional
            If given, only inject signals with the given simulation IDs.
        inj_filter_rejector: InjFilterRejector instance; optional, default=None
            If given send each injected waveform to the InjFilterRejector
            instance so that it can store a reduced representation of that
            injection if necessary.
        injection_sample_rate: float, optional
            The sample rate to generate the signal before injection

        Returns
        -------
        None

        Raises
        ------
        TypeError
            For invalid types of `strain`.
        """
        if strain.dtype not in (float32, float64):
            raise TypeError("Strain dtype must be float32 or float64, not " \
                    + str(strain.dtype))

        lalstrain = strain.lal()
        earth_travel_time = lal.REARTH_SI / lal.C_SI
        t0 = float(strain.start_time) - earth_travel_time
        t1 = float(strain.end_time) + earth_travel_time

        # pick lalsimulation injection function
        add_injection = injection_func_map[strain.dtype]

        delta_t = strain.delta_t
        if injection_sample_rate is not None:
            delta_t = 1.0 / injection_sample_rate

        injections = self.table
        if simulation_ids:
            injections = [inj for inj in injections \
                          if inj.simulation_id in simulation_ids]
        injection_parameters = []
        for inj in injections:
            f_l = inj.f_lower if f_lower is None else f_lower
            # roughly estimate if the injection may overlap with the segment
            # Add 2s to end_time to account for ringdown and light-travel delay
            end_time = inj.time_geocent + 2
            inj_length = tau0_from_mass1_mass2(inj.mass1, inj.mass2, f_l)
            # Start time is taken as twice approx waveform length with a 1s
            # safety buffer
            start_time = inj.time_geocent - 2 * (inj_length + 1)
            if end_time < t0 or start_time > t1:
                continue
            signal = self.make_strain_from_inj_object(inj, delta_t,
                    detector_name, f_lower=f_l, distance_scale=distance_scale)
            signal = resample_to_delta_t(signal, strain.delta_t, method='ldas')
            if float(signal.start_time) > t1:
                continue

            signal = signal.astype(strain.dtype)
            signal_lal = signal.lal()
            add_injection(lalstrain, signal_lal, None)
            injection_parameters.append(inj)
            if inj_filter_rejector is not None:
                sid = inj.simulation_id
                inj_filter_rejector.generate_short_inj_from_inj(signal, sid)

        strain.data[:] = lalstrain.data.data[:]

        injected = copy.copy(self)
        injected.table = lsctables.SimInspiralTable()
        injected.table += injection_parameters
        if inj_filter_rejector is not None:
            inj_filter_rejector.injection_params = injected
        return injected

    def make_strain_from_inj_object(self, inj, delta_t, detector_name,
                                    f_lower=None, distance_scale=1):
        """Make a h(t) strain time-series from an injection object as read from
        a sim_inspiral table, for example.

        Parameters
        -----------
        inj : injection object
            The injection object to turn into a strain h(t).
        delta_t : float
            Sample rate to make injection at.
        detector_name : string
            Name of the detector used for projecting injections.
        f_lower : {None, float}, optional
            Low-frequency cutoff for injected signals. If None, use value
            provided by each injection.
        distance_scale: {1, float}, optional
            Factor to scale the distance of an injection with. The default is
            no scaling.

        Returns
        --------
        signal : float
            h(t) corresponding to the injection.
        """
        f_l = inj.f_lower if f_lower is None else f_lower

        name, phase_order = legacy_approximant_name(inj.waveform)

        # compute the waveform time series
        hp, hc = get_td_waveform(
            inj, approximant=name, delta_t=delta_t,
            phase_order=phase_order,
            f_lower=f_l, distance=inj.distance,
            **self.extra_args)
        return projector(detector_name,
                         inj, hp, hc, distance_scale=distance_scale)

    def end_times(self):
        """Return the end times of all injections"""
        return [inj.time_geocent for inj in self.table]

    @staticmethod
    def write(filename, samples, write_params=None, static_args=None):
        """Writes the injection samples to the given xml.

        Parameters
        ----------
        filename : str
            The name of the file to write to.
        samples : io.FieldArray
            FieldArray of parameters.
        write_params : list, optional
            Only write the given parameter names. All given names must be keys
            in ``samples``. Default is to write all parameters in ``samples``.
        static_args : dict, optional
            Dictionary mapping static parameter names to values. These are
            written to the ``attrs``.
        """
        xmldoc = ligolw.Document()
        xmldoc.appendChild(ligolw.LIGO_LW())
        simtable = lsctables.New(lsctables.SimInspiralTable)
        xmldoc.childNodes[0].appendChild(simtable)
        if static_args is None:
            static_args = {}
        if write_params is None:
            write_params = samples.fieldnames
        for ii in range(samples.size):
            sim = lsctables.SimInspiral()
            # initialize all elements to None
            for col in sim.__slots__:
                setattr(sim, col, None)
            for field in write_params:
                data = samples[ii][field]
                set_sim_data(sim, field, data)
            # set any static args
            for (field, value) in static_args.items():
                set_sim_data(sim, field, value)
            simtable.append(sim)
        ligolw_utils.write_filename(xmldoc, filename, compress='auto')


# -----------------------------------------------------------------------------


class _HDFInjectionSet(metaclass=ABCMeta):
    """Manages sets of injections: reads injections from hdf files
    and injects them into time series.

    Parameters
    ----------
    sim_file : string
        Path to an hdf file containing injections.
    \**kwds :
        The rest of the keyword arguments are passed to the waveform generation
        function when generating injections.

    Attributes
    ----------
    filehandler
    table
    static_args
    extra_args
    required_params : tuple
        Parameter names that must exist in the injection HDF file in order to
        create an injection of that type.
    """

    _tableclass = pycbc.io.FieldArray
    injtype = None
    required_params = ()

    def __init__(self, sim_file, hdf_group=None, **kwds):
        # open the file
        fp = h5py.File(sim_file, 'r')
        group = fp if hdf_group is None else fp[hdf_group]
        self.filehandler = fp
        # get parameters
        parameters = list(group.keys())
        # get all injection parameter values
        injvals = {param: group[param][()] for param in parameters}
        # make sure Numpy S strings are loaded as strings and not bytestrings
        # (which could mess with approximant names, for example)
        for k in injvals:
            if injvals[k].dtype.kind == 'S':
                injvals[k] = injvals[k].astype('U')
        # if there were no variable args, then we only have a single injection
        if len(parameters) == 0:
            numinj = 1
        else:
            numinj = tuple(injvals.values())[0].size
        # add any static args in the file
        try:
            # ensure parameter names are string types
            self.static_args = group.attrs['static_args'].astype('U')
        except KeyError:
            self.static_args = []
        parameters.extend(self.static_args)
        # we'll expand the static args to be arrays with the same size as
        # the other values
        for param in self.static_args:
            val = group.attrs[param]
            # if val is a list or numpy array, we need to store it as an
            # object; otherwise, we'll get a shape mismatch between fields
            if isinstance(val, (np.ndarray, list, tuple)):
                arr = np.empty(numinj, dtype=object)
                for ii in range(numinj):
                    arr[ii] = val
            else:
                # otherwise, we can just repeat the value the needed number of
                # times
                arr = np.repeat(val, numinj)
            # make sure any byte strings are stored as strings instead
            if arr.dtype.char == 'S':
                arr = arr.astype('U')
            injvals[param] = arr
        # make sure required parameters are provided
        missing = set(self.required_params) - set(injvals.keys())
        if missing:
            raise ValueError("required parameter(s) {} not found in the given "
                             "injection file".format(', '.join(missing)))
        # initialize the table
        self.table = self._tableclass.from_kwargs(**injvals)
        # save the extra arguments
        self.extra_args = kwds

    @abstractmethod
    def apply(self, strain, detector_name, distance_scale=1,
              simulation_ids=None, inj_filter_rejector=None,
              **kwargs):
        """Adds injections to a detector's time series."""
        pass

    @abstractmethod
    def make_strain_from_inj_object(self, inj, delta_t, detector_name,
                                    distance_scale=1, **kwargs):
        """Make a h(t) strain time-series from an injection object.
        """
        pass

    @abstractmethod
    def end_times(self):
        """Return the end times of all injections"""
        pass

    @abstractmethod
    def supported_approximants(self):
        """Return a list of the supported approximants."""
        pass

    @classmethod
    def write(cls, filename, samples, write_params=None, static_args=None,
              **metadata):
        """Writes the injection samples to the given hdf file.

        Parameters
        ----------
        filename : str
            The name of the file to write to.
        samples : io.FieldArray
            FieldArray of parameters.
        write_params : list, optional
            Only write the given parameter names. All given names must be keys
            in ``samples``. Default is to write all parameters in ``samples``.
        static_args : dict, optional
            Dictionary mapping static parameter names to values. These are
            written to the ``attrs``.
        \**metadata :
            All other keyword arguments will be written to the file's attrs.
        """
        with h5py.File(filename, 'w') as fp:
            # write metadata
            if static_args is None:
                static_args = {}
            fp.attrs["static_args"] = list(map(str, static_args.keys()))
            fp.attrs['injtype'] = cls.injtype
            for key, val in metadata.items():
                fp.attrs[key] = val
            if write_params is None:
                write_params = samples.fieldnames
            for arg, val in static_args.items():
                try:
                    fp.attrs[arg] = val
                except TypeError:
                    # can get this in python 3 if the val was numpy.str_ type
                    # try decoding it and writing
                    fp.attrs[arg] = str(val)
            for field in write_params:
                try:
                    fp[field] = samples[field]
                except TypeError as e:
                    # can get this in python 3 if the val was a numpy.str_ type
                    # we'll try again as a string type
                    if samples[field].dtype.char == 'U':
                        fp[field] = samples[field].astype('S')
                    else:
                        raise e


class CBCHDFInjectionSet(_HDFInjectionSet):
    """Manages CBC injections.
    """
    _tableclass = pycbc.io.WaveformArray
    injtype = 'cbc'
    required_params = ('tc',)

    def apply(self, strain, detector_name, f_lower=None, distance_scale=1,
              simulation_ids=None,
              inj_filter_rejector=None,
              injection_sample_rate=None,):
        """Add injections (as seen by a particular detector) to a time series.

        Parameters
        ----------
        strain : TimeSeries
            Time series to inject signals into, of type float32 or float64.
        detector_name : string
            Name of the detector used for projecting injections.
        f_lower : {None, float}, optional
            Low-frequency cutoff for injected signals. If None, use value
            provided by each injection.
        distance_scale: {1, float}, optional
            Factor to scale the distance of an injection with. The default is
            no scaling.
        simulation_ids: iterable, optional
            If given, only inject signals with the given simulation IDs.
        inj_filter_rejector: InjFilterRejector instance; optional, default=None
            If given send each injected waveform to the InjFilterRejector
            instance so that it can store a reduced representation of that
            injection if necessary.
        injection_sample_rate: float, optional
            The sample rate to generate the signal before injection

        Returns
        -------
        None

        Raises
        ------
        TypeError
            For invalid types of `strain`.
        """
        if strain.dtype not in (float32, float64):
            raise TypeError("Strain dtype must be float32 or float64, not " \
                    + str(strain.dtype))

        lalstrain = strain.lal()
        if self.table[0]['approximant'] in fd_det:
            t0 = float(strain.start_time)
            t1 = float(strain.end_time)
        else:
            earth_travel_time = lal.REARTH_SI / lal.C_SI
            t0 = float(strain.start_time) - earth_travel_time
            t1 = float(strain.end_time) + earth_travel_time

        # pick lalsimulation injection function
        add_injection = injection_func_map[strain.dtype]

        delta_t = strain.delta_t
        if injection_sample_rate is not None:
            delta_t = 1.0 / injection_sample_rate

        injections = self.table
        if simulation_ids:
            injections = injections[list(simulation_ids)]

        injected_ids = []
        for ii, inj in enumerate(injections):
            f_l = inj.f_lower if f_lower is None else f_lower
            # roughly estimate if the injection may overlap with the segment
            # Add 2s to end_time to account for ringdown and light-travel delay
            end_time = inj.tc + 2
            inj_length = tau0_from_mass1_mass2(inj.mass1, inj.mass2, f_l)
            # Start time is taken as twice approx waveform length with a 1s
            # safety buffer
            start_time = inj.tc - 2 * (inj_length + 1)
            if end_time < t0 or start_time > t1:
                continue
            signal = self.make_strain_from_inj_object(inj, delta_t,
                     detector_name, f_lower=f_l,
                     distance_scale=distance_scale)
            signal = resample_to_delta_t(signal, strain.delta_t, method='ldas')
            if float(signal.start_time) > t1:
                continue

            signal = signal.astype(strain.dtype)
            signal_lal = signal.lal()
            add_injection(lalstrain, signal_lal, None)
            injected_ids.append(ii)
            if inj_filter_rejector is not None:
                inj_filter_rejector.generate_short_inj_from_inj(signal, ii)

        strain.data[:] = lalstrain.data.data[:]

        injected = copy.copy(self)
        injected.table = injections[np.array(injected_ids).astype(int)]
        if inj_filter_rejector is not None:
            if hasattr(inj_filter_rejector, 'injected'):
                prev_p = inj_filter_rejector.injection_params
                prev_id = inj_filter_rejector.injection_ids
                injected = np.concatenate([prev_p, injected])
                injected_ids = np.concatenate([prev_id, injected_ids])

            inj_filter_rejector.injection_params = injected
            inj_filter_rejector.injection_ids = injected_ids
        return injected

    def make_strain_from_inj_object(self, inj, delta_t, detector_name,
                                    f_lower=None, distance_scale=1):
        """Make a h(t) strain time-series from an injection object.

        Parameters
        -----------
        inj : injection object
            The injection object to turn into a strain h(t). Can be any
            object which has waveform parameters as attributes, such as an
            element in a ``WaveformArray``.
        delta_t : float
            Sample rate to make injection at.
        detector_name : string
            Name of the detector used for projecting injections.
        f_lower : {None, float}, optional
            Low-frequency cutoff for injected signals. If None, use value
            provided by each injection.
        distance_scale: {1, float}, optional
            Factor to scale the distance of an injection with. The default is
            no scaling.

        Returns
        --------
        signal : float
            h(t) corresponding to the injection.
        """
        if f_lower is None:
            f_l = inj.f_lower
        else:
            f_l = f_lower

        if inj['approximant'] in fd_det:
            strain = get_td_det_waveform_from_fd_det(
                        inj, delta_t=delta_t, f_lower=f_l,
                        ifos=detector_name, **self.extra_args)[detector_name]
            strain /= distance_scale
        else:
            # compute the waveform time series
            hp, hc = get_td_waveform(inj, delta_t=delta_t, f_lower=f_l,
                                     **self.extra_args)
            strain = projector(detector_name,
                               inj, hp, hc, distance_scale=distance_scale)
        return strain

    def end_times(self):
        """Return the end times of all injections"""
        return self.table.tc

    @staticmethod
    def supported_approximants():
        all_apprxs = []
        for d in [waveform.waveform.td_wav, waveform.waveform.fd_wav]:
            for key in d:
                all_apprxs.extend(d[key])
        all_apprxs.extend(waveform.waveform.fd_det)
        return list(set(all_apprxs))


class RingdownHDFInjectionSet(_HDFInjectionSet):
    """Manages a ringdown injection: reads injection from hdf file
    and injects it into time series.
    """
    injtype = 'ringdown'
    required_params = ('tc',)

    def apply(self, strain, detector_name, distance_scale=1,
              simulation_ids=None, inj_filter_rejector=None,
              injection_sample_rate=None):
        """Add injection (as seen by a particular detector) to a time series.

        Parameters
        ----------
        strain : TimeSeries
            Time series to inject signals into, of type float32 or float64.
        detector_name : string
            Name of the detector used for projecting injections.
        distance_scale: float, optional
            Factor to scale the distance of an injection with. The default (=1)
            is no scaling.
        simulation_ids: iterable, optional
            If given, only inject signals with the given simulation IDs.
        inj_filter_rejector: InjFilterRejector instance, optional
            Not implemented. If not ``None``, a ``NotImplementedError`` will
            be raised.
        injection_sample_rate: float, optional
            The sample rate to generate the signal before injection

        Returns
        -------
        None

        Raises
        ------
        NotImplementedError
            If an ``inj_filter_rejector`` is provided.
        TypeError
            For invalid types of `strain`.
        """
        if inj_filter_rejector is not None:
            raise NotImplementedError("Ringdown injections do not support "
                                      "inj_filter_rejector")
        if strain.dtype not in (float32, float64):
            raise TypeError("Strain dtype must be float32 or float64, not " \
                    + str(strain.dtype))

        lalstrain = strain.lal()

        # pick lalsimulation injection function
        add_injection = injection_func_map[strain.dtype]

        delta_t = strain.delta_t
        if injection_sample_rate is not None:
            delta_t = 1.0 / injection_sample_rate

        injections = self.table
        if simulation_ids:
            injections = injections[list(simulation_ids)]
        for ii in range(injections.size):
            injection = injections[ii]
            signal = self.make_strain_from_inj_object(
                injection, delta_t, detector_name,
                distance_scale=distance_scale)
            signal = resample_to_delta_t(signal, strain.delta_t, method='ldas')
            signal = signal.astype(strain.dtype)
            signal_lal = signal.lal()
            add_injection(lalstrain, signal_lal, None)

            strain.data[:] = lalstrain.data.data[:]

    def make_strain_from_inj_object(self, inj, delta_t, detector_name,
                                    distance_scale=1):
        """Make a h(t) strain time-series from an injection object as read from
        an hdf file.

        Parameters
        -----------
        inj : injection object
            The injection object to turn into a strain h(t).
        delta_t : float
            Sample rate to make injection at.
        detector_name : string
            Name of the detector used for projecting injections.
        distance_scale: float, optional
            Factor to scale the distance of an injection with. The default (=1)
            is no scaling.

        Returns
        --------
        signal : float
            h(t) corresponding to the injection.
        """
        # compute the waveform time series
        hp, hc = ringdown_td_approximants[inj['approximant']](
            inj, delta_t=delta_t, **self.extra_args)
        return projector(detector_name,
                         inj, hp, hc, distance_scale=distance_scale)

    def end_times(self):
        """Return the approximate end times of all injections.

        Currently, this just assumes all ringdowns are 2 seconds long.
        """
        # the start times are the tcs
        tcs = self.table.tc
        # FIXME: this could be figured out using the ringdown module
        return tcs + 2

    @staticmethod
    def supported_approximants():
        return list(waveform.ringdown_td_approximants.keys())


class IncoherentFromFileHDFInjectionSet(_HDFInjectionSet):
    """Manages injecting an arbitrary time series loaded from a file.

    The injections must have the following attributes set:

    * ``filename``: (str) the name of the file to load containing the time
      series. The file type and format can be a frame file or anything
      understood by :py:func:`pycbc.types.timeseries.load_timeseries`. If a
      frame file (ends in ``.gwf``) is specified, a ``channel`` attribute must
      also be set.

    * ``DETECTOR_gps_time``: (float) The GPS time at which the time series
      should be added to the ``DETECTOR`` data, where ``DETECTOR`` is the name
      of the instrument to inject into (e.g., ``h1_gps_time``). **The time
      series will only be injected into a detector if a GPS time is given for
      that detector.** Set to -inf, nan, or do not provide a GPS time for a
      particular detector if you do not want to inject into that detector.

    * ``ref_point``: (str or float) What to use as the reference time of the
      injected time series. The time series will be injected into the detector
      such that the ``ref_point`` in the time series occurs at the specifed
      ``DETECTOR_gps_time``. Options are: ``'start'``, ``'end'``, ``'center'``,
      ``'absmax'``, or a float giving the number of seconds into the time
      series.

    In addition, the following attributes may optionally be provided:

    * ``channel``: (str): If the filename points to a frame file, the channel
      to load in that file. Must be provided for frame files.

    * ``DETECTOR_phase_shift``: (float) Apply a phase shift to the time series
      before adding it to detector ``DETECTOR``.

    * ``DETECTOR_amp_scale``: (float) Scale the amplitude by the given amount
      before adding it to detector ``DETECTOR``.

    * ``slice_start``: (float) Slice the time series starting at
      ``ref_point + slice_start`` before injecting into the data. Measured in
      seconds.

    * ``slice_end``: (float) Slice the time series ending at
      ``ref_point + slice_end`` before injecting into the data. Measured in
      seconds.

    * ``left_taper_width``: (float) Taper the start of the time series (after
      slicing) using half a kaiser window over the given number of seconds.
      See `:py:func:waveform.utils.td_taper` for more details.

    * ``right_taper_width``: (float) Taper the end of the time series (after
      slicing) using half a kaiser window over the given number of seconds.
      See `:py:func:waveform.utils.td_taper` for more details.

    The signal will be resampled to the same sample rate as the data it is
    being injected into.

    In order to use with ``pycbc_create_injections``, set the ``approximant``
    name to ``'incoherent_from_file'``.
    """
    injtype = 'incoherent_from_file'
    required_params = ('filename', 'ref_point')
    _buffersize = 10
    _buffer = None
    _rtbuffer = None

    def end_times(self):
        raise NotImplementedError("IncoherentFromFile times cannot be "
                                  "determined without loading time series")

    @staticmethod
    def supported_approximants():
        return ['incoherent_from_file']

    def loadts(self, inj):
        """Loads an injection time series.

        After the first time a time series is loaded it will be added to an
        internal buffer for faster in case another injection uses the same
        series.
        """
        if self._buffer is None:
            # create the buffer
            self._buffer = LimitedSizeDict(size_limit=self._buffersize)
        try:
            return self._buffer[inj.filename]
        except KeyError:
            pass
        # not in buffer, so load
        if inj.filename.endswith('.gwf'):
            try:
                channel = inj.channel
            except AttributeError as _err:
                # Py3.XX: uncomment the "from _err" when we drop 2.7
                raise ValueError("Must provide a channel for "
                                 "frame files") #from _err
            ts = frame.read_frame(inj.filename, channel)
        else:
            ts = load_timeseries(inj.filename)
        # cache
        self._buffer[inj.filename] = ts
        return ts

    def set_ref_time(self, inj, ts):
        """Sets t=0 of the given time series based on what the given
        injection's ``ref_point`` is.
        """
        try:
            ref_point = inj.ref_point
        except AttributeError as _err:
            # Py3.XX: uncomment the "from _err" when we drop 2.7
            raise ValueError("Must provide a ref_point for {} injections"
                             .format(self.injtype))  #from _err
        # try to get from buffer
        if self._rtbuffer is None:
            self._rtbuffer = LimitedSizeDict(size_limit=self._buffersize)
        try:
            reftime = self._rtbuffer[inj.filename, ref_point]
        except KeyError:
            if ref_point == "start":
                reftime = 0.
            elif ref_point == "end":
                reftime = -len(ts)*ts.delta_t
            elif ref_point == "center":
                reftime = -len(ts)*ts.delta_t/2.
            elif ref_point == "absmax":
                reftime = -ts.abs_arg_max()*ts.delta_t
            elif isinstance(ref_point, (float, int)):
                reftime = -float(ref_point)
            else:
                raise ValueError("Unrecognized ref_point {} provided"
                                 .format(ref_point))
            self._rtbuffer[inj.filename, ref_point] = reftime
        ts._epoch = reftime

    @staticmethod
    def slice_and_taper(inj, ts):
        """Slices and tapers a timeseries based on the injection settings.

        This assumes that ``set_ref_time`` has been applied to the timeseries
        first. A copy of the time series will be returned even if no slicing
        or tapering is done.
        """
        try:
            tstart = inj.slice_start
        except AttributeError:
            tstart = ts.start_time
        try:
            tend = inj.slice_end
        except AttributeError:
            tend = ts.end_time
        ts = ts.time_slice(tstart, tend).copy()
        # now taper
        try:
            twidth = inj.left_taper_width
        except AttributeError:
            twidth = 0
        if twidth:
            ts = wfutils.td_taper(ts, ts.start_time, ts.start_time+twidth,
                                  side='left')
        try:
            twidth = inj.right_taper_width
        except AttributeError:
            twidth = 0
        if twidth:
            ts = wfutils.td_taper(ts, ts.end_time-twidth, ts.end_time,
                                  side='right')
        return ts

    def apply(self, strain, detector_name, distance_scale=1,
              injection_sample_rate=None, inj_filter_rejector=None):
        if inj_filter_rejector is not None:
            raise NotImplementedError("IncoherentFromFile injections do not "
                                      "support inj_filter_rejector")
        if injection_sample_rate is not None:
            delta_t = 1./injection_sample_rate
        else:
            delta_t = strain.delta_t
        injections = self.table
        for inj in injections:
            # Check if we should inject or not...
            # loading the time series like this is a bit brute-force, since
            # we really only need to know the delta_t and length of the
            # timeseries if the ref_point is anything but absmax, but that
            # would require adding logic to figure out how to get that metadata
            # based on the filetype and ref_point
            ts = self.loadts(inj)
            # set the ref time
            self.set_ref_time(inj, ts)
            # determine if we inject or not based on the times
            try:
                injtime = inj['{}_gps_time'.format(detector_name).lower()]
            except ValueError:
                injtime = -np.inf
            if np.isnan(injtime):
                # nan means don't inject
                injtime = -np.inf
            start_time = injtime + ts.start_time
            end_time = injtime + ts.end_time
            inject = (start_time < strain.end_time and
                      end_time > strain.start_time)
            if inject:
                ts = self.make_strain_from_inj_object(
                    inj, delta_t, detector_name,
                    distance_scale=distance_scale, ts=ts)
                if ts.delta_t != strain.delta_t:
                    ts = resample_to_delta_t(ts, strain.delta_t, method='ldas')
                strain.inject(ts, copy=False)

    def make_strain_from_inj_object(self, inj, delta_t, detector_name,
                                    distance_scale=1, ts=None):
        if ts is None:
            ts = load_timeseries(inj.filename)
            self.set_ref_time(inj, ts)
        # slice and taper
        ts = self.slice_and_taper(inj, ts)
        # shift reference to the detector time
        ts._epoch += inj['{}_gps_time'.format(detector_name).lower()]
        # resample
        ts = resample_to_delta_t(ts, delta_t, method='ldas')
        # apply any phase shift
        try:
            phase_shift = inj[
                '{}_phase_shift'.format(detector_name).lower()]
        except ValueError:
            phase_shift = 0
        if phase_shift:
            fs = ts.to_frequencyseries()
            fs *= np.exp(1j*phase_shift)
            ts = fs.to_timeseries()
        # apply any scaling
        try:
            amp_scale = inj[
                '{}_amp_scale'.format(detector_name).lower()]
        except ValueError:
            amp_scale = 1.
        amp_scale /= distance_scale
        ts *= amp_scale
        return ts


hdfinjtypes = {
    CBCHDFInjectionSet.injtype: CBCHDFInjectionSet,
    RingdownHDFInjectionSet.injtype: RingdownHDFInjectionSet,
    IncoherentFromFileHDFInjectionSet.injtype:
    IncoherentFromFileHDFInjectionSet,
}


def get_hdf_injtype(sim_file):
    """Gets the HDFInjectionSet class to use with the given file.

    This looks for the ``injtype`` in the given file's top level ``attrs``. If
    that attribute isn't set, will default to :py:class:`CBCHDFInjectionSet`.

    Parameters
    ----------
    sim_file : str
        Name of the file. The file must already exist.

    Returns
    -------
    HDFInjectionSet :
        The type of HDFInjectionSet to use.
    """
    with h5py.File(sim_file, 'r') as fp:
        try:
            ftype = fp.attrs['injtype']
        except KeyError:
            ftype = CBCHDFInjectionSet.injtype
    try:
        return hdfinjtypes[ftype]
    except KeyError:
        # may get a key error if the file type was stored as unicode instead
        # of string; if so, try decoding it
        try:
            ftype = str(ftype.decode())
        except AttributeError:
            # not actually a byte error; passing will reraise the KeyError
            pass
        return hdfinjtypes[ftype]


def hdf_injtype_from_approximant(approximant):
    """Gets the HDFInjectionSet class to use with the given approximant.

    Parameters
    ----------
    approximant : str
        Name of the approximant.

    Returns
    -------
    HDFInjectionSet :
        The type of HDFInjectionSet to use.
    """
    retcls = None
    for cls in hdfinjtypes.values():
        if approximant in cls.supported_approximants():
            retcls = cls
    if retcls is None:
        # none were found, raise an error
        raise ValueError("Injection file type unknown for approximant {}"
                         .format(approximant))
    return retcls


class InjectionSet(object):
    """Manages sets of injections and injects them into time series.

    Injections are read from either LIGOLW XML files or HDF files.

    Parameters
    ----------
    sim_file : string
        Path to an hdf file or a LIGOLW XML file that contains a
        SimInspiralTable.
    \**kwds :
        The rest of the keyword arguments are passed to the waveform generation
        function when generating injections.

    Attributes
    ----------
    table
    """

    def __init__(self, sim_file, **kwds):
        ext = os.path.basename(sim_file)
        if ext.endswith(('.xml', '.xml.gz', '.xmlgz')):
            self._injhandler = _XMLInjectionSet(sim_file, **kwds)
            self.indoc = self._injhandler.indoc
        else:
            # assume hdf file
            self._injhandler = get_hdf_injtype(sim_file)(sim_file, **kwds)
        self.table = self._injhandler.table
        self.extra_args = self._injhandler.extra_args
        self.apply = self._injhandler.apply
        self.make_strain_from_inj_object = \
            self._injhandler.make_strain_from_inj_object
        self.end_times = self._injhandler.end_times

    @staticmethod
    def write(filename, samples, write_params=None, static_args=None,
              injtype=None, **metadata):
        """Writes the injection samples to the given hdf file.

        Parameters
        ----------
        filename : str
            The name of the file to write to.
        samples : io.FieldArray
            FieldArray of parameters.
        write_params : list, optional
            Only write the given parameter names. All given names must be keys
            in ``samples``. Default is to write all parameters in ``samples``.
        static_args : dict, optional
            Dictionary mapping static parameter names to values. These are
            written to the ``attrs``.
        injtype : str, optional
            Specify which `HDFInjectionSet` class to use for writing. If not
            provided, will try to determine it by looking for an approximant in
            the ``static_args``, followed by the ``samples``.
        \**metadata :
            All other keyword arguments will be written to the file's attrs.
        """
        # DELETE the following "if" once xml is dropped
        ext = os.path.basename(filename)
        if ext.endswith(('.xml', '.xml.gz', '.xmlgz')):
            _XMLInjectionSet.write(filename, samples, write_params,
                                   static_args)
        else:
            # try determine the injtype if it isn't given
            if injtype is None:
                if static_args is not None and 'approximant' in static_args:
                    injcls = hdf_injtype_from_approximant(
                        static_args['approximant'])
                elif 'approximant' in samples.fieldnames:
                    apprxs = np.unique(samples['approximant'])
                    # make sure they all correspond to the same injection type
                    injcls = [hdf_injtype_from_approximant(a) for a in apprxs]
                    if not all(c == injcls[0] for c in injcls):
                        raise ValueError("injections must all be of the same "
                                         "type")
                    injcls = injcls[0]
                else:
                    raise ValueError("Could not find an approximant in the "
                                     "static args or samples to determine the "
                                     "injection type. Please specify an "
                                     "injtype instead.")
            else:
                injcls = hdfinjtypes[injtype]
            injcls.write(filename, samples, write_params, static_args,
                         **metadata)

    @staticmethod
    def from_cli(opt):
        """Return an instance of InjectionSet configured as specified
        on the command line.
        """
        if opt.injection_file is None:
            return None

        kwa = {}
        if opt.injection_f_ref is not None:
            kwa['f_ref'] = opt.injection_f_ref
        if opt.injection_f_final is not None:
            kwa['f_final'] = opt.injection_f_final
        return InjectionSet(opt.injection_file, **kwa)


class SGBurstInjectionSet(object):
    """Manages sets of sine-Gaussian burst injections: reads injections
    from LIGOLW XML files and injects them into time series.

    Parameters
    ----------
    sim_file : string
        Path to a LIGOLW XML file containing a SimBurstTable
        with injection definitions.

    Attributes
    ----------
    indoc
    table
    """

    def __init__(self, sim_file, **kwds):
        self.indoc = ligolw_utils.load_filename(
            sim_file, False, contenthandler=LIGOLWContentHandler)
        self.table = lsctables.SimBurstTable.get_table(self.indoc)
        self.extra_args = kwds

    def apply(self, strain, detector_name, f_lower=None, distance_scale=1):
        """Add injections (as seen by a particular detector) to a time series.

        Parameters
        ----------
        strain : TimeSeries
            Time series to inject signals into, of type float32 or float64.
        detector_name : string
            Name of the detector used for projecting injections.
        f_lower : {None, float}, optional
            Low-frequency cutoff for injected signals. If None, use value
            provided by each injection.
        distance_scale: {1, foat}, optional
            Factor to scale the distance of an injection with. The default is
            no scaling.

        Returns
        -------
        None

        Raises
        ------
        TypeError
            For invalid types of `strain`.
        """
        if strain.dtype not in (float32, float64):
            raise TypeError("Strain dtype must be float32 or float64, not " \
                    + str(strain.dtype))

        lalstrain = strain.lal()
        #detector = Detector(detector_name)
        earth_travel_time = lal.REARTH_SI / lal.C_SI
        t0 = float(strain.start_time) - earth_travel_time
        t1 = float(strain.end_time) + earth_travel_time

        # pick lalsimulation injection function
        add_injection = injection_func_map[strain.dtype]

        for inj in self.table:
            # roughly estimate if the injection may overlap with the segment
            end_time = inj.time_geocent
            #CHECK: This is a hack (10.0s); replace with an accurate estimate
            inj_length = 10.0
            eccentricity = 0.0
            polarization = 0.0
            start_time = end_time - 2 * inj_length
            if end_time < t0 or start_time > t1:
                continue

            # compute the waveform time series
            hp, hc = sim.SimBurstSineGaussian(float(inj.q),
                float(inj.frequency),float(inj.hrss),float(eccentricity),
                float(polarization),float(strain.delta_t))
            hp = TimeSeries(hp.data.data[:], delta_t=hp.deltaT, epoch=hp.epoch)
            hc = TimeSeries(hc.data.data[:], delta_t=hc.deltaT, epoch=hc.epoch)
            hp._epoch += float(end_time)
            hc._epoch += float(end_time)
            if float(hp.start_time) > t1:
                continue

            # compute the detector response, taper it if requested
            # and add it to the strain
            strain = wfutils.taper_timeseries(strain, inj.taper)
            signal_lal = hp.astype(strain.dtype).lal()
            add_injection(lalstrain, signal_lal, None)

        strain.data[:] = lalstrain.data.data[:]
