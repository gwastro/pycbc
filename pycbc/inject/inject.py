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
from abc import ABCMeta, abstractmethod
import lalsimulation as sim
import h5py
from pycbc import waveform
from pycbc.waveform import get_td_waveform, utils as wfutils
from pycbc.waveform import ringdown_td_approximants
from pycbc.types import float64, float32, TimeSeries
from pycbc.detector import Detector
import pycbc.io

injection_func_map = {
    np.dtype(float32): sim.SimAddInjectionREAL4TimeSeries,
    np.dtype(float64): sim.SimAddInjectionREAL8TimeSeries
}


#
# Remove everything between the dashed lines once we get rid of xml
# -----------------------------------------------------------------------------
#
from glue.ligolw import utils as ligolw_utils
from glue.ligolw import ligolw, table, lsctables

# dummy class needed for loading LIGOLW files
class LIGOLWContentHandler(ligolw.LIGOLWContentHandler):
    pass

lsctables.use_in(LIGOLWContentHandler)

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
    else:
        setattr(inj, sim_field, data)


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
        self.table = table.get_table(self.indoc, lsctables.SimInspiralTable.tableName)
        self.extra_args = kwds

    def apply(self, strain, detector_name, f_lower=None, distance_scale=1,
              simulation_ids=None, inj_filter_rejector=None):
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

        injections = self.table
        if simulation_ids:
            injections = [inj for inj in injections \
                          if inj.simulation_id in simulation_ids]
        injection_parameters = []
        for inj in injections:
            if f_lower is None:
                f_l = inj.f_lower
            else:
                f_l = f_lower
            # roughly estimate if the injection may overlap with the segment
            # Add 2s to end_time to account for ringdown and light-travel delay
            end_time = inj.get_time_geocent() + 2
            inj_length = sim.SimInspiralTaylorLength(
                strain.delta_t, inj.mass1 * lal.MSUN_SI,
                inj.mass2 * lal.MSUN_SI, f_l, 0)
            # Start time is taken as twice approx waveform length with a 1s
            # safety buffer
            start_time = inj.get_time_geocent() - 2 * (inj_length+1)
            if end_time < t0 or start_time > t1:
                continue
            signal = self.make_strain_from_inj_object(inj, strain.delta_t,
                     detector_name, f_lower=f_l, distance_scale=distance_scale)
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
        detector = Detector(detector_name)
        if f_lower is None:
            f_l = inj.f_lower
        else:
            f_l = f_lower

        name, phase_order = legacy_approximant_name(inj.waveform)

        # compute the waveform time series
        hp, hc = get_td_waveform(
            inj, approximant=name, delta_t=delta_t,
            phase_order=phase_order,
            f_lower=f_l, distance=inj.distance,
            **self.extra_args)

        hp /= distance_scale
        hc /= distance_scale

        hp._epoch += inj.get_time_geocent()
        hc._epoch += inj.get_time_geocent()

        # taper the polarizations
        hp_tapered = wfutils.taper_timeseries(hp, inj.taper)
        hc_tapered = wfutils.taper_timeseries(hc, inj.taper)

        # compute the detector response and add it to the strain
        signal = detector.project_wave(hp_tapered, hc_tapered,
                             inj.longitude, inj.latitude, inj.polarization)

        return signal

    def end_times(self):
        """Return the end times of all injections"""
        return [inj.get_time_geocent() for inj in self.table]

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
        ligolw_utils.write_filename(xmldoc, filename,
                                    gz=filename.endswith('gz'))


# -----------------------------------------------------------------------------


class _HDFInjectionSet(object):
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
    """
    __metaclass__ = ABCMeta
    _tableclass = pycbc.io.FieldArray
    injtype = None

    def __init__(self, sim_file, hdf_group=None, **kwds):
        # open the file
        fp = h5py.File(sim_file, 'r')
        group = fp if hdf_group is None else fp[hdf_group]
        self.filehandler = fp
        # get parameters
        parameters = group.keys()
        # get all injection parameter values
        injvals = {param: group[param][()] for param in parameters}
        # if there were no variable args, then we only have a single injection
        if len(parameters) == 0:
            numinj = 1
        else:
            numinj = injvals.values()[0].size
        # add any static args in the file
        try:
            self.static_args = group.attrs['static_args']
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
            injvals[param] = arr
        # make sure a coalescence time is specified for injections
        if 'tc' not in injvals:
            raise ValueError("no tc found in the given injection file; "
                             "this is needed to determine where to place the "
                             "injection")
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
            fp.attrs["static_args"] = static_args.keys()
            fp.attrs['injtype'] = cls.injtype
            for key, val in metadata.items():
                fp.attrs[key] = val
            if write_params is None:
                write_params = samples.fieldnames
            for arg, val in static_args.items():
                fp.attrs[arg] = val
            for field in write_params:
                fp[field] = samples[field]


class CBCHDFInjectionSet(_HDFInjectionSet):
    """Manages CBC injections.
    """
    _tableclass = pycbc.io.WaveformArray
    injtype = 'cbc'

    def apply(self, strain, detector_name, f_lower=None, distance_scale=1,
              simulation_ids=None, inj_filter_rejector=None):
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

        injections = self.table
        if simulation_ids:
            injections = injections[list(simulation_ids)]
        for ii in range(injections.size):
            inj = injections[ii]
            if f_lower is None:
                f_l = inj.f_lower
            else:
                f_l = f_lower
            # roughly estimate if the injection may overlap with the segment
            # Add 2s to end_time to account for ringdown and light-travel delay
            end_time = inj.tc + 2
            inj_length = sim.SimInspiralTaylorLength(
                strain.delta_t, inj.mass1 * lal.MSUN_SI,
                inj.mass2 * lal.MSUN_SI, f_l, 0)
            # Start time is taken as twice approx waveform length with a 1s
            # safety buffer
            start_time = inj.tc - 2 * (inj_length+1)
            if end_time < t0 or start_time > t1:
                continue
            signal = self.make_strain_from_inj_object(inj, strain.delta_t,
                     detector_name, f_lower=f_l, distance_scale=distance_scale)
            if float(signal.start_time) > t1:
                continue

            signal = signal.astype(strain.dtype)
            signal_lal = signal.lal()
            add_injection(lalstrain, signal_lal, None)
            if inj_filter_rejector is not None:
                inj_filter_rejector.generate_short_inj_from_inj(signal, ii)

        strain.data[:] = lalstrain.data.data[:]

        injected = copy.copy(self)
        injected.table = injections
        if inj_filter_rejector is not None:
            inj_filter_rejector.injection_params = injected
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
        detector = Detector(detector_name)
        if f_lower is None:
            f_l = inj.f_lower
        else:
            f_l = f_lower

        # compute the waveform time series
        hp, hc = get_td_waveform(inj, delta_t=delta_t, f_lower=f_l,
                                 **self.extra_args)

        hp /= distance_scale
        hc /= distance_scale

        hp._epoch += inj.tc
        hc._epoch += inj.tc

        # taper the polarizations
        try:
            hp_tapered = wfutils.taper_timeseries(hp, inj.taper)
            hc_tapered = wfutils.taper_timeseries(hc, inj.taper)
        except AttributeError:
            hp_tapered = hp
            hc_tapered = hc

        # compute the detector response and add it to the strain
        signal = detector.project_wave(hp_tapered, hc_tapered,
                             inj.ra, inj.dec, inj.polarization)

        return signal

    def end_times(self):
        """Return the end times of all injections"""
        return self.table.tc

    @staticmethod
    def supported_approximants():
        all_apprxs = []
        for d in [waveform.waveform.td_wav, waveform.waveform.fd_wav]:
            for key in d:
                all_apprxs.extend(d[key])
        return list(set(all_apprxs))


class RingdownHDFInjectionSet(_HDFInjectionSet):
    """Manages a ringdown injection: reads injection from hdf file
    and injects it into time series.
    """
    injtype = 'ringdown'

    def apply(self, strain, detector_name, distance_scale=1,
              simulation_ids=None, inj_filter_rejector=None):
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

        injections = self.table
        if simulation_ids:
            injections = injections[list(simulation_ids)]
        for ii in range(injections.size):
            injection = injections[ii]
            signal = self.make_strain_from_inj_object(
                injection, strain.delta_t, detector_name,
                distance_scale=distance_scale)
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
        detector = Detector(detector_name)

        # compute the waveform time series
        hp, hc = ringdown_td_approximants[inj['approximant']](
            inj, delta_t=delta_t, **self.extra_args)

        hp._epoch += inj['tc']
        hc._epoch += inj['tc']

        if distance_scale != 1:
            hp /= distance_scale
            hc /= distance_scale

        # compute the detector response and add it to the strain
        signal = detector.project_wave(hp, hc,
                             inj['ra'], inj['dec'], inj['polarization'])

        return signal

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


hdfinjtypes = {
    CBCHDFInjectionSet.injtype: CBCHDFInjectionSet,
    RingdownHDFInjectionSet.injtype: RingdownHDFInjectionSet,
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
        self.table = table.get_table(
            self.indoc, lsctables.SimBurstTable.tableName)
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
            end_time = inj.get_time_geocent()
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
