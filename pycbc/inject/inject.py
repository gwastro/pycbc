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

import numpy as np
import lal
import copy
import lalsimulation as sim
import h5py
from pycbc.waveform import get_td_waveform, utils as wfutils
from pycbc.waveform import ringdown_td_approximants
from pycbc_glue.ligolw import utils as ligolw_utils
from pycbc_glue.ligolw import ligolw, table, lsctables
from pycbc.types import float64, float32, TimeSeries
from pycbc.detector import Detector

injection_func_map = {
    np.dtype(float32): sim.SimAddInjectionREAL4TimeSeries,
    np.dtype(float64): sim.SimAddInjectionREAL8TimeSeries
}

# dummy class needed for loading LIGOLW files
class LIGOLWContentHandler(ligolw.LIGOLWContentHandler):
    pass

lsctables.use_in(LIGOLWContentHandler)

def legacy_approximant_name(apx):
    """Convert the old style xml approximant name to a name
    and phase_order. Alex: I hate this function. Please delet this when we
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
    

class InjectionSet(object):

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

    def getswigrow(self, glue_row):
        """Translates glue row from the table to libmetaio row"""
        import lalmetaio as lmt
        swigrow = lmt.SimInspiralTable()
        for simattr in lsctables.SimInspiralTable.validcolumns.keys():
            if simattr in ["waveform", "source", "numrel_data", "taper"]:
                setattr( swigrow, simattr, str(getattr(glue_row, simattr)) )
            else:
                setattr( swigrow, simattr, getattr(glue_row, simattr) )
        swigrow.geocent_end_time.gpsNanoSeconds = glue_row.geocent_end_time_ns
        return swigrow

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

    def getswigrow(self, glue_row):
        """Translates glue row from the table to libmetaio row"""
        import lalmetaio as lmt
        swigrow = lmt.SimBurstTable()
        for simattr in lsctables.SimBurstTable.validcolumns.keys():
            if simattr in ["waveform"]:
                setattr( swigrow, simattr, str(getattr(glue_row, simattr)) )
            else:
                setattr( swigrow, simattr, getattr(glue_row, simattr) )
        swigrow.geocent_end_time.gpsNanoSeconds = glue_row.geocent_end_time_ns
        return swigrow

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

class RingdownInjectionSet(object):

    """Manages a ringdown injection: reads injection from hdf file
    and injects it into time series.

    Parameters
    ----------
    hdf_file : string
        Path to hdf file containing injection definitions.

    Attributes
    ----------
    table
    """

    def __init__(self, hdf_file):
        self.name = hdf_file

        injfile = h5py.File(hdf_file,'r')
        pnames = map(str, injfile.keys())
        self.table = dict([(p, injfile[p].value) for p in pnames])
        injfile.close()

    def apply(self, strain, detector_name):
        """Add injection (as seen by a particular detector) to a time series.

        Parameters
        ----------
        strain : TimeSeries
            Time series to inject signals into, of type float32 or float64.
        detector_name : string
            Name of the detector used for projecting injections.

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

        # pick lalsimulation injection function
        add_injection = injection_func_map[strain.dtype]

        # Read injection parameters (for now, only one injection in file)
        injection = self.table

        signal = self.make_strain_from_inj_object(injection, strain.delta_t,
                     detector_name)
        signal = signal.astype(strain.dtype)
        signal_lal = signal.lal()
        add_injection(lalstrain, signal_lal, None)

        strain.data[:] = lalstrain.data.data[:]

    def make_strain_from_inj_object(self, inj, delta_t, detector_name):
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

        Returns
        --------
        signal : float
            h(t) corresponding to the injection.
        """
        detector = Detector(detector_name)

        # compute the waveform time series
        hp, hc = ringdown_td_approximants[inj['approximant']](
            delta_t=delta_t, **inj)

        hp._epoch += inj['tc']
        hc._epoch += inj['tc']

        # compute the detector response and add it to the strain
        signal = detector.project_wave(hp, hc,
                             inj['ra'], inj['dec'], inj['polarization'])

        return signal
