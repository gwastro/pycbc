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
"""This module provides utilities for injecting signals into data
"""

import numpy as np
import lal
import lalsimulation as sim
from pycbc.waveform import get_td_waveform
from glue.ligolw import utils as ligolw_utils
from glue.ligolw import table, lsctables
from pycbc.types import float64, float32, TimeSeries
from pycbc.detector import Detector

# map between tapering string in sim_inspiral
# table and lalsimulation constants
taper_map = {
    'TAPER_NONE': None,
    'TAPER_START': sim.LAL_SIM_INSPIRAL_TAPER_START,
    'TAPER_END': sim.LAL_SIM_INSPIRAL_TAPER_END,
    'TAPER_STARTEND': sim.LAL_SIM_INSPIRAL_TAPER_STARTEND
}

taper_func_map = {
    np.dtype(float32): sim.SimInspiralREAL4WaveTaper,
    np.dtype(float64): sim.SimInspiralREAL8WaveTaper
}

injection_func_map = {
    np.dtype(float32): sim.SimAddInjectionREAL4TimeSeries,
    np.dtype(float64): sim.SimAddInjectionREAL8TimeSeries
}

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
        self.indoc = ligolw_utils.load_filename(sim_file, False)
        self.table = table.get_table(self.indoc, lsctables.SimInspiralTable.tableName)
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

        if not strain.dtype in (float32, float64):
            raise TypeError("Strain dtype must be float32 or float64, not " \
                    + str(strain.dtype))

        lalstrain = strain.lal()    
        detector = Detector(detector_name)
        earth_travel_time = lal.LAL_REARTH_SI / lal.LAL_C_SI
        t0 = float(strain.start_time) - earth_travel_time
        t1 = float(strain.end_time) + earth_travel_time

        # pick lalsimulation tapering function
        taper = taper_func_map[strain.dtype]

        # pick lalsimulation injection function
        add_injection = injection_func_map[strain.dtype]

        for inj in self.table:
            if f_lower is None:
                f_l = inj.f_lower
            else:
                f_l = f_lower

            # roughly estimate if the injection may overlap with the segment
            end_time = inj.get_time_geocent()
            inj_length = sim.SimInspiralTaylorLength(
                    strain.delta_t, inj.mass1 * lal.LAL_MSUN_SI,
                    inj.mass2 * lal.LAL_MSUN_SI, f_l, 0)
            start_time = end_time - 2 * inj_length
            if end_time < t0 or start_time > t1:
                continue

            # compute the waveform time series
            hp, hc = get_td_waveform(
                    inj, approximant=inj.waveform, delta_t=strain.delta_t,
                    f_lower=f_l, distance=inj.distance * distance_scale, **self.extra_args)
            hp._epoch += float(end_time)
            hc._epoch += float(end_time)
            if float(hp.start_time) > t1:
                continue

            # compute the detector response, taper it if requested
            # and add it to the strain
            signal = detector.project_wave(
                    hp, hc, inj.longitude, inj.latitude, inj.polarization)
            signal_lal = signal.astype(strain.dtype).lal()
            if taper_map[inj.taper] is not None:
                taper(signal_lal.data, taper_map[inj.taper])
            add_injection(lalstrain, signal_lal, None)

        strain.data[:] = lalstrain.data.data[:]
        
    def maximum_snrs(self, psd):
        """ Calculate the maximum SNR (without noise), for each signal in the
        injection set, using the given PSD.        

        Parameters
        ----------


        Returns
        -------
        None

        """
        raise NotImplementedError("This is not yet implemented")
        
    def expected_snrs(self, template, psd):
        """ Calculate the expected SNR (without noise), for each signal in the
        injection set, using the given PSD, when matched with the given template.       

        Parameters
        ----------


        Returns
        -------
        None

        """
        raise NotImplementedError("This is not yet implemented")
    
    

