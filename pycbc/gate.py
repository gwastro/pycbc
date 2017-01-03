# Copyright (C) 2016 Collin Capano
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
""" Functions for applying gates to data.
"""

import logging
from pycbc import strain

def _gates_from_cli(opts, gate_opt):
    """Parses the given `gate_opt` into something understandable by
    `strain.gate_data`.
    """
    gates = {}
    if getattr(opts, gate_opt) is None:
        return gates
    for gate in getattr(opts, gate_opt):
        try:
            ifo, central_time, half_dur, taper_dur = gate.split(':')
            central_time = float(central_time)
            half_dur = float(half_dur)
            taper_dur = float(taper_dur)
        except ValueError:
            raise ValueError("--gate {} not formatted correctly; ".format(
                gate) + "see help")
        try:
            gates[ifo].append((central_time, half_dur, taper_dur))
        except KeyError:
            gates[ifo] = [(central_time, half_dur, taper_dur)]
    return gates


def gates_from_cli(opts):
    """Parses the --gate option into something understandable by
    `strain.gate_data`.
    """
    return _gates_from_cli(opts, 'gate')


def psd_gates_from_cli(opts):
    """Parses the --psd-gate option into something understandable by
    `strain.gate_data`.
    """
    return _gates_from_cli(opts, 'psd_gate')


def apply_gates_to_td(strain_dict, gates):
    """Applies the given dictionary of gates to the given dictionary of
    strain.

    Parameters
    ----------
    strain_dict : dict
        Dictionary of time-domain strain, keyed by the ifos.
    gates : dict
        Dictionary of gates. Keys should be the ifo to apply the data to,
        values are a tuple giving the central time of the gate, the half
        duration, and the taper duration.

    Returns
    -------
    dict
        Dictionary of time-domain strain with the gates applied.
    """
    # copy data to new dictionary
    outdict = dict(strain_dict.items())
    for ifo in gates:
        logging.info("Gating {} strain".format(ifo))
        outdict[ifo] = strain.gate_data(outdict[ifo], gates[ifo])
    return outdict


def apply_gates_to_fd(stilde_dict, gates):
    """Applies the given dictionary of gates to the given dictionary of
    strain in the frequency domain.

    Gates are applied by IFFT-ing the strain data to the time domain, applying
    the gate, then FFT-ing back to the frequency domain.

    Parameters
    ----------
    stilde_dict : dict
        Dictionary of frequency-domain strain, keyed by the ifos.
    gates : dict
        Dictionary of gates. Keys should be the ifo to apply the data to,
        values are a tuple giving the central time of the gate, the half
        duration, and the taper duration.

    Returns
    -------
    dict
        Dictionary of frequency-domain strain with the gates applied.
    """
    # copy data to new dictionary
    outdict = dict(stilde_dict.items())
    # create a time-domin strain dictionary to apply the gates to
    strain_dict = dict([[ifo, outdict[ifo].to_timeseries()] for ifo in gates])
    # apply gates and fft back to the frequency domain
    for ifo,d in apply_gates_to_td(strain_dict, gates).items():
        outdict[ifo] = d.to_frequencyseries()
    return outdict


def add_gate_option_group(parser):
    """Adds the options needed to apply gates to data.

    Parameters
    ----------
    parser : object
        ArgumentParser instance.
    """
    gate_group = parser.add_argument_group("Options for gating data.")

    gate_group.add_argument("--gate", nargs="+", type=str,
                            metavar="IFO:CENTRALTIME:HALFDUR:TAPERDUR",
                            help="Apply one or more gates to the data before "
                                 "filtering.")
    gate_group.add_argument("--gate-overwhitened", action="store_true",
                            help="Overwhiten data first, then apply the "
                                 "gates specified in --gate. Overwhitening "
                                 "allows for sharper tapers to be used, "
                                 "since lines are not blurred.")
    gate_group.add_argument("--psd-gate", nargs="+", type=str,
                            metavar="IFO:CENTRALTIME:HALFDUR:TAPERDUR",
                            help="Apply one or more gates to the data used "
                                 "for computing the PSD. Gates are applied "
                                 "prior to FFT-ing the data for PSD "
                                 "estimation.")
    return gate_group
