import logging

import numpy
from igwn_ligolw import ligolw, lsctables
from igwn_ligolw import utils as ligolw_utils

from pycbc import pnutils
from pycbc.constants import GAMMA, MTSUN_SI, PI, TWOPI
from pycbc.io.hdf import HFile
from pycbc.io.ligolw import (
    create_process_table,
    return_empty_sngl,
    return_search_summary,
)
from pycbc.tmpltbank.lambda_mapping import ethinca_order_from_string
from pycbc.waveform import get_waveform_filter_length_in_time as gwflit

logger = logging.getLogger(__name__)


def convert_to_sngl_inspiral_table(params, proc_id):
    """
    Convert a list of m1,m2,spin1z,spin2z values into a basic sngl_inspiral
    table with mass and spin parameters populated and event IDs assigned

    Parameters
    ----------
    params : iterable
        Each entry in the params iterable should be a sequence of
        [mass1, mass2, spin1z, spin2z] in that order
    proc_id : int
        Process ID to add to each row of the sngl_inspiral table

    Returns
    -------
    SnglInspiralTable
        Bank of templates in SnglInspiralTable format

    """
    sngl_inspiral_table = lsctables.SnglInspiralTable.new()
    col_names = ["mass1", "mass2", "spin1z", "spin2z"]

    for values in params:
        tmplt = return_empty_sngl()

        tmplt.process_id = proc_id
        for colname, value in zip(col_names, values):
            setattr(tmplt, colname, value)
        tmplt.mtotal, tmplt.eta = pnutils.mass1_mass2_to_mtotal_eta(
            tmplt.mass1, tmplt.mass2
        )
        tmplt.mchirp, _ = pnutils.mass1_mass2_to_mchirp_eta(tmplt.mass1, tmplt.mass2)
        tmplt.template_duration = 0  # FIXME
        tmplt.event_id = sngl_inspiral_table.get_next_id()
        sngl_inspiral_table.append(tmplt)

    return sngl_inspiral_table


def calculate_ethinca_metric_comps(
    metric_params,
    ethinca_params,
    mass1,
    mass2,
    spin1z=0.0,
    spin2z=0.0,
    full_ethinca=True,
):
    r"""
    Calculate the Gamma components needed to use the ethinca metric.
    At present this outputs the standard TaylorF2 metric over the end time
    and chirp times \tau_0 and \tau_3.
    A desirable upgrade might be to use the \chi coordinates [defined WHERE?]
    for metric distance instead of \tau_0 and \tau_3.
    The lower frequency cutoff is currently hard-coded to be the same as the
    bank layout options f_low and f0 (which must be the same as each other).

    Parameters
    ----------
    metric_params : MetricParameters instance
        Structure holding all the options for construction of the metric
        and the eigenvalues, eigenvectors and covariance matrix
        needed to manipulate the space.
    ethinca_params : EthincaParameters instance
        Structure holding options relevant to the ethinca metric computation.
    mass1 : float
        Mass of the heavier body in the considered template.
    mass2 : float
        Mass of the lighter body in the considered template.
    spin1z : float (optional, default=0)
        Spin of the heavier body in the considered template.
    spin2z : float (optional, default=0)
        Spin of the lighter body in the considered template.
    full_ethinca : boolean (optional, default=True)
        If True calculate the ethinca components in all 3 directions (mass1,
        mass2 and time). If False calculate only the time component (which is
        stored in Gamma0).

    Returns
    -------
    f_max_theor : float
        Value of the upper frequency cutoff given by the template parameters
        and the cutoff formula requested.

    gamma_vals : numpy_array
        Array holding 6 independent metric components in
        (end_time, tau_0, tau_3) coordinates to be stored in the Gamma0-5
        slots of a SnglInspiral object.

    """
    if (float(spin1z) != 0.0 or float(spin2z) != 0.0) and full_ethinca:
        raise NotImplementedError(
            "Ethinca cannot at present be calculated for nonzero component spins!"
        )
    f0 = metric_params.f0
    if f0 != metric_params.f_low:
        raise ValueError(
            "If calculating ethinca the bank f0 value must be equal to f-low!"
        )
    if ethinca_params.f_low is not None and (
        ethinca_params.f_low != metric_params.f_low
    ):
        raise NotImplementedError(
            "An ethinca metric f-low different from the"
            " bank metric f-low is not supported!"
        )

    twice_pn_order = ethinca_order_from_string(ethinca_params.pn_order)

    pi_fl = PI * f0
    total_mass, eta = pnutils.mass1_mass2_to_mtotal_eta(mass1, mass2)
    total_mass = total_mass * MTSUN_SI
    v0cube = total_mass * pi_fl
    v0 = v0cube ** (1.0 / 3.0)

    # Get theoretical cutoff frequency and work out the closest
    # frequency for which moments were calculated
    f_max_theor = pnutils.frequency_cutoff_from_name(
        ethinca_params.cutoff, mass1, mass2, spin1z, spin2z
    )
    f_maxes = list(metric_params.moments["J4"].keys())
    f_max_idx = abs(numpy.array(f_maxes, dtype=float) - f_max_theor).argmin()
    f_max = f_maxes[f_max_idx]

    # Set the appropriate moments
    Js = numpy.zeros([18, 3], dtype=float)
    for i in range(18):
        Js[i, 0] = metric_params.moments["J%d" % (i)][f_max]
        Js[i, 1] = metric_params.moments["log%d" % (i)][f_max]
        Js[i, 2] = metric_params.moments["loglog%d" % (i)][f_max]

    # Compute the time-dependent metric term.
    two_pi_flower_sq = TWOPI * f0 * TWOPI * f0
    gamma_vals = numpy.zeros([6], dtype=float)
    gamma_vals[0] = 0.5 * two_pi_flower_sq * (Js[(1, 0)] - (Js[(4, 0)] * Js[(4, 0)]))

    # If mass terms not required stop here
    if not full_ethinca:
        return f_max_theor, gamma_vals

    # 3pN is a mess, so split it into pieces
    a0 = 11583231236531 / 200286535680 - 5 * PI * PI - 107 * GAMMA / 14
    a1 = (-15737765635 / 130056192 + 2255 * PI * PI / 512) * eta
    a2 = (76055 / 73728) * eta * eta
    a3 = (-127825 / 55296) * eta * eta * eta
    alog = numpy.log(4 * v0)  # Log terms are tricky - be careful

    # Get the Psi coefficients
    Psi = [{}, {}]  # Psi = numpy.zeros([2,8,2],dtype=float)
    Psi[0][0, 0] = 3 / 5
    Psi[0][2, 0] = (743 / 756 + 11 * eta / 3) * v0 * v0
    Psi[0][3, 0] = 0.0
    Psi[0][4, 0] = (
        (-3058673 / 508032 + 5429 * eta / 504 + 617 * eta * eta / 24) * v0cube * v0
    )
    Psi[0][5, 1] = (-7729 * PI / 126) * v0cube * v0 * v0 / 3
    Psi[0][6, 0] = (
        (128 / 15)
        * (-3 * a0 - a1 + a2 + 3 * a3 + 107 * (1 + 3 * alog) / 14)
        * v0cube
        * v0cube
    )
    Psi[0][6, 1] = (6848 / 35) * v0cube * v0cube / 3
    Psi[0][7, 0] = (-15419335 / 63504 - 75703 * eta / 756) * PI * v0cube * v0cube * v0

    Psi[1][0, 0] = 0.0
    Psi[1][2, 0] = (3715 / 12096 - 55 * eta / 96) / PI / v0
    Psi[1][3, 0] = -3 / 2
    Psi[1][4, 0] = (
        (15293365 / 4064256 - 27145 * eta / 16128 - 3085 * eta * eta / 384) * v0 / PI
    )
    Psi[1][5, 1] = (193225 / 8064) * v0 * v0 / 3
    Psi[1][6, 0] = (
        (4 / PI)
        * (2 * a0 + a1 / 3 - 4 * a2 / 3 - 3 * a3 - 107 * (1 + 6 * alog) / 42)
        * v0cube
    )
    Psi[1][6, 1] = (-428 / PI / 7) * v0cube / 3
    Psi[1][7, 0] = (
        (77096675 / 1161216 + 378515 * eta / 24192 + 74045 * eta * eta / 8064)
        * v0cube
        * v0
    )

    # Set the appropriate moments
    Js = numpy.zeros([18, 3], dtype=float)
    for i in range(18):
        Js[i, 0] = metric_params.moments["J%d" % (i)][f_max]
        Js[i, 1] = metric_params.moments["log%d" % (i)][f_max]
        Js[i, 2] = metric_params.moments["loglog%d" % (i)][f_max]

    # Calculate the g matrix
    pn_terms = [(0, 0), (2, 0), (3, 0), (4, 0), (5, 1), (6, 0), (6, 1), (7, 0)]
    pn_terms = [term for term in pn_terms if term[0] <= twice_pn_order]

    # Now can compute the mass-dependent gamma values
    for m in [0, 1]:
        for k in pn_terms:
            gamma_vals[1 + m] += (
                0.5
                * two_pi_flower_sq
                * Psi[m][k]
                * (Js[(9 - k[0], k[1])] - Js[(12 - k[0], k[1])] * Js[(4, 0)])
            )

    g = numpy.zeros([2, 2], dtype=float)
    for m, n in [(0, 0), (0, 1), (1, 1)]:
        for k in pn_terms:
            for l in pn_terms:
                g[m, n] += (
                    Psi[m][k]
                    * Psi[n][l]
                    * (
                        Js[(17 - k[0] - l[0], k[1] + l[1])]
                        - Js[(12 - k[0], k[1])] * Js[(12 - l[0], l[1])]
                    )
                )
        g[m, n] = 0.5 * two_pi_flower_sq * g[m, n]
        g[n, m] = g[m, n]

    gamma_vals[3] = g[0, 0]
    gamma_vals[4] = g[0, 1]
    gamma_vals[5] = g[1, 1]

    return f_max_theor, gamma_vals


def output_sngl_inspiral_table(
    output_file, temp_bank, program_name="", opt_dict=None, outdoc=None, **kwargs
):  # pylint:disable=unused-argument
    """
    Function that converts the information produced by the various PyCBC bank
    generation codes into a valid LIGOLW XML file containing a sngl_inspiral
    table and outputs to file.

    Parameters
    ----------
    output_file : string
        Name of the file that the bank will be written to
    temp_bank : iterable
        Each entry in the temp_bank iterable should be a sequence of
        [mass1,mass2,spin1z,spin2z] in that order.
    program_name (key-word-argument) : string
        Name of the executable that has been run
    opt_dict (key-word argument) : dictionary
        Dictionary of the command line arguments passed to the program
    outdoc (key-word argument) : ligolw xml document
        If given add template bank to this representation of a xml document and
        write to disk. If not given create a new document.
    kwargs : optional key-word arguments
        Allows unused options to be passed to this function (for modularity)

    """
    if opt_dict is None:
        opt_dict = {}
    if outdoc is None:
        outdoc = ligolw.Document()
        outdoc.appendChild(ligolw.LIGO_LW())

    # get IFO to put in search summary table
    ifos = []
    if "channel_name" in opt_dict.keys():
        if opt_dict["channel_name"] is not None:
            ifos = [opt_dict["channel_name"][0:2]]

    proc = create_process_table(
        outdoc, program_name=program_name, detectors=ifos, options=opt_dict
    )
    proc_id = proc.process_id
    sngl_inspiral_table = convert_to_sngl_inspiral_table(temp_bank, proc_id)

    # set per-template low-frequency cutoff
    if (
        "f_low_column" in opt_dict
        and "f_low" in opt_dict
        and opt_dict["f_low_column"] is not None
    ):
        for sngl in sngl_inspiral_table:
            setattr(sngl, opt_dict["f_low_column"], opt_dict["f_low"])

    outdoc.childNodes[0].appendChild(sngl_inspiral_table)

    # get times to put in search summary table
    start_time = 0
    end_time = 0
    if "gps_start_time" in opt_dict.keys() and "gps_end_time" in opt_dict.keys():
        start_time = opt_dict["gps_start_time"]
        end_time = opt_dict["gps_end_time"]

    # make search summary table
    search_summary_table = lsctables.SearchSummaryTable.new()
    search_summary = return_search_summary(
        start_time, end_time, len(sngl_inspiral_table), ifos
    )
    search_summary_table.append(search_summary)
    outdoc.childNodes[0].appendChild(search_summary_table)

    # write the xml doc to disk
    ligolw_utils.write_filename(outdoc, output_file)


def output_bank_to_hdf(
    output_file,
    temp_bank,
    opt_dict=None,
    program_name="",
    approximant=None,
    output_duration=False,
    **kwargs,
):  # pylint:disable=unused-argument
    """
    Function that converts the information produced by the various PyCBC bank
    generation codes into a hdf5 file.

    Parameters
    ----------
    output_file : string
        Name of the file that the bank will be written to
    temp_bank : iterable
        Each entry in the temp_bank iterable should be a sequence of
        [mass1,mass2,spin1z,spin2z] in that order.
    program_name (key-word-argument) : string
        Name of the executable that has been run
    opt_dict (key-word argument) : dictionary
        Dictionary of the command line arguments passed to the program
    approximant : string
        The approximant to be outputted to the file,
        if output_duration is True, this is also used for that calculation.
    output_duration : boolean
        Output the duration of the template, calculated using
        get_waveform_filter_length_in_time, to the file.
    kwargs : optional key-word arguments
        Allows unused options to be passed to this function (for modularity)

    """
    bank_dict = {}
    mass1, mass2, spin1z, spin2z = list(zip(*temp_bank))
    bank_dict["mass1"] = mass1
    bank_dict["mass2"] = mass2
    bank_dict["spin1z"] = spin1z
    bank_dict["spin2z"] = spin2z

    # Add other values to the bank dictionary as appropriate
    if opt_dict is not None:
        bank_dict["f_lower"] = numpy.ones_like(mass1) * opt_dict["f_low"]
        argument_string = [f"{k}:{v}" for k, v in opt_dict.items()]

    if opt_dict is not None and opt_dict["output_f_final"]:
        bank_dict["f_final"] = numpy.ones_like(mass1) * opt_dict["f_upper"]

    if approximant:
        if not isinstance(approximant, bytes):
            appx = approximant.encode()
        bank_dict["approximant"] = numpy.repeat(appx, len(mass1))

    if output_duration:
        appx = approximant or "SPAtmplt"
        tmplt_durations = numpy.zeros_like(mass1)
        for i in range(len(mass1)):
            wfrm_length = gwflit(
                appx,
                mass1=mass1[i],
                mass2=mass2[i],
                f_lower=opt_dict["f_low"],
                phase_order=7,
            )
            tmplt_durations[i] = wfrm_length
        bank_dict["template_duration"] = tmplt_durations

    with HFile(output_file, "w") as bankf_out:
        bankf_out.attrs["program"] = program_name
        if opt_dict is not None:
            bankf_out.attrs["arguments"] = argument_string
        for k, v in bank_dict.items():
            bankf_out[k] = v


def output_bank_to_file(output_file, temp_bank, **kwargs):
    if output_file.endswith((".xml", ".xml.gz", ".xmlgz")):
        output_sngl_inspiral_table(output_file, temp_bank, **kwargs)
    elif output_file.endswith((".h5", ".hdf", ".hdf5")):
        output_bank_to_hdf(output_file, temp_bank, **kwargs)
    else:
        err_msg = f"Unrecognized extension for file {output_file}."
        raise ValueError(err_msg)
