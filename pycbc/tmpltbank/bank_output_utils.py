import logging

import numpy
from igwn_ligolw import ligolw, lsctables
from igwn_ligolw import utils as ligolw_utils

from pycbc import pnutils
from pycbc.io.hdf import HFile
from pycbc.io.ligolw import (
    create_process_table,
    return_empty_sngl,
    return_search_summary,
)
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
