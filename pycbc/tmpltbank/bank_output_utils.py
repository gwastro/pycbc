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
    outputFile, tempBank, programName="", optDict=None, outdoc=None, **kwargs
):  # pylint:disable=unused-argument
    """
    Function that converts the information produced by the various PyCBC bank
    generation codes into a valid LIGOLW XML file containing a sngl_inspiral
    table and outputs to file.

    Parameters
    ----------
    outputFile : string
        Name of the file that the bank will be written to
    tempBank : iterable
        Each entry in the tempBank iterable should be a sequence of
        [mass1,mass2,spin1z,spin2z] in that order.
    programName (key-word-argument) : string
        Name of the executable that has been run
    optDict (key-word argument) : dictionary
        Dictionary of the command line arguments passed to the program
    outdoc (key-word argument) : ligolw xml document
        If given add template bank to this representation of a xml document and
        write to disk. If not given create a new document.
    kwargs : optional key-word arguments
        Allows unused options to be passed to this function (for modularity)

    """
    if optDict is None:
        optDict = {}
    if outdoc is None:
        outdoc = ligolw.Document()
        outdoc.appendChild(ligolw.LIGO_LW())

    # get IFO to put in search summary table
    ifos = []
    if "channel_name" in optDict.keys():
        if optDict["channel_name"] is not None:
            ifos = [optDict["channel_name"][0:2]]

    proc = create_process_table(
        outdoc, program_name=programName, detectors=ifos, options=optDict
    )
    proc_id = proc.process_id
    sngl_inspiral_table = convert_to_sngl_inspiral_table(tempBank, proc_id)

    # set per-template low-frequency cutoff
    if (
        "f_low_column" in optDict
        and "f_low" in optDict
        and optDict["f_low_column"] is not None
    ):
        for sngl in sngl_inspiral_table:
            setattr(sngl, optDict["f_low_column"], optDict["f_low"])

    outdoc.childNodes[0].appendChild(sngl_inspiral_table)

    # get times to put in search summary table
    start_time = 0
    end_time = 0
    if "gps_start_time" in optDict.keys() and "gps_end_time" in optDict.keys():
        start_time = optDict["gps_start_time"]
        end_time = optDict["gps_end_time"]

    # make search summary table
    search_summary_table = lsctables.SearchSummaryTable.new()
    search_summary = return_search_summary(
        start_time, end_time, len(sngl_inspiral_table), ifos
    )
    search_summary_table.append(search_summary)
    outdoc.childNodes[0].appendChild(search_summary_table)

    # write the xml doc to disk
    ligolw_utils.write_filename(outdoc, outputFile)


def output_bank_to_hdf(
    outputFile,
    tempBank,
    optDict=None,
    programName="",
    approximant=None,
    output_duration=False,
    **kwargs,
):  # pylint:disable=unused-argument
    """
    Function that converts the information produced by the various PyCBC bank
    generation codes into a hdf5 file.

    Parameters
    ----------
    outputFile : string
        Name of the file that the bank will be written to
    tempBank : iterable
        Each entry in the tempBank iterable should be a sequence of
        [mass1,mass2,spin1z,spin2z] in that order.
    programName (key-word-argument) : string
        Name of the executable that has been run
    optDict (key-word argument) : dictionary
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
    mass1, mass2, spin1z, spin2z = list(zip(*tempBank))
    bank_dict["mass1"] = mass1
    bank_dict["mass2"] = mass2
    bank_dict["spin1z"] = spin1z
    bank_dict["spin2z"] = spin2z

    # Add other values to the bank dictionary as appropriate
    if optDict is not None:
        bank_dict["f_lower"] = numpy.ones_like(mass1) * optDict["f_low"]
        argument_string = [f"{k}:{v}" for k, v in optDict.items()]

    if optDict is not None and optDict["output_f_final"]:
        bank_dict["f_final"] = numpy.ones_like(mass1) * optDict["f_upper"]

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
                f_lower=optDict["f_low"],
                phase_order=7,
            )
            tmplt_durations[i] = wfrm_length
        bank_dict["template_duration"] = tmplt_durations

    with HFile(outputFile, "w") as bankf_out:
        bankf_out.attrs["program"] = programName
        if optDict is not None:
            bankf_out.attrs["arguments"] = argument_string
        for k, v in bank_dict.items():
            bankf_out[k] = v


def output_bank_to_file(outputFile, tempBank, **kwargs):
    if outputFile.endswith((".xml", ".xml.gz", ".xmlgz")):
        output_sngl_inspiral_table(outputFile, tempBank, **kwargs)
    elif outputFile.endswith((".h5", ".hdf", ".hdf5")):
        output_bank_to_hdf(outputFile, tempBank, **kwargs)
    else:
        err_msg = f"Unrecognized extension for file {outputFile}."
        raise ValueError(err_msg)
