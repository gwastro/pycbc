import logging
import os
import pathlib
import datetime
import numpy

from lal import gpstime as lalgps

logger = logging.getLogger('pycbc.io.live')


def maximum_string(numbers):
    """
    Find the maximum possible length string to match
    all values between two numbers

    Parameters
    ----------
    numbers : list of integers
        A list of integers from which to determine the longest
        common string prefix. E.g. '12345', '12346', '12356'
        returns '123'
    """
    # The max length of the number will be the integer above log10
    # of the biggest number
    maxlen = int(numpy.ceil(numpy.log10(max(numbers))))
    # Convert the numbers to (possibly leading zero-padded) strings
    strings = [f"{{n:0{maxlen:d}d}}".format(n=n) for n in numbers]
    # Count how many digits are the same:
    n_digits = 0
    for str_digit in zip(*strings):
        if len(numpy.unique(str_digit)) == 1:
            # This digit is the same for all numbers
            n_digits += 1
        else:
            break
    return strings[0][:n_digits]


def filter_file(filename, start_time, end_time):
    """
    Indicate whether the filename indicates that the file is within the
    start and end times
    Parameters
    ----------
    filename : string
        Filename which matches the format
        {id_string}-{start_time}-{duration}.hdf
    start_time : float
        Start of search window, i.e. GPS time of when the
        file cannot end before
    end_time : float
        End of search window, i.e. GPS time of when the
        file cannot start after

    Returns
    -------
    boolean
        Does any of the file lie within the start/end times
    """
    # FIX ME eventually - this uses the gps time and duration from the filename
    # Is there a better way? (i.e. trigger gps times in the file or
    # add an attribute)
    fend = filename.split('-')[-2:]
    file_start = float(fend[0])
    duration = float(fend[1][:-4])

    return ((file_start + duration) >= start_time) and (file_start <= end_time)


def add_live_trigger_selection_options(parser):
    """
    Add options required for obtaining the right set of PyCBC live triggers
    into an argument parser
    """
    finding_group = parser.add_argument_group('Trigger Finding')
    finding_group.add_argument(
        "--trigger-directory",
        metavar="PATH",
        required=True,
        help="Directory containing trigger files, directory "
             "can contain subdirectories. Required."
    )
    finding_group.add_argument(
        "--gps-start-time",
        type=int,
        required=True,
        help="Start time of the analysis. Integer, required"
    )
    finding_group.add_argument(
        "--gps-end-time",
        type=int,
        required=True,
        help="End time of the analysis. Integer, required"
    )
    finding_group.add_argument(
        "--date-directories",
        action="store_true",
        help="Indicate if the trigger files are stored in "
             "directories by date."
    )
    default_dd_format = "%Y_%m_%d"
    finding_group.add_argument(
        "--date-directory-format",
        default=default_dd_format,
        help="Format of date, see datetime strftime "
             "documentation for details. Default: "
             "%%Y_%%m_%%d"
    )
    finding_group.add_argument(
        "--file-identifier",
        default="H1L1V1-Live",
        help="String required in filename to be considered for "
             "analysis. Default: 'H1L1V1-Live'."
    )


def find_trigger_files(directory, gps_start_time, gps_end_time,
                       id_string='*', date_directories=False,
                       date_directory_format="%Y_%m_%d"):
    """
    Find a list of PyCBC live trigger files which are between the gps
    start and end times given
    """

    # Find the string at the start of the gps time which will match all
    # files in this range - this helps to cut which ones we need to
    # compare later
    num_match = maximum_string([gps_start_time, gps_end_time])

    # ** means recursive, so for large directories, this is expensive.
    # It is not too bad if date_directories is set, as we don't waste time
    # in directories where there cant be any files.
    glob_string = f'**/*{id_string}*{num_match}*.hdf'
    if date_directories:
        # convert the GPS times into dates, and only use the directories
        # of those dates to search
        # Add a day on either side to ensure we get files which straddle
        # the boundary
        one_day = datetime.timedelta(days=1)
        date_check = lalgps.gps_to_utc(gps_start_time) - one_day
        date_end = lalgps.gps_to_utc(gps_end_time) + one_day
        matching_files = []
        while date_check < date_end:
            date_dir = date_check.strftime(date_directory_format)
            subdir = os.path.join(directory, date_dir)
            matching_files_gen = pathlib.Path(subdir).glob(glob_string)
            matching_files += [f.as_posix() for f in matching_files_gen]
            date_check += one_day
    else:
        # Grab all hdf files in the directory
        matching_files_gen = pathlib.Path(directory).glob(glob_string)
        matching_files = [f.as_posix() for f in matching_files_gen]

    # Is the file in the time window?
    matching_files = [f for f in matching_files
                      if filter_file(f, gps_start_time, gps_end_time)]

    return sorted(matching_files)


def find_trigger_files_from_cli(args):
    """
    Wrapper around the find_trigger_files function to use when called using
    options from the add_live_trigger_selection_options function
    """
    return find_trigger_files(
        args.trigger_directory,
        args.gps_start_time,
        args.gps_end_time,
        id_string=args.file_identifier,
        date_directories=args.date_directories,
        date_directory_format=args.date_directory_format
    )


__all__ = [
    'add_live_trigger_selection_options',
    'find_trigger_files',
    'find_trigger_files_from_cli',
]
