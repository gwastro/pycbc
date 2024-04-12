# Copyright (C) 2023 Arthur Tolley, Gareth Cabourn Davies
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
""" This module contains functions for optimizing the signal-to-noise ratio
of triggers produced by PyCBC Live. Also contained within this module are the
command line arguments required and options group for the SNR optimization.
This module is primarily used in the pycbc_optimize_snr program.
"""

import logging
import subprocess
from datetime import datetime
from dateutil.relativedelta import relativedelta
import time

logger = logging.getLogger('pycbc.live.supervision')

def symlink(target, link_name):
    """
    Create a symbolic link replacing the destination and checking for
    errors.
    """
    cp = subprocess.run([
        'ln', '-sf', target, link_name
    ])
    if cp.returncode:
        raise subprocess.SubprocessError(
            f"Could not link {target} to {link_name}"
        )


def dict_to_args(opts_dict):
    """
    Convert an option dictionary into a list to be used by subprocess.run
    """
    dargs = []
    for option in opts_dict.keys():
        dargs.append('--' + option.strip())
        value = opts_dict[option]
        if len(value.split()) > 1:
            # value is a list, append individually
            for v in value.split():
                dargs.append(v)
        elif not value:
            # option is a flag, do nothing
            continue
        else:
            # Single value option - easy enough
            dargs.append(value)
    return dargs


def mail_volunteers_error(controls, mail_body_lines, subject):
    """
    Email a list of people, defined by mail-volunteers-file
    To be used for errors or unusual occurences
    """
    with open(controls['mail-volunteers-file'], 'r') as mail_volunteers_file:
        volunteers = [volunteer.strip() for volunteer in
                      mail_volunteers_file.readlines()]
    logger.info("Emailing %s with warnings", ' '.join(volunteers))
    mail_command = [
        'mail',
        '-s',
        subject
    ]
    mail_command += volunteers
    mail_body = '\n'.join(mail_body_lines)
    subprocess.run(mail_command, input=mail_body, text=True)

def run_and_error(command_arguments):
    """
    Wrapper around subprocess.run to catch errors and send emails if required
    """
    logger.info("Running " + " ".join(command_arguments))
    command_output = subprocess.run(command_arguments, capture_output=True)
    if command_output.returncode:
        error_contents = [' '.join(command_arguments),
                          command_output.stderr.decode()]
        mail_volunteers_error(controls, error_contents,
            f"PyCBC live could not run {command_arguments[0]}")
        err_msg = f"Could not run {command_arguments[0]}"
        raise subprocess.SubprocessError(err_msg)

def wait_for_utc_time(target_str):
    """Wait until the UTC time is as given by `target_str`, in HH:MM:SS format.
    """
    target_hour, target_minute, target_second = map(int, target_str.split(':'))
    now = datetime.utcnow()
    # for today's target, take now and replace the time
    target_today = now + relativedelta(
        hour=target_hour, minute=target_minute, second=target_second
    )
    # for tomorrow's target, take now, add one day, and replace the time
    target_tomorrow = now + relativedelta(
        days=1, hour=target_hour, minute=target_minute, second=target_second
    )
    next_target = target_today if now <= target_today else target_tomorrow
    sleep_seconds = (next_target - now).total_seconds()
    logger.info('Waiting %.0f s', sleep_seconds)
    time.sleep(sleep_seconds)
