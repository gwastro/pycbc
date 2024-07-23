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

"""
This module contains functions for supervising codes to run regularly
during pycbc_live production, taking input from the search and returning
files which can be used in the search.
This module is primarily used in the pycbc_live_supervise_* programs.
"""

import logging
import subprocess
import time
import os
from datetime import datetime
from dateutil.relativedelta import relativedelta

import pycbc

logger = logging.getLogger('pycbc.live.supervision')


def symlink(target, link_name):
    """
    Create a symbolic link replacing the destination and checking for
    errors.
    """
    # Ensure that the target and link name are absolute paths
    target = os.path.abspath(target)
    link_name = os.path.abspath(link_name)
    logger.info("Linking %s to %s", target, link_name)
    try:
        subprocess.run(['ln', '-sf', target, link_name], check=True)
    except subprocess.CalledProcessError as sub_err:
        logging.error("Could not link %s to %s", target, link_name)
        raise sub_err


def dict_to_args(opts_dict):
    """
    Convert an option dictionary into a list to be used by subprocess.run
    """
    dargs = []
    for option, value in opts_dict.items():
        dargs.append('--' + option.strip())
        if value == '':
            # option is a flag, do nothing
            continue
        if len(value.split()) > 1:
            # value is a list, append individually
            for v in value.split():
                dargs.append(v)
        else:
            # Single value option - append once
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
    try:
        subprocess.run(mail_command, input=mail_body, text=True, check=True)
    except subprocess.CalledProcessError as sub_err:
        logging.error("Could not send mail on error")
        raise sub_err


def run_and_error(command_arguments, controls):
    """
    Wrapper around subprocess.run to catch errors and send emails if required
    """
    logger.info("Running %s", " ".join(command_arguments))
    command_output = subprocess.run(
        command_arguments,
        capture_output=True
    )

    if command_output.returncode:
        error_contents = [' '.join(command_arguments), '\n',
                          command_output.stderr.decode()]
        if 'mail-volunteers-file' in controls:
            mail_volunteers_error(
                controls,
                error_contents,
                f"PyCBC live could not run {command_arguments[0]}"
            )
        err_msg = f"Could not run {command_arguments[0]}:\n"
        err_msg += ' '.join(error_contents)
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


def ensure_directories(control_values, day_str):
    """
    Ensure that the required directories exist
    """
    output_dir = os.path.join(
        control_values['output-directory'],
        day_str
    )
    pycbc.makedir(output_dir)
    if 'public-dir' in control_values:
        # The public directory wil be in subdirectories for the year, month,
        # day, e.g. 2024_04_12 will be in 2024/04/12.
        public_dir = os.path.join(
            control_values['public-dir'],
            *day_str.split('_')
        )
        pycbc.makedir(public_dir)
