# Copyright (C) 2017 Christopher M. Biwer
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# self.option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
""" This modules defines functions for reading and samples that the
inference samplers generate and are stored in an ASCII TXT file.
"""

import numpy


class InferenceTXTFile(object):
    """ A class that has extra functions for handling reading the samples
    from posterior-only TXT files.

    Parameters
    -----------
    path : str
        The path to the TXT file.
    mode : {None, str}
        The mode to open the file. Only accepts "r" or "rb" for reading.
    delimiter : str
        Delimiter to use for TXT file. Default is space-delimited.
    """
    name = "txt"
    comments = ""
    delimiter = " "

    def __init__(self, path, mode=None, delimiter=None):
        self.path = path
        self.delimiter = delimiter if delimiter is not None else self.delimiter
        if mode in ["r", "rb"]:
            self.mode = mode
        else:
            raise ValueError("Mode for InferenceTXTFile must be 'r' or 'rb'.")

    @classmethod
    def write(cls, output_file, samples, labels, delimiter=None):
        """ Writes a text file with samples.

        Parameters
        -----------
        output_file : str
            The path of the file to write.
        samples : FieldArray
            Samples to write to file.
        labels : list
            A list of strings to include as header in TXT file.
        delimiter : str
            Delimiter to use in TXT file.
        """
        delimiter = delimiter if delimiter is not None else cls.delimiter
        header = delimiter.join(labels)
        numpy.savetxt(output_file, samples,
                      comments=cls.comments, header=header,
                      delimiter=delimiter)
