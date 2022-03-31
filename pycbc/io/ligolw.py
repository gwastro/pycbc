# Copyright (C) 2020 Leo Singer, 2021 Tito Dal Canton
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""Tools for dealing with LIGOLW XML files."""

import os
import sys
from ligo.lw import lsctables
from ligo.lw.ligolw import Param, LIGOLWContentHandler \
    as OrigLIGOLWContentHandler
from ligo.lw.lsctables import TableByName
from ligo.lw.table import Column, TableStream
from ligo.lw.types import FormatFunc, FromPyType, ToPyType
from ligo.lw.utils import process as ligolw_process
import pycbc.version as pycbc_version


__all__ = ('default_null_value',
           'return_empty_sngl',
           'return_search_summary',
           'legacy_row_id_converter',
           'LIGOLWContentHandler')

ROWID_PYTYPE = int
ROWID_TYPE = FromPyType[ROWID_PYTYPE]
ROWID_FORMATFUNC = FormatFunc[ROWID_TYPE]
IDTypes = set([u"ilwd:char", u"ilwd:char_u"])


def default_null_value(col_name, col_type):
    """
    Associate a sensible "null" default value to a given LIGOLW column type.
    """
    if col_type in ['real_4', 'real_8']:
        return 0.
    if col_type in ['int_4s', 'int_8s']:
        # this case includes row IDs
        return 0
    if col_type == 'lstring':
        return ''
    raise NotImplementedError(('Do not know how to initialize column '
                               '{} of type {}').format(col_name, col_type))

def return_empty_sngl(nones=False):
    """
    Function to create a SnglInspiral object where all columns are populated
    but all are set to values that test False (ie. strings to '', floats/ints
    to 0, ...). This avoids errors when you try to create a table containing
    columns you don't care about, but which still need populating. NOTE: This
    will also produce a process_id and event_id with 0 values. For most
    applications these should be set to their correct values.

    Parameters
    ----------
    nones : bool (False)
        If True, just set all columns to None.

    Returns
    --------
    lsctables.SnglInspiral
        The "empty" SnglInspiral object.
    """

    sngl = lsctables.SnglInspiral()
    cols = lsctables.SnglInspiralTable.validcolumns
    for entry in cols:
        col_name = Column.ColumnName(entry)
        value = None if nones else default_null_value(col_name, cols[entry])
        setattr(sngl, col_name, value)
    return sngl

def return_search_summary(start_time=0, end_time=0, nevents=0, ifos=None):
    """
    Function to create a SearchSummary object where all columns are populated
    but all are set to values that test False (ie. strings to '', floats/ints
    to 0, ...). This avoids errors when you try to create a table containing
    columns you don't care about, but which still need populating. NOTE: This
    will also produce a process_id with 0 values. For most applications these
    should be set to their correct values.

    It then populates columns if given them as options.

    Returns
    --------
    lsctables.SeachSummary
        The "empty" SearchSummary object.
    """
    if ifos is None:
        ifos = []

    # create an empty search summary
    search_summary = lsctables.SearchSummary()
    cols = lsctables.SearchSummaryTable.validcolumns
    for entry in cols:
        col_name = Column.ColumnName(entry)
        value = default_null_value(col_name, cols[entry])
        setattr(search_summary, col_name, value)

    # fill in columns
    if ifos:
        search_summary.instruments = ifos
    if nevents:
        search_summary.nevents = nevents
    if start_time and end_time:
        search_summary.in_start_time = int(start_time)
        search_summary.in_start_time_ns = int(start_time % 1 * 1e9)
        search_summary.in_end_time = int(end_time)
        search_summary.in_end_time_ns = int(end_time % 1 * 1e9)
        search_summary.out_start_time = int(start_time)
        search_summary.out_start_time_ns = int(start_time % 1 * 1e9)
        search_summary.out_end_time = int(end_time)
        search_summary.out_end_time_ns = int(end_time % 1 * 1e9)

    return search_summary

def create_process_table(document, program_name=None, detectors=None,
                         comment=None):
    """Create a LIGOLW process table with sane defaults, add it to a LIGOLW
    document, and return it.
    """
    if program_name is None:
        program_name = os.path.basename(sys.argv[0])

    # ligo.lw does not like `cvs_entry_time` being an empty string
    cvs_entry_time = pycbc_version.date or None

    process = ligolw_process.register_to_xmldoc(
            document, program_name, {}, version=pycbc_version.version,
            cvs_repository='pycbc/'+pycbc_version.git_branch,
            cvs_entry_time=cvs_entry_time, instruments=detectors,
            comment=comment)
    return process

def legacy_row_id_converter(ContentHandler):
    """Convert from old-style to new-style row IDs on the fly.

    This is loosely adapted from :func:`ligo.lw.utils.ilwd.strip_ilwdchar`.

    Notes
    -----
    When building a ContentHandler, this must be the _outermost_ decorator,
    outside of :func:`ligo.lw.lsctables.use_in`, :func:`ligo.lw.param.use_in`,
    or :func:`ligo.lw.table.use_in`.
    """

    def endElementNS(self, uri_localname, qname,
                     __orig_endElementNS=ContentHandler.endElementNS):
        """Convert values of <Param> elements from ilwdchar to int."""
        if isinstance(self.current, Param) and self.current.Type in IDTypes:
            old_type = ToPyType[self.current.Type]
            old_val = str(old_type(self.current.pcdata))
            new_value = ROWID_PYTYPE(old_val.split(":")[-1])
            self.current.Type = ROWID_TYPE
            self.current.pcdata = ROWID_FORMATFUNC(new_value)
        __orig_endElementNS(self, uri_localname, qname)

    remapped = {}

    def startColumn(self, parent, attrs,
                    __orig_startColumn=ContentHandler.startColumn):
        """Convert types in <Column> elements from ilwdchar to int.

        Notes
        -----
        This method is adapted from
        :func:`ligo.lw.utils.ilwd.strip_ilwdchar`.

        """
        result = __orig_startColumn(self, parent, attrs)

        # If this is an ilwdchar column, then create a function to convert its
        # rows' values for use in the startStream method below.
        if result.Type in IDTypes:
            old_type = ToPyType[result.Type]

            def converter(old_value):
                return ROWID_PYTYPE(str(old_type(old_value)).split(":")[-1])

            remapped[(id(parent), result.Name)] = converter
            result.Type = ROWID_TYPE

        # If this is an ilwdchar column, then normalize the column name.
        if parent.Name in TableByName:
            validcolumns = TableByName[parent.Name].validcolumns
            if result.Name not in validcolumns:
                stripped_column_to_valid_column = {
                    Column.ColumnName(name): name for name in validcolumns}
                if result.Name in stripped_column_to_valid_column:
                    result.setAttribute(
                        'Name', stripped_column_to_valid_column[result.Name])

        return result

    def startStream(self, parent, attrs,
                    __orig_startStream=ContentHandler.startStream):
        """Convert values in table <Stream> elements from ilwdchar to int.

        Notes
        -----
        This method is adapted from
        :meth:`ligo.lw.table.TableStream.config`.

        """
        result = __orig_startStream(self, parent, attrs)
        if isinstance(result, TableStream):
            loadcolumns = set(parent.columnnames)
            if parent.loadcolumns is not None:
                # FIXME:  convert loadcolumns attributes to sets to
                # avoid the conversion.
                loadcolumns &= set(parent.loadcolumns)
            result._tokenizer.set_types([
                (remapped.pop((id(parent), colname), pytype)
                 if colname in loadcolumns else None)
                for pytype, colname
                in zip(parent.columnpytypes, parent.columnnames)])
        return result

    ContentHandler.endElementNS = endElementNS
    ContentHandler.startColumn = startColumn
    ContentHandler.startStream = startStream

    return ContentHandler


@legacy_row_id_converter
@lsctables.use_in
class LIGOLWContentHandler(OrigLIGOLWContentHandler):
    "Dummy class needed for loading LIGOLW files"
