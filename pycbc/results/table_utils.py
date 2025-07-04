# Copyright (C) 2014 Alex Nitz
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
""" This module provides functions to generate sortable html tables
"""
import mako.template
import uuid
import copy
import numpy

google_table_template = mako.template.Template("""
    <script type='text/javascript' src='https://www.google.com/jsapi'></script>
    <script type='text/javascript'>
      google.load('visualization', '1', {packages:['table']});
      google.setOnLoadCallback(drawTable);
      function drawTable() {
        var data = new google.visualization.DataTable();
        % for type, name in column_descriptions:
            data.addColumn('${str(type)}', '${str(name)}');
        % endfor
        data.addRows(${data});

        % if format_strings is not None:
            % for i, format_string in enumerate(format_strings):
                % if format_string is not None:
                    var formatter = new google.visualization.NumberFormat({pattern:'${format_string}'});
                    formatter.format(data, ${i});
                % endif
            % endfor
        % endif
        var table = new google.visualization.Table(document.getElementById('${div_id}'));
        table.draw(data, {showRowNumber: 'true',
                          page: '${page_enable}',
                          allowHtml: 'true',
                          pageSize: ${page_size}});
      }
    </script>
    <div id='${div_id}'></div>
""")

def html_table(columns, names, page_size=None, format_strings=None):
    """ Return an html table of this data

    Parameters
    ----------
    columns : list of numpy arrays
    names : list of strings
        The list of columns names
    page_size : {int, None}, optional
        The number of items to show on each page of the table
    format_strings : {lists of strings, None}, optional
        The ICU format string for this column, None for no formatting. All
    columns must have a format string if provided.

    Returns
    -------
    html_table : str
        A str containing the html code to display a table of this data
    """
    if page_size is None:
        page = 'disable'
    else:
        page = 'enable'

    div_id = uuid.uuid4()

    column_descriptions = []
    for column, name in zip(columns, names):
        if column.dtype.kind in 'iuf':
            # signed and unsigned integers and floats
            ctype = 'number'
        else:
            # this comprises strings, bools, complex, void, etc
            # but we will convert all those to str in a moment
            ctype = 'string'
        column_descriptions.append((ctype, name))

    data = []
    for row in zip(*columns):
        row2 = []
        # the explicit conversions here are to make sure the JS code
        # sees proper numbers and not something like 'np.float64(12)'
        for item, column in zip(row, columns):
            if column.dtype.kind == 'f':
                row2.append(float(item))
            elif column.dtype.kind in 'iu':
                row2.append(int(item))
            else:
                row2.append(str(item))
        data.append(row2)

    return google_table_template.render(
        div_id=div_id,
        page_enable=page,
        column_descriptions=column_descriptions,
        page_size=page_size,
        data=data,
        format_strings=format_strings,
    )

static_table_template = mako.template.Template("""
    <table class="table">
        % for row in range(n_rows):
            % if titles is not None:
                <tr>
                % if row_labels is not None:
                    <td>
                    </td>
                % endif
                % for i in range(n_columns):
                    <th>
                        ${titles[row * n_columns + i]}
                    </th>
                % endfor
                </tr>
            % endif

            % for i in range(len(data)):
                <tr>
                % if row_labels is not None:
                    <td>
                        ${row_labels[i]}
                    </td>
                % endif
                % for j in range(n_columns):
                    <td>
                        ${data[i][row * n_columns + j]}
                    </td>
                % endfor
                </tr>
            % endfor
        % endfor
    </table>
""")

def static_table(data, titles=None, columns_max=None, row_labels=None):
    """ Return an html table of this data

    Parameters
    ----------
    data : two-dimensional string array
        Array containing the cell values
    titles : numpy array
        Vector str of titles, must be the same length as data
    columns_max : integer or None
        If given, will restrict the number of columns in the table
    row_labels : list of strings
        Optional list of row labels to be given as the first cell in
        each data row. Does not count towards columns_max

    Returns
    -------
    html_table : str
        A string containing the html table.
    """
    data = copy.deepcopy(data)
    titles = copy.deepcopy(titles)
    row_labels = copy.deepcopy(row_labels)
    drows, dcols = numpy.array(data).shape
    if titles is not None and not len(titles) == dcols:
        raise ValueError("titles and data lengths do not match")

    if row_labels is not None and not len(row_labels) == drows:
        raise ValueError(
            "row_labels must be the same number of rows supplied to data"
        )

    if columns_max is not None:
        n_rows = int(numpy.ceil(len(data[0]) / columns_max))
        n_columns = min(columns_max, len(data[0]))
        if len(data[0]) < n_rows * n_columns:
            # Pad the data and titles with empty strings
            n_missing = int(n_rows * n_columns - len(data[0]))
            data = numpy.hstack((data, numpy.zeros((len(data), n_missing), dtype='U1')))
            if titles is not None:
                titles += [' '] * n_missing
    else:
        n_rows = 1
        n_columns = len(data[0])

    return static_table_template.render(
        data=data,
        titles=titles,
        n_columns=n_columns,
        n_rows=n_rows,
        row_labels=row_labels,
    )

