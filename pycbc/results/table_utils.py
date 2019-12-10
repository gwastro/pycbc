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
import mako.template, uuid

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
        if column.dtype.kind == 'S':
            ctype = 'string'
        else:
            ctype = 'number'
        column_descriptions.append((ctype, name))

    data = []
    for item in zip(*columns):
        data.append(list(item))

    return google_table_template.render(div_id=div_id,
                                page_enable=page,
                                column_descriptions = column_descriptions,
                                page_size=page_size,
                                data=data,
                                format_strings=format_strings,
                               )

static_table_template = mako.template.Template("""
    <table class="table">
        % if titles is not None:
            <tr>
            % for i in range(len(titles)):
                <th>
                    ${titles[i]}
                </th>
            % endfor
            </tr>
        % endif

        % for i in range(len(data)):
            <tr>
            % for j in range(len(data[i])):
                <td>
                    ${data[i][j]}
                </td>
            % endfor
            </tr>
        % endfor
    </table>
""")

def static_table(data, titles=None):
    """ Return an html tableo of this data

    Parameters
    ----------
    data : two-dimensional numpy string array
        Array containing the cell values
    titles : numpy array
        Vector str of titles

    Returns
    -------
    html_table : str
        A string containing the html table.
    """
    return static_table_template.render(data=data, titles=titles)

