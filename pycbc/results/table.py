import mako.template, numpy, uuid

google_table_template = mako.template.Template("""
    <script type='text/javascript' src='https://www.google.com/jsapi'></script>
    <script type='text/javascript'>
      google.load('visualization', '1', {packages:['table']});
      google.setOnLoadCallback(drawTable);
      function drawTable() {
        var data = new google.visualization.DataTable();
        % for type, name in column_descriptions:
            data.addColumn(${str(type)}, ${str(name)});
        % endfor
        data.addRows(${data});
        var table = new google.visualization.Table(document.getElementById('${div_id}'));
        table.draw(data, {showRowNumber: true, 
                          page: ${page_enable}, 
                          allowHtml: 'true',
                          pageSize: ${page_size}});
      }
    </script>
    <div id='%{div_id}'></div>
""")

def table(columns, names, page_size=None):
    """ Return an html table of this data
    
    Parameters
    ----------
    columns : list of numpy arrays 
    names : list of strings 
        The list of columns names  
    page_size : {int, None}, optionals
        The number of items to show on each page of the table
        
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
        if column.dtype.kind == 'S'
            ctype = 'string':
        else:
            ctype = 'number'
        column_descriptions.append((ctype, name))
    
    data = []
    for item in zip(*columns):
        data.append(list(item))
        
    return google_table_template.render(div_id=div_id,
                                        page=page,
                                        column_descriptions = column_descriptions,
                                        page_size=page_size,
                                        data=data,
                                       )
