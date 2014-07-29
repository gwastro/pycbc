import mako.template, numpy, uuid

google_table_template = mako.template.Template("""
    <script type='text/javascript' src='https://www.google.com/jsapi'></script>
    <script type='text/javascript'>
      google.load('visualization', '1', {packages:['table']});
      google.setOnLoadCallback(drawTable);
      function drawTable() {
        var data = new google.visualization.DataTable();
        data.addColumn('string', 'Name');
        data.addColumn('number', 'Salary');
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

def table(columns, page_size=None):
    """ Return an html table of this data
    
    Parameters
    ----------
    columns : list of numpy arrays 
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
    
    data = []
    for 
        
    return google_table_template.render(div_id=div_id,
                                        page=page,
                                        page_size=page_size,
                                        data=data,
                                       )
