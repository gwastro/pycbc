from flask import Flask
from flask import render_template_string

source_html = u'''
<tbody>
{% for x in the_best %}
    <tr>
    </tr>
{% endfor %}
</tbody>
'''

filename = "test.html"

app = Flask('pycbc.results')
app.test_request_context('/home/cbiwer/projects/pycbc_plotting/make_static_html')

context = {'the_best' : [1,2,3,4,5,6]}

with open(filename, "w") as outfh:
    full_html = render_template_string(source_html, context)
    outfh.write(full_html)
