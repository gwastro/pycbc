"""
This Module contains generic utility functions for creating plots within
PyCBC.
"""
import os.path, pycbc.version
import configparser as ConfigParser
from html.parser import HTMLParser
from xml.sax.saxutils import escape, unescape

escape_table = {
                '"': "&quot;",
                "'": "&apos;",
                "@": "&#64;",
                }
unescape_table = {
                  "&#64;" : "@",
                 }
for k, v in escape_table.items():
    unescape_table[v] = k

def html_escape(text):
    """ Sanitize text for html parsing """
    return escape(text, escape_table)

class MetaParser(HTMLParser):
    def __init__(self):
        self.metadata = {}
        HTMLParser.__init__(self)

    def handle_data(self, data):
        pass

    def handle_starttag(self, tag, attrs):
        attr= {}
        for key, value in attrs:
            attr[key] = value
        if tag == 'div' and 'class' in attr and attr['class'] == 'pycbc-meta':
            self.metadata[attr['key']] = unescape(attr['value'], unescape_table)


def save_html_with_metadata(fig, filename, fig_kwds, kwds):
    """ Save a html output to file with metadata """
    if isinstance(fig, str):
        text = fig
    else:
        from mpld3 import fig_to_html
        text = fig_to_html(fig, **fig_kwds)

    f = open(filename, 'w')
    for key, value in kwds.items():
        value = escape(value, escape_table)
        line = "<div class=pycbc-meta key=\"%s\" value=\"%s\"></div>" % (str(key), value)
        f.write(line)

    f.write(text)

def load_html_metadata(filename):
    """ Get metadata from html file """
    parser = MetaParser()
    data = open(filename, 'r').read()

    if 'pycbc-meta' in data:
        print("LOADING HTML FILE %s" % filename)
    parser.feed(data)
    cp = ConfigParser.ConfigParser(parser.metadata)
    cp.add_section(os.path.basename(filename))
    return cp

def save_png_with_metadata(fig, filename, fig_kwds, kwds):
    """ Save a matplotlib figure to a png with metadata
    """
    from PIL import Image, PngImagePlugin
    fig.savefig(filename, **fig_kwds)

    im = Image.open(filename)
    meta = PngImagePlugin.PngInfo()

    for key in kwds:
        meta.add_text(str(key), str(kwds[key]))

    im.save(filename, "png", pnginfo=meta)

def save_pdf_with_metadata(fig, filename, fig_kwds, kwds):
    """Save a matplotlib figure to a PDF file with metadata.
    """
    # https://stackoverflow.com/a/17462125
    from matplotlib.backends.backend_pdf import PdfPages

    with PdfPages(filename) as pdfp:
        fig.savefig(pdfp, format='pdf', **fig_kwds)
        metadata = pdfp.infodict()
        for key in kwds:
            if str(key).lower() == 'title':
                # map the title to the official PDF keyword (capitalized)
                metadata['Title'] = str(kwds[key])
            else:
                metadata[str(key)] = str(kwds[key])

def load_png_metadata(filename):
    from PIL import Image
    data = Image.open(filename).info
    cp = ConfigParser.ConfigParser(data)
    cp.add_section(os.path.basename(filename))
    return cp

_metadata_saver = {'.png': save_png_with_metadata,
                   '.html': save_html_with_metadata,
                   '.pdf': save_pdf_with_metadata,
                  }
_metadata_loader = {'.png': load_png_metadata,
                    '.html': load_html_metadata,
                   }

def save_fig_with_metadata(fig, filename, fig_kwds=None, **kwds):
    """ Save plot to file with metadata included. Kewords translate to metadata
    that is stored directly in the plot file. Limited format types available.

    Parameters
    ----------
    fig: matplotlib figure
        The matplotlib figure to save to the file
    filename: str
        Name of file to store the plot.
    """
    if fig_kwds is None:
        fig_kwds = {}
    try:
        extension = os.path.splitext(filename)[1]
        kwds['version'] = pycbc.version.git_verbose_msg
        _metadata_saver[extension](fig, filename, fig_kwds, kwds)
    except KeyError:
        raise TypeError('Cannot save file %s with metadata, extension %s not '
                        'supported at this time' % (filename, extension))

def load_metadata_from_file(filename):
    """ Load the plot related metadata saved in a file

    Parameters
    ----------
    filename: str
        Name of file load metadata from.

    Returns
    -------
    cp: ConfigParser
        A configparser object containing the metadata
    """
    try:
        extension = os.path.splitext(filename)[1]
        return _metadata_loader[extension](filename)
    except KeyError:
        raise TypeError('Cannot read metadata from file %s, extension %s not '
                        'supported at this time' % (filename, extension))
