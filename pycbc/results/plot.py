"""
This Module contains generic utility functions for creating plots within 
PyCBC. 
"""
import os.path
import ConfigParser

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

def load_png_metadata(filename):
    from PIL import Image, PngImagePlugin
    data = Image.open(filename).info
    cp = ConfigParser.ConfigParser(data)
    cp.add_section(os.path.basename(filename))
    return cp
    
_metadata_saver = {'.png':save_png_with_metadata,
                  }
_metadata_loader = {'.png':load_png_metadata,
                   }

def save_fig_with_metadata(fig, filename, fig_kwds={}, **kwds):
    """ Save plot to file with metadata included. Kewords translate to metadata
    that is stored directly in the plot file. Limited format types available.
    
    Parameters
    ----------
    fig: matplotlib figure
        The matplotlib figure to save to the file
    filename: str
        Name of file to store the plot.
    """
    try:
        extension = os.path.splitext(filename)[1]
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
