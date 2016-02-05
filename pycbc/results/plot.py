""" Plotting utilities and premade plot configurations
"""

def hist_overflow(val, val_max, **kwds):
    """ Make a histogram with an overflow bar above val_max """
    import pylab, numpy

    overflow = len(val[val>=val_max])
    pylab.hist(val[val<val_max], **kwds)
        
    if 'color' in kwds:
        color = kwds['color']
    else:
        color = None    
 
    if overflow > 0:
        rect = pylab.bar(val_max+0.05, overflow, .5, color=color)[0]
        pylab.text(rect.get_x(), 
                   1.10*rect.get_height(), '%s+' % val_max)
    
