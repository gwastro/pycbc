""" Plotting utilities and premade plot configurations
"""

def hist_overflow(val, val_max, **kwds):
    """ Make a histogram with an overflow bar above val_max """
    import pylab

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


def add_style_opt_to_parser(parser, default=None):
    """Adds an option to set the matplotlib style to a parser.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        The parser to add the option to.
    default : str, optional
        The default style to use. Default, None, will result in the default
        matplotlib style to be used.
    """
    from matplotlib import pyplot
    parser.add_argument('--mpl-style', default=default,
                        choices=['default']+pyplot.style.available+['xkcd'],
                        help='Set the matplotlib style to use.')


def set_style_from_cli(opts):
    """Uses the mpl-style option to set the style for plots.

    Note: This will change the global rcParams.
    """
    from matplotlib import pyplot
    if opts.mpl_style == 'xkcd':
        # this is treated differently for some reason
        pyplot.xkcd()
    elif opts.mpl_style is not None:
        pyplot.style.use(opts.mpl_style)
