""" Utilities for managing matplotlib colors and mapping ifos to color
"""

_ifo_color_map = {'H1':'red', 'L1':'green', 'V1':'magenta'}

def ifo_color(ifo):
    return _ifo_color_map[ifo]
