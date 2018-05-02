""" Utilities for managing matplotlib colors and mapping ifos to color
"""

_ifo_color_map = {
    'G1': '#222222',  # dark gray
    'K1': '#ffb200',  # yellow/orange
    'H1': '#ee0000',  # red
    'I1': '#b0dd8b',  # light green
    'L1': '#4ba6ff',  # blue
    'V1': '#9b59b6',  # magenta/purple
}


def ifo_color(ifo):
    return _ifo_color_map[ifo]
