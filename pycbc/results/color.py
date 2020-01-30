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

_source_color_map = {
    'BNS': '#A2C8F5',   # light blue
    'NSBH': '#FFB482',  # light orange
    'BBH': '#FE9F9B',   # light red
    'Mass Gap': '#8EE5A1',  # light green
    'GNS': '#98D6CB',   # turquoise
    'GG': '#79BB87',    # green
    'BHG': '#C6C29E'    # dark khaki
}


def ifo_color(ifo):
    return _ifo_color_map[ifo]


def source_color(source):
    return _source_color_map[source]
