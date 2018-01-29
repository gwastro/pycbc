""" Utilities for managing matplotlib colors and mapping ifos to color
"""

# try and import color scheme from GWpy, see
# https://gwpy.github.io/docs/stable/plotter/colors.html for details
try:
    from gwpy.plotter.colors import GW_OBSERVATORY_COLORS as _ifo_color_map
except ImportError:
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
