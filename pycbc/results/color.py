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

_coinc_color_map = {
    'H1L1': '#caff00',  # lime green
    'H1V1': '#088e49',  # dark green
    'L1V1': '#00ffe3',  # cyan
    'H1L1V1': '#ca9210',  # brown
    'H1K1': '#0cd75d',  # bright green
    'K1L1': '#df9aff',  # light purple
    'H1K1L1': '#ff866d',  # coral
    'K1V1': '#6d61ff',  # blue violet
    'H1K1V1': '#b65910',  # sienna
    'K1L1V1': '#00c2ce',  # turquoise
    'H1K1L1V1': '#82a200',  # olive
    'H1I1': '#20869e',  # teal
    'I1L1': '#b2baff',  # light blue
    'H1I1L1': '#4db68e',  # sea green
    'I1V1': '#8a7520',  # olive drab
    'H1I1V1': '#92e7ff',  # sky blue
    'I1L1V1': '#928ad2',  # medium purple
    'H1I1L1V1': '#ffe77d',  # light yellow
    'I1K1': '#c2c600',  # chartreuse
    'H1I1K1': '#10ff00',  # neon green
    'I1K1L1': '#087ddb',  # dodger blue
    'H1I1K1L1': '#5dffa6',  # aquamarine
    'I1K1V1': '#f36100',  # orange red
    'H1I1K1V1': '#ce5dff',  # medium orchid
    'I1K1L1V1': '#00927d',  # dark cyan
    'H1I1K1L1V1': '#00ae31',  # forest green
}

def ifo_color(ifo):
    """ Return a color for the IFO
    """
    return _ifo_color_map[ifo]


def source_color(source):
    """ Return a color to indicate the source type """
    return _source_color_map[source]

def coinc_color(coinc):
    """ Return a color for the coincidence type

    Parameters
    coinc : string
        A string for the IFOs in the coincidence.
        This will be in alphabetical order, i.e. H1L1V1.

    Returns
    string : The RGB color for the corresponding coinc
    """
    if len(coinc) == 1:
        if coinc in _ifo_color_map:
            return _ifo_color_map[coinc]
        raise KeyError(f"Unknown IFO coincidence '{coinc}'")
    try:
        return _coinc_color_map[coinc]
    except KeyError:
        raise KeyError(f"Unknown coincidence type '{coinc}'") from None
