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
    'H1L1': '#caff00',
    'H1V1': '#088e49',
    'L1V1': '#00ffe3',
    'H1L1V1': '#ca9210',
    'H1K1': '#0cd75d',
    'K1L1': '#df9aff',
    'H1K1L1': '#ff866d',
    'K1V1': '#6d61ff',
    'H1K1V1': '#b65910',
    'K1L1V1': '#00c2ce',
    'H1K1L1V1': '#82a200',
    'H1I1': '#20869e',
    'I1L1': '#b2baff',
    'H1I1L1': '#4db68e',
    'I1V1': '#8a7520',
    'H1I1V1': '#92e7ff',
    'I1L1V1': '#928ad2',
    'H1I1L1V1': '#ffe77d',
    'I1K1': '#c2c600',
    'H1I1K1': '#10ff00',
    'I1K1L1': '#087ddb',
    'H1I1K1L1': '#5dffa6',
    'I1K1V1': '#f36100',
    'H1I1K1V1': '#ce5dff',
    'I1K1L1V1': '#00927d',
    'H1I1K1L1V1': '#00ae31',
    'G1H1': '#0000aa',
    'G1L1': '#650400',
    'G1H1L1': '#550075',
    'G1V1': '#004500',
    'G1H1V1': '#615d61',
    'G1L1V1': '#a69a96',
    'G1H1L1V1': '#ebe3eb',
    'G1K1': '#593d00',
    'G1H1K1': '#00495d',
    'G1K1L1': '#4d498e',
    'G1H1K1L1': '#28614d',
    'G1K1V1': '#1035ff',
    'G1H1K1V1': '#efbe9e',
    'G1K1L1V1': '#824539',
    'G1H1K1L1V1': '#a61400',
    'G1I1': '#201c49',
    'G1H1I1': '#aac6be',
    'G1I1L1': '#8200e7',
    'G1H1I1L1': '#496900',
    'G1I1V1': '#718269',
    'G1H1I1V1': '#86a2b6',
    'G1I1L1V1': '#493951',
    'G1H1I1L1V1': '#a2aa6d',
    'G1I1K1': '#391804',
    'G1H1I1K1': '#aa7165',
    'G1I1K1L1': '#657196',
    'G1H1I1K1L1': '#00392d',
    'G1I1K1V1': '#c6b2ce',
    'G1H1I1K1V1': '#c2f7db',
    'G1I1K1L1V1': '#414535',
    'G1H1I1K1L1V1': '#242800',
}

def ifo_color(ifo):
    """ Return a color defining the IFO of interest
    """
    return _ifo_color_map[ifo]


def source_color(source):
    """ Return a standard color to indicate the source type """
    return _source_color_map[source]

def coinc_color(coinc):
    """ Return a color for the coincidence type 
    Parameters:
    coinc: tuple
        A tuple of strings for the IFOs in the coincidence.
        
    Returns:
    string
        The color from the corresponding color map
    """
    if len(coinc) == 1 and coinc in _ifo_color_map:
        return _ifo_color_map[coinc]
    else:
        return _coinc_color_map[coinc]
