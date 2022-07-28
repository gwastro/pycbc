import logging
from . import live_pastro as livepa


def insert_live_pastro_option_group(parser):
    """ Add low-latency p astro options to the argparser object.

    Parameters
    ----------
    parser : object
        ArgumentParser instance.

    Returns
    -------
    live_pastro_group :
        Argument group object
    """

    live_pastro_group = parser.add_argument_group('Options for live p_astro '
                                                  'calculation')
    live_pastro_group.add_argument('--p-astro-spec',
        help='File containing information to set up p_astro calculation')

    return live_pastro_group


# Only one choice so far, allow for more later
_check_spec = {
      'template_param_bins': livepa.check_template_param_bin_data
}

_read_bank = {
      'template_param_bins': livepa.read_template_bank_param
}

_do_calc = {
      'template_param_bins': livepa.template_param_bin_calc
}


class PAstroData():
    """ Class for managing live p_astro calculation persistent info """
    def __init__(self, specfile, bank):
        """
        Read in spec file and extract relevant info from bank

        Parameters
        ----------
        specfile: str
            Path to file giving method and static data used in calculation
        bank: str
            Path to hdf template bank file
        """
        if specfile is None:
            self.do = False
        else:
            self.do = True
              
            with open(specfile) as specf:
                self.spec_json = json.load(specf)
            try:
                self.method = self.spec_json['method']
            except:
                raise ValueError("Can't find 'method' in p_astro spec file!")
            logging.info('Setting up p_astro data with method %s', self.method)
            self.spec = _check_spec[self.method](self.spec_json)
            self.bank = _read_bank[self.method](self.spec, bank)

    def do_pastro_calc(self, trigger_data, horizons):
        """ No-op, or call the despatch dictionary to evaluate p_astro """
        if not self.do:
            return None, None

        logging.info('Computing p_astro')
        p_astro, p_terr = _do_calc[self.method](self, trigger_data, horizons)
        return p_astro, p_terr


__all__ = [
    "insert_live_pastro_option_group",
    "verify_live_pastro_options",
    "PAstroData"
]
