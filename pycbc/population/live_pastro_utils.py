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

    # Only one choice so far, allow for more later
    live_pastro_group.add_argument('--p-astro-method',
                                   choices=['template_param_bins'])
    live_pastro_group.add_argument('--p-astro-spec-data',
        help='File containing information to set up p_astro calculation')

    return live_pastro_group


def verify_live_pastro_options(args):
    """ Input check for live p_astro options

    Parameters
    ----------
    args : argparse Namespace object
        Populated namespace output by parse_args()

    Returns
    -------
    args : argparse Namespace object
    """
    if args.p_astro_method is not None and args.p_astro_spec_data is None:
        raise RuntimeError('Need a p astro data spec file for method %s! ' %
                           args.p_astro_method)
    if args.p_astro_method is None and args.p_astro_spec_data is not None:
        raise RuntimeError('Need to specify p_astro method to do a p_astro '
                           'calculation!')
    return args


_read_spec = {
      'template_param_bins': livepa.read_template_param_bin_data
}

_read_bank = {
      'template_param_bins': livepa.read_template_bank_param
}

_do_calc = {
      'template_param_bins': livepa.template_param_bin_calc
}


class PAstroData():
    """ Class for managing live p_astro calculation persistent info """
    def __init__(self, pastro_method, specfile, bank):
        """
        Read in spec file and extract relevant info from bank

        Parameters
        ----------
        pastro_method: str
            Name of calculation method to be used
        specfile: str
            Path to file containing static data used in calculation
        bank: str
            Path to hdf template bank file
        """
        if pastro_method is None:
            self.do = False
        else:
            self.do = True
            self.method = pastro_method
            logging.info('Setting up p_astro data with method %s', method)
            self.spec = _read_spec[self.method](specfile)
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
