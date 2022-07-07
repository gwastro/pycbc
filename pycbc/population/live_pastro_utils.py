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
    live_pastro_group.add_argument('--p-astro-method', choices=\
          ['template_param_bins'])
    live_pastro_group.add_argument('--p-astro-spec-data', help='File '
          'containing information to set up p_astro calculation')

    return live_pastro_group


def verify_live_pastro_options(args):
    """ Input check for live p astro options

    Parameters
    ----------
    args : argparse Namespace object
        Populated namespace output by parse_args()

    Returns
    -------
    args : argparse Namespace object
    """

    # Convenience attr do_p_astro
    if args.p_astro_method is None:
        args.do_p_astro = False
    else:
        args.do_p_astro = True

    if args.do_p_astro and args.p_astro_spec_data is None:
        raise RuntimeError('Need a p astro data spec file for method %s! ' %
                           args.p_astro_method)
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
    def __init__(self, args, bank):
        """ Describe """
        if args.do_p_astro:
            self.do = True
            self.method = args.p_astro_method
            self.spec = _read_spec[self.method](args.p_astro_spec_data)
            self.bank = _read_bank[self.method](self.spec, bank)
        else:
            self.do = False

    def do_pastro_calc(self, trigger_data, horizons):
        """ Describe """
        if not self.do:
            return None

        p_astro = _do_calc[self.method](self, trigger_data, horizons)
        return p_astro
