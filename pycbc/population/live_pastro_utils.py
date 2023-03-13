import logging
import json
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

    live_pastro_group = parser.add_argument_group('Options for live p_astro')
    live_pastro_group.add_argument('--p-astro-spec',
                                   help='File containing information to set '
                                   'up p_astro calculation')

    return live_pastro_group


# Choices of p astro calc method
_check_spec = {
    'template_param_bins': livepa.check_template_param_bin_data,
    'template_param_bins_types': livepa.check_template_param_bin_data,
    'template_param_bins_types_farlim':
        livepa.check_template_param_bin_farlim_data
}

_read_bank = {
    'template_param_bins': livepa.read_template_bank_param,
    'template_param_bins_types': livepa.read_template_bank_param,
    'template_param_bins_types_farlim': livepa.read_template_bank_param
}

_do_calc = {
    'template_param_bins': livepa.template_param_bin_pa,
    'template_param_bins_types': livepa.template_param_bin_types_pa,
    'template_param_bins_types_farlim':
        livepa.template_param_bin_types_farlim_pa
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
            except KeyError as ke:
                raise ValueError("Can't find 'method' in p_astro spec file!") \
                    from ke
            logging.info('Setting up p_astro data with method %s', self.method)
            self.spec = _check_spec[self.method](self.spec_json)
            self.bank = _read_bank[self.method](self.spec, bank)

    def apply_significance_limits(self, trigger_data):
        """
        If the network SNR and FAR indicate saturation of the FAR estimate,
        set them to the fixed values given in the specification.
        """
        # This only happens for double or triple events
        if len(trigger_data['triggered']) == 1:
            return trigger_data

        if len(trigger_data['triggered']) > 1:
            farlim = self.spec['limit_far']
            snrlim = self.spec['limit_snr']
            # Only do anything if FAR and SNR are beyond given limits
            if trigger_data['far'] > farlim or \
                    trigger_data['network_snr'] < snrlim:
                return trigger_data

            logging.debug('Truncating FAR and SNR from %f, %f to %f, %f',
                          trigger_data['far'], trigger_data['network_snr'],
                          farlim, snrlim)
            trigger_data['network_snr'] = snrlim
            trigger_data['far'] = farlim
            return trigger_data

        raise RuntimeError('Number of triggered ifos must be >0 !')

    def do_pastro_calc(self, trigger_data, horizons):
        """ No-op, or call the despatch dictionary to evaluate p_astro """
        if not self.do:
            return None, None

        logging.info('Computing p_astro')
        p_astro, p_terr = _do_calc[self.method](self, trigger_data, horizons)
        return p_astro, p_terr


__all__ = [
    "insert_live_pastro_option_group",
    "PAstroData"
]
