import logging
import h5py
import numpy
from lal import YRJUL_SI as lal_s_per_yr
from pycbc.tmpltbank import bank_conversions as bankconv
from pycbc.events import triggers
from . import fgmc_functions as fgmcfun


def check_template_param_bin_data(spec_json):
    """
    Parameters
    ----------
    spec_json: JSON dictionary-like object
        Result of parsing json file containing static data

    Returns
    -------
    spec_json: dictionary
    """
    # Check the necessary data are present
    assert 'param' in spec_json
    assert 'bin_edges' in spec_json  # should be a list of floats
    assert 'sig_per_yr_binned' in spec_json  # signal rate per bin (per year)
    # Do the lengths of bin arrays match?
    assert len(spec_json['bin_edges']) == \
           len(spec_json['sig_per_yr_binned']) + 1
    assert 'ref_bns_horizon' in spec_json  # float
    assert 'netsnr_thresh' in spec_json  # float

    return spec_json


def read_template_bank_param(spec_d, bankf):
    """
    Parameters
    ----------
    spec_d: dictionary
        Prerequisite data for p astro calc
    bankf: string
        Path to HDF5 template bank file

    Returns
    -------
    bank_data: dictionary
        Template counts binned over specified param
    """
    with h5py.File(bankf, 'r') as bank:
        # All the templates
        tids = numpy.arange(len(bank['mass1']))
        # Get param vals
        logging.info('Getting %s values from bank', spec_d['param'])
        parvals = bankconv.get_bank_property(spec_d['param'], bank, tids)
    counts, edges = numpy.histogram(parvals, bins=spec_d['bin_edges'])
    bank_data = {'bin_edges': edges, 'tcounts': counts, 'num_t': counts.sum()}
    logging.info('Binned template counts: %s', counts)

    return bank_data


def noise_density_from_far(far, exp_fac):
    """
    Exponential model of noise rate density per time per (reweighted) SNR
    b(rho) ~ k exp(-alpha * rho),
    where alpha is the 'index', yields the relation
    b(rho) = alpha * FAR(rho),
    where FAR is the integral of b(rho) from rho to infinity.
    E.g. fits to single-ifo noise typically yield alpha ~ 6
    """
    return exp_fac * far


def trials_type(ntriggered, nactive):
    """
    The trials factor previously applied to an individual event type FAR
    For single triggers, the factor is the number of active ifos
    For coincs, the factor is either 1 (in double time) or 6 (in triple time)
    6 accounts for both the trials over coinc type and pvalue (non-)followup
    %%% NOTE - ONLY VALID FOR 2- OR 3-IFO SEARCH %%%
    """
    if ntriggered == 1:
        return nactive
    if ntriggered == 2 and nactive == 2:
        return 1
    if ntriggered == 2 and nactive == 3:
        return 6
    # All valid inputs are exhausted, throw an error
    raise ValueError(f"I don't know what to do with {ntriggered} triggered and"
                     f" {nactive} active ifos!")


def signal_pdf_from_snr(netsnr, thresh):
    """ FGMC approximate signal distribution ~ SNR ** -4
    """
    return numpy.exp(fgmcfun.log_rho_fg_analytic(netsnr, thresh))


def signal_rate_rescale(horizons, ref_dhor):
    """
    Compute a factor proportional to the rate of signals with given network SNR
    to account for network sensitivity variation relative to a reference state
    """
    # Combine sensitivities over ifos in a way analogous to network SNR
    net_horizon = sum(hor ** 2. for hor in horizons.values()) ** 0.5
    # signal rate is proportional to horizon distance cubed
    return net_horizon ** 3. / ref_dhor ** 3.


def signal_rate_trig_type(horizons, sens_ifos, trig_ifos):
    """
    Compute a factor accounting for the fraction of signals seen as a given
    trigger type
    """
    # Single-ifo time
    if len(sens_ifos) == 1:
        assert len(trig_ifos) == 1
        return 1.
    # Single trigger in multi-ifo time
    if len(trig_ifos) == 1:
        # Sensitive volume scales with horizon^3
        # Suppress horizon by sqrt(2) wrt coincs
        return (horizons[trig_ifos[0]] / 2**0.5) ** 3. /\
            sum([horizons[i] ** 3. for i in sens_ifos])
    # Double coinc : volume determined by less sensitive ifo
    # Compare to 2nd most sensitive ifo over the observing network
    return sorted([horizons[i] for i in trig_ifos])[0] ** 3. /\
        sorted([horizons[i] for i in sens_ifos])[-2] ** 3.


def template_param_bin_pa(padata, trdata, horizons):
    """
    Parameters
    ----------
    padata: PAstroData instance
        Static information on p astro calculation
    trdata: dictionary
        Trigger properties
    horizons: dictionary
        BNS horizon distances keyed on ifo

    Returns
    -------
    p_astro, p_terr: tuple of floats
    """
    massspin = (trdata['mass1'], trdata['mass2'],
                trdata['spin1z'], trdata['spin2z'])
    trig_param = triggers.get_param(padata.spec['param'], None, *massspin)
    # NB digitize gives '1' for first bin, '2' for second etc.
    bind = numpy.digitize(trig_param, padata.bank['bin_edges']) - 1
    logging.debug('Trigger %s is in bin %i', padata.spec['param'], bind)

    # Get noise rate density
    if 'bg_fac' not in padata.spec:
        expfac = 6.
    else:
        expfac = padata.spec['bg_fac']

    # FAR is in Hz, therefore convert to rate per year (per SNR)
    dnoise = noise_density_from_far(trdata['far'], expfac) * lal_s_per_yr
    logging.debug('FAR %.3g, noise density per yr per SNR %.3g',
                  trdata['far'], dnoise)
    # Scale by fraction of templates in bin
    dnoise *= padata.bank['tcounts'][bind] / padata.bank['num_t']
    logging.debug('Noise density in bin %.3g', dnoise)

    # Get signal rate density per year at given SNR
    dsig = signal_pdf_from_snr(trdata['network_snr'],
                               padata.spec['netsnr_thresh'])
    logging.debug('SNR %.3g, signal pdf %.3g', trdata['network_snr'], dsig)
    dsig *= padata.spec['sig_per_yr_binned'][bind]
    logging.debug('Signal density per yr per SNR in bin %.3g', dsig)
    # Scale by network sensitivity accounting for BNS horizon distances
    dsig *= signal_rate_rescale(horizons, padata.spec['ref_bns_horizon'])
    logging.debug('After horizon rescaling %.3g', dsig)

    p_astro = dsig / (dsig + dnoise)
    logging.debug('p_astro %.4g', p_astro)
    return p_astro, 1 - p_astro


def template_param_bin_types_pa(padata, trdata, horizons):
    """
    Parameters
    ----------
    padata: PAstroData instance
        Static information on p astro calculation
    trdata: dictionary
        Trigger properties
    horizons: dictionary
        BNS horizon distances keyed on ifo

    Returns
    -------
    p_astro, p_terr: tuple of floats
    """
    massspin = (trdata['mass1'], trdata['mass2'],
                trdata['spin1z'], trdata['spin2z'])
    trig_param = triggers.get_param(padata.spec['param'], None, *massspin)
    # NB digitize gives '1' for first bin, '2' for second etc.
    bind = numpy.digitize(trig_param, padata.bank['bin_edges']) - 1
    logging.debug('Trigger %s is in bin %i', padata.spec['param'], bind)

    # Get noise rate density
    if 'bg_fac' not in padata.spec:
        expfac = 6.
    else:
        expfac = padata.spec['bg_fac']

    # Ifos over trigger threshold
    tr_ifos = trdata['triggered']

    # FAR is in Hz, therefore convert to rate per year (per SNR)
    dnoise = noise_density_from_far(trdata['far'], expfac) * lal_s_per_yr
    logging.debug('FAR %.3g, noise density per yr per SNR %.3g',
                  trdata['far'], dnoise)
    # Scale by fraction of templates in bin
    dnoise *= padata.bank['tcounts'][bind] / padata.bank['num_t']
    logging.debug('Noise density in bin %.3g', dnoise)
    # Back out trials factor to give noise density for triggered event type
    dnoise /= float(trials_type(len(tr_ifos), len(trdata['sensitive'])))
    logging.debug('Divide by previously applied trials factor: %.3g', dnoise)

    # Get signal rate density per year at given SNR
    dsig = signal_pdf_from_snr(trdata['network_snr'],
                               padata.spec['netsnr_thresh'])
    logging.debug('SNR %.3g, signal pdf %.3g', trdata['network_snr'], dsig)
    dsig *= padata.spec['sig_per_yr_binned'][bind]
    logging.debug('Total signal density per yr per SNR in bin %.3g', dsig)
    # Scale by network sensitivity accounting for BNS horizons
    dsig *= signal_rate_rescale(horizons, padata.spec['ref_bns_horizon'])
    logging.debug('After network horizon rescaling %.3g', dsig)
    # Scale by relative signal rate in triggered ifos
    dsig *= signal_rate_trig_type(horizons, trdata['sensitive'], tr_ifos)
    logging.debug('After triggered ifo rate rescaling %.3g', dsig)

    p_astro = dsig / (dsig + dnoise)
    logging.debug('p_astro %.4g', p_astro)
    return p_astro, 1 - p_astro


__all__ = [
    "check_template_param_bin_data",
    "read_template_bank_param",
    "noise_density_from_far",
    "signal_pdf_from_snr",
    "signal_rate_rescale",
    "signal_rate_trig_type",
    "template_param_bin_pa",
    "template_param_bin_types_pa",
]
