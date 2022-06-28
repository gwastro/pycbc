import json
import h5py
import numpy
from lal import YRJUL_SI as lal_s_per_yr
#from pycbc import conversions
from pycbc.tmpltbank import bank_conversions as bankconv


def read_template_param_bin_data(spec_file):
    """
    Parameters

    spec_file: json file with various pieces of data

    Returns

    pa_spec: dictionary giving prerequisite data for p astro calc
    """
    pa_spec = json.load(spec_file)
    # check that the file has the right contents
    assert 'param' in pa_spec
    assert 'bin_edges' in pa_spec  # should be a list of floats
    assert 'sig_per_yr_binned' in pa_spec  # signal rate per bin (per year)
    # do the lengths of bin arrays match?
    assert len(pa_spec['bin_edges']) == len(pa_spec['sig_per_yr_binned']) + 1
    assert 'ref_bns_range' in pa_spec  # float

    return pa_spec


def read_template_bank_param(spec_data, bankf):
    """
    Parameters

    spec_data: dictionary giving prerequisite data for p astro calc

    bankf: HDF5 template bank file

    Returns

    bank_data: dictionary giving bank information binned over specified param
    """
    bank = h5py.File(bankf, 'r')
    # All the templates
    tids = numpy.arange(len(bank['mass1']))
    # Get param vals
    parvals = bankconv.get_bank_property(spec_data['param'], bank, tids)
    counts, edges = numpy.histogram(parvals, bins=spec_data['bin_edges'])
    bank_data = {'bin_edges': edges, 'tcounts': counts, 'num_t': counts.sum()}

    return bank_data


def template_param_bin_calc(padata, trigger_data):
    """
    Parameters

    padata: PAstroData instance storing static information on p astro calculation

    trigger_data: dictionary with trigger properties

    Returns

    p_astro: float
    """
    trig_param = get_param(trigger_data, padata.spec['param'])
    # Which bin is the trigger in?
    bind = numpy.digitize(trig_param, padata.bank['bin_edges'])

    # Get noise rate density and scale by fraction of templates in bin
    # FAR is in Hz, therefore convert to rate per year
    dnoise = noise_density_from_far(trigger_data['far']) * lal_s_per_yr
    dnoise *= padata.bank['tcounts'][bind] / padata.bank['num_t']

    # Get signal rate density at given SNR
    dsig = signal_pdf_from_snr(trigger_data['network_snr'], padata.spec['netsnr_thresh'])\
                                       * padata.spec['sig_per_yr_binned'][bind]
    # Scale by network sensitivity accounting for BNS ranges
    dsig *= signal_rate_from_ranges(trigger_data, padata.spec['ref_bns_range'])

    # Placeholder p astro result
    p_astro = dsig / (dsig + dnoise)
    return p_astro
