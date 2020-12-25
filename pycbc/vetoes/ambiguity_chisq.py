# Copyright (C) 2020 Bhooshan Gadre
## Ambiguity chisq based on sec V of https://arxiv.org/abs/1708.03605

import numpy as np
import itertools
import copy
import logging

import pycbc.pnutils as pnu
import pycbc.noise as pn
import pycbc.filter as pf
import pycbc.psd as pp
import pycbc.noise as pn
import pycbc.waveform as pw
import pycbc.types as pt
from pycbc.conversions import chi_eff


def segment_snrs(filters, stilde, psd, flow, fhigh):
    """ This functions calculates the snr of each bank veto template against
    the segment
    Parameters
    ----------
    filters: list of FrequencySeries
        The list of bank veto templates filters.
    stilde: FrequencySeries
        The current segment of data.
    psd: FrequencySeries
    low_frequency_cutoff: float
    Returns
    -------
    snr (list): List of snr time series.
    """
    snrs = []

    for template in filters:
        snr = pf.matched_filter(template, stilde, None, flow, fhigh, pf.sigmasq(template, psd, flow, fhigh))
        snrs.append(snr)
    return snrs

def simple_inner(htilde, stilde, psd, flow, fhigh):
    kmin, kmax = pf.get_cutoff_indices(flow, fhigh, psd.delta_f, 2*len(psd)+1)
    indices = slice(kmin, kmax)
    norm = pf.sigma(htilde, psd, flow, fhigh) * pf.sigma(stilde, psd, flow, fhigh)
    return (np.conjugate(htilde.data[indices])*stilde.data[indices]/psd.data[indices]).sum()*4.0*psd.delta_f / norm

def get_cov_gg(filters, psd, flow, fhigh):
    len_filters = len(filters)
    cov = np.zeros((len_filters, len_filters))
    for i, j in itertools.combinations_with_replacement(range(len_filters), 2):
        cov[i, j] = simple_inner(filters[i], filters[j], psd, flow, fhigh).real
        cov[j,i] = simple_inner(filters[j], filters[i], psd, flow, fhigh).real
        # cov[j,i] = cov[i, j]
    return cov

def get_cov_gh(filters, htilde, psd, flow, fhigh):
    return [simple_inner(htilde, stilde, psd, flow, fhigh) for stilde in filters]

def get_cov_matrix(snr, cov_gg, cov_gh):
    cov = np.zeros_like(cov_gg)
    for i, j in itertools.combinations_with_replacement(range(len(cov_gh)), 2):
        cov[i, j] = cov_gg[i, j] - (snr.real*cov_gh[i].real + snr.imag*cov_gh[i].imag)*(snr.real*cov_gh[j].real + snr.imag*cov_gh[j].imag)
        cov[j,i] = cov[i, j]
    return cov

def get_eval_rot_mat(cov):
    eig, evec = np.linalg.eig(cov)
    sort = eig.argsort()
    return eig[sort], evec[sort, :]

def get_vector(snr, snr_id, seg_snrs, cov_gh):
    vec = np.zeros(len(cov_gh))
    for i in range(len(cov_gh)):
        vec[i] = seg_snrs[i][snr_id].real - snr.real*cov_gh[i].real - snr.imag*cov_gh[i].imag
    return vec

def get_chisq(vec, eig, rot_mat, threshold=0.05):
    ## Variaous conditions to impose sane chisq and DoF
    rel_eig = eig/eig[-1] > threshold
    dof = rel_eig.sum()
    rot_vec = np.dot(vec, rot_mat)
    chi =(rot_vec[rel_eig]**2.0/eig[rel_eig]).sum()

    logging.info('Eigenvalues min, max, ratio: {}, {}, {}'.format(eig[rel_eig][0], eig[rel_eig][-1], eig[rel_eig][-1]/eig[rel_eig][0]))

    return chi, dof

def compute_chisq(snrs, snr_ids, seg_snrs, cov_gg, cov_gh, threshold=1e-3):
    chis, dofs = [], []
    for snr, snr_id in zip(snrs, snr_ids):
        cov = get_cov_matrix(snr, cov_gg, cov_gh)
        eig, rot_mat = get_eval_rot_mat(cov)
        vec = get_vector(snr, snr_id, seg_snrs, cov_gh)
        chi, dof = get_chisq(vec, eig, rot_mat, threshold)
        chis.append(chi)
        dofs.append(dof)

    return np.array(chis), np.array(dofs)


_snr = None
def inner(vec1, vec2, psd=None,
          low_frequency_cutoff=None, high_frequency_cutoff=None,
          v1_norm=None, v2_norm=None):

    htilde = pf.make_frequency_series(vec1)
    stilde = pf.make_frequency_series(vec2)

    N = (len(htilde)-1) * 2

    global _snr
    _snr = None
    if _snr is None or _snr.dtype != htilde.dtype or len(_snr) != N:
        _snr = pt.zeros(N,dtype=pt.complex_same_precision_as(vec1))
        snr, corr, snr_norm = pf.matched_filter_core(htilde,stilde,psd,low_frequency_cutoff,
                                                  high_frequency_cutoff, v1_norm, out=_snr)
    if v2_norm is None:
        v2_norm = pf.sigmasq(stilde, psd, low_frequency_cutoff, high_frequency_cutoff)

    snr.data = snr.data * snr_norm / np.sqrt(v2_norm)

    return snr

class SingleDetAmbiguityChisq(object):
    returns = {'ambiguity_chisq': np.float32, 'ambiguity_chisq_dof' : np.int}
    def __init__(self, status, bank, snr_threshold, flow, fmax,
            mc_bound, eta_bound, chie_bound, min_filters, max_filters,
            time_indices=[0], condition_threshold=1.0e-1, per_template=False, use_full_bank=True):
        if status:
            self.do = True

            self.column_name = "ambiguity_chisq"
            self.table_dof_name = "ambiguity_chisq_dof"
            self.snr_threshold = snr_threshold
            self.condition_threshold = condition_threshold
            self.min_filters = min_filters
            self.max_filters = max_filters
            self.flow = flow
            self.fmax = fmax

            if ('mchirp' not in bank.table.fieldnames) or  ('eta' not in bank.table.fieldnames):
                mchirp, eta = pnu.mass1_mass2_to_mchirp_eta(bank.table['mass1'], bank.table['mass2'])
                if 'mchirp' not in bank.table.fieldnames:
                    bank.table.add_fields(mchirp, 'mchirp')
                if 'eta' not in bank.table.fieldnames:
                    bank.table.add_fields(eta, 'eta')
            if 'chi_eff' not in bank.table.fieldnames:
                chie = chi_eff(bank.table['mass1'], bank.table['mass2'], bank.table['spin1z'], bank.table['spin2z'])
                bank.table.add_fields(chie, 'chi_eff')
            self.bank = bank
            self.use_full_bank = use_full_bank
            if use_full_bank:
                self.get_filters = self.get_full_bank
            else:
                if per_template:
                    # self.choose_filters = filters_for_template(self.bank, min_filters, max_filters)
                    self.choose_filters = FiltersForTemplate(self.bank, mc_bound=mc_bound, eta_bound=eta_bound, chie_bound=chie_bound ,min_filters=min_filters, max_filters=max_filters)
                    self.get_filters = self.choose_filters.get_relevant_filters
                else:
                    self.choose_filters = FiltersRegions(self.bank, mc_bound=mc_bound, eta_bound=eta_bound, chie_bound=chie_bound ,min_filters=min_filters, max_filters=max_filters)
                    self.get_filters = self.choose_filters.get_relevant_filters

            self.time_idices = np.array(time_indices) # Currently only trigger time = 0 indices are supported
            self._cache_seg_snrs = {} # seg_snrs
            self._cache_filter_matches = {} # (g, g') matches between filters # Can be time series in future due to 'time_indices'
            self._cache_template_filter_matches = {} # (h, g) matches between template and filter # Can be time series in future 'time_indices'

        else:
            self.do = False

    @staticmethod
    def insert_option_group(parser):
        group = parser.add_argument_group("ambiguity chi-square")
        group.add_argument("--ambi-status", action="store_true", default=True)
        group.add_argument("--ambi-snr-threshold", type=float,
            help="Minimum SNR threshold to use ambi chisq")
        group.add_argument("--ambi-condition-threshold", type=float, default=0.05,
            help="Minimum SNR threshold to use ambi chisq")
        group.add_argument("--ambi-veto-bank-file", type=str, help="bank file for ambiguity chisq")
        group.add_argument("--ambi-veto-use-full-bank", action="store_true", default=True)
        group.add_argument("--ambi-min-filters", type=int, default=10,
                help="maximum filters to be used for ambiguity chisq from the veto bank")
        group.add_argument("--ambi-max-filters", type=int, default=30,
                help="maximum filters to be used for ambiguity chisq from the veto bank")
        group.add_argument('--ambi-mc-bins', nargs='*', type=float, help="chirp mass bin boundaries")
        group.add_argument('--ambi-eta-bins', nargs='*', type=float, help="symmetric mass ratio bin boundaries: eta range [0.25, 0)")
        group.add_argument('--ambi-chi-eff-bins', nargs='*', type=float, help="chi eff bin boundaries: chi eff range (-1, 1)")
        group.add_argument("--ambi-per-template", action="store_true", default=False)

    @classmethod
    def from_cli(cls, args, bank): # , match_class):
        flow = args.low_frequency_cutoff
        fmax = args.sample_rate/2.
        mc_bound = args.ambi_mc_bins if args.ambi_mc_bins else []
        eta_bound = args.ambi_eta_bins if args.ambi_eta_bins else [1e-4, 0.25]
        chie_bound = args.ambi_chi_eff_bins if args.ambi_chi_eff_bins else [-1, 1]
        mc_bound.sort()
        eta_bound.sort()
        chie_bound.sort()
        return cls(args.ambi_status, bank, args.ambi_snr_threshold, flow, fmax,
                mc_bound, eta_bound, chie_bound, args.ambi_min_filters, args.ambi_max_filters,
                per_template=args.ambi_per_template, use_full_bank=args.ambi_veto_use_full_bank,
                condition_threshold=args.ambi_condition_threshold)

    def get_full_bank(self, htilde):
        mchirp, eta = pnu.mass1_mass2_to_mchirp_eta(htilde.params.mass1, htilde.params.mass2)
        chie = chi_eff(htilde.params.mass1, htilde.params.mass2, htilde.params.spin1z, htilde.params.spin2z)
        h_params = {'mchirp': mchirp, 'eta': eta, 'chi_eff': chie}
        # Just to make sure we are not having the trigger template or a very
        # close template
        idx = np.abs(self.bank.table['mchirp'] - h_params['mchirp'])
        idx += np.abs(self.bank.table['eta'] - h_params['eta'])
        idx += np.abs(self.bank.table['chi_eff'] - h_params['chi_eff'])
        idx = (idx > 1.0e-3)
        filter_list = copy.copy(self.bank)
        filter_list.table = filter_list.table[idx]
        return list(filter_list)

    def cache_seg_snrs(self, htilde, stilde, psd):
        key = (id(htilde), stilde.end_time, id(psd))
        if key not in self._cache_seg_snrs:
            filters = self.get_filters(htilde)
            self._cache_seg_snrs[key] = segment_snrs(filters, stilde, psd, self.flow, self.fmax)  ## Assumed data is overwhitened
        return self._cache_seg_snrs[key]

    ## This can be made better by actually caching match between commom filters
    ### But we are using 'id', filters may not live long enough if memory problem. Need to check what happens
    def cache_filter_matches(self, htilde, psd):
        key = (id(htilde), id(psd))
        if key not in self._cache_seg_snrs:
            filters = self.get_filters(htilde)
            self._cache_filter_matches[key] = get_cov_gg(filters, psd, self.flow, self.fmax)
        return self._cache_filter_matches[key]

    def cache_template_filter_matches(self, htilde, psd):
        # key = (id(htilde), id(psd))
        key = (id(htilde), id(psd))
        if key not in self._cache_template_filter_matches:
            filters = self.get_filters(htilde)
            self._cache_template_filter_matches[key] = get_cov_gh(filters, htilde, psd, self.flow, self.fmax)
        return self._cache_template_filter_matches[key]

    def values(self, snrs, snr_ids, htilde, stilde, psd, snr_threshold=None, condition_threshold=None):
        if self.do:

            filters = self.get_filters(htilde)
            if len(filters) > self.min_filters:
                logging.info('Ambiguity chi-squares are on! ...')

                thre = snr_threshold if snr_threshold else self.snr_threshold
                logging.info('Gathering relevant triggers')
                if thre:
                    rel_ids = np.abs(snrs) > thre
                else:
                    rel_ids = np.repeat(True, len(snrs))

                seg_snrs = self.cache_seg_snrs(htilde, stilde, psd)
                cov_gg = self.cache_filter_matches(htilde, psd)
                cov_gh = self.cache_template_filter_matches(htilde, psd)

                cond_thre = condition_threshold if condition_threshold else self.condition_threshold
                logging.info('Computing ambiguity chi-squares ...')
                chisq, dof = compute_chisq(snrs[rel_ids], snr_ids[rel_ids], seg_snrs, cov_gg, cov_gh, cond_thre)
                logging.info(chisq/dof, dof)
                return chisq/dof, dof
            else:
                return None, None
        else:
            return None, None


class FiltersRegions(object):
    def __init__(self, bank, mc_bound=[], eta_bound=[1e-4, 0.25], chie_bound=[-1, 1], min_filters=10, max_filters=20):

        logging.info("Using for bins")
        if ('mchirp' not in bank.table.fieldnames) or  ('eta' not in bank.table.fieldnames):
            mchirp, eta = pnu.mass1_mass2_to_mchirp_eta(bank.table['mass1'], bank.table['mass2'])
            if 'mchirp' not in bank.table.fieldnames:
                bank.table.add_fields(mchirp, 'mchirp')
            if 'eta' not in bank.table.fieldnames:
                bank.table.add_fields(eta, 'eta')
        if 'chi_eff' not in bank.table.fieldnames:
            chie = chi_eff(bank.table['mass1'], bank.table['mass2'], bank.table['spin1z'], bank.table['spin2z'])
            bank.table.add_fields(chie, 'chi_eff')

        self.bank = bank
        self.min_filters = min_filters
        self.max_filters = max_filters
        self.ranges = {'mchirp': mc_bound if mc_bound else np.array([bank.table['mchirp'].min()*0.95,bank.table['mchirp'].max()*1.05]),
                'eta': np.array(eta_bound), 'chi_eff': np.array(chie_bound)} ## make sure to have min and max of the mf banks covered
        self._cache_relevant_filters = {}

    def idx_from_ranges(self,  htilde):
            mchirp, eta = pnu.mass1_mass2_to_mchirp_eta(htilde.params.mass1, htilde.params.mass2)
            chie = chi_eff(htilde.params.mass1, htilde.params.mass2, htilde.params.spin1z, htilde.params.spin2z)
            h_params = {'mchirp': mchirp, 'eta': eta, 'chi_eff': chie}

            idx = np.abs(self.bank.table['mchirp'] - h_params['mchirp']) > 1e-5 # Just to make sure we are not having the trigger template
            logging.info("not me left over are {}".format(idx.sum()))
            for k, v in self.ranges.items():
                bin = np.searchsorted(v, h_params[k])
                if bin < 1 or bin >= len(v): # If out of boundaries given, we are not calculating anything
                    idx *= np.zeros_like(idx)
                    logging.info("after {}, out of the bound: left over are {}".format(k, idx.sum()))
                    logging.info("bins are {} and template is {} ".format(v, h_params[k]))
                    break
                else:
                    idx = idx * (self.bank.table[k] > v[bin-1]) * (self.bank.table[k] < v[bin])
                    logging.info("after {} left over are {}".format(k, idx.sum()))
                    logging.info("bins are {} and template is {} ".format(v, h_params[k]))

            return idx

    def get_relevant_filters(self, htilde):
        logging.info("Getting relevant templates")
        key = id(htilde)
        if key not in self._cache_relevant_filters:
            idx = self.idx_from_ranges(htilde)
            filter_list = copy.copy(self.bank)
            filter_list.table = filter_list.table[idx]
            if idx.sum() > self.max_filters:
                mchirp, eta = pnu.mass1_mass2_to_mchirp_eta(htilde.params.mass1, htilde.params.mass2)
                chie = chi_eff(htilde.params.mass1, htilde.params.mass2, htilde.params.spin1z, htilde.params.spin2z)
                h_params = {'mchirp': mchirp, 'eta': eta, 'chi_eff': chie}
                filter_list.table = filter_list.table[np.argsort(np.abs(filter_list.table['mchirp'] - h_params['mchirp']))[:self.max_filters]]
            self._cache_relevant_filters[key] = list(filter_list) # list of filters to be used in chisq
        logging.info("Getting relevant templates... Done!")
        return self._cache_relevant_filters[key]

    def clear_filters(self, htilde):
        key = id(htilde)
        if key in self._cache_relevant_filters:
            del self._cache_relevant_filters[key]


class FiltersForTemplate(object):
    def __init__(self, bank, mc_bound=[], eta_bound=[1e-4, 0.25], chie_bound=[-1, 1], min_filters=10, max_filters=20):

        logging.info("Using for template")
        if ('mchirp' not in bank.table.fieldnames) or  ('eta' not in bank.table.fieldnames):
            mchirp, eta = pnu.mass1_mass2_to_mchirp_eta(bank.table['mass1'], bank.table['mass2'])
            if 'mchirp' not in bank.table.fieldnames:
                bank.table.add_fields(mchirp, 'mchirp')
            if 'eta' not in bank.table.fieldnames:
                bank.table.add_fields(eta, 'eta')
        if 'chi_eff' not in bank.table.fieldnames:
            chie = chi_eff(bank.table['mass1'], bank.table['mass2'], bank.table['spin1z'], bank.table['spin2z'])
            bank.table.add_fields(chie, 'chi_eff')

        self.bank = bank
        self.min_filters = min_filters
        self.max_filters = max_filters
        self.ranges = {'mchirp': mc_bound if mc_bound else np.array([bank.table['mchirp'].min()*0.95,bank.table['mchirp'].max()*1.05]),
                'eta': np.array(eta_bound), 'chi_eff': np.array(chie_bound)} ## make sure to have min and max of the mf banks covered
        self._cache_relevant_filters = {}

    def idx_from_ranges(self,  htilde):
            mchirp, eta = pnu.mass1_mass2_to_mchirp_eta(htilde.params.mass1, htilde.params.mass2)
            chie = chi_eff(htilde.params.mass1, htilde.params.mass2, htilde.params.spin1z, htilde.params.spin2z)
            h_params = {'mchirp': mchirp, 'eta': eta, 'chi_eff': chie}

            idx = np.abs(self.bank.table['mchirp'] - h_params['mchirp']) > 1e-5 # Just to make sure we are not having the trigger template
            for k, v in self.ranges.items():
                bin = np.searchsorted(v, h_params[k])
                if bin < 1 or bin >= len(v): # If out of boundaries given, we are not calculating anything
                    idx *= np.zeros_like(idx)
                else:
                    idx = idx * (self.bank.table[k] > v[bin-1]) * (self.bank.table[k] < v[bin])

            return idx

    def get_relevant_filters(self, htilde):
        logging.info("Getting relevant templates")
        key = id(htilde)
        if key not in self._cache_relevant_filters:
            idx = self.idx_from_ranges(htilde)
            filter_list = copy.copy(self.bank)
            filter_list.table = filter_list.table[idx]
            if idx.sum() > self.max_filters:
                mchirp, eta = pnu.mass1_mass2_to_mchirp_eta(htilde.params.mass1, htilde.params.mass2)
                chie = chi_eff(htilde.params.mass1, htilde.params.mass2, htilde.params.spin1z, htilde.params.spin2z)
                h_params = {'mchirp': mchirp, 'eta': eta, 'chi_eff': chie}
                filter_list.table = filter_list.table[np.argsort(np.abs(filter_list.table['mchirp'] - h_params['mchirp']))[:self.max_filters]]
            self._cache_relevant_filters[key] = list(filter_list) # list of filters to be used in chisq
        logging.info("Getting relevant templates... Done!")
        return self._cache_relevant_filters[key]

    def clear_filters(self, htilde):
        key = id(htilde)
        if key in self._cache_relevant_filters:
            del self._cache_relevant_filters[key]

