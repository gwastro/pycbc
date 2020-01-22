""" utilities for assigning FAR to single detector triggers
"""
import h5py
import numpy as np
from pycbc.events import ranking, trigger_fits as fits
from pycbc.types import MultiDetOptionAction
from pycbc import conversions as conv
from pycbc import bin_utils

class LiveSingle(object):
    def __init__(self, ifo,
                 newsnr_threshold=10.0,
                 reduced_chisq_threshold=5,
                 duration_threshold=0,
                 fit_file=None,
                 sngl_ifar_est_dist=None,
                 fixed_ifar=None):
        self.ifo = ifo
        if sngl_ifar_est_dist in ["fixed", None]:
            fit_bins = None
            fit_rates = None
            fit_thresh = None
            fit_coeffs = None
        else:
            with h5py.File(fit_file, 'r') as fit_file:
                bin_edges = fit_file['bins_edges'][:]
                fit_bins = bin_utils.IrregularBins(bin_edges)
                dist_grp = fit_file[ifo + '/' + sngl_ifar_est_dist]
                live_time = fit_file[ifo].attrs['live_time']
                fit_rates = dist_grp['counts'][:] / live_time
                fit_coeffs = dist_grp['fit_coeff'][:]
                fit_thresh = fit_file.attrs['fit_threshold']
        self.thresholds = {
            "newsnr": newsnr_threshold,
            "reduced_chisq": reduced_chisq_threshold,
            "duration": duration_threshold}
        self.fit_info = {
            "fixed_ifar": fixed_ifar,
            "bins": fit_bins,
            "rates": fit_rates,
            "coeffs": fit_coeffs,
            "thresh": fit_thresh}

    @staticmethod
    def insert_args(parser):
        parser.add_argument('--single-newsnr-threshold', nargs='+',
                            type=float, action=MultiDetOptionAction,
                            help='Newsnr min threshold for single triggers. '
                                 'Can be given as a single value or as '
                                 'detector-value pairs, e.g. H1:6 L1:7 V1:6.5')
        parser.add_argument('--single-reduced-chisq-threshold', nargs='+',
                            type=float, action=MultiDetOptionAction,
                            help='Maximum reduced chi-squared threshold for '
                                 'single triggers. Can be given as a single '
                                 'value or as detector-value pairs, e.g. '
                                 'H1:2 L1:2 V1:3')
        parser.add_argument('--single-duration-threshold', nargs='+',
                            type=float, action=MultiDetOptionAction,
                            help='Minimum duration threshold for single '
                                 'triggers. Can be given as a single value '
                                 'or as detector-value pairs, e.g. H1:6 L1:6 '
                                 'V1:8')
        parser.add_argument('--single-fixed-ifar', nargs='+',
                            type=float, action=MultiDetOptionAction,
                            help='A fixed value for IFAR, still uses cuts '
                                 'defined by command line. Can be given as '
                                 'a single value or as detector-value pairs, '
                                 'e.g. H1:0.001 L1:0.001 V1:0.0005')
        parser.add_argument('--single-fit-file',
                            help='File which contains definitons of fit '
                                 'coefficients and counts for specific '
                                 'single trigger IFAR fitting.')
        parser.add_argument('--sngl-ifar-est-dist', nargs='+',
                            default='conservative',
                            choices=['conservative', 'mean'],
                            action=MultiDetOptionAction,
                            help='Which trigger distribution to use when '
                                 'calculating IFAR of single triggers. '
                                 'Default conservative. Can be given as '
                                 'a single value or as detector-value pairs, '
                                 'e.g. H1:mean L1:mean V1:conservative')

    @classmethod
    def from_cli(cls, args, ifo):
        sngl_opts_required = all([args.single_fit_file,
                                  args.single_reduced_chisq_threshold,
                                  args.single_duration_threshold,
                                  args.single_newsnr_threshold])
        if args.enable_single_detector_background and not sngl_opts_required:
            raise RuntimeError("Single detector trigger options "
                               "(--single-fit-file, "
                               "--single-reduced-chisq-threshold, "
                               "--single-duration-threshold, "
                               "--single-newsnr-threshold) "
                               "must ALL be given if single detector "
                               "background is enabled")
        return cls(
           ifo, newsnr_threshold=args.single_newsnr_threshold[ifo],
           reduced_chisq_threshold=args.single_reduced_chisq_threshold[ifo],
           duration_threshold=args.single_duration_threshold[ifo],
           fixed_ifar=args.single_fixed_ifar,
           fit_file=args.single_fit_file,
           sngl_ifar_est_dist=args.sngl_ifar_est_dist[ifo]
           )

    def check(self, trigs, data_reader):
        """ Look for a single detector trigger that passes the thresholds in
        the current data.
        """

        if len(trigs['snr']) == 0:
            return None

        # Apply cuts to trigs before clustering
        # Cut on snr so that triggers which could not reach newsnr
        # threshold do not have newsnr calculated
        valid_idx = (trigs['template_duration'] >
                     self.thresholds['duration']) & \
                    (trigs['chisq'] <
                     self.thresholds['reduced_chisq']) & \
                    (trigs['snr'] >
                     self.thresholds['newsnr'])
        if not np.any(valid_idx):
            return None

        cutdurchi_trigs = {k: trigs[k][valid_idx] for k in trigs}

        # This uses the pycbc live convention of chisq always meaning the
        # reduced chisq.
        nsnr_all = ranking.newsnr(cutdurchi_trigs['snr'],
                                  cutdurchi_trigs['chisq'])
        nsnr_idx = nsnr_all > self.thresholds['newsnr']
        if not np.any(nsnr_idx):
            return None

        cutall_trigs = {k: cutdurchi_trigs[k][nsnr_idx]
                        for k in trigs}

        # 'cluster' by taking the maximal newsnr value over the trigger set
        i = nsnr_all[nsnr_idx].argmax()

        nsnr = nsnr_all[nsnr_idx][i]
        dur = cutall_trigs['template_duration'][i]

        # create the coincidence
        fake_coinc = {'foreground/%s/%s' % (self.ifo, k):
                      cutall_trigs[k][i] for k in trigs}
        fake_coinc['foreground/stat'] = nsnr
        fake_coinc['foreground/ifar'] = self.calculate_ifar(nsnr, dur)
        fake_coinc['HWINJ'] = data_reader.near_hwinj()

        return fake_coinc

    def calculate_ifar(self, newsnr, duration):
        if self.fit_info['fixed_ifar']:
            return self.fit_info['fixed_ifar']
        dur_bin = self.fit_info['bins'][duration]
        rate = self.fit_info['rates'][dur_bin]
        coeff = self.fit_info['coeffs'][dur_bin]
        rate_louder = rate * fits.cum_fit('exponential', [newsnr],
                                          coeff, self.fit_info['thresh'])[0]
        # apply a trials factor of the number of duration bins
        rate_louder *= len(self.fit_info['rates'])
        return conv.sec_to_year(1. / rate_louder)
