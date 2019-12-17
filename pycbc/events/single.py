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
                 fit_bins=None,
                 fit_coeffs=None,
                 fit_rates=None,
                 fixed_ifar=None,
                 fit_thresh=None):
        self.ifo = ifo
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
                            type=float, action=MultiDetOptionAction)
        parser.add_argument('--single-reduced-chisq-threshold', nargs='+',
                            type=float, action=MultiDetOptionAction)
        parser.add_argument('--single-fixed-ifar', nargs='+', default=None,
                            type=float, action=MultiDetOptionAction)
        parser.add_argument('--single-fit-file',
                            help="File which contains definitons of fit "
                                 "coefficients and counts for specific "
                                 "single trigger IFAR fitting.")
        parser.add_argument('--sngl-ifar-est-dist', nargs='+',
                            default='conservative',
                            action=MultiDetOptionAction)
        parser.add_argument('--single-duration-threshold', nargs='+',
                            type=float, action=MultiDetOptionAction)

    @classmethod
    def from_cli(cls, args, ifo):
        if args.sngl_ifar_est_dist[ifo] not in \
                ['conservative', 'mean', 'fixed']:
            raise NotImplementedError("Single fitting distribution must be "
                                      "conservative, mean or fixed")
        if args.sngl_ifar_est_dist[ifo] == "fixed":
            fit_bins = None
            fit_rates = None
            fit_thresh = None
            fit_coeffs = None
        else:
            with h5py.File(args.single_fit_file, 'r') as fit_file:
                bin_edges = fit_file['bins_edges'][:]
                fit_bins = bin_utils.IrregularBins(bin_edges)
                dist_grp = fit_file[ifo + '/' + args.sngl_ifar_est_dist[ifo]]
                live_time = fit_file[ifo].attrs['live_time']
                fit_rates = dist_grp['counts'][:] / live_time
                fit_coeffs = dist_grp['fit_coeff'][:]
                fit_thresh = fit_file.attrs['fit_threshold']
        return cls(
           ifo, newsnr_threshold=args.single_newsnr_threshold[ifo],
           reduced_chisq_threshold=args.single_reduced_chisq_threshold[ifo],
           duration_threshold=args.single_duration_threshold[ifo],
           fixed_ifar=args.single_fixed_ifar,
           fit_bins=fit_bins,
           fit_coeffs=fit_coeffs,
           fit_rates=fit_rates,
           fit_thresh=fit_thresh
           )

    def check(self, trigs, data_reader):
        """ Look for a single detector trigger that passes the thresholds in
        the current data.
        """
        if len(trigs['snr']) == 0:
            return None

        # Apply cuts to trigs before clustering
        valid_idx = (trigs['template_duration'] >
                     self.thresholds['duration']) & \
                    (trigs['chisq'] < self.thresholds['reduced_chisq'])
        if not np.any(valid_idx):
            return None
        cutdurchi_trigs = {k: trigs[k][valid_idx] for k in trigs}
        nsnr_all = ranking.newsnr(cutdurchi_trigs['snr'],
                                  cutdurchi_trigs['chisq'])
        nsnr_idx = nsnr_all > self.thresholds['newsnr']
        if not np.any(nsnr_idx):
            return None
        cutall_trigs = {k: cutdurchi_trigs[k][nsnr_idx]
                        for k in trigs}

        # 'cluster' by taking the maximal newsnr value over the trigger set
        i = nsnr_all[nsnr_idx].argmax()

        # This uses the pycbc live convention of chisq always meaning the
        # reduced chisq.
        rchisq = cutall_trigs['chisq'][i]
        nsnr = ranking.newsnr(cutall_trigs['snr'][i], rchisq)
        dur = cutall_trigs['template_duration'][i]

        if nsnr > self.thresholds['newsnr'] and \
                dur > self.thresholds['duration'] and \
                rchisq < self.thresholds['reduced_chisq']:

            fake_coinc = {'foreground/%s/%s' % (self.ifo, k):
                          cutall_trigs[k][i] for k in trigs}
            fake_coinc['foreground/stat'] = nsnr
            fake_coinc['foreground/ifar'] = self.calculate_ifar(nsnr, dur)
            fake_coinc['HWINJ'] = data_reader.near_hwinj()
            print("Single found: {}".format(fake_coinc))
            exit()
            return fake_coinc
        return None

    def calculate_ifar(self, newsnr, duration):
        if self.fit_info['fixed_ifar']:
            return self.fit_info['fixed_ifar']
        dur_bin = self.fit_info['bins'][duration]
        rate = self.fit_info['rates'][dur_bin]
        coeff = self.fit_info['coeffs'][dur_bin]
        rate_louder = rate * fits.cum_fit('exponential', [newsnr],
                                          coeff, self.fit_info['thresh'])[0]
        # apply a trials factor of the number of duration bins
        rate_louder /= len(self.fit_info['rates'])
        return conv.sec_to_year(1. / rate_louder)
