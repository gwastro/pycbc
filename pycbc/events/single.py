""" utilities for assigning FAR to single detector triggers
"""
import sys
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
                 sngl_ifo_est_dist='conservative',
                 fixed_ifar=None):
        if sngl_ifo_est_dist not in ['conservative', 'mean', 'fixed']:
            raise NotImplementedError("Single fitting distribution must be "
                                      "conservative, mean or fixed")
        self.ifo = ifo
        self.newsnr_threshold = newsnr_threshold
        self.reduced_chisq_threshold = reduced_chisq_threshold
        self.duration_threshold = duration_threshold
        self.sngl_ifo_est_dist = sngl_ifo_est_dist
        self.fit_file = fit_file
        self.fixed_ifar = fixed_ifar

    @staticmethod
    def insert_args(parser):
        parser.add_argument('--single-newsnr-threshold', nargs='+',
                            type=float, action=MultiDetOptionAction)
        parser.add_argument('--single-reduced-chisq-threshold', nargs='+',
                            type=float, action=MultiDetOptionAction)
        parser.add_argument('--single-fixed-ifar', nargs='+', default=None,
                            type=float, action=MultiDetOptionAction)
        parser.add_argument('--fit-definition-file',
                            help="File which contains definitons of fit "
                                 "coefficients and counts for specific "
                                 "single trigger IFAR fitting.")
        if '--single-fixed-ifar' in sys.argv:
            parser.add_argument('--sngl-ifar-est-dist', nargs='+',
                                default='fixed', action=MultiDetOptionAction)
        else:
            parser.add_argument('--sngl-ifar-est-dist', nargs='+',
                                default='conservative',
                                action=MultiDetOptionAction)
        parser.add_argument('--single-duration-threshold', nargs='+',
                            type=float, action=MultiDetOptionAction)

    @classmethod
    def from_cli(cls, args, ifo):
        return cls(
           ifo, newsnr_threshold=args.single_newsnr_threshold[ifo],
           reduced_chisq_threshold=args.single_reduced_chisq_threshold[ifo],
           duration_threshold=args.single_duration_threshold[ifo],
           fit_file=h5py.File(args.fit_definition_file),
           sngl_ifo_est_dist=args.sngl_ifar_est_dist[ifo],
           fixed_ifar=args.single_fixed_ifar
           )

    def check(self, trigs, data_reader):
        """ Look for a single detector trigger that passes the thresholds in
        the current data.
        """
        if len(trigs['snr']) == 0:
            return None

        # Apply cuts to trigs before clustering
        valid_idx = (trigs['template_duration'] > self.duration_threshold) & \
                    (trigs['chisq'] < self.reduced_chisq_threshold)
        if not np.count_nonzero(valid_idx):
            return None
        cutdurchi_trigs = {k: trigs[k][valid_idx] for k in trigs}
        nsnr_all = ranking.newsnr(cutdurchi_trigs['snr'],
                                  cutdurchi_trigs['chisq'])
        nsnr_idx = nsnr_all > self.newsnr_threshold
        if not np.count_nonzero(nsnr_idx):
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

        if nsnr > self.newsnr_threshold and \
                dur > self.duration_threshold and \
                rchisq < self.reduced_chisq_threshold:

            fake_coinc = {'foreground/%s/%s' % (self.ifo, k):
                          cutall_trigs[k][i] for k in trigs}
            fake_coinc['foreground/stat'] = nsnr
            fake_coinc['foreground/ifar'] = self.calculate_ifar(nsnr, dur)
            fake_coinc['HWINJ'] = data_reader.near_hwinj()
            return fake_coinc
        return None

    def calculate_ifar(self, newsnr, duration):
        if self.sngl_ifo_est_dist == 'fixed':
            return self.fixed_ifar
        dur_bin_edges = self.fit_file['bins_edges'][:]
        duration_bins = bin_utils.IrregularBins(dur_bin_edges)
        dur_bin = duration_bins[duration]
        dist_grp = self.fit_file[self.ifo + '/' + self.sngl_ifo_est_dist]
        count = dist_grp['counts'][dur_bin]
        coeff = dist_grp['fit_coeff'][dur_bin]
        fit_live_time = self.fit_file[self.ifo].attrs['live_time']
        fit_thresh = self.fit_file.attrs['fit_threshold']
        n_louder = count * fits.cum_fit('exponential', [newsnr],
                                        coeff, fit_thresh)[0]
        return conv.sec_to_year(fit_live_time / n_louder)
