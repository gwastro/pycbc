""" utilities for assigning FAR to single detector triggers
"""
import logging
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
        self.fit_group = fit_file[ifo + '/' + sngl_ifo_est_dist]
        self.fit_live_time = fit_file[ifo].attrs['live_time']
        self.fit_threshold = fit_file.attrs['fit_threshold']
        self.fixed_ifar = fixed_ifar
        dur_bin_edges = fit_file['bins_edges'][:]
        self.duration_bins = bin_utils.IrregularBins(dur_bin_edges)

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

    def check(self, triggers, data_reader):
        """ Look for a single detector trigger that passes the thresholds in
        the current data.
        """
        if len(triggers['snr']) == 0:
            return None

        # Apply cuts to triggers before clustering
        valid_idx = (triggers['template_duration'] > self.duration_threshold) & \
                        (triggers['chisq'] < self.reduced_chisq_threshold)
        if not np.count_nonzero(valid_idx):
            return None
        cutdurchi_triggers = {k: triggers[k][valid_idx] for k in triggers}
        nsnr_all = ranking.newsnr(cutdurchi_triggers['snr'],
                                  cutdurchi_triggers['chisq'])
        nsnr_idx = nsnr_all > self.newsnr_threshold
        if not np.count_nonzero(nsnr_idx):
            return None
        cutall_triggers = {k: cutdurchi_triggers[k][nsnr_idx] for k in triggers}

        #'cluster' by taking the maximal newsnr value over the trigger set
        i = nsnr_all[nsnr_idx].argmax()

        # This uses the pycbc live convention of chisq always meaning the
        # reduced chisq.
        rchisq = cutall_triggers['chisq'][i]
        nsnr = ranking.newsnr(cutall_triggers['snr'][i], rchisq)
        dur = cutall_triggers['template_duration'][i]

        if nsnr > self.newsnr_threshold and \
                dur > self.duration_threshold and \
                rchisq < self.reduced_chisq_threshold:
            
            fake_coinc = {'foreground/%s/%s' % (self.ifo, k): cutall_triggers[k][i]
                          for k in triggers}
            fake_coinc['foreground/stat'] = nsnr
            fake_coinc['foreground/ifar'] = self.calculate_ifar(nsnr, dur)
            fake_coinc['HWINJ'] = data_reader.near_hwinj()
            return fake_coinc
        return None

    def calculate_ifar(self, newsnr, duration):
        if self.sngl_ifo_est_dist == 'fixed':
            return self.fixed_ifar
        dur_bin = self.duration_bins[duration]
        count = self.fit_group['counts'][dur_bin]
        coeff = self.fit_group['fit_coeff'][dur_bin]
        n_louder = count * fits.cum_fit('exponential', [newsnr],
                                        coeff, self.fit_threshold)[0]
        return conv.sec_to_year(self.fit_live_time / n_louder)
