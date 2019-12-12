""" utilities for assigning FAR to single detector triggers
"""
import logging
import sys
import numpy as np
from pycbc.events import ranking, trigger_fits as fits
from pycbc.types import MultiDetOptionAction
from pycbc import conversions as conv

class LiveSingle(object):
    def __init__(self, ifo,
                 newsnr_threshold=10.0,
                 reduced_chisq_threshold=5,
                 duration_threshold=0,
                 sngl_ifo_est_dist='conservative',
                 fixed_ifar=None):
        self.ifo = ifo
        self.newsnr_threshold = newsnr_threshold
        self.reduced_chisq_threshold = reduced_chisq_threshold
        self.duration_threshold = duration_threshold
        self.sngl_ifo_est_dist = sngl_ifo_est_dist
        self.fixed_ifar = fixed_ifar

    @staticmethod
    def insert_args(parser):
        parser.add_argument('--single-newsnr-threshold', nargs='+',
                            type=float, action=MultiDetOptionAction)
        parser.add_argument('--single-reduced-chisq-threshold', nargs='+',
                            type=float, action=MultiDetOptionAction)
        parser.add_argument('--single-fixed-ifar', nargs='+', default=None,
                            type=float, action=MultiDetOptionAction)
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
        dur_idx = triggers['template_duration'] > self.duration_threshold
        if not np.count_nonzero(dur_idx): return None
        cutdur_triggers = {k: triggers[k][dur_idx] for k in triggers}
        chisq_idx = cutdur_triggers['chisq'] < self.reduced_chisq_threshold
        if not np.count_nonzero(chisq_idx): return None
        cutchisq_triggers = {k: cutdur_triggers[k][chisq_idx] for k in triggers}
        nsnr_all = ranking.newsnr(cutchisq_triggers['snr'], cutchisq_triggers['chisq'])
        nsnr_idx = nsnr_all > self.newsnr_threshold
        if not np.count_nonzero(nsnr_idx): return None
        cutall_triggers = {k: cutchisq_triggers[k][nsnr_idx] for k in triggers}

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
            fake_coinc['foreground/ifar'] = self.calculate_ifar(nsnr)
            fake_coinc['HWINJ'] = data_reader.near_hwinj()
            return fake_coinc
        return None

    def calculate_ifar(self, newsnr, fit_threshold=6.5):
        if self.sngl_ifo_est_dist == 'fixed':
            return self.fixed_ifar[self.ifo]
        if self.sngl_ifo_est_dist == 'conservative':
            cct = SINGLE_FITS_COEFF_COUNT_TIME['conservative']
        elif self.sngl_ifo_est_dist == 'individual':
            cct = SINGLE_FITS_COEFF_COUNT_TIME[self.ifo]
        print('using cct values {}'.format(cct))
        if newsnr < fit_threshold:
            n_louder = cct[1] + 1
            logging.info('newsnr is below fit threshold, this is an upper '
                         'bound to IFAR')
        else:
            n_louder = cct[1] * fits.cum_fit('exponential', [newsnr], cct[0],
                                             fit_threshold)[0]
        return conv.sec_to_year(cct[2] / (n_louder + 1))


SINGLE_FITS_COEFF_COUNT_TIME = {
    'H1': (5.69015, 107893373, 10184192),
    'L1': (5.72713, 95714610, 10908800),
    'V1': (5.58611, 112008750, 11468656),
    'conservative': (5.03662, 112008750, 10184192)
}
