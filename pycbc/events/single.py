""" utilities for assigning FAR to single detector triggers
"""
from pycbc.events import ranking, trigger_fits as fits
from pycbc.types import MultiDetOptionAction
from pycbc import conversions as conv
import logging

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
        parser.add_argument('--sngl-ifar-est-dist', nargs='+',
                            action=MultiDetOptionAction)
        parser.add_argument('--single-duration-threshold', nargs='+',
                            type=float, action=MultiDetOptionAction)

    @classmethod
    def from_cli(cls, args, ifo):
        if args.single_fixed_ifar:
            self.sngl_ifo_est_dist = 'fixed'
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

        i = triggers['snr'].argmax()
        # This uses the pycbc live convention of chisq always meaning the
        # reduced chisq.
        rchisq = triggers['chisq'][i]
        nsnr = ranking.newsnr(triggers['snr'][i], rchisq)
        dur = triggers['template_duration'][i]

        if nsnr >  self.newsnr_threshold: # and \
#                dur > self.duration_threshold: #and \
#                rchisq < self.reduced_chisq_threshold:
            fake_coinc = {'foreground/%s/%s' % (self.ifo, k): triggers[k][i]
                          for k in triggers}
            fake_coinc['foreground/stat'] = nsnr
            fake_coinc['foreground/ifar'] = self.calculate_ifar(nsnr)
            fake_coinc['HWINJ'] = data_reader.near_hwinj()
            print("Newsnr(=stat): {}\nifar: {}\nduration: {}\nrchisq: {}\nsnr: {}".format(nsnr, fake_coinc['foreground/ifar'], dur, rchisq, triggers['snr'][i]))
            return fake_coinc
        return None

    def calculate_ifar(self, newsnr, fit_threshold=6.0):
        if self.sngl_ifo_est_dist == 'fixed':
            return self.fixed_ifar[self.ifo]
        if self.sngl_ifo_est_dist == 'conservative':
            cct = single_fits_coeff_count_time['conservative']
        elif self.sngl_ifo_est_dist == 'individual':
            cct = single_fits_coeff_count_time[self.ifo]
        print('using cct values {}'.format(cct))
        if newsnr < fit_threshold:
            n_louder = cct[1] + 1
            logging.WARN('newsnr is below fit threshold, this is an upper '
                         'bound to IFAR')
        else:
            n_louder = cct[1] * fits.cum_fit('exponential', [newsnr], cct[0],
                                             fit_threshold)
        return conv.sec_to_year(cct[2] / (n_louder + 1))[0]


single_fits_coeff_count_time = {
    'H1': (5.69015, 107893373, 10184192),
    'L1': (5.72713, 95714610, 10908800),
    'V1': (5.58611, 112008750, 11468656),
    'conservative': (5.03662, 112008750, 10184192)
}
