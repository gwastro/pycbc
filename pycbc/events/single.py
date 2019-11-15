""" utilities for assigning FAR to single detector triggers
"""
from pycbc.events import ranking, trigger_fits as fits
from pycbc.types import MultiDetOptionAction
from pycbc import bin_utils, conversions as conv

class LiveSingle(object):
    def __init__(self, ifo,
                 newsnr_threshold=10.0,
                 reduced_chisq_threshold=5,
                 duration_threshold=0):
        self.ifo = ifo
        self.newsnr_threshold = newsnr_threshold
        self.reduced_chisq_threshold = reduced_chisq_threshold
        self.duration_threshold = duration_threshold

    @staticmethod
    def insert_args(parser):
        parser.add_argument('--single-newsnr-threshold', nargs='+',
                            type=float, action=MultiDetOptionAction)
        parser.add_argument('--single-reduced-chisq-threshold', nargs='+',
                            type=float, action=MultiDetOptionAction)
        parser.add_argument('--single-fixed-ifar', nargs='+',
                            type=float, action=MultiDetOptionAction)
        parser.add_argument('--single-duration-threshold', nargs='+',
                            type=float, action=MultiDetOptionAction)

    @classmethod
    def from_cli(cls, args, ifo):
        return cls(
           ifo, newsnr_threshold=args.single_newsnr_threshold[ifo],
           reduced_chisq_threshold=args.single_reduced_chisq_threshold[ifo],
           fixed_ifar=args.single_fixed_ifar[ifo],
           duration_threshold=args.single_duration_threshold[ifo],
           )

    def check(self, triggers, data_reader, conservative=False):
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

        if nsnr > self.newsnr_threshold and \
                rchisq < self.reduced_chisq_threshold and \
                dur > self.duration_threshold:
            fake_coinc = {'foreground/%s/%s' % (self.ifo, k): triggers[k][i]
                          for k in triggers}
            fake_coinc['foreground/stat'] = nsnr
            fake_coinc['foreground/ifar'] = calculate_ifar(triggers, nsnr, ifo,
                                                fixed_ifar=self.fixed_ifar,
                                                conservative)
            fake_coinc['HWINJ'] = data_reader.near_hwinj()
            return fake_coinc
        return None

class LiveSingleFarThreshold(LiveSingle):
    def __init__(self, ifo,
                 newsnr_threshold=10.0,
                 reduced_chisq_threshold=5,
                 duration_threshold=0,
                 fixed_ifar=0):
        super(LiveSingleFarThreshold).__init__(self, ifo,
                 newsnr_threshold=10.0,
                 reduced_chisq_threshold=5,
                 duration_threshold=0)
        self.fixed_ifar = fixed_ifar

class LiveSingleCalculateIfar(LiveSingle):
    def __init__(self, ifo,
                 newsnr_threshold=10.0,
                 reduced_chisq_threshold=5,
                 duration_threshold=0):
        super(LiveSingleFarThreshold).__init__(self, ifo,
                 newsnr_threshold=10.0,
                 reduced_chisq_threshold=5,
                 duration_threshold=0)
        self.fixed_ifar = False

single_fits_coeff_count_time = {
    'H1' : (5.69015, 107893373, 10184192),
    'L1' : (5.72713,  95714610, 10908800),
    'V1' : (5.58611, 112008750, 11468656),
    'conservative' : (5.03662, 28909432, 10184192)
}

def calculate_ifar(triggers, newsnr, ifo, fit_threshold=6.0,
                   fixed_ifar=False, conservative=False):
    if fixed_ifar: return fixed_ifar
    # Fit coefficients for each ifo based on O3a trigger fits files
    # 1 / (total counts=1.6e7 * exp(-alpha=5.7 * (SNR=9.776 - thresh=6) / number_templates=200000 / live_time=1.6e6) and convert to years
    if conservative:
        cct_to_use = 'conservative'
    else:
        cct_to_use = ifo
    cct = single_fits_coeff_count_time[cct_to_use]
    n_louder = cct[1] * fits.cum_fit('exponential', newsnr, cct[0], fit_threshold)
    return conv.sec_to_year(cct[2] / (n_louder))
