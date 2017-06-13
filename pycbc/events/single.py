""" utilities for assigning FAR to single detector triggers
"""
from pycbc.events import newsnr
from pycbc.io.live import SingleForGraceDB

class LiveSingleFarThreshold(object):
    def __init__(self, ifo, 
                 newsnr_threshold=10.0,
                 reduced_chisq_threshold=5,
                 duration_threshold=0,
                 fixed_ifar=0):
        self.ifo = ifo
        self.newsnr_threshold = newsnr_threshold
        self.reduced_chisq_threshold = reduced_chisq_threshold
        self.fixed_ifar = fixed_ifar
        self.duration_threshold = duration_threshold

    @staticmethod
    def insert_args(parser):
        parser.add_argument('--single-newsnr-threshold', type=float)
        parser.add_argument('--single-reduced-chisq-threshold', type=float)
        parser.add_argument('--single-fixed-ifar', type=float)
        parser.add_argument('--single-duration-threshold', type=float)

    @classmethod
    def from_cli(cls, args, ifo):
        return cls(ifo, newsnr_threshold=args.single_newsnr_threshold,
                   reduced_chisq_threshold=args.single_reduced_chisq_threshold,
                   fixed_ifar=args.single_fixed_ifar,
                   duration_threshold=args.single_duration_threshold,
                   )

    def check(self, triggers, data_reader, **kwds):
        """ Look for a single detector trigger that passes the thresholds in 
        the current data.
        """
        if len(triggers['snr']) == 0:
            return None 

        i = triggers['snr'].argmax()
        # This uses the pycbc live convention of chisq always meaning the 
        # reduced chisq. 
        rchisq = triggers['chisq'][i]
        nsnr = newsnr(triggers['snr'][i], rchisq)
        dur = triggers['template_duration'][i]

        if (nsnr > self.newsnr_threshold and
           rchisq < self.reduced_chisq_threshold and
            dur > self.duration_threshold):
            d = {key: triggers[key][i] for key in triggers}
            d['stat'] = nsnr
            d['ifar'] = self.fixed_ifar
            return SingleForGraceDB(self.ifo, d, 
                                    hardware_injection=data_reader.near_hwinj(), 
                                    **kwds)
        else:
            return None        
