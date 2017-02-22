""" Chisq based on sine-gaussian tiles """

import numpy
import logging

from pycbc.waveform.utils import apply_fseries_time_shift
from pycbc.filter import sigma
from pycbc.waveform import sinegauss
from pycbc.vetoes.chisq import SingleDetPowerChisq
from pycbc.events import newsnr

class SingleDetSGChisq(SingleDetPowerChisq):
    """Class that handles precomputation and memory management for efficiently
    running the sine-Gaussian chisq
    """
    returns = {'sg_chisq': numpy.float32}
    
    def __init__(self, bank, num_bins=0,
                       snr_threshold=None,
                       chisq_locations=None):
        """ Create sine-Gaussian Chisq Calculator

        Parameters
        ----------
        bank: pycbc.waveform.TemplateBank
            The template bank that will be processed.
        num_bins: str
            The string determining the number of power chisq bins
        snr_threshold: float
            The threshold to calculate the sine-Gaussian chisq
        chisq_locations: list of strs
            List of strings which detail where to place a sine-Gaussian.
            The format is 'region-boolean:q1-offset1,q2-offset2'.
            The offset is relative to the end frequency of the approximant.
            The region is a boolean expresion such as 'mtotal>40' and indicates
            where to apply this set of sine-Gaussians.
        """
        if snr_threshold is not None:
            self.do = True
            self.num_bins = num_bins
            self.snr_threshold = snr_threshold
            self.params = {}
            for descr in chisq_locations:
                region, values = descr.split(":")
                mask = bank.table.parse_boolargs([(1, region), (0, 'else')])[0]
                hashes = bank.table['template_hash'][mask.astype(bool)]
                for h in hashes:
                    self.params[h] = values
        else:
            self.do = False
            
    @staticmethod
    def insert_option_group(parser):
        group = parser.add_argument_group("Sine-Gaussian Chisq")
        group.add_argument("--sgchisq-snr-threshold", type=float,
            help="Minimum SNR threshold to use SG chisq")
        group.add_argument("--sgchisq-locations", type=str, nargs='+',
            help="The frequencies and quality factors of the sine-gaussians"
                 " to use. The format is 'region-boolean:q1-offset1,q2-offset2'."
                 "The offset is relative to the end frequency of the approximant."
                 "The region is a boolean expresion such as 'mtotal>40' and indicates "
                 "where to apply this set of sine-Gaussians.")

    @classmethod
    def from_cli(cls, args, bank, chisq_bins):
        return cls(bank, chisq_bins,
                   args.sgchisq_snr_threshold,
                   args.sgchisq_locations)

    def values(self, stilde, template, psd, snrv, snr_norm,
                     bchisq, bchisq_dof, indices):
        """ Calculate sine-Gaussian chisq

        Parameters
        ----------
        stilde: pycbc.types.Frequencyseries
            The overwhitened strain
        template: pycbc.types.Frequencyseries
            The waveform template being analyzed
        psd: pycbc.types.Frequencyseries
            The power spectral density of the data
        snrv: numpy.ndarray
            The peak unnormalized complex SNR values
        snr_norm: float
            The normalization factor for the snr
        bchisq: numpy.ndarray
            The Bruce Allen power chisq values for these triggers
        bchisq_dof: numpy.ndarray
            The degrees of freedom of the Bruce chisq
        indics: numpy.ndarray
            The indices of the snr peaks.

        Returns
        -------
        chisq: Array
            Chisq values, one for each sample index
        """
        if not self.do:
            return None

        if template.params.template_hash not in self.params:
            return numpy.ones(len(snrv))
        values = self.params[template.params.template_hash].split(',')
        
        # Get the chisq bins to use as the frequency reference point
        bins = self.cached_chisq_bins(template, psd)

        # This is implemented slowly, so let's not call it often, OK?
        chisq = numpy.ones(len(snrv))
        for i, snrvi in enumerate(snrv):
            #Skip if newsnr too low
            snr = abs(snrvi * snr_norm)
            nsnr = newsnr(snr, bchisq[i] / bchisq_dof[i])
            if nsnr < self.snr_threshold:
                continue

            N = (len(template) - 1) * 2
            dt = 1.0 / (N * template.delta_f)
            kmin = int(template.f_lower / psd.delta_f)
            time = float(template.epoch) + dt * indices[i]
            # Shift the time of interest to be centered on 0
            stilde_shift = apply_fseries_time_shift(stilde, -time)

            # Only apply the sine-Gaussian in a +-50 Hz range around the 
            # central frequency
            qwindow = 50
            chisq[i] = 0
            
            # Estimate the maximum frequency up to which the waveform has
            # power by approximating power per frequency
            # as constant over the last 2 chisq bins. We cannot use the final
            # chisq bin edge as it does not have to be where the waveform
            # terminates.
            fstep = (bins[-2] - bins[-3])
            fpeak = (bins[-2] + fstep) * template.delta_f
            
            # This is 90% of the Nyquist frequency of the data
            # This allows us to avoid issues near Nyquist due to resample
            # Filtering
            fstop = len(stilde) * stilde.delta_f * 0.9
            
            dof = 0
            # Calculate the sum of SNR^2 for the sine-Gaussians specified
            for descr in values:
                # Get the q and frequency offset from the descriptor
                q, offset = descr.split('-')
                q, offset = float(q), float(offset)
                fcen = fpeak + offset
                flow = max(kmin * template.delta_f, fcen - qwindow)
                fhigh = fcen + qwindow

                # If any sine-gaussian tile has an upper frequency near 
                # nyquist return 1 instead.
                if fhigh > fstop:
                    return numpy.ones(len(snrv))

                kmin = int(flow / template.delta_f)
                kmax = int(fhigh / template.delta_f)

                #Calculate sine-gaussian tile
                gtem = sinegauss.fd_sine_gaussian(1.0, q, fcen, flow,
                                      len(template) * template.delta_f,
                                      template.delta_f).astype(numpy.complex64)
                gsigma = sigma(gtem, psd=psd,
                                     low_frequency_cutoff=flow, 
                                     high_frequency_cutoff=fhigh)
                #Calculate the SNR of the tile
                gsnr = (gtem[kmin:kmax] * stilde_shift[kmin:kmax]).sum()
                gsnr *= 4.0 * gtem.delta_f / gsigma
                chisq[i] += abs(gsnr)**2.0
                dof += 2
            if dof == 0:
                chisq[i] = 1
            else:
                chisq[i] /= dof
            logging.info('Found chisq %s', chisq[i])
        return chisq

