"""
Chisq based on sine-gaussian tiles.
See https://arxiv.org/abs/1709.08974 for a discussion.
"""

import numpy

from pycbc.events import ranking
from pycbc.filter import sigma
from pycbc.opt import LimitedSizeDict
from pycbc.types import complex_same_precision_as
from pycbc.vetoes.chisq import SingleDetPowerChisq
from pycbc.waveform import sinegauss
from pycbc.waveform.utils import apply_fseries_time_shift


class SingleDetSGChisq(SingleDetPowerChisq):
    """
    Class that handles precomputation and memory management for efficiently
    running the sine-Gaussian chisq
    """

    name = "sg"
    columns = {"sg_chisq": numpy.float32}
    default_values = {"sg_chisq": 1}
    returns = {"sg_chisq": numpy.float32}

    def __init__(self, bank, num_bins=0, snr_threshold=None, chisq_locations=None):
        """
        Create sine-Gaussian Chisq Calculator

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
            The region is a boolean expression such as 'mtotal>40' indicating
            which templates to apply this set of sine-Gaussians to.

        """
        self.snr_threshold = snr_threshold
        self.num_bins = num_bins
        self._tile_cache = LimitedSizeDict(size_limit=256)
        if chisq_locations is not None:
            self.do = True
            self.params = {}
            for descr in chisq_locations:
                region, values = descr.split(":")
                mask = bank.table.parse_boolargs([(1, region), (0, "else")])[0]
                hashes = bank.table["template_hash"][mask.astype(bool)]
                for h in hashes:
                    self.params[h] = values
        else:
            self.do = False

    @staticmethod
    def insert_option_group(parser):
        group = parser.add_argument_group("Sine-Gaussian Chisq")
        group.add_argument(
            "--sgchisq-snr-threshold",
            type=float,
            help="Minimum SNR threshold to use SG chisq",
        )
        group.add_argument(
            "--sgchisq-locations",
            type=str,
            nargs="+",
            help="Frequency offsets and quality factors of the sine-Gaussians"
            " to use, format 'region-boolean:q1-offset1,q2-offset2'. "
            "Offset is relative to the end frequency of the approximant."
            " Region is a boolean expression selecting templates to "
            "apply the sine-Gaussians to, ex. 'mtotal>40'",
        )

    @classmethod
    def from_cli(cls, args, bank=None, chisq_bins=None, **kwargs):
        if chisq_bins is None:
            chisq_bins = str(args.chisq_bins)
        return cls(bank, chisq_bins, args.sgchisq_snr_threshold, args.sgchisq_locations)

    def evaluate_batch(self, batch, j):
        """Evaluate sine-Gaussian chisq for trigger j in a batch."""
        stilde = batch.stilde
        overwhitened = getattr(stilde, "overwhitened", None)
        template = batch.template(j)
        val = self.values(
            overwhitened if overwhitened is not None else stilde,
            template,
            stilde.psd,
            numpy.array([batch.snr[j]]),
            1,
            indices=numpy.array([batch.time_index[j]]),
            overwhitened=overwhitened is not None,
        )
        return {"sg_chisq": val[0]}

    def values(
        self,
        stilde,
        template,
        psd,
        snrv,
        snr_norm,
        bchisq=None,
        bchisq_dof=None,
        indices=None,
        overwhitened=True,
        subtract_signal=True,
    ):
        """
        Calculate sine-Gaussian chisq

        Parameters
        ----------
        stilde: pycbc.types.FrequencySeries
            The data, divided by its PSD (overwhitened) unless overwhitened is
            False
        template: pycbc.types.FrequencySeries
            The waveform template being analyzed
        psd: pycbc.types.FrequencySeries
            The power spectral density of the data
        snrv: numpy.ndarray
            The peak unnormalized complex SNR values
        snr_norm: float
            The normalization factor for the snr
        bchisq: numpy.ndarray, optional
            The Bruce Allen power chisq values for these triggers
        bchisq_dof: numpy.ndarray, optional
            The degrees of freedom of the Bruce chisq
        indices: numpy.ndarray, optional
            The indices of the snr peaks.
        overwhitened: bool, optional
            Whether stilde is already divided by the PSD.
        subtract_signal: bool, optional
            Whether to subtract the expected signal overlap from each tile.

        Returns
        -------
        chisq: Array
            Chisq values, one for each sample index

        """
        if not self.do:
            return None

        chisq = numpy.ones(len(snrv))
        if template.params.template_hash not in self.params:
            return chisq

        snrs = numpy.array(snrv) * snr_norm
        if bchisq is not None and self.snr_threshold is not None:
            # Skip if newsnr too low
            rchisq = numpy.array(bchisq) / numpy.array(bchisq_dof)
            nsnr = ranking.newsnr(abs(snrs), rchisq)
            above = numpy.flatnonzero(nsnr >= self.snr_threshold)
            if not len(above):
                return chisq
        else:
            above = range(len(snrv))

        tiles = self.cached_tiles(
            stilde, psd, template, subtract_signal=subtract_signal
        )
        if tiles is None:
            return chisq

        klo = min(tile[1] for tile in tiles)
        khi = max(tile[2] for tile in tiles)
        band = stilde[:khi] if overwhitened else stilde[:khi] / psd[:khi]

        tlen = (len(stilde) - 1) * 2
        dt = 1.0 / (tlen * stilde.delta_f)
        for i in above:
            time = float(template.epoch) + dt * indices[i]
            # Shift the time of interest to be centered on 0
            stilde_shift = apply_fseries_time_shift(band, -time, kmin=klo)
            chisq[i] = 0
            # Calculate the sum of SNR^2 for the sine-Gaussians specified
            for gtem, kmin, kmax, gsigma, overlap in tiles:
                # Calculate the SNR of the tile
                gsnr = (gtem[kmin:kmax] * stilde_shift[kmin:kmax]).sum()
                gsnr *= 4.0 * gtem.delta_f / gsigma
                if subtract_signal:
                    gsnr -= snrs[i] * overlap
                    chisq[i] += abs(gsnr) ** 2.0 / (1.0 - abs(overlap) ** 2.0)
                else:
                    chisq[i] += abs(gsnr) ** 2.0
            chisq[i] /= 2 * len(tiles)
        return chisq

    def peak_frequency(self, template, psd=None):
        """Estimate the maximum frequency up to which the waveform has
        power by approximating power per frequency as constant over the
        last 2 chisq bins. We cannot use the final chisq bin edge as it
        does not have to be where the waveform terminates.
        """
        bins = self.cached_chisq_bins(template, getattr(template, "psd", psd))
        fstep = bins[-2] - bins[-3]
        return float((bins[-2] + fstep) * template.delta_f)

    def cached_tiles(self, stilde, psd, template, subtract_signal=True):
        template_hash = template.params.template_hash
        fpeak = self.peak_frequency(template, psd)
        key = (
            template_hash,
            id(psd),
            fpeak,
            template.f_lower,
            len(stilde),
            stilde.delta_f,
            subtract_signal,
        )
        if key not in self._tile_cache:
            # Only apply the sine-Gaussian in a +-50 Hz range around the
            # central frequency
            qwindow = 50
            # This is 90% of the Nyquist frequency of the data
            # This allows us to avoid issues near Nyquist due to resample
            # Filtering
            fstop = len(stilde) * stilde.delta_f * 0.9
            delta_f = stilde.delta_f
            fmin = int(template.f_lower / psd.delta_f) * delta_f
            if subtract_signal:
                step = round(template.delta_f / delta_f)
                t_psd = getattr(template, "psd", psd[::step][: len(template)])
                hsigma = sigma(
                    template, psd=t_psd, low_frequency_cutoff=template.f_lower
                )
            tiles = []
            for descr in self.params[template_hash].split(","):
                # Get the q and frequency offset from the descriptor
                q, offset = descr.split("-", 1)
                q, offset = float(q), float(offset)
                fcen = fpeak + offset
                flow = max(fmin, fcen - qwindow)
                fhigh = fcen + qwindow
                # If any sine-gaussian tile has an upper frequency near
                # nyquist return 1 instead.
                if fhigh > fstop:
                    self._tile_cache[key] = None
                    return None
                kmin = int(flow / delta_f)
                kmax = int(fhigh / delta_f)
                # Calculate sine-gaussian tile
                gtem = sinegauss.fd_sine_gaussian(
                    1.0, q, fcen, flow, (kmax + 1) * delta_f, delta_f
                ).astype(complex_same_precision_as(psd))
                gsigma = sigma(
                    gtem,
                    psd=psd,
                    low_frequency_cutoff=flow,
                    high_frequency_cutoff=fhigh,
                )
                if subtract_signal:
                    kl, kh = kmin // step, kmax // step
                    overlap = (
                        gtem[::step][kl:kh] * template[kl:kh] / t_psd[kl:kh]
                    ).sum()
                    overlap *= 4.0 * template.delta_f / (gsigma * hsigma)
                else:
                    overlap = 0.0
                tiles.append((gtem, kmin, kmax, gsigma, overlap))
            self._tile_cache[key] = tiles
        return self._tile_cache[key]
