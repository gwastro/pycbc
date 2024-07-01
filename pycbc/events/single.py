""" utilities for assigning FAR to single detector triggers
"""
import logging
import copy
import h5py
import numpy as np

from pycbc.events import trigger_fits as fits, stat
from pycbc.types import MultiDetOptionAction
from pycbc import conversions as conv
from pycbc import bin_utils

logger = logging.getLogger('pycbc.events.single')


class LiveSingle(object):
    def __init__(self, ifo,
                 ranking_threshold=10.0,
                 reduced_chisq_threshold=5,
                 duration_threshold=0,
                 fit_file=None,
                 sngl_ifar_est_dist=None,
                 fixed_ifar=None,
                 maximum_ifar=None,
                 statistic=None,
                 sngl_ranking=None,
                 stat_files=None,
                 **kwargs):
        self.ifo = ifo
        self.fit_file = fit_file
        self.sngl_ifar_est_dist = sngl_ifar_est_dist
        self.fixed_ifar = fixed_ifar
        self.maximum_ifar = maximum_ifar

        stat_class = stat.get_statistic(statistic)
        self.stat_calculator = stat_class(
            sngl_ranking,
            stat_files,
            ifos=[ifo],
            **kwargs
        )

        self.thresholds = {
            "ranking": ranking_threshold,
            "reduced_chisq": reduced_chisq_threshold,
            "duration": duration_threshold}

    @staticmethod
    def insert_args(parser):
        parser.add_argument('--single-ranking-threshold', nargs='+',
                            type=float, action=MultiDetOptionAction,
                            help='Single ranking threshold for '
                                 'single-detector events. Can be given '
                                 'as a single value or as detector-value '
                                 'pairs, e.g. H1:6 L1:7 V1:6.5')
        parser.add_argument('--single-reduced-chisq-threshold', nargs='+',
                            type=float, action=MultiDetOptionAction,
                            help='Maximum reduced chi-squared threshold for '
                                 'single triggers. Calcuated after any PSD '
                                 'variation reweighting is applied. Can be '
                                 'given as a single value or as '
                                 'detector-value pairs, e.g. H1:2 L1:2 V1:3')
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
        parser.add_argument('--single-maximum-ifar', nargs='+',
                            type=float, action=MultiDetOptionAction,
                            help='A maximum possible value for IFAR for '
                                 'single-detector events. Can be given as '
                                 'a single value or as detector-value pairs, '
                                 'e.g. H1:100 L1:1000 V1:50')
        parser.add_argument('--single-fit-file',
                            help='File which contains definitons of fit '
                                 'coefficients and counts for specific '
                                 'single trigger IFAR fitting.')
        parser.add_argument('--sngl-ifar-est-dist', nargs='+',
                            action=MultiDetOptionAction,
                            help='Which trigger distribution to use when '
                                 'calculating IFAR of single triggers. '
                                 'Can be given as a single value or as '
                                 'detector-value pairs, e.g. H1:mean '
                                 'L1:mean V1:conservative')

    @staticmethod
    def verify_args(args, parser, ifos):
        sngl_opts = [args.single_reduced_chisq_threshold,
                     args.single_duration_threshold,
                     args.single_ranking_threshold,
                     args.sngl_ifar_est_dist]

        sngl_opts_str = ("--single-reduced-chisq-threshold, "
                         "--single-duration-threshold, "
                         "--single-ranking-threshold, "
                         "--sngl-ifar-est-dist")

        if any(sngl_opts) and not all(sngl_opts):
            parser.error(f"Single detector trigger options ({sngl_opts_str}) "
                         "must either all be given or none.")

        if args.enable_single_detector_upload \
                and not args.enable_gracedb_upload:
            parser.error("--enable-single-detector-upload requires "
                         "--enable-gracedb-upload to be set.")

        sngl_optional_opts = [args.single_fixed_ifar,
                              args.single_fit_file,
                              args.single_maximum_ifar]
        sngl_optional_opts_str = ("--single-fixed-ifar, "
                                  "--single-fit-file,"
                                  "--single-maximum-ifar")
        if any(sngl_optional_opts) and not all(sngl_opts):
            parser.error("Optional singles options "
                         f"({sngl_optional_opts_str}) given but not all "
                         f"required options ({sngl_opts_str}) are.")

        for ifo in ifos:
            # Check which option(s) are needed for each IFO and if they exist:

            # Notes for the logic here:
            # args.sngl_ifar_est_dist.default_set is True if single value has
            # been set to be the same for all values
            # bool(args.sngl_ifar_est_dist) is True if option is given
            if args.sngl_ifar_est_dist and \
                    not args.sngl_ifar_est_dist.default_set \
                    and not args.sngl_ifar_est_dist[ifo]:
                # Option has been given, different for each IFO,
                # and this one is not present
                parser.error("All IFOs required in --single-ifar-est-dist "
                             "if IFO-specific options are given.")

            if args.sngl_ifar_est_dist[ifo] is None:
                # Default - no singles being used
                continue

            if not args.sngl_ifar_est_dist[ifo] == 'fixed':
                if not args.single_fit_file:
                    # Fixed IFAR option doesnt need the fits file
                    parser.error(f"Single detector trigger fits file must be "
                                 "given if --single-ifar-est-dist is not "
                                 f"fixed for all ifos (at least {ifo} has "
                                 f"option {args.sngl_ifar_est_dist[ifo]}).")
                if ifo in args.single_fixed_ifar:
                    parser.error(f"Value {args.single_fixed_ifar[ifo]} given "
                                 f"for {ifo} in --single-fixed-ifar, but "
                                 f"--single-ifar-est-dist for {ifo} "
                                 f"is {args.sngl_ifar_est_dist[ifo]}, not "
                                 "fixed.")
            else:
                # Check that the fixed IFAR value has actually been
                # given if using this instead of a distribution
                if not args.single_fixed_ifar[ifo]:
                    parser.error(f"--single-fixed-ifar must be "
                                 "given if --single-ifar-est-dist is fixed. "
                                 f"This is true for at least {ifo}.")

        # Return value is a boolean whether we are analysing singles or not
        # The checks already performed mean that all(sngl_opts) is okay
        return all(sngl_opts)

    @classmethod
    def from_cli(cls, args, ifo):
        # Allow None inputs
        stat_files = args.statistic_files or []
        stat_keywords = args.statistic_keywords or []

        # flatten the list of lists of filenames to a single list
        # (may be empty)
        stat_files = sum(stat_files, [])

        kwargs = stat.parse_statistic_keywords_opt(stat_keywords)
        return cls(
           ifo, ranking_threshold=args.single_ranking_threshold[ifo],
           reduced_chisq_threshold=args.single_reduced_chisq_threshold[ifo],
           duration_threshold=args.single_duration_threshold[ifo],
           fixed_ifar=args.single_fixed_ifar,
           maximum_ifar=args.single_maximum_ifar[ifo],
           fit_file=args.single_fit_file,
           sngl_ifar_est_dist=args.sngl_ifar_est_dist[ifo],
           statistic=args.ranking_statistic,
           sngl_ranking=args.sngl_ranking,
           stat_files=stat_files,
           **kwargs
           )

    def check(self, trigs, data_reader):
        """ Look for a single detector trigger that passes the thresholds in
        the current data.
        """

        # Apply cuts to trigs before clustering
        # Cut on snr so that triggers which could not reach the ranking
        # threshold do not have ranking calculated
        if 'psd_var_val' in trigs:
            # We should apply the PSD variation rescaling, as this can
            # re-weight the SNR to be above SNR
            trig_chisq = trigs['chisq'] / trigs['psd_var_val']
            trig_snr = trigs['snr'] / (trigs['psd_var_val'] ** 0.5)
        else:
            trig_chisq = trigs['chisq']
            trig_snr = trigs['snr']

        valid_idx = (trigs['template_duration'] >
                     self.thresholds['duration']) & \
                    (trig_chisq <
                     self.thresholds['reduced_chisq']) & \
                    (trig_snr >
                     self.thresholds['ranking'])
        if not np.any(valid_idx):
            return None

        cut_trigs = {k: trigs[k][valid_idx] for k in trigs}

        # Convert back from the pycbc live convention of chisq always
        # meaning the reduced chisq.
        trigsc = copy.copy(cut_trigs)
        trigsc['ifo'] = self.ifo
        trigsc['chisq'] = cut_trigs['chisq'] * cut_trigs['chisq_dof']
        trigsc['chisq_dof'] = (cut_trigs['chisq_dof'] + 2) / 2

        # Calculate the ranking reweighted SNR for cutting
        single_rank = self.stat_calculator.get_sngl_ranking(trigsc)
        sngl_idx = single_rank > self.thresholds['ranking']
        if not np.any(sngl_idx):
            return None

        cutall_trigs = {k: trigsc[k][sngl_idx]
                        for k in trigs}

        # Calculate the ranking statistic
        sngl_stat = self.stat_calculator.single(cutall_trigs)
        rank = self.stat_calculator.rank_stat_single((self.ifo, sngl_stat))

        # 'cluster' by taking the maximal statistic value over the trigger set
        i = rank.argmax()

        # calculate the (inverse) false-alarm rate
        ifar = self.calculate_ifar(
            rank[i],
            trigsc['template_duration'][i]
        )
        if ifar is None:
            return None

        # fill in a new candidate event
        candidate = {
            f'foreground/{self.ifo}/{k}': cutall_trigs[k][i] for k in trigs
        }
        candidate['foreground/stat'] = rank[i]
        candidate['foreground/ifar'] = ifar
        candidate['HWINJ'] = data_reader.near_hwinj()
        return candidate

    def calculate_ifar(self, sngl_ranking, duration):
        logger.info("Calculating IFAR")
        if self.fixed_ifar and self.ifo in self.fixed_ifar:
            return self.fixed_ifar[self.ifo]

        try:
            with h5py.File(self.fit_file, 'r') as fit_file:
                bin_edges = fit_file['bins_edges'][:]
                live_time = fit_file[self.ifo].attrs['live_time']
                thresh = fit_file.attrs['fit_threshold']
                dist_grp = fit_file[self.ifo][self.sngl_ifar_est_dist]
                rates = dist_grp['counts'][:] / live_time
                coeffs = dist_grp['fit_coeff'][:]
        except FileNotFoundError:
            logger.error(
                'Single fit file %s not found; '
                'dropping a potential single-detector candidate!',
                self.fit_file
            )
            return None

        bins = bin_utils.IrregularBins(bin_edges)
        dur_bin = bins[duration]

        rate = rates[dur_bin]
        coeff = coeffs[dur_bin]
        if np.isnan(coeff) or np.isnan(rate):
            logger.warning(
                "Single trigger fits are not valid - singles "
                "cannot be assessed for this detector at this time."
            )
            return None

        rate_louder = rate * fits.cum_fit(
            'exponential',
            [sngl_ranking],
            coeff,
            thresh
        )[0]

        # apply a trials factor of the number of duration bins
        rate_louder *= len(rates)

        return min(conv.sec_to_year(1. / rate_louder), self.maximum_ifar)
