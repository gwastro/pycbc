""" utilities for assigning FAR to single detector triggers
"""
import logging
import copy
import threading
import time
import numpy as np

from pycbc.events import trigger_fits as fits, stat
from pycbc.types import MultiDetOptionAction
from pycbc import conversions as conv
from pycbc.io.hdf import HFile
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
                 stat_features=None,
                 stat_keywords=None,
                 sngl_ranking=None,
                 stat_files=None,
                 statistic_refresh_rate=None,
                 **kwargs):
        """
        Parameters
        ----------
        ifo: str
            Name of the ifo that is being analyzed
        newsnr_threshold: float
            Minimum value for the reweighted SNR of the event under
            consideration. Which reweighted SNR is defined by sngl_ranking
        reduced_chisq_threshold: float
            Maximum value for the reduced chisquared of the event under
            consideration
        duration_threshold: float
            Minimum value for the duration of the template which found the
            event under consideration
        fit_file: str or path
            (optional) the file containing information about the
            single-detector event significance distribution fits
        sngl_ifar_est_dist: str
            Which trigger distribution to use when calculating IFAR of
            single-detector events
        fixed_ifar: float
            (optional) give a fixed IFAR value to any event which passes the
            threshold criteria
        statistic: str
            The name of the statistic to rank events.
        stat_features: list of str
            The names of features to use for the statistic
        stat_keywords: list of key:value strings
            argument for statistic keywords
        sngl_ranking: str
            The single detector ranking to use with the background statistic
        stat_files: list of strs
            List of filenames that contain information used to construct
            various coincident statistics.
        maximum_ifar: float
            The largest inverse false alarm rate in years that we would like to
            calculate.
        statistic_refresh_rate: float
            How regularly to run the update_files method on the statistic
            class (in seconds), default not do do this
        kwargs: dict
            Additional options for the statistic to use. See stat.py
            for more details on statistic options.
        """
        self.ifo = ifo
        self.fit_file = fit_file
        self.sngl_ifar_est_dist = sngl_ifar_est_dist
        self.fixed_ifar = fixed_ifar
        self.maximum_ifar = maximum_ifar

        self.time_stat_refreshed = time.time()
        self.stat_calculator_lock = threading.Lock()
        self.statistic_refresh_rate = statistic_refresh_rate

        stat_class = stat.get_statistic(statistic)
        stat_extras = stat.parse_statistic_feature_options(
            stat_features,
            stat_keywords,
        )
        self.stat_calculator = stat_class(
            sngl_ranking,
            stat_files,
            ifos=[ifo],
            **stat_extras
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
           stat_features=args.statistic_features,
           stat_keywords=args.statistic_keywords,
           stat_files=stat_files,
           statistic_refresh_rate=args.statistic_refresh_rate,
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
        with self.stat_calculator_lock:
            single_rank = self.stat_calculator.get_sngl_ranking(trigsc)

        sngl_idx = single_rank > self.thresholds['ranking']
        if not np.any(sngl_idx):
            return None

        cutall_trigs = {k: trigsc[k][sngl_idx]
                        for k in trigs}

        # Calculate the ranking statistic
        with self.stat_calculator_lock:
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
            f'foreground/{self.ifo}/{k}': cut_trigs[k][sngl_idx][i]
            for k in trigs
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
            with HFile(self.fit_file, 'r') as fit_file:
                bin_edges = fit_file['bins_edges'][:]
                live_time = fit_file[self.ifo].attrs['live_time']
                thresh = fit_file[self.ifo].attrs['fit_threshold']
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

    def start_refresh_thread(self):
        """
        Start a thread managing whether the stat_calculator will be updated
        """
        if self.statistic_refresh_rate is None:
            logger.info("Statistic refresh disabled for %s", self.ifo)
            return
        thread = threading.Thread(
            target=self.refresh_statistic,
            daemon=True,
            name="Stat refresh " + self.ifo
        )
        logger.info("Starting %s statistic refresh thread", self.ifo)
        thread.start()

    def refresh_statistic(self):
        """
        Function to refresh the stat_calculator at regular intervals
        """
        while True:
            # How long since the statistic was last updated?
            since_stat_refresh = time.time() - self.time_stat_refreshed
            if since_stat_refresh > self.statistic_refresh_rate:
                self.time_stat_refreshed = time.time()
                logger.info(
                    "Checking %s statistic for updated files", self.ifo
                )
                with self.stat_calculator_lock:
                    self.stat_calculator.check_update_files()
            # Sleep one second for safety
            time.sleep(1)
            # Now use the time it took the check / update the statistic
            since_stat_refresh = time.time() - self.time_stat_refreshed
            logger.debug(
                "%s statistic: Waiting %.3fs for next refresh",
                self.ifo,
                self.statistic_refresh_rate - since_stat_refresh
            )
            time.sleep(self.statistic_refresh_rate - since_stat_refresh)
