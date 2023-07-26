""" utilities for assigning FAR to single detector triggers
"""
import logging
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
                 sngl_ifar_est_dist=None,
                 fixed_ifar=None):
        self.ifo = ifo
        self.fit_file = fit_file
        self.sngl_ifar_est_dist = sngl_ifar_est_dist
        self.fixed_ifar = fixed_ifar

        self.thresholds = {
            "newsnr": newsnr_threshold,
            "reduced_chisq": reduced_chisq_threshold,
            "duration": duration_threshold}

    @staticmethod
    def insert_args(parser):
        parser.add_argument('--single-newsnr-threshold', nargs='+',
                            type=float, action=MultiDetOptionAction,
                            help='Newsnr min threshold for single triggers. '
                                 'Can be given as a single value or as '
                                 'detector-value pairs, e.g. H1:6 L1:7 V1:6.5')
        parser.add_argument('--single-reduced-chisq-threshold', nargs='+',
                            type=float, action=MultiDetOptionAction,
                            help='Maximum reduced chi-squared threshold for '
                                 'single triggers. Can be given as a single '
                                 'value or as detector-value pairs, e.g. '
                                 'H1:2 L1:2 V1:3')
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
                     args.single_newsnr_threshold,
                     args.sngl_ifar_est_dist]
        sngl_opts_str = ("--single-reduced-chisq-threshold, "
                         "--single-duration-threshold, "
                         "--single-newsnr-threshold, "
                         "--sngl-ifar-est-dist")

        if any(sngl_opts) and not all(sngl_opts):
            parser.error(f"Single detector trigger options ({sngl_opts_str}) "
                         "must either all be given or none.")

        if args.enable_single_detector_upload \
                and not args.enable_gracedb_upload:
            parser.error("--enable-single-detector-upload requires "
                         "--enable-gracedb-upload to be set.")

        sngl_optional_opts = [args.single_fixed_ifar,
                              args.single_fit_file]
        sngl_optional_opts_str = ("--single-fixed-ifar, "
                                  "--single-fit-file")
        if any(sngl_optional_opts) and not all(sngl_opts):
            parser.error("Optional singles options "
                         f"({sngl_optional_opts_str}) given but no "
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
        return cls(
           ifo, newsnr_threshold=args.single_newsnr_threshold[ifo],
           reduced_chisq_threshold=args.single_reduced_chisq_threshold[ifo],
           duration_threshold=args.single_duration_threshold[ifo],
           fixed_ifar=args.single_fixed_ifar,
           fit_file=args.single_fit_file,
           sngl_ifar_est_dist=args.sngl_ifar_est_dist[ifo]
           )

    def check(self, trigs, data_reader):
        """ Look for a single detector trigger that passes the thresholds in
        the current data.
        """

        # Apply cuts to trigs before clustering
        # Cut on snr so that triggers which could not reach newsnr
        # threshold do not have newsnr calculated
        valid_idx = (trigs['template_duration'] >
                     self.thresholds['duration']) & \
                    (trigs['chisq'] <
                     self.thresholds['reduced_chisq']) & \
                    (trigs['snr'] >
                     self.thresholds['newsnr'])
        if not np.any(valid_idx):
            return None

        cutdurchi_trigs = {k: trigs[k][valid_idx] for k in trigs}

        # This uses the pycbc live convention of chisq always meaning the
        # reduced chisq.
        nsnr_all = ranking.newsnr(cutdurchi_trigs['snr'],
                                  cutdurchi_trigs['chisq'])
        nsnr_idx = nsnr_all > self.thresholds['newsnr']
        if not np.any(nsnr_idx):
            return None

        cutall_trigs = {k: cutdurchi_trigs[k][nsnr_idx]
                        for k in trigs}

        # 'cluster' by taking the maximal newsnr value over the trigger set
        i = nsnr_all[nsnr_idx].argmax()

        # calculate the (inverse) false-alarm rate
        nsnr = nsnr_all[nsnr_idx][i]
        dur = cutall_trigs['template_duration'][i]
        ifar = self.calculate_ifar(nsnr, dur)
        if ifar is None:
            return None

        # fill in a new candidate event
        candidate = {
            f'foreground/{self.ifo}/{k}': cutall_trigs[k][i] for k in trigs
        }
        candidate['foreground/stat'] = nsnr
        candidate['foreground/ifar'] = ifar
        candidate['HWINJ'] = data_reader.near_hwinj()
        return candidate

    def calculate_ifar(self, sngl_ranking, duration):
        logging.info("Calculating IFAR")
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
            logging.error(
                'Single fit file %s not found; '
                'dropping a potential single-detector candidate!',
                self.fit_file
            )
            return None

        bins = bin_utils.IrregularBins(bin_edges)
        dur_bin = bins[duration]

        rate = rates[dur_bin]
        coeff = coeffs[dur_bin]
        rate_louder = rate * fits.cum_fit('exponential', [sngl_ranking],
                                          coeff, thresh)[0]
        # apply a trials factor of the number of duration bins
        rate_louder *= len(rates)
        return conv.sec_to_year(1. / rate_louder)
