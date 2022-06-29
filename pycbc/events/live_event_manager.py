"""
Define the class used to manage events in pycbc_live
"""

import numpy, os, itertools, logging, lal, h5py, sys, subprocess
from mpi4py import MPI as mpi
from shutil import which
import pycbc
from pycbc import makedir, version
from pycbc.detector import ppdets
from pycbc.events.ranking import newsnr
from pycbc.filter import compute_followup_snr_series
from pycbc.filter import followup_event_significance
from pycbc.io.hdf import recursively_save_dict_contents_to_group
from pycbc.io.live import SingleCoincForGraceDB
from pycbc.waveform.waveform import props

def ptof(p, ft):
    """Convert p-value to FAR via foreground time `ft`.
    """
    return numpy.log1p(-p) / -ft

def ftop(f, ft):
    """Convert FAR to p-value via foreground time `ft`.
    """
    return 1 - numpy.exp(-ft * f)

def combine_ifar_pvalue(ifar, pvalue, livetime):
    """Convert original IFAR to p-value and combine with followup p-value.
    """
    from scipy.stats import combine_pvalues
    # NB units of ifar and livetime must be the same
    _, pnew = combine_pvalues([ftop(1. / ifar, livetime), pvalue])
    nifar = 1. / ptof(pnew, livetime)
    # take max of original IFAR and combined IFAR & apply trials factor
    return numpy.maximum(ifar, nifar) / 2.

class LiveEventManager(object):
    def __init__(self, output_path, mc_area_args, ifos,
                 ranking_statistic=None,
                 snr_opt_timeout=None,
                 channel_name=None,
                 processing_scheme=None,
                 use_date_prefix=False,
                 ifar_upload_threshold=None,
                 ifar_double_followup_threshold=None,
                 enable_single_detector_upload=None,
                 pval_livetime=None,
                 enable_gracedb_upload=False,
                 gracedb_testing=True,
                 gracedb_server=None,
                 gracedb_search=None,
                 run_snr_optimization=False):

        self.path = output_path
        self.ifos = ifos
        self.snr_opt_timeout = snr_opt_timeout

        self.ranking_statistic = ranking_statistic
        self.mc_area_args = mc_area_args
        self.channel_name = channel_name

        # Figure out what we are supposed to process within the pool of MPI processes
        self.comm = mpi.COMM_WORLD
        self.size = self.comm.Get_size()
        self.rank = self.comm.Get_rank()

        self.processing_scheme = processing_scheme

        self.use_date_prefix = use_date_prefix
        self.ifar_upload_threshold = ifar_upload_threshold
        self.ifar_double_followup_threshold = ifar_upload_threshold
        self.pvalue_livetime = pval_livetime

        self.gracedb_testing = gracedb_testing
        self.gracedb_server = gracedb_server
        self.gracedb_search = gracedb_search
        self.enable_gracedb_upload = enable_gracedb_upload
        self.enable_single_detector_upload = enable_single_detector_upload

        self.run_snr_optimization = run_snr_optimization

        # Keep track of which events have been uploaded
        self.last_few_coincs_uploaded = [0.0] * 10

        if self.run_snr_optimization and self.rank == 0:
            # preestimate the number of CPU cores that we can afford giving to
            # followup processes without slowing down the main search
            available_cores = len(os.sched_getaffinity(0))
            bg_cores = len(tuple(itertools.combinations(ifos, 2)))
            analysis_cores = 1 + bg_cores
            self.fu_cores = available_cores - analysis_cores
            if self.fu_cores <= 0:
                logging.warning('Insufficient number of CPU cores (%d) to '
                                'run search and trigger followups. Uploaded '
                                'triggers will momentarily increase the lag',
                                available_cores)
                self.fu_cores = 1


    def commit_results(self, results):
        logging.info('Committing triggers')
        self.comm.gather(results, root=0)

    def barrier(self):
        self.comm.Barrier()

    def barrier_status(self, status):
        return self.comm.allreduce(status, op=mpi.LAND)

    def gather_results(self):
        """ Collect results from the mpi subprocesses and collate them into
        contiguous sets of arrays.
        """

        if self.rank != 0:
            raise RuntimeError('Not root process')

        logging.info('Gathering triggers')
        all_results = self.comm.gather(None, root=0)
        data_ends = [a[1] for a in all_results if a is not None]
        results = [a[0] for a in all_results if a is not None]

        combined = {}
        for ifo in results[0]:
            # check if any of the results returned invalid
            try:
                for r in results:
                    if r[ifo] is False:
                        raise ValueError
            except ValueError:
                continue
            combined[ifo] = {}
            for key in results[0][ifo]:
                combined[ifo][key] = numpy.concatenate([r[ifo][key] for r in results])

        return combined, data_ends[0]

    def compute_followup_data(self, ifos, triggers, data_readers, bank,
                              followup_ifos=None, recalculate_ifar=False):
        """Figure out which of the followup detectors are usable, and compute
        SNR time series for all the available detectors.
        """
        out = {}
        followup_ifos = [] if followup_ifos is None else followup_ifos

        template_id = triggers['foreground/' + ifos[0] + '/template_id']
        htilde = bank[template_id]

        coinc_times = {ifo: triggers['foreground/' + ifo + '/end_time'] for ifo in ifos}

        # Get the SNR series for the ifos that made the initial coinc
        for ifo in ifos:
            # NOTE we only check the state/DQ of followup IFOs here.
            # IFOs producing the coincidence are assumed to also
            # produce valid SNR series.
            snr_series = compute_followup_snr_series(
                    data_readers[ifo], htilde, coinc_times[ifo],
                    check_state=False)

            if snr_series is not None:
                out[ifo] = {'snr_series': snr_series}

        # Determine if the other ifos can contribute to the coincident event
        for ifo in followup_ifos:
            snr_series, ptime, pvalue, sigma2 = followup_event_significance(
                    ifo, data_readers[ifo], bank, template_id, coinc_times)
            if snr_series is not None:
                out[ifo] = {'snr_series': snr_series}
                self.get_followup_info(ifos[0], ifo, triggers, snr_series,
                                       ptime, pvalue, sigma2,
                                       recalculate_ifar=recalculate_ifar)

        # the SNR time series sample rate can vary slightly due to
        # rounding errors, so force all of them to be identical
        fix_delta_t = None
        for ifo in out:
            if 'snr_series' not in out[ifo]:
                continue
            if fix_delta_t is None:
                fix_delta_t = out[ifo]['snr_series']._delta_t
            else:
                out[ifo]['snr_series']._delta_t = fix_delta_t

        return out

    def get_followup_info(self, coinc_ifo, ifo, triggers, snr_series, ptime,
                          pvalue, sigma2, recalculate_ifar=False):
        # Copy the common fields from the other detector
        # ignore fields that contain detector-specific data
        fields_to_ignore = set(['end_time', 'snr', 'stat', 'coa_phase',
                                'chisq', 'chisq_dof', 'sg_chisq', 'sigmasq'])
        for key in set(triggers):
            if 'foreground/{}/'.format(coinc_ifo) in key:
                _, _, name = key.split('/')
                if name in fields_to_ignore:
                    continue
                triggers['foreground/{}/{}'.format(ifo, name)] = triggers[key]

        # Set the detector-specific fields for which we have data
        snr_series_peak = snr_series.at_time(ptime, nearest_sample=True)
        base = 'foreground/{}/'.format(ifo)
        triggers[base + 'end_time'] = float(ptime)
        triggers[base + 'snr'] = triggers[base + 'stat'] = abs(snr_series_peak)
        triggers[base + 'coa_phase'] = numpy.angle(snr_series_peak)
        triggers[base + 'sigmasq'] = sigma2
        if recalculate_ifar:
            # Calculate new ifar
            triggers['foreground/ifar'] = combine_ifar_pvalue(
                    triggers['foreground/ifar'], pvalue, self.pvalue_livetime)
            triggers['foreground/pvalue_{}'.format(ifo)] = pvalue

    def check_coincs(self, ifos, coinc_results, psds, f_low,
                     data_readers, bank):
        """ Perform any followup and save zerolag triggers to a coinc xml file
        """
        #check for hardware injection
        for ifo in ifos:
            if data_readers[ifo].near_hwinj():
                coinc_results['HWINJ'] = True
                break

        if 'foreground/ifar' in coinc_results:
            logging.info('computing followup data for coinc')

            coinc_ifos = coinc_results['foreground/type'].split('-')
            followup_ifos = list(set(ifos) - set(coinc_ifos))

            double_ifar = coinc_results['foreground/ifar']
            if double_ifar < self.ifar_double_followup_threshold:
                coinc_results['foreground/NO_FOLLOWUP'] = True
                return

            fud = self.compute_followup_data(coinc_ifos, coinc_results,
                                             data_readers, bank,
                                             followup_ifos=followup_ifos,
                                             recalculate_ifar=True)

            live_ifos = [ifo for ifo in fud if 'snr_series' in fud[ifo]]

            event = SingleCoincForGraceDB(live_ifos, coinc_results, bank=bank,
                                          psds=psds, followup_data=fud,
                                          low_frequency_cutoff=f_low,
                                          channel_names=self.channel_name,
                                          mc_area_args=self.mc_area_args,
                                          )

            fname = 'coinc-{}-{}.xml.gz'.format(event.merger_time,
                                                pycbc.random_string(6))
            fname = os.path.join(self.path, fname)
            logging.info('Coincident candidate! Saving as %s', fname)

            # verbally explain some details not obvious from the other info
            comment = ('Trigger produced as a {} coincidence. '
                       'FAR is based on all listed detectors.<br />'
                       'Two-detector ranking statistic: {}<br />'
                       'Followup detectors: {}')
            comment = comment.format(ppdets(coinc_ifos),
                                     self.ranking_statistic,
                                     ppdets(followup_ifos))

            ifar = coinc_results['foreground/ifar']
            if self.enable_gracedb_upload and self.ifar_upload_threshold < ifar:
                self.gid = event.upload(fname, gracedb_server=self.gracedb_server,
                                   testing=self.gracedb_testing,
                                   extra_strings=[comment],
                                   search=self.gracedb_search)
                # Keep track of the last few coincs uploaded in order to
                # prevent singles being uploaded as well for coinc events
                self.last_few_coincs_uploaded.append(event.merger_time)
                # Only need to keep a few (10) events
                self.last_few_coincs_uploaded = \
                    self.last_few_coincs_uploaded[-10:]
            else:
                self.gid = None
                event.save(fname)

            if self.run_snr_optimization and self.ifar_upload_threshold < ifar:
                self.setup_snr_optimization(coinc_ifos, coinc_results, bank, fname, data_readers, upload=True)


    def check_singles(self, results, bank, data_reader, psds,
                      sngl_estimator, f_low):
        active = [k for k in results if results[k] != None]
        for ifo in active:
            single = sngl_estimator[ifo].check(results[ifo], data_reader[ifo])

            if single is None:
                continue

            followup_ifos = [i for i in active if i is not ifo]
            fud = self.compute_followup_data([ifo], single, data_reader,
                                             bank, followup_ifos=followup_ifos,
                                             recalculate_ifar=False)
            ifar = single['foreground/ifar']
            # apply a trials factor of the number of active detectors
            ifar /= len(active)
            single['foreground/ifar'] = ifar

            live_ifos = [i for i in fud if 'snr_series' in fud[i]]

            event = SingleCoincForGraceDB(live_ifos, single, bank=bank,
                                          psds=psds, followup_data=fud,
                                          low_frequency_cutoff=f_low,
                                          channel_names=self.channel_name,
                                          mc_area_args=self.mc_area_args,
                                          )

            fname = 'single-%s-%s.xml.gz' % (ifo, event.merger_time)
            fname = os.path.join(self.path, fname)
            logging.info('Single-detector candidate! Saving as %s', fname)

            # verbally explain some details not obvious from the other info
            comment = ('Trigger produced as a {0} single. '
                       'FAR is based on {0} only.<br />'
                       'Followup detectors: {1}')
            comment = comment.format(ifo, ppdets(followup_ifos))

            # Has a coinc event at this time been uploaded recently?
            # If so, skip upload
            coinc_tdiffs =  abs(event.merger_time -
                                self.last_few_coincs_uploaded)
            nearby_coincs = coinc_tdiffs < 0.1
            if any(nearby_coincs):
                logging.info("Single detector event at time %.2f not "
                             "uploaded due to coinc event at %.2f.",
                             event.merger_time,
                             self.last_few_coincs_uploaded[nearby_coincs])

            if self.enable_single_detector_upload \
                    and self.ifar_upload_threshold < ifar \
                    and not any(nearby_coincs):
                event.upload(fname, gracedb_server=self.gracedb_server,
                             testing=self.gracedb_testing,
                             extra_strings=[comment],
                             search=self.gracedb_search)
                gdb_upload_opt = True
            else:
                event.save(fname)
                gdb_upload_opt = False

            if self.run_snr_optimization \
                    and self.ifar_upload_threshold < ifar:
                self.setup_snr_optimization([ifo], single, bank, fname,
                                            data_reader, upload=gdb_upload_opt)


    def setup_snr_optimization(self, ifos, results, bank, fname,
                               data_readers, upload=False):
        """
        Set up and run pycbc_optimize_snr for a found result

        Parameters
        ----------
        ifos: list
            The ifos involved in the event to be used for the optimization

        results: dictionary
            Dictionary with keys describing the parameters of the event

        bank: LiveFilterBank
            Template bank as used in the search

        fname: string
            The name of the xml file containing event parameters. Must end
            in 'xml.gz', this will be the basis of the output file names.

        data_readers: dictionary of StrainBuffer
            IFO-keyed dictionary containing the StrainBuffer for each detector

        upload: Boolean, default False
            Flag whether to supply the --enable-gracedb-upload option or not
            to pycbc_optimize_snr
        """
        template_id = \
            results['foreground/{}/template_id'.format(ifos[0])]
        p = props(bank.table[template_id])
        apr = p.pop('approximant')
        min_buffer = bank.minimum_buffer + 0.5
        buff_size = \
            pycbc.waveform.get_waveform_filter_length_in_time(apr, **p)

        tlen = bank.round_up((buff_size + min_buffer) \
           * bank.sample_rate)
        flen = int(tlen / 2 + 1)
        delta_f = bank.sample_rate / float(tlen)
        cmd = 'timeout {} '.format(self.snr_opt_timeout)
        exepath = which('pycbc_optimize_snr')
        cmd += exepath + ' '
        data_fils_str = '--data-files '
        psd_fils_str = '--psd-files '
        for ifo in ifos:
            curr_fname = \
                fname.replace('.xml.gz',
                              '_{}_data_overwhitened.hdf'.format(ifo))
            curr_data = data_readers[ifo].overwhitened_data(delta_f)
            curr_data.save(curr_fname)
            data_fils_str += '{}:{} ' .format(ifo, curr_fname)
            curr_fname = fname.replace('.xml.gz',
                                       '_{}_psd.hdf'.format(ifo))
            curr_psd = curr_data.psd
            curr_psd.save(curr_fname)
            psd_fils_str += '{}:{} ' .format(ifo, curr_fname)
        cmd += data_fils_str
        cmd += psd_fils_str

        curr_fname = fname.replace('.xml.gz', '_attributes.hdf')
        with h5py.File(curr_fname, 'w') as hdfp:
            for ifo in ifos:
                curr_time = \
                    results['foreground/{}/end_time'.format(ifo)]
                hdfp['coinc_times/{}'.format(ifo)] = curr_time
            f_end = bank.end_frequency(template_id)
            if f_end is None or f_end >= (flen * delta_f):
                f_end = (flen-1) * delta_f
            hdfp['flen'] = flen
            hdfp['delta_f'] = delta_f
            hdfp['f_end'] = f_end
            hdfp['sample_rate'] = bank.sample_rate
            hdfp['flow'] = bank.table[template_id].f_lower
            hdfp['mass1'] = bank.table[template_id].mass1
            hdfp['mass2'] = bank.table[template_id].mass2
            hdfp['spin1z'] = bank.table[template_id].spin1z
            hdfp['spin2z'] = bank.table[template_id].spin2z
            hdfp['template_duration'] = \
                bank.table[template_id].template_duration
            hdfp['ifar'] = results['foreground/ifar']
            if self.gid is not None:
                hdfp['gid'] = self.gid

            for ifo in self.channel_name:
                hdfp['channel_names/{}'.format(ifo)] = \
                    self.channel_name[ifo]

            recursively_save_dict_contents_to_group(hdfp,
                                                    'mc_area_args/',
                                                    self.mc_area_args)

        cmd += '--params-file {} '.format(curr_fname)
        cmd += '--approximant {} '.format(apr)
        cmd += '--gracedb-server {} '.format(self.gracedb_server)
        cmd += '--gracedb-search {} '.format(self.gracedb_search)
        if not self.gracedb_testing:
            cmd += '--production '
        cmd += '--verbose '
        cmd += '--output-path {} '.format(self.path)
        if self.enable_gracedb_upload and upload:
            cmd += '--enable-gracedb-upload '

        cmd += '--cores {} '.format(self.fu_cores)
        if self.processing_scheme:
            # we will use the cores for multiple workers of the
            # optimization routine, so we force the processing scheme
            # to a single core here.  This may be enforcing some
            # assumptions about the optimal way to do FFTs on the
            # machine. However, the dominant cost of pycbc_optimize_snr
            # is expected to be in waveform generation, which is
            # unlikely to benefit from a processing scheme with more
            # than 1 thread anyway.
            opt_scheme = self.processing_scheme.split(':')[0]
            cmd += '--processing-scheme {}:1 '.format(opt_scheme)

        log_fname = fname.replace('.xml.gz', '_optimize_snr.log')

        logging.info('Running %s with log to %s', cmd, log_fname)
        with open(log_fname, "w") as logfile:
            subprocess.Popen(
                cmd, shell=True, stdout=logfile, stderr=logfile
            )


    def dump(self, results, name, store_psd=False, time_index=None,
             store_loudest_index=False, raw_results=None, gates=None):
        """Save the results from this time block to an hdf output file.
        """
        if self.use_date_prefix:
            tm = lal.GPSToUTC(int(time_index))
            subdir = '{:04d}_{:02d}_{:02d}'.format(tm[0], tm[1], tm[2])
            makedir(os.path.join(self.path, subdir))
            fname = os.path.join(self.path, subdir, name) + '.hdf'
        else:
            makedir(self.path)
            fname = os.path.join(self.path, name) + '.hdf'

        with h5py.File(fname, 'w') as f:
            f.attrs['pycbc_version'] = version.git_verbose_msg
            f.attrs['command_line'] = sys.argv
            f.attrs['num_live_detectors'] = len(self.live_detectors)

            def h5py_unicode_workaround(stuff):
                # workaround for the fact that h5py cannot handle Unicode
                # numpy arrays that show up with Python 3
                if hasattr(stuff, 'dtype') and stuff.dtype.kind == 'U':
                    return [s.encode() for s in stuff]
                return stuff

            for ifo in results:
                for k in results[ifo]:
                    f['%s/%s' % (ifo, k)] = \
                            h5py_unicode_workaround(results[ifo][k])

            for key in raw_results:
                f[key] = h5py_unicode_workaround(raw_results[key])

            if store_loudest_index:
                for ifo in results:
                    if 'snr' in results[ifo]:
                        s = numpy.array(results[ifo]['snr'], ndmin=1)
                        c = numpy.array(results[ifo]['chisq'], ndmin=1)
                        nsnr = numpy.array(newsnr(s, c), ndmin=1) if len(s) > 0 else []

                        # loudest by newsnr
                        nloudest = numpy.argsort(nsnr)[::-1][0:store_loudest_index]

                        # loudest by snr
                        sloudest = numpy.argsort(s)[::-1][0:store_loudest_index]
                        f[ifo]['loudest'] = numpy.union1d(nloudest, sloudest)

            for ifo in (gates or {}):
                gate_dtype = [('center_time', float),
                              ('zero_half_width', float),
                              ('taper_width', float)]
                f['{}/gates'.format(ifo)] = \
                        numpy.array(gates[ifo], dtype=gate_dtype)

        for ifo in (store_psd or {}):
            if store_psd[ifo] is not None:
                store_psd[ifo].save(fname, group='%s/psd' % ifo)
