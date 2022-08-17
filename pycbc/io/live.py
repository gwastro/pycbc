import logging
import os
import pycbc
import numpy
import lal
import json
from ligo.lw import ligolw
from ligo.lw import lsctables
from ligo.lw import utils as ligolw_utils
from pycbc import version as pycbc_version
from pycbc import pnutils
from pycbc.io.ligolw import (
    return_empty_sngl, create_process_table, make_psd_xmldoc, snr_series_to_xml
)
from pycbc.results import generate_asd_plot
from pycbc.results import ifo_color
from pycbc.results import source_color
from pycbc.mchirp_area import calc_probabilities


class CandidateForGraceDB(object):
    """This class provides an interface for uploading candidates to GraceDB.
    """

    def __init__(self, ifos, coinc_results, **kwargs):
        """Initialize a representation of a zerolag candidate for upload to
        GraceDB.

        Parameters
        ----------
        ifos: list of strs
            A list of the ifos participating in this candidate.
        coinc_results: dict of values
            A dictionary of values. The format is defined in
            `pycbc/events/coinc.py` and matches the on-disk representation in
            the hdf file for this time.
        psds: dict of FrequencySeries
            Dictionary providing PSD estimates for all involved detectors.
        low_frequency_cutoff: float
            Minimum valid frequency for the PSD estimates.
        high_frequency_cutoff: float, optional
            Maximum frequency considered for the PSD estimates. Default None.
        followup_data: dict of dicts, optional
            Dictionary providing SNR time series for each detector, to be used
            in sky localization with BAYESTAR. The format should be
            `followup_data['H1']['snr_series']`. More detectors can be present
            than given in `ifos`. If so, the extra detectors will only be used
            for sky localization.
        channel_names: dict of strings, optional
            Strain channel names for each detector. Will be recorded in the
            `sngl_inspiral` table.
        padata: PAstroData instance
            Organizes info relevant to the astrophysical probability of the
            candidate.
        mc_area_args: dict of dicts, optional
            Dictionary providing arguments to be used in source probability
            estimation with `pycbc/mchirp_area.py`.
        """
        self.template_id = coinc_results[f'foreground/{ifos[0]}/template_id']
        self.coinc_results = coinc_results
        self.ifos = ifos
        self.basename = None

        # remember if this should be marked as HWINJ
        self.is_hardware_injection = ('HWINJ' in coinc_results
                                      and coinc_results['HWINJ'])

        # Check if we need to apply a time offset (this may be permerger)
        self.time_offset = 0
        rtoff = f'foreground/{ifos[0]}/time_offset'
        if rtoff in coinc_results:
            self.time_offset = coinc_results[rtoff]

        if 'followup_data' in kwargs:
            fud = kwargs['followup_data']
            assert len({fud[ifo]['snr_series'].delta_t for ifo in fud}) == 1, \
                    "delta_t for all ifos do not match"
            self.snr_series = {ifo: fud[ifo]['snr_series'] for ifo in fud}
            usable_ifos = fud.keys()
            followup_ifos = list(set(usable_ifos) - set(ifos))

            for ifo in self.snr_series:
                self.snr_series[ifo].start_time += self.time_offset
        else:
            self.snr_series = None
            usable_ifos = ifos
            followup_ifos = []

        # Set up the bare structure of the xml document
        outdoc = ligolw.Document()
        outdoc.appendChild(ligolw.LIGO_LW())

        # FIXME is it safe (in terms of downstream operations) to let
        # `program_name` default to the actual script name?
        proc_id = create_process_table(outdoc, program_name='pycbc',
                                       detectors=usable_ifos).process_id

        # Set up coinc_definer table
        coinc_def_table = lsctables.New(lsctables.CoincDefTable)
        coinc_def_id = lsctables.CoincDefID(0)
        coinc_def_row = lsctables.CoincDef()
        coinc_def_row.search = "inspiral"
        coinc_def_row.description = "sngl_inspiral<-->sngl_inspiral coincs"
        coinc_def_row.coinc_def_id = coinc_def_id
        coinc_def_row.search_coinc_type = 0
        coinc_def_table.append(coinc_def_row)
        outdoc.childNodes[0].appendChild(coinc_def_table)

        # Set up coinc inspiral and coinc event tables
        coinc_id = lsctables.CoincID(0)
        coinc_event_table = lsctables.New(lsctables.CoincTable)
        coinc_event_row = lsctables.Coinc()
        coinc_event_row.coinc_def_id = coinc_def_id
        coinc_event_row.nevents = len(usable_ifos)
        coinc_event_row.instruments = ','.join(usable_ifos)
        coinc_event_row.time_slide_id = lsctables.TimeSlideID(0)
        coinc_event_row.process_id = proc_id
        coinc_event_row.coinc_event_id = coinc_id
        coinc_event_row.likelihood = 0.
        coinc_event_table.append(coinc_event_row)
        outdoc.childNodes[0].appendChild(coinc_event_table)

        # Set up sngls
        sngl_inspiral_table = lsctables.New(lsctables.SnglInspiralTable)
        coinc_event_map_table = lsctables.New(lsctables.CoincMapTable)

        sngl_populated = None
        network_snrsq = 0
        for sngl_id, ifo in enumerate(usable_ifos):
            sngl = return_empty_sngl(nones=True)
            sngl.event_id = lsctables.SnglInspiralID(sngl_id)
            sngl.process_id = proc_id
            sngl.ifo = ifo
            names = [n.split('/')[-1] for n in coinc_results
                     if f'foreground/{ifo}' in n]
            for name in names:
                val = coinc_results[f'foreground/{ifo}/{name}']
                if name == 'end_time':
                    val += self.time_offset
                    sngl.end = lal.LIGOTimeGPS(val)
                else:
                    try:
                        setattr(sngl, name, val)
                    except AttributeError:
                        pass
            if sngl.mass1 and sngl.mass2:
                sngl.mtotal, sngl.eta = pnutils.mass1_mass2_to_mtotal_eta(
                        sngl.mass1, sngl.mass2)
                sngl.mchirp, _ = pnutils.mass1_mass2_to_mchirp_eta(
                        sngl.mass1, sngl.mass2)
                sngl_populated = sngl
            if sngl.snr:
                sngl.eff_distance = sngl.sigmasq ** 0.5 / sngl.snr
                network_snrsq += sngl.snr ** 2.0
            if 'channel_names' in kwargs and ifo in kwargs['channel_names']:
                sngl.channel = kwargs['channel_names'][ifo]
            sngl_inspiral_table.append(sngl)

            # Set up coinc_map entry
            coinc_map_row = lsctables.CoincMap()
            coinc_map_row.table_name = 'sngl_inspiral'
            coinc_map_row.coinc_event_id = coinc_id
            coinc_map_row.event_id = sngl.event_id
            coinc_event_map_table.append(coinc_map_row)

            if self.snr_series is not None:
                snr_series_to_xml(self.snr_series[ifo], outdoc, sngl.event_id)

        # Set merger time to the average of the ifo peaks
        self.merger_time = numpy.mean(
                    [coinc_results[f'foreground/{ifo}/end_time']
                     for ifo in ifos]) + self.time_offset

        # For subthreshold detectors, respect BAYESTAR's assumptions and checks
        bayestar_check_fields = ('mass1 mass2 mtotal mchirp eta spin1x '
                                 'spin1y spin1z spin2x spin2y spin2z').split()
        for sngl in sngl_inspiral_table:
            if sngl.ifo in followup_ifos:
                for bcf in bayestar_check_fields:
                    setattr(sngl, bcf, getattr(sngl_populated, bcf))
                sngl.end = lal.LIGOTimeGPS(self.merger_time)

        outdoc.childNodes[0].appendChild(coinc_event_map_table)
        outdoc.childNodes[0].appendChild(sngl_inspiral_table)

        # Set up the coinc inspiral table
        coinc_inspiral_table = lsctables.New(lsctables.CoincInspiralTable)
        coinc_inspiral_row = lsctables.CoincInspiral()
        # This seems to be used as FAP, which should not be in gracedb
        coinc_inspiral_row.false_alarm_rate = 0
        coinc_inspiral_row.minimum_duration = 0.
        coinc_inspiral_row.instruments = tuple(usable_ifos)
        coinc_inspiral_row.coinc_event_id = coinc_id
        coinc_inspiral_row.mchirp = sngl_populated.mchirp
        coinc_inspiral_row.mass = sngl_populated.mtotal
        coinc_inspiral_row.end_time = sngl_populated.end_time
        coinc_inspiral_row.end_time_ns = sngl_populated.end_time_ns
        coinc_inspiral_row.snr = network_snrsq ** 0.5
        far = 1.0 / (lal.YRJUL_SI * coinc_results['foreground/ifar'])
        coinc_inspiral_row.combined_far = far
        coinc_inspiral_table.append(coinc_inspiral_row)
        outdoc.childNodes[0].appendChild(coinc_inspiral_table)

        # Append the PSDs
        self.psds = kwargs['psds']
        psds_lal = {}
        for ifo in self.psds:
            psd = self.psds[ifo]
            kmin = int(kwargs['low_frequency_cutoff'] / psd.delta_f)
            fseries = lal.CreateREAL8FrequencySeries(
                "psd", psd.epoch, kwargs['low_frequency_cutoff'], psd.delta_f,
                lal.StrainUnit**2 / lal.HertzUnit, len(psd) - kmin)
            fseries.data.data = psd.numpy()[kmin:] / pycbc.DYN_RANGE_FAC ** 2.0
            psds_lal[ifo] = fseries
        make_psd_xmldoc(psds_lal, outdoc)

        # P astro calculation
        if 'padata' in kwargs:
            padata = kwargs['padata']
            trigger_data = {
                'mass1': sngl_populated.mass1,
                'mass2': sngl_populated.mass2,
                'spin1z': sngl_populated.spin1z,
                'spin2z': sngl_populated.spin2z,
                'network_snr': network_snrsq ** 0.5,
                'far': far}
            horizons = {ifo: self.psds[ifo].dist for ifo in self.psds}
            self.p_astro, self.p_terr = \
                                  padata.do_pastro_calc(trigger_data, horizons)
        else:
            self.p_astro, self.p_terr = None, None

        # Source probabilities estimation
        if 'mc_area_args' in kwargs:
            eff_distances = [sngl.eff_distance for sngl in sngl_inspiral_table]
            self.probabilities = calc_probabilities(coinc_inspiral_row.mchirp,
                                                    coinc_inspiral_row.snr,
                                                    min(eff_distances),
                                                    kwargs['mc_area_args'])
        else:
            self.probabilities = None

        # Combine p astro and source probs
        if self.p_astro is not None and self.probabilities is not None:
            self.astro_probs = {cl: pr * self.p_astro for
                                cl, pr in self.probabilities.items()}
            self.astro_probs['p_terr'] = self.p_terr
        else:
            self.astro_probs = None

        self.outdoc = outdoc
        self.time = sngl_populated.end

    def save(self, fname):
        """Write a file representing this candidate in a LIGOLW XML format
        compatible with GraceDB.

        Parameters
        ----------
        fname: str
            Name of file to write to disk.
        """
        ligolw_utils.write_filename(self.outdoc, fname, compress='auto')

        if self.basename is None:
            # here assume compression
            self.basename = fname.replace('.xml.gz', '')

        # Save multi-cpt p astro as json
        if self.astro_probs is not None:
            self.multipa_file = self.basename + '_p_astro.json'
            with open(self.multipa_file, 'w') as multipaf:
                json.dump(self.astro_probs, multipaf)
            logging.info('Multi p_astro file saved as %s', self.multipa_file)
            # Don't save any other files!
            return

        # Save source probabilities in a json file
        if self.probabilities is not None:
            self.prob_file = self.basename + '_probs.json'
            with open(self.prob_file, 'w') as probf:
                json.dump(self.probabilities, probf)
            logging.info('Source probabilities file saved as %s', self.prob_file)
            return

        # Save p astro / p terr as json
        if self.p_astro is not None:
            self.pastro_file = self.basename + '_pa_pterr.json'
            with open(self.pastro_file, 'w') as pastrof:
                json.dump({'p_astro': self.p_astro, 'p_terr': self.p_terr},
                          pastrof)
            logging.info('P_astro file saved as %s', self.pastro_file)
        return

    def upload(self, fname, gracedb_server=None, testing=True,
               extra_strings=None, search='AllSky'):
        """Upload this candidate to GraceDB, and annotate it with a few useful
        plots and comments.

        Parameters
        ----------
        fname: str
            The name to give the xml file associated with this trigger
        gracedb_server: string, optional
            URL to the GraceDB web API service for uploading the event.
            If omitted, the default will be used.
        testing: bool
            Switch to determine if the upload should be sent to gracedb as a
            test trigger (True) or a production trigger (False).
        search: str
            String going into the "search" field of the GraceDB event.
        """
        from ligo.gracedb.rest import GraceDb
        import matplotlib
        matplotlib.use('Agg')
        import pylab as pl

        if fname.endswith('.xml.gz'):
            self.basename = fname.replace('.xml.gz', '')
        elif fname.endswith('.xml'):
            self.basename = fname.replace('.xml', '')
        else:
            raise ValueError("Upload filename must end in .xml or .xml.gz, got"
                             " %s" % fname)

        # First make sure the event is saved on disk
        # as GraceDB operations can fail later
        self.save(fname)

        gid = None
        try:
            # try connecting to GraceDB
            gracedb = GraceDb(gracedb_server) \
                    if gracedb_server is not None else GraceDb()

            # create GraceDB event
            group = 'Test' if testing else 'CBC'
            r = gracedb.createEvent(group, "pycbc", fname, search).json()
            gid = r["graceid"]
            logging.info("Uploaded event %s", gid)

            if self.is_hardware_injection:
                gracedb.writeLabel(gid, 'INJ')
                logging.info("Tagging event %s as an injection", gid)

            # add info for tracking code version
            gracedb_tag_with_version(gracedb, gid)

            extra_strings = [] if extra_strings is None else extra_strings
            for text in extra_strings:
                gracedb.writeLog(gid, text, tag_name=['analyst_comments'])
        except Exception as exc:
            logging.error('Something failed during the upload/annotation of '
                          'event %s on GraceDB. The event may not have been '
                          'uploaded!', fname)
            logging.error(str(exc))

        # plot the SNR timeseries and noise PSDs
        if self.snr_series is not None:
            snr_series_fname = self.basename + '.hdf'
            snr_series_plot_fname = self.basename + '_snr.png'
            asd_series_plot_fname = self.basename + '_asd.png'

            pl.figure()
            ref_time = int(self.merger_time)
            for ifo in sorted(self.snr_series):
                curr_snrs = self.snr_series[ifo]
                curr_snrs.save(snr_series_fname, group='%s/snr' % ifo)
                pl.plot(curr_snrs.sample_times - ref_time, abs(curr_snrs),
                        c=ifo_color(ifo), label=ifo)
                if ifo in self.ifos:
                    base = 'foreground/{}/'.format(ifo)
                    snr = self.coinc_results[base + 'snr']
                    mt = (self.coinc_results[base + 'end_time']
                          + self.time_offset)
                    pl.plot([mt - ref_time], [snr], c=ifo_color(ifo),
                            marker='x')
            pl.legend()
            pl.xlabel('GPS time from {:d} (s)'.format(ref_time))
            pl.ylabel('SNR')
            pl.savefig(snr_series_plot_fname)
            pl.close()

            generate_asd_plot(self.psds, asd_series_plot_fname)

            # Additionally save the PSDs into the snr_series file
            for ifo in sorted(self.psds):
                # Undo dynamic range factor
                curr_psd = self.psds[ifo].astype(numpy.float64)
                curr_psd /= pycbc.DYN_RANGE_FAC ** 2.0
                curr_psd.save(snr_series_fname, group='%s/psd' % ifo)

        # Only make the pie plot if no multi cpt p_astro
        if self.probabilities is not None and self.astro_probs is None:
            self.prob_plotf = self.prob_file.replace('.json', '.png')
            # Don't try to plot zero probabilities
            prob_plot = {k: v for (k, v) in self.probabilities.items()
                         if v != 0.0}
            labels, sizes = zip(*prob_plot.items())
            colors = [source_color(label) for label in labels]
            fig, ax = pl.subplots()
            ax.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%',
                   textprops={'fontsize': 15})
            ax.axis('equal')
            fig.savefig(self.prob_plotf)
            pl.close()

        if gid is None:
            # Don't try to do anything else!
            return gid

        # Upload SNR series in HDF format and plots
        if self.snr_series is not None:
            try:
                gracedb.writeLog(
                    gid, 'SNR timeseries HDF file upload',
                    filename=snr_series_fname
                )
                gracedb.writeLog(
                    gid, 'SNR timeseries plot upload',
                    filename=snr_series_plot_fname,
                    tag_name=['background'],
                    displayName=['SNR timeseries']
                )
                gracedb.writeLog(
                    gid, 'ASD plot upload',
                    filename=asd_series_plot_fname,
                    tag_name=['psd'], displayName=['ASDs']
                )
            except Exception as exc:
                logging.error('Failed to upload SNR timeseries and ASD for %s',
                               gid)
                logging.error(str(exc))

        # Upload multi-cpt p_astro JSON
        if self.astro_probs is not None:
            try:
                gracedb.writeLog(
                    gid, 'Multi-component p_astro JSON file upload',
                    filename=self.multipa_file,
                    tag_name=['em_follow']
                )
                logging.info('Uploaded multi p_astro for %s', gid)
            except Exception as exc:
                logging.error('Failed to upload multi p_astro file for %s', gid)
                logging.error(str(exc))
            # Don't do anything else!
            return gid

        # If there is no multi p_astro, upload source probabilities in JSON
        # format and plot
        if self.probabilities is not None:
            try:
                gracedb.writeLog(
                    gid, 'Source probabilities JSON file upload',
                    filename=self.prob_file,
                    tag_name=['pe']
                )
                logging.info('Uploaded source probabilities for %s', gid)
                gracedb.writeLog(
                    gid, 'Source probabilities plot upload',
                    filename=self.prob_plotf,
                    tag_name=['pe']
                )
                logging.info('Uploaded source probabilities pie chart for %s',
                             gid)
            except Exception as exc:
                logging.error(
                    'Failed to upload source probability results for %s', gid)
                logging.error(str(exc))

        # If there is p_astro but no probabilities, upload p_astro JSON
        if self.p_astro is not None:
            try:
                gracedb.writeLog(
                    gid, '2-component p_astro JSON file upload',
                    filename=self.pastro_file,
                    tag_name=['sig_info']
                )
                logging.info('Uploaded p_astro for %s', gid)
            except Exception as exc:
                logging.error('Failed to upload p_astro file for %s', gid)
                logging.error(str(exc))

        return gid


def gracedb_tag_with_version(gracedb, event_id):
    """Add a GraceDB log entry reporting PyCBC's version and install location.
    """
    version_str = 'Using PyCBC version {}{} at {}'
    version_str = version_str.format(
            pycbc_version.version,
            ' (release)' if pycbc_version.release else '',
            os.path.dirname(pycbc.__file__))
    gracedb.writeLog(event_id, version_str)


__all__ = ['CandidateForGraceDB', 'gracedb_tag_with_version']
