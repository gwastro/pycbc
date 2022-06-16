import logging
import os
import pycbc
import numpy
import lal
import json
from ligo.lw import ligolw
from ligo.lw import lsctables
from ligo.lw import utils as ligolw_utils
from ligo.lw.param import Param as LIGOLWParam
from ligo.lw.array import Array as LIGOLWArray
from pycbc import version as pycbc_version
from pycbc import pnutils
from pycbc.io.ligolw import return_empty_sngl, create_process_table
from pycbc.results import generate_asd_plot
from pycbc.results import ifo_color
from pycbc.results import source_color
from pycbc.mchirp_area import calc_probabilities


#FIXME Legacy build PSD xml helpers, delete me when we move away entirely from
# xml formats
def _build_series(series, dim_names, comment, delta_name, delta_unit):
    Attributes = ligolw.sax.xmlreader.AttributesImpl
    elem = ligolw.LIGO_LW(
            Attributes({'Name': str(series.__class__.__name__)}))
    if comment is not None:
        elem.appendChild(ligolw.Comment()).pcdata = comment
    elem.appendChild(ligolw.Time.from_gps(series.epoch, 'epoch'))
    elem.appendChild(LIGOLWParam.from_pyvalue('f0', series.f0, unit='s^-1'))
    delta = getattr(series, delta_name)
    if numpy.iscomplexobj(series.data.data):
        data = numpy.row_stack((numpy.arange(len(series.data.data)) * delta,
                             series.data.data.real, series.data.data.imag))
    else:
        data = numpy.row_stack((numpy.arange(len(series.data.data)) * delta,
                                series.data.data))
    a = LIGOLWArray.build(series.name, data, dim_names=dim_names)
    a.Unit = str(series.sampleUnits)
    dim0 = a.getElementsByTagName(ligolw.Dim.tagName)[0]
    dim0.Unit = delta_unit
    dim0.Start = series.f0
    dim0.Scale = delta
    elem.appendChild(a)
    return elem

def snr_series_to_xml(snr_series, document, sngl_inspiral_id):
    """Save an SNR time series into an XML document, in a format compatible
    with BAYESTAR.
    """
    snr_lal = snr_series.lal()
    snr_lal.name = 'snr'
    snr_lal.sampleUnits = ''
    snr_xml = _build_series(snr_lal, ('Time', 'Time,Real,Imaginary'), None,
                            'deltaT', 's')
    snr_node = document.childNodes[-1].appendChild(snr_xml)
    eid_param = LIGOLWParam.from_pyvalue('event_id', sngl_inspiral_id)
    snr_node.appendChild(eid_param)

def make_psd_xmldoc(psddict, xmldoc=None):
    """Add a set of PSDs to a LIGOLW XML document. If the document is not
    given, a new one is created first.
    """
    xmldoc = ligolw.Document() if xmldoc is None else xmldoc.childNodes[0]

    # the PSDs must be children of a LIGO_LW with name "psd"
    root_name = 'psd'
    Attributes = ligolw.sax.xmlreader.AttributesImpl
    lw = xmldoc.appendChild(
        ligolw.LIGO_LW(Attributes({'Name': root_name})))

    for instrument, psd in psddict.items():
        xmlseries = _build_series(psd, ('Frequency,Real', 'Frequency'),
                                  None, 'deltaF', 's^-1')
        fs = lw.appendChild(xmlseries)
        fs.appendChild(LIGOLWParam.from_pyvalue('instrument', instrument))
    return xmldoc

class SingleCoincForGraceDB(object):
    """Create xml files and submit them to gracedb from PyCBC Live"""
    def __init__(self, ifos, coinc_results, **kwargs):
        """Initialize a ligolw xml representation of a zerolag trigger
        for upload from pycbc live to gracedb.

        Parameters
        ----------
        ifos: list of strs
            A list of the ifos participating in this trigger.
        coinc_results: dict of values
            A dictionary of values. The format is defined in
            pycbc/events/coinc.py and matches the on disk representation
            in the hdf file for this time.
        psds: dict of FrequencySeries
            Dictionary providing PSD estimates for all involved detectors.
        low_frequency_cutoff: float
            Minimum valid frequency for the PSD estimates.
        high_frequency_cutoff: float, optional
            Maximum frequency considered for the PSD estimates. Default None.
        followup_data: dict of dicts, optional
            Dictionary providing SNR time series for each detector,
            to be used in sky localization with BAYESTAR. The format should
            be `followup_data['H1']['snr_series']`. More detectors can be
            present than given in `ifos`. If so, the extra detectors will only
            be used for sky localization.
        channel_names: dict of strings, optional
            Strain channel names for each detector.
            Will be recorded in the sngl_inspiral table.
        mc_area_args: dict of dicts, optional
            Dictionary providing arguments to be used in source probability
            estimation with pycbc/mchirp_area.py
        """
        self.template_id = coinc_results['foreground/%s/template_id' % ifos[0]]
        self.coinc_results = coinc_results
        self.ifos = ifos

        # remember if this should be marked as HWINJ
        self.is_hardware_injection = ('HWINJ' in coinc_results
                                      and coinc_results['HWINJ'])

        # Check if we need to apply a time offset (this may be permerger)
        self.time_offset = 0
        rtoff = 'foreground/{}/time_offset'.format(ifos[0])
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
                     if 'foreground/%s' % ifo in n]
            for name in names:
                val = coinc_results['foreground/%s/%s' % (ifo, name)]
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
                sngl.eff_distance = (sngl.sigmasq)**0.5 / sngl.snr
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

        # set merger time to the average of the ifo peaks
        self.merger_time = numpy.mean(
                    [coinc_results['foreground/{}/end_time'.format(ifo)]
                     for ifo in ifos]) + self.time_offset

        # for subthreshold detectors, respect BAYESTAR's assumptions and checks
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

        # append the PSDs
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

        # source probabilities estimation
        if 'mc_area_args' in kwargs:
            eff_distances = [sngl.eff_distance for sngl in sngl_inspiral_table]
            probabilities = calc_probabilities(coinc_inspiral_row.mchirp,
                                               coinc_inspiral_row.snr,
                                               min(eff_distances),
                                               kwargs['mc_area_args'])
            self.probabilities = probabilities
        else:
            self.probabilities = None

        self.outdoc = outdoc
        self.time = sngl_populated.end

    def save(self, filename):
        """Write this trigger to gracedb compatible xml format

        Parameters
        ----------
        filename: str
            Name of file to write to disk.
        """
        ligolw_utils.write_filename(self.outdoc, filename, compress='auto')

        # save source probabilities in a json file
        if self.probabilities is not None:
            prob_fname = filename.replace('.xml.gz', '_probs.json')
            with open(prob_fname, 'w') as prob_outfile:
                json.dump(self.probabilities, prob_outfile)
            logging.info('Source probabilities file saved as %s', prob_fname)

    def upload(self, fname, gracedb_server=None, testing=True,
               extra_strings=None, search='AllSky'):
        """Upload this trigger to gracedb

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

        # first of all, make sure the event is saved on disk
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
            if fname.endswith('.xml.gz'):
                snr_series_fname = fname.replace('.xml.gz', '.hdf')
            else:
                snr_series_fname = fname.replace('.xml', '.hdf')
            snr_series_plot_fname = snr_series_fname.replace('.hdf',
                                                             '_snr.png')
            asd_series_plot_fname = snr_series_fname.replace('.hdf',
                                                             '_asd.png')
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

        if self.probabilities is not None:
            prob_fname = fname.replace('.xml.gz', '_probs.json')
            prob_plot_fname = prob_fname.replace('.json', '.png')

            prob_plot = {k: v for (k, v) in self.probabilities.items()
                         if v != 0.0}
            labels, sizes = zip(*prob_plot.items())
            colors = [source_color(label) for label in labels]
            fig, ax = pl.subplots()
            ax.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%',
                   textprops={'fontsize': 15})
            ax.axis('equal')
            fig.savefig(prob_plot_fname)
            pl.close()

        # upload SNR series in HDF format and plots
        if gid is not None and self.snr_series is not None:
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
                logging.error('Failed to upload plots for %s', gid)
                logging.error(str(exc))

        # upload source probabilities in JSON format and plot
        if gid is not None and self.probabilities is not None:
            try:
                gracedb.writeLog(
                    gid, 'Source probabilities JSON file upload',
                    filename=prob_fname, tag_name=['em_follow']
                )
                logging.info('Uploaded source probabilities for event %s', gid)
                gracedb.writeLog(
                    gid, 'Source probabilities plot upload',
                    filename=prob_plot_fname,
                    tag_name=['em_follow']
                )
                logging.info('Uploaded source probabilities pie chart for '
                             'event %s', gid)
            except Exception as exc:
                logging.error(
                    'Failed to upload source probability results for %s',
                    gid
                )
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

__all__ = ['SingleCoincForGraceDB', 'make_psd_xmldoc', 'snr_series_to_xml',
           'gracedb_tag_with_version']
