from lal import LIGOTimeGPS, YRJUL_SI
from glue.ligolw import ligolw
from glue.ligolw import lsctables
from glue.ligolw import utils as ligolw_utils
from glue.ligolw.utils import process as ligolw_process
from pycbc import version as pycbc_version
from pycbc import pnutils
from pycbc.tmpltbank import return_empty_sngl
import logging
import pycbc


class SingleCoincForGraceDB(object):
    """Create xml files and submit them to gracedb from PyCBC Live"""
    def __init__(self, ifos, coinc_results):
        """Initialize a ligolw xml representation of a zerolag trigger
        for upload from pycbc live to gracedb.

        Parameters
        ----------
        ifos: list of strs
            A list of the ifos pariticipating in this trigger
        coinc_results: dict of values
            A dictionary of values. The format is define
            in pycbc/events/coinc.py
        and matches the on disk representation in the hdf file for this time.
        """
        # remember if this should be marked as HWINJ
        self.is_hardware_injection = False
        if 'HWINJ' in coinc_results:
            self.is_hardware_injection = True

        # Set up the bare structure of the xml document
        outdoc = ligolw.Document()
        outdoc.appendChild(ligolw.LIGO_LW())

        proc_id = ligolw_process.register_to_xmldoc(
            outdoc, 'pycbc',
            {}, ifos=ifos, comment='', version=pycbc_version.git_hash,
            cvs_repository='pycbc/'+pycbc_version.git_branch,
            cvs_entry_time=pycbc_version.date).process_id

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
        coinc_event_row.nevents = len(ifos)
        coinc_event_row.instruments = ','.join(ifos)
        coinc_event_row.time_slide_id = lsctables.TimeSlideID(0)
        coinc_event_row.process_id = proc_id
        coinc_event_row.coinc_event_id = coinc_id
        coinc_event_row.likelihood = 0.
        coinc_event_table.append(coinc_event_row)
        outdoc.childNodes[0].appendChild(coinc_event_table)

        # Set up sngls
        sngl_inspiral_table = lsctables.New(lsctables.SnglInspiralTable)
        coinc_event_map_table = lsctables.New(lsctables.CoincMapTable)

        sngl_id = 0
        for ifo in ifos:
            names = [n.split('/')[-1] for n in coinc_results
                     if 'foreground/%s' % ifo in n]
            sngl_id += 1
            sngl = return_empty_sngl()
            sngl.event_id = lsctables.SnglInspiralID(sngl_id)
            sngl.ifo = ifo
            for name in names:
                val = coinc_results['foreground/%s/%s' % (ifo, name)]
                if name == 'end_time':
                    sngl.set_end(LIGOTimeGPS(val))
                else:
                    try:
                        setattr(sngl, name, val)
                    except AttributeError:
                        pass
            sngl.mtotal, sngl.eta = pnutils.mass1_mass2_to_mtotal_eta(
                    sngl.mass1, sngl.mass2)
            sngl.mchirp, _ = pnutils.mass1_mass2_to_mchirp_eta(
                    sngl.mass1, sngl.mass2)
            sngl.eff_distance = (sngl.sigmasq)**0.5 / sngl.snr
            sngl_inspiral_table.append(sngl)

            # Set up coinc_map entry
            coinc_map_row = lsctables.CoincMap()
            coinc_map_row.table_name = 'sngl_inspiral'
            coinc_map_row.coinc_event_id = coinc_id
            coinc_map_row.event_id = sngl.event_id
            coinc_event_map_table.append(coinc_map_row)

        outdoc.childNodes[0].appendChild(coinc_event_map_table)
        outdoc.childNodes[0].appendChild(sngl_inspiral_table)

        # Set up the coinc inspiral table
        coinc_inspiral_table = lsctables.New(lsctables.CoincInspiralTable)
        coinc_inspiral_row = lsctables.CoincInspiral()
        # This seems to be used as FAP, which should not be in gracedb
        coinc_inspiral_row.false_alarm_rate = 0
        coinc_inspiral_row.minimum_duration = 0.
        coinc_inspiral_row.set_ifos(ifos)
        coinc_inspiral_row.coinc_event_id = coinc_id
        coinc_inspiral_row.mchirp = sngl.mchirp
        coinc_inspiral_row.mass = sngl.mtotal
        coinc_inspiral_row.end_time = sngl.end_time
        coinc_inspiral_row.end_time_ns = sngl.end_time_ns
        coinc_inspiral_row.snr = coinc_results['foreground/stat']
        far = 1.0 / (YRJUL_SI * coinc_results['foreground/ifar'])
        coinc_inspiral_row.combined_far = far
        coinc_inspiral_table.append(coinc_inspiral_row)
        outdoc.childNodes[0].appendChild(coinc_inspiral_table)
        self.outdoc = outdoc

    def save(self, filename):
        """Write this trigger to gracedb compatible xml format

        Parameters
        ----------
        filename: str
            Name of file to write to disk.
        """
        ligolw_utils.write_filename(self.outdoc, filename)

    def upload(self, fname, psds, low_frequency_cutoff, testing=True):
        """Upload this trigger to gracedb

        Parameters
        ----------
        fname: str
            The name to give the xml file associated with this trigger
        pds: dict of pybc.types.FrequencySeries
            A ifo keyed dictionary of psds to be uploaded in association
        with this trigger.
        low_frequency_cutoff: float
            The low frequency cutoff of the psds.
        testing: bool
            Switch to determine if the upload should be sent to gracedb as a
        test trigger (True) or a production trigger (False)
        """
        from ligo.gracedb.rest import GraceDb
        import lal
        import lal.series

        self.save(fname)
        if testing:
            group = 'Test'
        else:
            group = 'CBC'

        gracedb = GraceDb()
        r = gracedb.createEvent(group, "pycbc", fname, "AllSky").json()
        logging.info("Uploaded event %s.", r["graceid"])

        if self.is_hardware_injection:
            gracedb.writeLabel(r['graceid'], 'INJ')
            logging.info("Tagging event %s as an injection", r["graceid"])

        # Convert our psds to the xml psd format.
        # FIXME: we should not use lal.series!!!
        psds_lal = {}
        for ifo in psds:
            psd = psds[ifo]
            kmin = int(low_frequency_cutoff / psd.delta_f)
            fseries = lal.CreateREAL8FrequencySeries(
                "psd", psd.epoch, low_frequency_cutoff, psd.delta_f,
                lal.StrainUnit**2 / lal.HertzUnit, len(psd) - kmin)
            fseries.data.data = psd.numpy()[kmin:] / pycbc.DYN_RANGE_FAC ** 2.0
            psds_lal[ifo] = fseries

        psd_xmldoc = lal.series.make_psd_xmldoc(psds_lal)
        ligolw_utils.write_filename(psd_xmldoc, "tmp_psd.xml.gz", gz=True)
        gracedb.writeLog(r["graceid"],
                         "PyCBC PSD estimate from the time of event",
                         "psd.xml.gz", open("tmp_psd.xml.gz", "rb").read(),
                         "psd").json()
        logging.info("Uploaded file psd.xml.gz to event %s.", r["graceid"])
