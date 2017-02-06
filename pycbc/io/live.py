import logging
import pycbc
import numpy
import lal
from glue.ligolw import ligolw
from glue.ligolw import lsctables
from glue.ligolw import utils as ligolw_utils
from glue.ligolw.utils import process as ligolw_process
from glue.ligolw import param as ligolw_param
from pycbc import version as pycbc_version
from pycbc import pnutils
from pycbc.tmpltbank import return_empty_sngl
from pycbc.types import TimeSeries, zeros
from pycbc.filter import correlate
from pycbc.fft import ifft

#FIXME Legacy build PSD xml helpers, delete me when we move away entirely from
# xml formats
def _build_series(series, dim_names, comment, delta_name, delta_unit):
    from glue.ligolw import array as ligolw_array
    Attributes = ligolw.sax.xmlreader.AttributesImpl
    elem = ligolw.LIGO_LW(Attributes({u"Name": unicode(series.__class__.__name__)}))
    if comment is not None:
        elem.appendChild(ligolw.Comment()).pcdata = comment
    elem.appendChild(ligolw.Time.from_gps(series.epoch, u"epoch"))
    elem.appendChild(ligolw_param.from_pyvalue(u"f0", series.f0, unit=u"s^-1"))
    delta = getattr(series, delta_name)
    if numpy.iscomplexobj(series.data.data):
        data = numpy.row_stack((numpy.arange(len(series.data.data)) * delta,
                             series.data.data.real, series.data.data.imag))
    else:
        data = numpy.row_stack((numpy.arange(len(series.data.data)) * delta, series.data.data))
    a = ligolw_array.from_array(series.name, data, dim_names=dim_names)
    a.Unit = str(series.sampleUnits)
    dim0 = a.getElementsByTagName(ligolw.Dim.tagName)[0]
    dim0.Unit = delta_unit
    dim0.Start = series.f0
    dim0.Scale = delta
    elem.appendChild(a)
    return elem

def make_psd_xmldoc(psddict):
    Attributes = ligolw.sax.xmlreader.AttributesImpl
    xmldoc = ligolw.Document()
    root_name = u"psd"
    lw = xmldoc.appendChild(ligolw.LIGO_LW(Attributes({u"Name": root_name})))
    for instrument, psd in psddict.items():
        xmlseries = _build_series(psd, (u"Frequency,Real", u"Frequency"),
                                  None, 'deltaF', 's^-1')
        fs = lw.appendChild(xmlseries)
        fs.appendChild(ligolw_param.from_pyvalue(u"instrument", instrument))
    return xmldoc

class SingleCoincForGraceDB(object):
    """Create xml files and submit them to gracedb from PyCBC Live"""
    def __init__(self, ifos, coinc_results, **kwargs):
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
        self.ifos = ifos
        self.template_id = coinc_results['foreground/%s/template_id' % self.ifos[0]]
        
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
        sngl_event_id_map = {}
        for ifo in ifos:
            names = [n.split('/')[-1] for n in coinc_results
                     if 'foreground/%s' % ifo in n]
            sngl_id += 1
            sngl = return_empty_sngl()
            sngl.event_id = lsctables.SnglInspiralID(sngl_id)
            sngl_event_id_map[ifo] = sngl.event_id
            sngl.ifo = ifo
            for name in names:
                val = coinc_results['foreground/%s/%s' % (ifo, name)]
                if name == 'end_time':
                    sngl.set_end(lal.LIGOTimeGPS(val))
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
        far = 1.0 / (lal.YRJUL_SI * coinc_results['foreground/ifar'])
        coinc_inspiral_row.combined_far = far
        coinc_inspiral_table.append(coinc_inspiral_row)
        outdoc.childNodes[0].appendChild(coinc_inspiral_table)
        self.outdoc = outdoc
        self.time = sngl.get_end()

        # compute SNR time series
        self.upload_snr_series = kwargs['upload_snr_series']
        if self.upload_snr_series:
            data_readers = kwargs['data_readers']
            bank = kwargs['bank']
            htilde = bank[self.template_id]
            self.snr_series = {}
            self.snr_series_psd = {}
            for ifo in self.ifos:
                stilde = data_readers[ifo].overwhitened_data(htilde.delta_f)
                norm = 4.0 * htilde.delta_f / (htilde.sigmasq(stilde.psd) ** 0.5)
                qtilde = zeros((len(htilde)-1)*2, dtype=htilde.dtype)
                correlate(htilde, stilde, qtilde)
                snr = qtilde * 0
                ifft(qtilde, snr)

                valid_end = int(len(qtilde) - data_readers[ifo].trim_padding)
                valid_start = int(valid_end - data_readers[ifo].blocksize * data_readers[ifo].sample_rate)
                seg = slice(valid_start, valid_end)
                snr = snr[seg]
                snr *= norm
                delta_t = 1.0 / data_readers[ifo].sample_rate
                start = data_readers[ifo].start_time
                snr = TimeSeries(snr, delta_t=delta_t, epoch=start)
                self.snr_series[ifo] = snr
                self.snr_series_psd[ifo] = stilde.psd

                # store the on-source slice of the series into the XML doc
                snr_onsource_time = coinc_results['foreground/%s/end_time' % ifo] - snr.start_time
                snr_onsource_dur = lal.REARTH_SI / lal.C_SI
                onsource_idx = round(snr_onsource_time * snr.sample_rate)
                onsource_start = onsource_idx - int(snr.sample_rate * snr_onsource_dur / 2)
                onsource_end = onsource_idx + int(snr.sample_rate * snr_onsource_dur / 2)
                onsource_slice = slice(onsource_start, onsource_end + 1)
                snr_lal = snr[onsource_slice].lal()
                snr_lal.name = 'snr'
                snr_lal.sampleUnits = ''
                snr_xml = _build_series(
                        snr_lal, (u'Time', u'Time,Real,Imaginary'), None,
                        'deltaT', 's')
                snr_node = outdoc.childNodes[-1].appendChild(snr_xml)
                eid_param = ligolw_param.new_param(
                        u'event_id', u'ilwd:char',
                        unicode(sngl_event_id_map[ifo]))
                snr_node.appendChild(eid_param)

    def save(self, filename):
        """Write this trigger to gracedb compatible xml format

        Parameters
        ----------
        filename: str
            Name of file to write to disk.
        """
        ligolw_utils.write_filename(self.outdoc, filename)

    def upload(self, fname, psds, low_frequency_cutoff,
               testing=True,
               extra_strings=None,
               ):
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

        self.save(fname)
        extra_strings = [] if extra_strings is None else extra_strings
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

        psds_lal = {}
        for ifo in psds:
            psd = psds[ifo]
            kmin = int(low_frequency_cutoff / psd.delta_f)
            fseries = lal.CreateREAL8FrequencySeries(
                "psd", psd.epoch, low_frequency_cutoff, psd.delta_f,
                lal.StrainUnit**2 / lal.HertzUnit, len(psd) - kmin)
            fseries.data.data = psd.numpy()[kmin:] / pycbc.DYN_RANGE_FAC ** 2.0
            psds_lal[ifo] = fseries
        psd_xmldoc = make_psd_xmldoc(psds_lal)

        ligolw_utils.write_filename(psd_xmldoc, "tmp_psd.xml.gz", gz=True)
        gracedb.writeLog(r["graceid"],
                         "PyCBC PSD estimate from the time of event",
                         "psd.xml.gz", open("tmp_psd.xml.gz", "rb").read(),
                         "psd").json()
        gracedb.writeLog(r["graceid"],
            "using pycbc code hash %s" % pycbc_version.git_hash).json()
        for text in extra_strings:
            gracedb.writeLog(r["graceid"], text).json()
        logging.info("Uploaded file psd.xml.gz to event %s.", r["graceid"])

        if self.upload_snr_series:
            snr_series_fname = fname + '.hdf'
            for ifo in self.ifos:
                self.snr_series[ifo].save(snr_series_fname,
                                          group='%s/snr' % ifo)
                self.snr_series_psd[ifo].save(snr_series_fname,
                                              group='%s/psd' % ifo)
            GraceDb().writeFile(r['graceid'], snr_series_fname)

        return r['graceid']

class SingleForGraceDB(SingleCoincForGraceDB):
    """Create xml files and submit them to gracedb from PyCBC Live"""
    def __init__(self, ifo, sngls_dict, hardware_injection=False):
        """Initialize a ligolw xml representation of this single trigger for
        upload to gracedb

        Parameters
        ----------
        ifo: str
            The IFO that the trigger came from.
        sngls_dict: dict
            Dictionary of singles parameters. Must include template parameters
        and both 'ifar' and 'stat' values.
        """
        fake_coinc = {}
        fake_coinc['foreground/stat'] = sngls_dict.pop('stat')
        fake_coinc['foreground/ifar'] = sngls_dict.pop('ifar')
        for key in sngls_dict:
            fake_coinc['foreground/%s/%s' % (ifo, key)] = sngls_dict[key]
        if hardware_injection:
            fake_coinc['HWINJ'] = True
        SingleCoincForGraceDB.__init__(self, [ifo], fake_coinc)
 
