from pycbc.calibration.recalibrate import Recalibrate

models = {
    Recalibrate.name : Recalibrate
}

def read_model_from_config(cp, strain, ifo, section="calibration"):
    """Returns an instance of the calibration model specified in the
    given configuration file.

    Parameters
    ----------
    cp : WorflowConfigParser
        An open config file to read.
    strain :
    ifo : string
        The detector wh
    section : {"calibration", string}
        Section name from which to retrieve the model.

    Returns
    -------
    Instance
        An instance of the calibration model class.
    """
    name = cp.get_opt_tag(section, "name")
    freq, fc0, c0, d0, a_tst0, a_pu0, fs0, qinv0 = Recalibrate.fromconfig(cp, ifo)
    recalibrator = models[name](strain, freq=freq, fc0=fc0, c0=c0, d0=d0,
                                a_tst0=a_tst0, a_pu0=a_pu0, fs0=fs0,
                                qinv0=qinv0)
    return recalibrator
