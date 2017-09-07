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
    strain : FrequencySeries
        The strain data that will be recalibrated.
    ifo : string
        The detector (H1, L1)  whose model will be loaded.
    section : {"calibration", string}
        Section name from which to retrieve the model.

    Returns
    -------
    instance
        An instance of the calibration model class.
    """
    name = cp.get_opt_tag(section, "name")
    model_args = models[name].from_config(cp, ifo)
    freq, fc0, c0, d0, a_tst0, a_pu0, fs0, qinv0 = model_args
    recalibrator = models[name](strain, freq=freq, fc0=fc0, c0=c0, d0=d0,
                                a_tst0=a_tst0, a_pu0=a_pu0, fs0=fs0,
                                qinv0=qinv0)
    return recalibrator
