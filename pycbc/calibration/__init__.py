from pycbc.calibration.recalibrate import Recalibrate

models = {
    Recalibrate.name : Recalibrate
}

def read_model_from_config(cp, ifo, section="calibration"):
    """Returns an instance of the calibration model specified in the
    given configuration file.

    Parameters
    ----------
    cp : WorflowConfigParser
        An open config file to read.
    ifo : string
        The detector (H1, L1)  whose model will be loaded.
    section : {"calibration", string}
        Section name from which to retrieve the model.

    Returns
    -------
    instance
        An instance of the calibration model class.
    """
    name = cp.get_opt_tag(section, "name", None)
    recalibrator = models[name].from_config(cp, ifo.lower(), section)

    return recalibrator
