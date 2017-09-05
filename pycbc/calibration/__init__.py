from pycbc.calibration.recalibrate import Recalibrate

models = {
    Recalibrate.name : Recalibrate
}

def read_model_from_config(cp, section="calibration"):
    """

    Parameters
    ----------
    cp : WorflowConfigParser
        An open config file to read.
    section : {"calibration", string}
        Section name from which to retrieve the model.

    Returns
    -------

    """
    pass
