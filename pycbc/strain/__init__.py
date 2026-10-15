from .gate import (
    add_gate_option_group,
    apply_gates_to_fd,
    apply_gates_to_td,
    gates_from_cli,
    psd_gates_from_cli,
)
from .recalibrate import CubicSpline, PhysicalModel
from .strain import (
    StrainBuffer,
    StrainSegments,
    detect_loud_glitches,
    from_cli,
    from_cli_multi_ifos,
    from_cli_single_ifo,
    gate_data,
    insert_strain_option_group,
    insert_strain_option_group_multi_ifo,
    verify_strain_options,
    verify_strain_options_multi_ifo,
)

models = {CubicSpline.name: CubicSpline, PhysicalModel.name: PhysicalModel}


def read_model_from_config(cp, ifo, section="calibration"):
    """
    Returns an instance of the calibration model specified in the
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
    model = cp.get_opt_tag(section, f"{ifo.lower()}_model", None)
    recalibrator = models[model].from_config(cp, ifo.lower(), section)

    return recalibrator
