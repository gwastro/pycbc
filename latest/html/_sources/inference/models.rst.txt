.. _models_detailed:

---------------------------------------------
Details of common Models in PyCBC Inference
---------------------------------------------

Commonly used likelihood models are compared below. PyCBC Inference
also has the capability to interface with external models. See the
`interactive tutorial <https://colab.research.google.com/github/gwastro/pycbc-tutorials/blob/master/tutorial/inference_9_AddingCustomModels.ipynb>`_ if you'd like to learn how to add you own
model for handling a new problem such as GRBs, Kilonova or whatever you can
imagine.

==============================================
Standard models with full waveform generation
==============================================

.. card:: Gaussian Noise

    ``'gaussian_noise'`` :py:class:`pycbc.inference.models.gaussian_noise.GaussianNoise`

    The `gaussian_noise` model is one of the most generic models in
    PyCBC Inference. It is designed for gravitational-wave data that is
    assumed to be Gaussian and stationary. Otherwise, there are no
    restricts on the types of waveform models that may be used.
    Waveform models supported by either `get_fd_waveform` or
    `get_td_waveform` can be used as waveform models. Use this model only if
    none other fits and is more optimized for your specific problem. Because
    this waveform model is generic is often slower than the other more
    specialized models.

    Supported Marginalizations: None
    +++
    Earth Rotation:❌ LISA:❌ Higher Modes:✅
    :ref:`Example <inference_example_bbh>`

.. card:: Marginalized Phase

    ``'marginalized_phase'`` :py:class:`pycbc.inference.models.marginalized_gaussian_noise.MarginalizedPhaseGaussianNoise`


    The `marginalized_phase` model is a straightforward extension to
    the `gaussian_noise` model that marginalizes the signal over an overall
    phase. This can account for the orbital phase of a gravitational-wave
    signal which is usually a nuissance parameter. However, this
    implementation is only correct for signals that only include the
    dominant gravitational-wave mode (or 22 mode). This is because each
    mode has a different functional dependence on the orbital phase.


    Supported Marginalizations: coa_phase (automatic, dominant-mode)
    +++
    Earth Rotation:❌ LISA:❌ Higher Modes:❌

.. card:: Marginalized Polarization

    ``'marginalized_polarization'`` :py:class:`pycbc.inference.models.marginalized_gaussian_noise.MarginalizedPolarization`


    The `marginalized_polarization` model is a straightforward extension to
    the `gaussian_noise` model that marginalizes the signal over
    the polarization angle of the gravitational-wave signal. This is
    done by filtering both the plus and cross polarizations separately
    and then numerically integrating the resulting inner products over
    a grid of polarization angle values. The density of the grid
    is selectable. Additional marginalizations can also be optionally enabled.
    This model is often used in perference to the `marginalized_phase` model
    because it can support waveform models with higher order modes.

    Supported Marginalizations: polarization (automatic), coa_phase (dominant-mode), distance
    +++
    Earth Rotation:❌ LISA:❌ Higher Modes:✅

.. card:: Marginalized Time

    ``'marginalized_time'`` :py:class:`pycbc.inference.models.marginalized_gaussian_noise.MarginalizedTime`


    The `marginalized_time` model calculates a time series of inner products
    between the signal model and data so that the likelihood can be
    numerically marginalized over. The marginalization uses a weighted
    monte-carlo which samples the most likely regions of the time series
    more densly. The time series is also sub-sample interpolated. Higher
    mode waveforms are also permitted as both the plus and cross polarizations
    are explicitly handled separately. The time marginalization can also be
    done in concert with sky and various intrinsic marginalizations.

    Supported Marginalizations: tc (time), distance, coa_phase (dominant mode), polarization, ra dec
    +++
    Earth Rotation:❌ LISA:❌ Higher Modes:✅
    :ref:`Example <marginalized_time_example>`

.. card:: Marginalized Higher Mode Phase

    ``'marginalized_hmpolphase'`` :py:class:`pycbc.inference.models.marginalized_gaussian_noise.MarginalizedHMPolPhase`


    The `marginalized_hmpolphase` model numerically marginalizes the likelihood
    over a grid of both polarization values and orbital phase values. It
    explicitly handles this for higher-order-mode waveforms by calculating
    all inner products between the data and signal on a mode-by-mode basis.
    The likelihoods are then assembled for each polarization and phase value
    and numerically integrated over. This model can only be used for
    waveform approximants which support the PyCBC higher mode interface
    as we need to be able to calculate the each mode separately.

    Supported Marginalizations: polarization (automatic), coa_phase (automatic)
    +++
    Earth Rotation:❌ LISA:❌ Higher Modes:✅

=========================================
Heterodyne / Relative Models
=========================================

.. card:: Relative Model

    ``'relative'`` :py:class:`pycbc.inference.models.relbin.Relative`

    The `relative` model uses a reference signal provided in its
    configuration to expand the likelihood in terms of differences between
    the reference signal and a target signal. If the reference is close to the
    peak in the likelihood, the reasonable portions of parameter space to
    explore will only have small phase deviations from the reference signal
    (i.e. several cycles). This allows us to represent the ration between
    target signal and our referencee using a piece-wise linear approximation.
    The target signal only needs to calculated at edges of this approximation.
    This model thus only supports waveforms which can efficiently generate
    a signal only at a given set of frequency values. Higher order modes
    are not recommended due their possible violation of the approximation
    that the ration between target and reference signals is a slowly varying
    smooth function. In this case use the `multi_signal` model and treat
    use a `relative` model for each mode. Where the model approximations hold,
    use this in preference
    to the models that need to generate a full waveform for every likelihood
    as these will usually be much faster.

    There is also support in this model for use with :ref:`LISA Sangria data analysis <inference_example_lisa_smbhb_ldc>` and :ref:`LISA injection data analysis <inference_example_lisa_smbhb_inj>`.

    Supported Marginalizations: distance, coa_phase (dominant mode), polarization
    +++
    Earth Rotation:✅ LISA:✅ Higher Modes:❌
    :ref:`Example <relative_example>`

.. card:: Relative Time

    ``'relative_time'`` :py:class:`pycbc.inference.models.relbin.RelativeTime`

    The `relative_time` model extends the `relative model` by using calculating
    the likelihood for a grid of possible merger times. A weighted monte-carlo
    is performed to integrate over the merger time. Sub-sample interpolation
    of the time series is performed. Sky location marginalization among others
    can be done in concert.

    Supported Marginalizations: distance, coa_phase (dominant mode), polarization, ra dec
    +++
    Earth Rotation:❌ LISA:❌ Higher Modes:❌

.. card:: Relative Time Dominant-mode Only

    ``'relative_time'`` :py:class:`pycbc.inference.models.relbin.RelativeTime`

    The `relative_time_dom` model further specializes the model to completely
    preclude use with higher order mode signals. For this restriction, it
    gains the ability to marginalize over inclination

    Supported Marginalizations: distance, coa_phase (dominant mode), polarization, inclination, ra dec
    +++
    Earth Rotation:❌ LISA:❌ Higher Modes:❌

.. card:: Brute force LISA sky modes

    ``'brute_lisa_sky_modes_marginalize'`` :py:class:`pycbc.inference.models.relbin.Relative`

    The models does a brute force marginalization over the LISA sky mode
    degeneracies. It is built upon the `relative` model

    Supported Marginalizations: distance, coa_phase (dominant mode)
    +++
    Earth Rotation:❌ LISA:✅ Higher Modes:❌

=========================================
Extrinsic Parameter Only Models
=========================================

.. card:: Single Template

    ``'single_template'`` :py:class:`pycbc.inference.models.single_template.SingleTemplate`

    The `single_template` model is only for extrinsic parameter estimation e.g.
    sky location, distance, inclination, etc. It speeds up parameter estimation
    by completely avoiding waveform generation while calculating likelihoods.
    A single reference signal is generated with fixed intrinsic (masses, spins, etc)
    parameters. The likelihoods can then be precalculated up to constant factors
    which vary with the extrinsic parameters. Only dominant-mode signals are supported
    as the plus and cross polarizations are assumed to be different only by
    a phase. With this model all supported parameters may be marginalized over
    or any subset.

    Supported Marginalizations: tc (time), distance, coa_phase (dominant mode), polarization, ra dec
    +++
    Earth Rotation:❌ LISA:❌ Higher Modes:❌
    :ref:`Example <single_template_examples>`

=================================================
Composite Models
=================================================

.. card:: Hierarchical

    :py:class:`pycbc.inference.models.hierarchical.HierarchicalModel`

    The hierachical model is a container for submodels. Each submodel
    makes an indepedent contribution to an overall likelihood function.

    See :ref:`the full page on the hierarchical model <hierachical_model>`.

.. card:: Multiple Signal

    ``'multi_signal'`` :py:class:`pycbc.inference.models.hierarchical.MultiSignalModel`

    This is a container for submodels where each model shares use of the same
    data and hence there may be cross terms between the several signals.
    This model requires support from the submodel to handle cross-term
    calculations. Some supported models include the `gaussian_noise`, `relative`
    and `single_template`.
