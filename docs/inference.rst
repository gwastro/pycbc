###################################################################
PyCBC inference documentation (``pycbc.inference``)
###################################################################

.. Following are useful substituions for classes and modules
.. GaussianNoise:
.. |GaussianNoise| replace:: :py:class:`GaussianNoise <pycbc.inference.models.gaussian_noise.GaussianNoise>`
.. MarginalizedPhaseGaussianNoise:
.. |MarginalizedPhaseGaussianNoise| replace:: :py:class:`MarginalizedPhaseGaussianNoise <pycbc.inference.models.marginalized_gaussian_noise.MarginalizedPhaseGaussianNoise>`

===================
Introduction
===================

This page gives details on how to use the various parameter estimation
executables and modules available in PyCBC. The ``pycbc.inference`` subpackage
contains classes and functions for evaluating probability distributions,
likelihoods, and running Bayesian samplers.

==================================================
Sampling the parameter space (``pycbc_inference``)
==================================================

--------
Overview
--------

The executable ``pycbc_inference`` is designed to sample the parameter space
and save the samples in an HDF file. A high-level description of the
``pycbc_inference`` algorithm is

#. Read priors from a configuration file.

#. Setup the model to use. If the model uses data, then:

    * Read gravitational-wave strain from a gravitational-wave model or use
      recolored fake strain.

    * Estimate a PSD.

#. Run a sampler to estimate the posterior distribution of the model.

#. Write the samples and metadata to an HDF file.

The model, data, sampler, parameters to vary and their priors are specified in
one or more configuration files, which are passed to the program using the
``--config-file`` option. Other command-line options determine what
parallelization settings to use. For a full listing of all options run
``pycbc_inference --help``. Below, we give details on how to set up a
configuration file and provide examples of how to run ``pycbc_inference``.

------------------------------------------------
Configuring the model, sampler, priors, and data
------------------------------------------------

The configuration file(s) uses :py:class:`WorkflowConfigParser
<pycbc.workflow.configuration.WorkflowConfigParser>` syntax.  The required
sections are: ``[model]``, ``[sampler]``, and ``[variable_params]``.  In
addition, multiple ``[prior]`` sections must be provided that define the prior
distribution to use for the parameters in ``[variable_params]``. If a model
uses data a ``[data]`` section must also be provided.

These sections may be split up over multiple files. In that case, all of the
files should be provided as space-separated arguments to the ``--config-file``.
Providing multiple files is equivalent to providing a single file with
everything across the files combined. If the same section is specified in
multiple files, the all of the options will be combined.

Configuration files allow for referencing values in other sections using the
syntax ``${section|option}``. See the examples below for an example of this.
When providing multiple configuration files, sections in other files may be
referenced, since in the multiple files are combined into a single file in
memory when the files are loaded.

^^^^^^^^^^^^^^^^^^^^^
Configuring the model
^^^^^^^^^^^^^^^^^^^^^

The ``[model]`` section sets up what model to use for the analysis. At minimum,
a ``name`` argument must be provided, specifying which model to use. For
example:

.. code-block:: ini

   [model]
   name = gaussian_noise

In this case, the |GaussianNoise| model would be used. (Examples of using this
model on a BBH injection and on GW150914 are given below.) Other arguments to
configure the model may also be set in this section. The recognized arguments
depend on the model. The currently available models are:

.. include:: _include/models-table.rst

Refer to the models' ``from_config`` method to see what configuration arguments
are available.

Any model name that starts with ``test_`` is an analytic test distribution that
requires no data or waveform generation. See the section below on running on an
analytic distribution for more details.


^^^^^^^^^^^^^^^^^^^^^^^
Configuring the sampler
^^^^^^^^^^^^^^^^^^^^^^^

The ``[sampler]`` section sets up what sampler to use for the analysis. As
with the ``[model]`` section, a ``name`` must be provided to specify which
sampler to use. The currently available samplers are:

.. include:: _include/samplers-table.rst

See :ref:`example of trying different samplers<inference_example_samplers>`

Configuration options for the sampler should also be specified in the
``[sampler]`` section. For example:

.. code-block:: ini

   [sampler]
   name = emcee
   nwalkers = 5000
   niterations = 1000
   checkpoint-interval = 100

This would tell ``pycbc_inference`` to run the
:py:class:`EmceeEnsembleSampler <pycbc.inference.sampler.emcee.EmceeEnsembleSampler>`
with 5000 walkers for 1000 iterations, checkpointing every 100th iteration.
Refer to the samplers' ``from_config`` method to see what configuration options
are available.

Burn-in tests may also be configured for MCMC samplers in the config file. The
options for the burn-in should be placed in ``[sampler-burn_in]``. At minimum,
a ``burn-in-test`` argument must be given in this section. This argument
specifies which test(s) to apply. Multiple tests may be combined using standard
python logic operators. For example:

.. code-block:: ini

   [sampler-burn_in]
   burn-in-test = nacl & max_posterior

In this case, the sampler would be considered to be burned in when both the
``nacl`` *and* ``max_posterior`` tests were satisfied. Setting this to ``nacl |
max_postrior`` would instead consider the sampler to be burned in when either
the ``nacl`` *or* ``max_posterior`` tests were satisfied. For more information
on what tests are available, see the :py:mod:`pycbc.inference.burn_in` module.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Thinning samples (MCMC only)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The default behavior for the MCMC samplers (``emcee``, ``emcee_pt``) is to save
every iteration of the Markov chains to the output file. This can quickly lead
to very large files. For example, a BBH analysis (~15 parameters) with 200
walkers, 20 temperatures may take ~50 000 iterations to acquire ~5000
independent samples. This will lead to a file that is ~ 50 000 iterations x 200
walkers x 20 temperatures x 15 parameters x 8 bytes ~ 20GB.  Quieter signals
can take an order of magnitude more iterations to converge, leading to O(100GB)
files. Clearly, since we only obtain 5000 independent samples from such a run,
the vast majority of these samples are of little interest.

To prevent large file size growth, samples may be thinned before they are
written to disk. Two thinning options are available, both of which are set in
the ``[sampler]`` section of the configuration file. They are:

 * ``thin-interval``: This will thin the samples by the given integer before
   writing the samples to disk. File sizes can still grow unbounded, but at
   a slower rate. The interval must be less than the checkpoint interval.
 * ``max-samples-per-chain``: This will cap the maximum number of samples per
   walker and per temperature to the given integer. This ensures that file
   sizes never exceed ~ ``max-samples-per-chain`` x ``nwalkers`` x ``ntemps``
   x ``nparameters`` x 8 bytes. Once the limit is reached,
   samples will be thinned on disk, and new samples will be thinned to match.
   The thinning interval will grow with longer runs as a result. To ensure
   that enough samples exist to determine burn in and to measure an
   autocorrelation length, ``max-samples-per-chain`` must be greater than
   or equal to 100.

The thinned interval that was used for thinning samples is saved to the output
file's ``thinned_by`` attribute (stored in the HDF file's ``.attrs``).  Note
that this is not the autocorrelation length (ACL), which is the amount that the
samples need to be further thinned to obtain independent samples.


.. note::

    In the output file creates by the MCMC samplers, we adopt the convention
    that "iteration" means iteration of the sampler, not index of the samples.
    For example, if a burn in test is used, ``burn_in_iteration`` will be
    stored to the ``sampler_info`` group in the output file. This gives the
    iteration of the sampler at which burn in occurred, not the sample on disk.
    To determine  which samples an iteration corresponds to in the file, divide
    iteration by ``thinned_by``.

    Likewise, we adopt the convention that autocorrelation **length** (ACL) is
    the autocorrelation length of the thinned samples (the number of samples on
    disk that you need to skip to get independent samples) whereas
    autocorrelation **time** (ACT) is the autocorrelation length in terms of
    iteration (it is the number of **iterations** that you need to skip to get
    independent samples); i.e., ``ACT = thinned_by x ACL``. The ACT is (up to
    measurement resolution) independent of the thinning used, and thus is
    useful for comparing the performance of the sampler.



^^^^^^^^^^^^^^^^^^^^^
Configuring the prior
^^^^^^^^^^^^^^^^^^^^^

What parameters to vary to obtain a posterior distribution are determined by
``[variable_params]`` section. For example:

.. code-block:: ini

   [variable_params]
   x =
   y =

This would tell ``pycbc_inference`` to sample a posterior over two parameters
called ``x`` and ``y``.

A prior must be provided for every parameter in ``[variable_params]``. This
is done by adding sections named ``[prior-{param}]`` where ``{param}`` is the
name of the parameter the prior is for. For example, to provide a prior for the
``x`` parameter in the above example, you would need to add a section called
``[prior-x]``. If the prior couples more than one parameter together in a joint
distribution, the parameters should be provided as a ``+`` separated list,
e.g., ``[prior-x+y+z]``.

The prior sections specify what distribution to use for the parameter's prior,
along with any settings for that distribution.  Similar to the ``model`` and
``sampler`` sections, each ``prior`` section must have a ``name`` argument that
identifies the distribution to use. Distributions are defined in the
:py:mod:`pycbc.distributions` module. The currently available distributions
are:

.. include:: _include/distributions-table.rst

^^^^^^^^^^^^^^^^^
Static parameters
^^^^^^^^^^^^^^^^^

A ``[static_params]`` section may be provided to list any parameters that
will remain fixed throughout the run. For example:

.. code-block:: ini

   [static_params]
   approximant = IMRPhenomPv2
   f_lower = 18

^^^^^^^^^^^^
Setting data
^^^^^^^^^^^^

Many models, such as the |GaussianNoise| model, require data to be provided. To
do so, a ``[data]`` section must be included that provides information about
what data to load, and how to condition it.

The type of data to be loaded depends on the model. For example, if you are
using the |GaussianNoise| or |MarginalizedPhaseGaussianNoise| models (the
typical case), one will need to load gravitational-wave data.  This is
accomplished using tools provided in the :py:mod:`pycbc.strain` module. The
full set of options are:

.. include:: _include/inference_data_opts-table.rst

As indicated in the table, the ``psd-model`` and ``fake-strain`` options can
accept an analytical PSD as an argument. The available PSD models are:

.. include:: _include/psd_models-table.rst

-------------------------------
Advanced configuration settings
-------------------------------

The following are additional settings that may be provided in the configuration
file, in order to do more sophisticated analyses.

^^^^^^^^^^^^^^^^^^^
Sampling transforms
^^^^^^^^^^^^^^^^^^^

One or more of the ``variable_params`` may be transformed to a different
parameter space for purposes of sampling. This is done by specifying a
``[sampling_params]`` section. This section specifies which
``variable_params`` to replace with which parameters for sampling. This must be
followed by one or more ``[sampling_transforms-{sampling_params}]`` sections
that provide the transform class to use. For example, the following would cause
the sampler to sample in chirp mass (``mchirp``) and mass ratio (``q``) instead
of ``mass1`` and ``mass2``:

.. code-block:: ini

   [sampling_params]
   mass1, mass2: mchirp, q

   [sampling_transforms-mchirp+q]
   name = mass1_mass2_to_mchirp_q

Transforms are provided by the :py:mod:`pycbc.transforms` module. The currently
available transforms are:

.. include:: _include/transforms-table.rst


.. note::
   Both a ``jacobian`` and ``inverse_jacobian`` must be defined in order to use
   a transform class for a sampling transform. Not all transform classes in
   :py:mod:`pycbc.transforms` have these defined. Check the class
   documentation to see if a Jacobian is defined.

^^^^^^^^^^^^^^^^^^^
Waveform transforms
^^^^^^^^^^^^^^^^^^^

There can be any number of ``variable_params`` with any name. No parameter name
is special (with the exception of parameters that start with ``calib_``; see
below).

However, when doing parameter estimation with CBC waveforms, certain parameter
names must be provided for waveform generation. The parameter names recognized
by the CBC waveform generators are:

.. include:: _include/waveform-parameters.rst

It is possible to specify a ``variable_param`` that is not one of these
parameters. To do so, you must provide one or more
``[waveforms_transforms-{param}]`` section(s) that define transform(s) from the
arbitrary ``variable_params`` to the needed waveform parameter(s) ``{param}``.
For example, in the following we provide a prior on ``chirp_distance``. Since
``distance``, not ``chirp_distance``, is recognized by the CBC waveforms
module, we provide a transform to go from ``chirp_distance`` to ``distance``:

.. code-block:: ini

   [variable_params]
   chirp_distance =

   [prior-chirp_distance]
   name = uniform
   min-chirp_distance = 1
   max-chirp_distance = 200

   [waveform_transforms-distance]
   name = chirp_distance_to_distance

A useful transform for these purposes is the
:py:class:`CustomTransform <pycbc.transforms.CustomTransform>`, which allows
for arbitrary transforms using any function in the :py:mod:`pycbc.conversions`,
:py:mod:`pycbc.coordinates`, or :py:mod:`pycbc.cosmology` modules, along with
numpy math functions. For example, the following would use the I-Love-Q
relationship :py:meth:`pycbc.conversions.dquadmon_from_lambda` to relate the
quadrupole moment of a neutron star ``dquad_mon1`` to its tidal deformation
``lambda1``:

.. code-block:: ini

   [variable_params]
   lambda1 =

   [waveform_transforms-dquad_mon1]
   name = custom
   inputs = lambda1
   dquad_mon1 = dquadmon_from_lambda(lambda1)

.. note::
   A Jacobian is not necessary for waveform transforms, since the transforms
   are only being used to convert a set of parameters into something that the
   waveform generator understands. This is why in the above example we are
   able to use a custom transform without needing to provide a Jacobian.

Some common transforms are pre-defined in the code. These are: the mass
parameters ``mass1`` and ``mass2`` can be substituted with ``mchirp`` and
``eta`` or ``mchirp`` and ``q``.  The component spin parameters ``spin1x``,
``spin1y``, and ``spin1z`` can be substituted for polar coordinates
``spin1_a``, ``spin1_azimuthal``, and ``spin1_polar`` (ditto for ``spin2``).

^^^^^^^^^^^^^^^^^^^^^^
Calibration parameters
^^^^^^^^^^^^^^^^^^^^^^

If any calibration parameters are used (prefix ``calib_``), a ``[calibration]``
section must be included. This section must have a ``name`` option that
identifies what calibration model to use. The models are described in
:py:mod:`pycbc.calibration`. The ``[calibration]`` section must also include
reference values ``fc0``, ``fs0``, and ``qinv0``, as well as paths to ASCII
transfer function files for the test mass actuation, penultimate mass
actuation, sensing function, and digital filter for each IFO being used in the
analysis. E.g. for an analysis using H1 only, the required options would be
``h1-fc0``, ``h1-fs0``, ``h1-qinv0``, ``h1-transfer-function-a-tst``,
``h1-transfer-function-a-pu``, ``h1-transfer-function-c``,
``h1-transfer-function-d``.

^^^^^^^^^^^
Constraints
^^^^^^^^^^^

One or more constraints may be applied to the parameters; these are
specified by the ``[constraint]`` section(s). Additional constraints may be
supplied by adding more ``[constraint-{tag}]`` sections. Any tag may be used; the
only requirement is that they be unique. If multiple constraint sections are
provided, the union of all constraints is applied. Alternatively, multiple
constraints may be joined in a single argument using numpy's logical operators.

The parameter that constraints are applied to may be any parameter in
``variable_params`` or any output parameter of the transforms. Functions may be
applied to these parameters to obtain constraints on derived parameters. Any
function in the conversions, coordinates, or cosmology module may be used,
along with any numpy ufunc. So, in the following example, the mass ratio (q) is
constrained to be <= 4 by using a function from the conversions module.

.. code-block:: ini

   [variable_params]
   mass1 =
   mass2 =

   [prior-mass1]
   name = uniform
   min-mass1 = 3
   max-mass1 = 12

   [prior-mass2]
   name = uniform
   min-mass2 = 1
   min-mass2 = 3

   [constraint-1]
   name = custom
   constraint_arg = q_from_mass1_mass2(mass1, mass2) <= 4

------------------------------
Checkpointing and output files
------------------------------

While ``pycbc_inference`` is running it will create a checkpoint file which
is named ``{output-file}.checkpoint``, where ``{output-file}`` was the name
of the file you specified with the ``--output-file`` command. When it
checkpoints it will dump results to this file; when finished, the file is
renamed to ``{output-file}``. A ``{output-file}.bkup`` is also created, which
is a copy of the checkpoint file. This is kept in case the checkpoint file gets
corrupted during writing. The ``.bkup`` file is deleted at the end of the run,
unless ``--save-backup`` is turned on.

When ``pycbc_inference`` starts, it checks if either
``{output-file}.checkpoint`` or ``{output-file}.bkup`` exist (in that order).
If at least one of them exists, ``pycbc_inference`` will attempt to load them
and continue to run from the last checkpoint state they were in.

The output/checkpoint file are HDF files. To peruse the structure of the file
you can use the `h5ls` command-line utility. More advanced utilities for
reading and writing from/to them are provided by the sampler IO classes in
:py:mod:`pycbc.inference.io`. To load one of these files in python do:

.. code-block:: python

   from pycbc.inference import io
   fp = io.loadfile(filename, "r")

Here, ``fp`` is an instance of a sampler IO class. Basically, this is an
instance of an :py:mod:`h5py.File <h5py:File>` handler, with additional
convenience functions added on top. For example, if you want all of the samples
of all of the variable parameters in the file, you can do:

.. code-block:: python

   samples = fp.read_samples(fp.variable_params)

This will return a :py:class:`FieldArray <pycbc.io.record.FieldArray>` of all
of the samples.

Each sampler has it's own sampler IO class that adds different convenience
functions, depending on the sampler that was used. For more details on these
classes, see the :py:mod:`pycbc.inference.io` module.

===============================================
Examples
===============================================

.. toctree::
    :maxdepth: 1

    inference/examples/sampler_platter.rst
    inference/examples/analytic.rst
    inference/examples/bbh.rst
    inference/examples/gw150914.rst
    inference/examples/single.rst
    inference/examples/relative.rst

===============================================
Visualizing the Posteriors
===============================================

.. toctree::
   :maxdepth: 1

   inference/viz.rst

===============================================
Workflows
===============================================

.. toctree::
   :maxdepth: 1

   workflow/pycbc_make_inference_workflow
   workflow/pycbc_make_inference_inj_workflow

===============================================
For Developers
===============================================

.. toctree::
    :maxdepth: 1

    inference/sampler_api.rst
    inference/io.rst
