###################################################################
PyCBC inference documentation (``pycbc.inference``)
###################################################################

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

#. Estimate a PSD from a model or data.

#. Read gravitational-wave strain from a gravitational-wave model or use recolored fake strain.

#. Read priors from configuration file.

#. Construct prior-weighted likelihood function from PSD, gravitational-wave strain, and priors.

#. Run sampler that walks around parameter space and calculates the prior-weighted likelihood function.

---------------------------------------------------
Options for samplers, likelihood models, and priors
---------------------------------------------------

For a full listing of all options run ``pycbc_inference --help``. In this subsection we reference documentation for Python classes that contain more information about choices for samplers, likelihood models, and priors.

The user specifies the sampler on the command line with the ``--sampler`` option.
A complete list of samplers is given in ``pycbc_inference --help``.
These samplers are described in :py:class:`pycbc.inference.sampler_kombine.KombineSampler`, :py:class:`pycbc.inference.sampler_emcee.EmceeEnsembleSampler`, and :py:class:`pycbc.inference.sampler_emcee.EmceePTSampler`.
In addition to ``--sampler`` the user will need to specify the number of walkers to use ``--nwalkers``, and for parallel-tempered samplers the number of temperatures ``--ntemps``. You also need to either specify the number of iterations to run for using ``--niterations`` **or** the number of independent samples to collect using ``--n-independent-samples``. For the latter, a burn-in function must be specified using ``--burn-in-function``. In this case, the program will run until the sampler has burned in, at which point the number of independent samples equals the number of walkers. If the number of independent samples desired is greater than the number of walkers, the program will continue to run until it has collected the specified number of independent samples (to do this, an autocorrelation length is computed at each checkpoint to determine how many iterations need to be skipped to obtain independent samples).

The user specifies the likelihood model on the command line with the ``--likelihood-evaluator`` option. Any choice that starts with ``test_`` is an analytic test distribution that requires no data or waveform generation; see the section below on running on an analytic distribution for more details. For running on data, use ``--likelihood-evaluator gaussian``; this uses :py:class:`pycbc.inference.likelihood.GaussianLikelihood` for evaluating posteriors. Examples of using this on a BBH injection and on GW150914 are given below.

The user specifies a configuration file that defines the priors with the ``--config-files`` option.
The syntax of the configuration file is described in the following subsection.

-------------------------
Configuration file syntax
-------------------------

Configuration files follow the ``ConfigParser`` syntax.
There are two required sections.
One is a ``[variable_args]`` section that contains a list of parameters that we will vary to obtain a posterior distribution and the other is ``[static_args]`` section that contains a list of parameters are held fixed through out the run.

Each parameter in ``[variable_args]`` must have a subsection in ``[prior]``.
To create a subsection use the ``-`` char, e.g. for one of the mass parameters do ``[prior-mass1]``.

Each prior subsection must have a ``name`` option that identifies what prior to use.
These distributions are described in :py:mod:`pycbc.distributions`.
A list of all distributions that can be used is found with

.. literalinclude:: ../examples/distributions/list_distributions.py
.. command-output:: python ../examples/distributions/list_distributions.py

One or more of the ``variable_args`` may be transformed to a different parameter space for purposes of sampling. This is done by specifying a ``[sampling_parameters]`` section. This section specifies which ``variable_args`` to replace with which parameters for sampling. This must be followed by one or more ``[sampling_transforms-{sampling_params}]`` sections that provide the transform class to use. For example, the following would cause the sampler to sample in chirp mass (``mchirp``) and mass ratio (``q``) instead of ``mass1`` and ``mass2``::

    [sampling_parameters]
    mass1, mass2: mchirp, q

    [sampling_transforms-mchirp+q]
    name = mass1_mass2_to_mchirp_q

For a list of all possible transforms see :py:mod:`pycbc.transforms`.

There can be any number of ``variable_args`` with any name. No parameter name is special (with the exception of parameters that start with ``calib_``; see below). However, in order to generate waveforms, certain parameters must be provided for waveform generation. If you would like to specify a ``variable_arg`` that is not one of these parameters, then you must provide a ``[waveforms_transforms-{param}]`` section that provides a transform from the arbitrary ``variable_args`` to the needed waveform parameter(s) ``{param}``. For example, in the following we provide a prior on ``chirp_distance``. Since ``distance``, not ``chirp_distance``, is recognized by the CBC waveforms module, we provide a transform to go from ``chirp_distance`` to ``distance``::

    [variable_args]
    chirp_distance =

    [prior-chirp_distance]
    name = uniform
    min-chirp_distance = 1
    max-chirp_distance = 200

    [waveform_transforms-distance]
    name = chirp_distance_to_distance

Any class in the transforms module may be used. A useful transform for these purposes is the :py:class:`pycbc.transforms.CustomTransform`, which allows for arbitrary transforms using any function in the :py:mod:`pycbc.conversions`, :py:mod:`pycbc.coordinates`, or :py:mod:`pycbc.cosmology` modules, along with numpy math functions. For example, the following would use the I-Love-Q relationship :py:meth:`pycbc.conversions.dquadmon_from_lambda` to relate the quadrupole moment of a neutron star ``dquad_mon1`` to its tidal deformation ``lambda1``::

    [variable_args]
    lambda1 =

    [waveform_transforms-dquad_mon1]
    name = custom
    inputs = lambda1
    dquad_mon1 = dquadmon_from_lambda(lambda1)

A list of all parameters that are understood by the CBC waveform generator can be found with:

.. literalinclude:: ../examples/inference/list_parameters.py
.. command-output:: python ../examples/inference/list_parameters.py

Some common transforms are pre-defined in the code. These are: the mass parameters ``mass1`` and ``mass2`` can be substituted with ``mchirp`` and ``eta`` or ``mchirp`` and ``q``.
The component spin parameters ``spin1x``, ``spin1y``, and ``spin1z`` can be substituted for polar coordinates ``spin1_a``, ``spin1_azimuthal``, and ``spin1_polar`` (ditto for ``spin2``).

If any calibration parameters are used (prefix ``calib_``), a ``[calibration]`` section must be included. This section must have a ``name`` option that identifies what calibration model to use. The models are described in :py:mod:`pycbc.calibration`. The ``[calibration]`` section must also include reference values ``fc0``, ``fs0``, and ``qinv0``, as well as paths to ASCII transfer function files for the test mass actuation, penultimate mass actuation, sensing function, and digital filter for each IFO being used in the analysis. E.g. for an analysis using H1 only, the required options would be ``h1-fc0``, ``h1-fs0``, ``h1-qinv0``, ``h1-transfer-function-a-tst``, ``h1-transfer-function-a-pu``, ``h1-transfer-function-c``, ``h1-transfer-function-d``.

Simple examples are given in the subsections below.

-----------------------------------
Running on an analytic distribution
-----------------------------------

Several analytic distributions are available to run tests on. These can be run quickly on a laptop to check that a sampler is working properly.

This example demonstrates how to sample a 2D normal distribution with the ``emcee`` sampler. First, create the following configuration file (named ``normal2d.ini``)::

    [variable_args]
    x =
    y =

    [prior-x]
    name = uniform
    min-x = -10
    max-x = 10

    [prior-y]
    name = uniform
    min-y = -10
    max-y = 10

Then run::

    pycbc_inference --verbose \
        --config-files normal2d.ini \
        --output-file normal2d.hdf \
        --sampler emcee \
        --niterations 100 \
        --nwalkers 5000 \
        --likelihood-evaluator test_normal

This will run the ``emcee`` sampler on the 2D analytic normal distribution with 5000 walkers for 100 iterations.

To plot the posterior distribution after the last iteration, run::

    pycbc_inference_plot_posterior --verbose \
        --input-file normal2d.hdf \
        --output-file posterior-normal2d.png \
        --plot-scatter \
        --plot-contours \
        --plot-marginal \
        --z-arg loglr \
        --iteration -1

This will plot each walker's position as a single point colored by the log likelihood ratio at that point, with the 50th and 90th percentile contours drawn. See below for more information about using ``pycbc_inference_plot_posterior``.

To make a movie showing how the walkers evolved, run::

    pycbc_inference_plot_movie --verbose \
        --input-file normal2d.hdf \
        --output-prefix frames-normal2d \
        --movie-file normal2d_mcmc_evolution.mp4 \
        --cleanup \
        --plot-scatter \
        --plot-contours \
        --plot-marginal \
        --z-arg loglr \
        --frame-step 1

**Note:** you need ``ffmpeg`` installed for the mp4 to be created. See below for more information on using ``pycbc_inference_plot_movie``.

The number of dimensions of the distribution is set by the number of ``variable_args`` in the configuration file. The names of the ``variable_args`` do not matter, just that the prior sections use the same names (in this example ``x`` and ``y`` were used, but ``foo`` and ``bar`` would be equally valid). A higher (or lower) dimensional distribution can be tested by simply adding more (or less) ``variable_args``.

Which analytic distribution is used is set by the ``--likelihood-evaluator`` option. By setting to ``test_normal`` we used :py:class:`pycbc.inference.likelihood.TestNormal`. To see the list of available likelihood classes run ``pycbc_inference --help``; any choice for ``--likelihood-evaluator`` that starts with ``test_`` is analytic. The other analytic distributions available are: :py:class:`pycbc.inference.likelihood.TestEggbox`, :py:class:`pycbc.inference.likelihood.TestRosenbrock`, and :py:class:`pycbc.inference.likelihood.TestVolcano`. As with ``test_normal``, the dimensionality of these test distributions is set by the number of ``variable_args`` in the configuration file. The ``test_volcano`` distribution must be two dimensional, but all of the other distributions can have any number of dimensions. The configuration file syntax for the other test distributions is the same as in this example. Indeed, with this configuration file one only needs to change the ``--likelihood-evaluator`` argument to try (2D versions of) the other distributions.

------------------------------
BBH software injection example
------------------------------

This example recovers the parameters of a precessing binary black-hole (BBH).

An example configuration file (named ``inference.ini``) is::

    [variable_args]
    ; waveform parameters that will vary in MCMC
    tc =
    mass1 =
    mass2 =
    spin1_a =
    spin1_azimuthal =
    spin1_polar =
    spin2_a =
    spin2_azimuthal =
    spin2_polar =
    distance =
    coa_phase =
    inclination =
    polarization =
    ra =
    dec =

    [static_args]
    ; waveform parameters that will not change in MCMC
    approximant = IMRPhenomPv2
    f_lower = 18
    f_ref = 20

    [prior-tc]
    ; coalescence time prior
    name = uniform
    min-tc = 1126259461.8
    max-tc = 1126259462.2

    [prior-mass1]
    name = uniform
    min-mass1 = 10.
    max-mass1 = 80.

    [prior-mass2]
    name = uniform
    min-mass2 = 10.
    max-mass2 = 80.

    [prior-spin1_a]
    name = uniform
    min-spin1_a = 0.0
    max-spin1_a = 0.99

    [prior-spin1_polar+spin1_azimuthal]
    name = uniform_solidangle
    polar-angle = spin1_polar
    azimuthal-angle = spin1_azimuthal

    [prior-spin2_a]
    name = uniform
    min-spin2_a = 0.0
    max-spin2_a = 0.99

    [prior-spin2_polar+spin2_azimuthal]
    name = uniform_solidangle
    polar-angle = spin2_polar
    azimuthal-angle = spin2_azimuthal

    [prior-distance]
    ; following gives a uniform volume prior
    name = uniform_radius
    min-distance = 10
    max-distance = 1000

    [prior-coa_phase]
    ; coalescence phase prior
    name = uniform_angle

    [prior-inclination]
    ; inclination prior
    name = sin_angle

    [prior-ra+dec]
    ; sky position prior
    name = uniform_sky

    [prior-polarization]
    ; polarization prior
    name = uniform_angle

    ;
    ;   Sampling transforms
    ;
    [sampling_parameters]
    ; parameters on the left will be sampled in
    ; parametes on the right
    mass1, mass2 : mchirp, q

    [sampling_transforms-mchirp+q]
    ; inputs mass1, mass2
    ; outputs mchirp, q
    name = mass1_mass2_to_mchirp_q

An example of generating an injection::

    # define waveform parameters
    TRIGGER_TIME=1126259462.0
    INJ_APPROX=SEOBNRv2threePointFivePN
    MASS1=37.
    MASS2=32.
    RA=2.21535724066
    DEC=-1.23649695537
    INC=2.5
    COA_PHASE=1.5
    POLARIZATION=1.75
    DISTANCE=100000 # in kpc
    INJ_F_MIN=18.
    TAPER="start"

    # path of injection file that will be created in the example
    INJ_PATH=injection.xml.gz

    # lalapps_inspinj requires degrees on the command line
    LONGITUDE=`python -c "import numpy; print ${RA} * 180/numpy.pi"`
    LATITUDE=`python -c "import numpy; print ${DEC} * 180/numpy.pi"`
    INC=`python -c "import numpy; print ${INC} * 180/numpy.pi"`
    POLARIZATION=`python -c "import numpy; print ${POLARIZATION} * 180/numpy.pi"`
    COA_PHASE=`python -c "import numpy; print ${COA_PHASE} * 180/numpy.pi"`

    # create injection file
    lalapps_inspinj \
        --output ${INJ_PATH} \
        --seed 1000 \
        --f-lower ${INJ_F_MIN} \
        --waveform ${INJ_APPROX} \
        --amp-order 7 \
        --gps-start-time ${TRIGGER_TIME} \
        --gps-end-time ${TRIGGER_TIME} \
        --time-step 1 \
        --t-distr fixed \
        --l-distr fixed \
        --longitude ${LONGITUDE} \
        --latitude ${LATITUDE} \
        --d-distr uniform \
        --min-distance ${DISTANCE} \
        --max-distance ${DISTANCE} \
        --i-distr fixed \
        --fixed-inc ${INC} \
        --coa-phase-distr fixed \
        --fixed-coa-phase ${COA_PHASE} \
        --polarization ${POLARIZATION} \
        --m-distr fixMasses \
        --fixed-mass1 ${MASS1} \
        --fixed-mass2 ${MASS2} \
        --taper-injection ${TAPER} \
        --disable-spin

An example of running ``pycbc_inference`` to analyze the injection in fake data::

    # injection parameters
    TRIGGER_TIME=1126259462.0
    INJ_PATH=injection.xml.gz

    # sampler parameters
    CONFIG_PATH=inference.ini
    OUTPUT_PATH=inference.hdf
    SEGLEN=8
    PSD_INVERSE_LENGTH=4
    IFOS="H1 L1"
    STRAIN="H1:aLIGOZeroDetHighPower L1:aLIGOZeroDetHighPower"
    SAMPLE_RATE=2048
    F_MIN=20
    N_UPDATE=500
    N_WALKERS=5000
    N_SAMPLES=5000
    N_CHECKPOINT=1000
    PROCESSING_SCHEME=cpu

    # the following sets the number of cores to use; adjust as needed to
    # your computer's capabilities
    NPROCS=12

    # get coalescence time as an integer
    TRIGGER_TIME_INT=${TRIGGER_TIME%.*}

    # start and end time of data to read in
    GPS_START_TIME=$((${TRIGGER_TIME_INT} - ${SEGLEN}))
    GPS_END_TIME=$((${TRIGGER_TIME_INT} + ${SEGLEN}))

    # run sampler
    # specifies the number of threads for OpenMP
    # Running with OMP_NUM_THREADS=1 stops lalsimulation
    # to spawn multiple jobs that would otherwise be used
    # by pycbc_inference and cause a reduced runtime.
    OMP_NUM_THREADS=1 \
    pycbc_inference --verbose \
        --seed 12 \
        --instruments ${IFOS} \
        --gps-start-time ${GPS_START_TIME} \
        --gps-end-time ${GPS_END_TIME} \
        --psd-model ${STRAIN} \
        --psd-inverse-length ${PSD_INVERSE_LENGTH} \
        --fake-strain ${STRAIN} \
        --fake-strain-seed 44 \
        --sample-rate ${SAMPLE_RATE} \
        --low-frequency-cutoff ${F_MIN} \
        --channel-name H1:FOOBAR L1:FOOBAR \
        --injection-file ${INJ_PATH} \
        --config-file ${CONFIG_PATH} \
        --output-file ${OUTPUT_PATH} \
        --processing-scheme ${PROCESSING_SCHEME} \
        --sampler kombine \
        --burn-in-function max_posterior \
        --update-interval ${N_UPDATE} \
        --likelihood-evaluator gaussian \
        --nwalkers ${N_WALKERS} \
        --n-independent-samples ${N_SAMPLES} \
        --checkpoint-interval ${N_CHECKPOINT} \
        --nprocesses ${NPROCS} \
        --save-strain \
        --save-psd \
        --save-stilde \
        --force

----------------
GW150914 example
----------------

With a minor change to the ``tc`` prior, you can reuse ``inference.ini`` from the previous example to analyze the data containing GW150914. Change the ``[prior-tc]`` section to::

    [prior-tc]
    ; coalescence time prior
    name = uniform
    min-tc = 1126259462.32
    max-tc = 1126259462.52

Next, you need to obtain the real LIGO data containing GW150914. Do one of
the following:

* **If you are a LIGO member and are running on a LIGO Data Grid cluster:**
  you can use the LIGO data server to automatically obtain the frame files.
  Simply set the following environment variables::

    FRAMES="--frame-type H1:H1_HOFT_C02 L1:L1_HOFT_C02"
    CHANNELS="H1:H1:DCS-CALIB_STRAIN_C02 L1:L1:DCS-CALIB_STRAIN_C02"

* **If you are not a LIGO member, or are not running on a LIGO Data Grid
  cluster:** you need to obtain the data from the
  `LIGO Open Science Center <https://losc.ligo.org>`_. First run the following
  commands to download the needed frame files to your working directory::

    wget https://losc.ligo.org/s/events/GW150914/H-H1_LOSC_4_V2-1126257414-4096.gwf
    wget https://losc.ligo.org/s/events/GW150914/L-L1_LOSC_4_V2-1126257414-4096.gwf

  Then set the following enviornment variables::

    FRAMES="--frame-files H1:H-H1_LOSC_4_V2-1126257414-4096.gwf L1:L-L1_LOSC_4_V2-1126257414-4096.gwf"
    CHANNELS="H1:LOSC-STRAIN L1:LOSC-STRAIN"

Now run::

    # trigger parameters
    TRIGGER_TIME=1126259462.42

    # data to use
    # the longest waveform covered by the prior must fit in these times
    SEARCH_BEFORE=6
    SEARCH_AFTER=2

    # use an extra number of seconds of data in addition to the data specified
    PAD_DATA=8

    # PSD estimation options
    PSD_ESTIMATION="H1:median L1:median"
    PSD_INVLEN=4
    PSD_SEG_LEN=16
    PSD_STRIDE=8
    PSD_DATA_LEN=1024

    # sampler parameters
    CONFIG_PATH=inference.ini
    OUTPUT_PATH=inference.hdf
    IFOS="H1 L1"
    SAMPLE_RATE=2048
    F_HIGHPASS=15
    F_MIN=20
    N_UPDATE=500
    N_WALKERS=5000
    N_SAMPLES=5000
    N_CHECKPOINT=1000
    PROCESSING_SCHEME=cpu

    # the following sets the number of cores to use; adjust as needed to
    # your computer's capabilities
    NPROCS=12

    # get coalescence time as an integer
    TRIGGER_TIME_INT=${TRIGGER_TIME%.*}

    # start and end time of data to read in
    GPS_START_TIME=$((${TRIGGER_TIME_INT} - ${SEARCH_BEFORE} - ${PSD_INVLEN}))
    GPS_END_TIME=$((${TRIGGER_TIME_INT} + ${SEARCH_AFTER} + ${PSD_INVLEN}))

    # start and end time of data to read in for PSD estimation
    PSD_START_TIME=$((${TRIGGER_TIME_INT} - ${PSD_DATA_LEN}/2))
    PSD_END_TIME=$((${TRIGGER_TIME_INT} + ${PSD_DATA_LEN}/2))

    # run sampler
    # specifies the number of threads for OpenMP
    # Running with OMP_NUM_THREADS=1 stops lalsimulation
    # to spawn multiple jobs that would otherwise be used
    # by pycbc_inference and cause a reduced runtime.
    OMP_NUM_THREADS=1 \
    pycbc_inference --verbose \
        --seed 12 \
        --instruments ${IFOS} \
        --gps-start-time ${GPS_START_TIME} \
        --gps-end-time ${GPS_END_TIME} \
        --channel-name ${CHANNELS} \
        ${FRAMES} \
        --strain-high-pass ${F_HIGHPASS} \
        --pad-data ${PAD_DATA} \
        --psd-estimation ${PSD_ESTIMATION} \
        --psd-start-time ${PSD_START_TIME} \
        --psd-end-time ${PSD_END_TIME} \
        --psd-segment-length ${PSD_SEG_LEN} \
        --psd-segment-stride ${PSD_STRIDE} \
        --psd-inverse-length ${PSD_INVLEN} \
        --sample-rate ${SAMPLE_RATE} \
        --low-frequency-cutoff ${F_MIN} \
        --config-file ${CONFIG_PATH} \
        --output-file ${OUTPUT_PATH} \
        --processing-scheme ${PROCESSING_SCHEME} \
        --sampler kombine \
        --burn-in-function max_posterior \
        --update-interval ${N_UPDATE} \
        --likelihood-evaluator gaussian \
        --nwalkers ${N_WALKERS} \
        --n-independent-samples ${N_SAMPLES} \
        --checkpoint-interval ${N_CHECKPOINT} \
        --nprocesses ${NPROCS} \
        --save-strain \
        --save-psd \
        --save-stilde \
        --force

----------------------------------------------------
HDF output file handler (``pycbc.io.InferenceFile``)
----------------------------------------------------

The executable ``pycbc_inference`` will write a HDF file with all the samples from each walker along with the PSDs and some meta-data about the sampler.
There is a handler class ``pycbc.io.InferenceFile`` that extends ``h5py.File``.
To read the output file you can do::

    from pycbc.io import InferenceFile
    fp = InferenceFile("cbc_example-n1e4.hdf", "r")

To get all samples for ``distance`` from the first walker you can do::

    samples = fp.read_samples("distance", walkers=0)
    print samples.distance

The function ``InferenceFile.read_samples`` includes the options to thin the samples.
By default the function will return samples beginning at the end of the burn-in to the last written sample, and will use the autocorrelation length (ACL) calculated by ``pycbc_inference`` to select the indepdedent samples.
You can supply ``thin_start``, ``thin_end``, and ``thin_interval`` to override this. To read all samples you would do::

    samples = fp.read_samples("distance", walkers=0, thin_start=0, thin_end=-1, thin_interval=1)
    print samples.distance

Some standard parameters that are derived from the variable arguments (listed via ``fp.variable_args``) can also be retrieved. For example, if ``fp.variable_args`` includes ``mass1`` and ``mass2``, then you can retrieve the chirp mass with::

   samples = fp.read_samples("mchirp")
   print samples.mchirp

In this case, ``fp.read_samples`` will retrieve ``mass1`` and ``mass2`` (since they are needed to compute chirp mass); ``samples.mchirp`` then returns an array of the chirp mass computed from ``mass1`` and ``mass2``.

For more information, including the list of predefined derived parameters, see :py:class:`pycbc.io.InferenceFile`.

===============================================
Visualizing the Posteriors
===============================================

.. toctree::
   :maxdepth: 1

   inference/viz.rst

===============================================
Workflows (``pycbc_make_inference_workflow``)
=============================================== 

.. toctree::
   :maxdepth: 1

   workflow/pycbc_make_inference_workflow
