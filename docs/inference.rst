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
In addition to ``--sampler`` the user will need to specify the number of walkers to use ``--nwalkers``, the number of iterations to go until ``--niterations``, and for parallel-tempered samplers the number of temperatures ``--ntemps``.
If the sampler has a built-in burn-in function it will be used by default, otherwise you can skill the burn-in with ``--skip-burn-in`` or set a minimum number of iterations for burn-in with ``--min-burn-in``.

The user specifies the likelihood model on the command line with the ``--likelihood`` option.
At the moment there is only a single choice ``--likelihood gaussian`` that is described in :py:class:`pycbc.inference.likelihood.GaussianLikelihood`.

The user specifies a configuration file that defines the priors with the ``--config-files`` option.
The syntax of the configuration file is described in the subsection below.

-------------------------
Configuration file syntax
-------------------------

Configuration files follow the ``ConfigParser`` syntax.
There are two required sections.
One is a ``[variable_args]`` section that contains a list of varying parameters and the other is ``[static_args]`` section that contains a list of parameters that do not vary.

A list of all parameters that can be used is found with

.. literalinclude:: ../examples/inference/list_parameters.py
.. command-output:: python ../examples/inference/list_parameters.py

The mass parameters ``mass1`` and ``mass2`` can be substituted for ``mchirp`` and ``eta``, or ``mchirp`` and ``q``.
The component spin parameters ``spin1x``, ``spin1y``, and ``spin1z`` can be substituted for polar coordinates ``spin1_a``, ``spin1_azimuthal``, and ``spin1_polar``.

Each parameter in ``[variable_args]`` must have a subsection in ``[prior]``.
To create a subsection use the ``-`` char, eg. for chirp mass do ``[prior-mchirp]``.

Each prior subsection must have a ``name`` option that identifies what prior to use.
These distributions are described in :py:mod:`pycbc.distributions`.
A list of all distributions that can be used is found with

.. literalinclude:: ../examples/distributions/list_distributions.py
.. command-output:: python ../examples/distributions/list_distributions.py

A simple example is given in the subsection below.

------------------------------
BBH software injection example
------------------------------

This example recovers the parameters of a precessing binary black-hole (BBH).

An example configuration file (named ``inference.ini``) is::

    [variable_args]
    ; waveform parameters that will vary in MCMC
    tc =
    mchirp =
    q =
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
    f_lower = 19.0

    [prior-tc]
    ; coalescence time prior
    name = uniform
    min-tc = 1126259461.8
    max-tc = 1126259462.2

    [prior-mchirp]
    name = uniform
    min-mchirp = 7.
    max-mchirp = 40.

    [prior-q]
    name = uniform
    min-q = 1.
    max-q = 5.

    [prior-spin1_a]
    name = uniform
    min-spin1_a = 0.0
    max-spin1_a = 0.9

    [prior-spin1_azimuthal]
    name = uniform
    min-spin1_azimuthal = 0.
    max-spin1_azimuthal = 6.283185307179586

    [prior-spin1_polar]
    name = sin_angle

    [prior-spin2_a]
    name = uniform
    min-spin2_a = 0.0
    max-spin2_a = 0.9

    [prior-spin2_azimuthal]
    name = uniform
    min-spin2_azimuthal = 0.
    max-spin2_azimuthal = 6.283185307179586

    [prior-spin2_polar]
    name = sin_angle

    [prior-distance]
    ; distance prior
    name = uniform
    min-distance = 10
    max-distance = 500

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
    INJ_F_MIN=28.
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
    F_MIN=30.
    N_UPDATE=500
    N_WALKERS=5000
    N_ITERATIONS=12000
    N_CHECKPOINT=1000
    PROCESSING_SCHEME=cpu
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
        --skip-burn-in \
        --update-interval ${N_UPDATE} \
        --likelihood-evaluator gaussian \
        --nwalkers ${N_WALKERS} \
        --niterations ${N_ITERATIONS} \
        --checkpoint-interval ${N_CHECKPOINT} \
        --checkpoint-fast \
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

Then run::

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
    PSD_SEG_LEN=8
    PSD_STRIDE=4
    PSD_DATA_LEN=1024

    # frame type and channel
    FRAMES="H1:H1_HOFT_C02 L1:L1_HOFT_C02"
    CHANNELS="H1:H1:DCS-CALIB_STRAIN_C02 L1:L1:DCS-CALIB_STRAIN_C02"

    # sampler parameters
    CONFIG_PATH=inference.ini
    OUTPUT_PATH=inference.hdf
    SEGLEN=8
    IFOS="H1 L1"
    SAMPLE_RATE=2048
    F_HIGHPASS=20
    F_MIN=30.
    N_UPDATE=500
    N_WALKERS=5000
    N_ITERATIONS=12000
    N_CHECKPOINT=1000
    PROCESSING_SCHEME=cpu
    NPROCS=12

    # get coalescence time as an integer
    TRIGGER_TIME_INT=${TRIGGER_TIME%.*}

    # start and end time of data to read in
    GPS_START_TIME=$((${TRIGGER_TIME_INT} - ${SEARCH_BEFORE} - ${PSD_INVLEN}))
    GPS_END_TIME=$((${TRIGGER_TIME_INT} + ${SEARCH_AFTER} + ${PSD_INVLEN}))

    # start and end time of data to read in for PSD estimation
    PSD_START_TIME=$((${GPS_START_TIME} - ${PSD_DATA_LEN}/2))
    PSD_END_TIME=$((${GPS_END_TIME} + ${PSD_DATA_LEN}/2))

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
        --frame-type ${FRAMES} \
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
        --skip-burn-in \
        --update-interval ${N_UPDATE} \
        --likelihood-evaluator gaussian \
        --nwalkers ${N_WALKERS} \
        --niterations ${N_ITERATIONS} \
        --checkpoint-interval ${N_CHECKPOINT} \
        --checkpoint-fast \
        --nprocesses ${NPROCS} \
        --save-strain \
        --save-psd \
        --save-stilde \
        --force

To get data we used ``--frame-type`` which will query the LIGO Data
Replicator (LDR) server to locate the frame files for us. You can also
use ``--frame-files`` or ``--frame-cache`` if you have a list or LAL cache
file of frames you wish to use.

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

=============================================================
Plotting the posteriors (``pycbc_inference_plot_posteriors``)
=============================================================

--------
Overview
--------

There is an executable that can plot the posteriors called
``pycbc_inference_plot_posteriors``. You can use ``--plot-scatter``
to plot a each sample as a point or ``--plot-density`` to plot a density map.

By default the plotting executables will plot all the parameters in the input
file. In order to specific a different set of variables to plot use the
``--parameters`` option. Examples how to use this option are shown below.

By default the plotting executables will plot samples beginning at the end of
the burn in. If the burn-in was skipped, then it starts from the first sample.
It will then use a sample every autocorrelation length along the chain.
Examples how to plot a specific iteration or change how the thinning is
performed is shown in the examples below.

You may plot a z-axis on the 2-D histograms using the ``--z-arg`` option.
For a list of options use ``pycbc_inference_plot_posterior --help``.

-----------------------------
Plotting a specific iteration
-----------------------------

An example of plotting the posteriors at a specific iteration::

    ITER=4999
    INPUT_FILE=inference.hdf
    OUTPUT_FILE=scatter.png
    pycbc_inference_plot_posterior \
        --iteration ${ITER} \
        --input-file ${INPUT_FILE} \
        --output-file ${OUTPUT_FILE} \
        --plot-scatter \
        --plot-marginal \
        --z-arg logplr \
        --parameters "ra*12/pi:$\alpha$ (h)" \
                     "dec*180/pi:$\delta$ (deg)" \
                     "polarization*180/pi:$\psi$ (deg)" \
                     mchirp q spin1_a spin1_azimuthal spin1_polar \
                     spin2_a spin2_azimuthal spin2_polar \
                     "inclination*180/pi:$\iota$ (deg)" distance \
                     "coa_phase*180/pi:$\phi_0$ (deg)" tc

-----------------------------------
Plotting a thinned chain of samples
-----------------------------------

There are also options for thinning the chains of samples from the command line, an example starting at the 6000-th iteration and taking every 2000-th iteration until the 12000-th iteration::

    THIN_START=5999
    THIN_INTERVAL=2000
    THIN_END=11999
    INPUT_FILE=inference.hdf
    OUTPUT_FILE=scatter.png
    pycbc_inference_plot_posterior \
        --input-file ${INPUT_FILE} \
        --output-file ${OUTPUT_FILE} \
        --plot-scatter \
        --thin-start ${THIN_START} \
        --thin-interval ${THIN_INTERVAL} \
        --thin-end ${THIN_END} \
        --plot-marginal \
        --z-arg logplr \
        --parameters "ra*12/pi:$\alpha$ (h)" \
                     "dec*180/pi:$\delta$ (deg)" \
                     "polarization*180/pi:$\psi$ (deg)" \
                     mchirp q spin1_a spin1_azimuthal spin1_polar \
                     spin2_a spin2_azimuthal spin2_polar \
                     "inclination*180/pi:$\iota$ (deg)" distance \
                     "coa_phase*180/pi:$\phi_0$ (deg)" tc
