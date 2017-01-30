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

=================================================
Sampling the parameter space: ``pycbc_inference``
=================================================

--------
Overview
--------

The executable ``pycbc_inference`` is designed to sample the parameter space
and save the samples in an HDF file. A high-level description of
``pycbc_inference`` is::

#. Estimate a PSD from a model or data.

#. Read gravitational-wave strain from a gravitational-wave model or use recolored fake strain.

#. Read priors from configuration file.

#. Run sampler

The user specifies the sampler on the command line with the ``--sampler``
option. A list of available samplers is::

    <code>

The user specifies the likelihood model on the command line with
the ``--likelihood`` option. At the moment there is only a single
choice ``--likelihood gaussian``.

The user specifies a configuration file that defines the priors with the
``--config-files`` option. A description of the configuration file is given
in the subsection below.

The user 

------------------------------
BBH software injection example
------------------------------

This example recovers the parameters of a non-spinning binary black-hole (BBH)
software injection in fake data.

.. code-block:: bash

    # define coalescence time, observed masses, and waveform parameters
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

    # sampler parameters
    OUTPUT=cbc_example-n1e4.hdf
    SEGLEN=8
    PSD_INVERSE_LENGTH=4
    IFOS="H1 L1"
    STRAIN="H1:aLIGOZeroDetHighPower L1:aLIGOZeroDetHighPower"
    SAMPLE_RATE=2048
    F_MIN=30.
    N_WALKERS=5000
    N_ITERATIONS=1000
    N_CHECKPOINT=100
    PROCESSING_SCHEME=cpu
    NPROCS=12
    CONFIG_PATH=inference.ini

    # get coalescence time as an integer
    TRIGGER_TIME_INT=${TRIGGER_TIME%.*}

    # start and end time of data to read in
    GPS_START_TIME=$((${TRIGGER_TIME_INT} - ${SEGLEN}))
    GPS_END_TIME=$((${TRIGGER_TIME_INT} + ${SEGLEN}))

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

    # run sampler
    # specifies the number of threads for OpenMP
    # Running with OMP_NUM_THREADS=1 stops lalsimulation
    # to spawn multiple jobs that would otherwise be used
    # by pycbc_inference and cause a reduced runtime.
    OMP_NUM_THREADS=1 \
    pycbc_inference --verbose \
        --instruments ${IFOS} \
        --gps-start-time ${GPS_START_TIME} \
        --gps-end-time ${GPS_END_TIME} \
        --psd-model ${STRAIN} \
        --psd-inverse-length ${PSD_INVERSE_LENGTH} \
        --fake-strain ${STRAIN} \
        --sample-rate ${SAMPLE_RATE} \
        --low-frequency-cutoff ${F_MIN} \
        --channel-name H1:FOOBAR L1:FOOBAR \
        --injection-file ${INJ_PATH} \
        --processing-scheme ${PROCESSING_SCHEME} \
        --sampler kombine \
        --likelihood-evaluator gaussian \
        --nwalkers ${N_WALKERS} \
        --niterations ${N_ITERATIONS} \
        --config-file ${CONFIG_PATH} \
        --output-file ${OUTPUT} \
        --checkpoint-interval ${N_CHECKPOINT} \
        --nprocesses ${NPROCS}

An example configuration file (named ``inference.ini`` above) is::

    [variable_args]
    ; waveform parameters that will vary in MCMC
    tc =
    mass1 =
    mass2 =
    distance =
    coa_phase =
    inclination =
    polarization =
    ra =
    dec =

    [static_args]
    ; waveform parameters that will not change in MCMC
    approximant = SEOBNRv2_ROM_DoubleSpin
    f_lower = 28.0

    [prior-tc]
    ; coalescence time prior
    name = uniform
    min-tc = 1126259461.8
    max-tc= 1126259462.2

    [prior-mass1]
    ; component mass prior
    name = uniform
    min-mass1 = 10.
    max-mass1 = 80.

    [prior-mass2]
    ; component mass prior
    name = uniform
    min-mass2 = 10.
    max-mass2 = 80.

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
    name = uniform_angle
    min-inclination = 0
    max-inclination = 1

    [prior-ra+dec]
    ; sky position prior
    name = uniform_sky

    [prior-polarization]
    ; polarization prior
    name = uniform_angle


---------------------------------------------------
HDF output file handler: ``pycbc.io.InferenceFile``
---------------------------------------------------

The executable ``pycbc_inference`` will write a HDF file with all the samples from each walker along with the PSDs and some meta-data about the sampler. There is a handler class ``pycbc.io.InferenceFile`` that extends ``h5py.File``. To read the output file you can do::

    from pycbc.io import InferenceFile
    fp = InferenceFile("cbc_example-n1e4.hdf.hdf", "r")

To get all samples for ``mass1`` from the first walker you can do::

    samples = fp.read_samples("mass1", walkers=0)
    print samples.mass1

The function ``InferenceFile.read_samples`` includes the options to thin the samples. By default the function will return samples beginning at the end of the burn-in to the last written sample, and will use the autocorrelation length (ACL) calcualted by ``pycbc_inference`` to select the indepdedent samples. You can supply ``thin_start``, ``thin_end``, and ``thin_interval`` to override this. To read all samples you would do::

    samples = fp.read_samples("mass1", walkers=0, thin_start=0, thin_end=-1, thin_interval=1)
    print samples.mass1

Some standard parameters that are derived from the variable arguments (listed via ``fp.variable_args``) can also be retrieved. For example, if ``fp.variable_args`` includes ``mass1`` and ``mass2``, then you can retrieve the chirp mass with::

   samples = samples = fp.read_samples("mchirp")
   print samples.mchirp

In this case, ``fp.read_samples`` will retrieve ``mass1`` and ``mass2`` (since they are needed to compute chirp mass); ``samples.mchirp`` then returns an array of the chirp mass computed from ``mass1`` and ``mass2``.

For more information, including the list of predefined derived parameters, see the docstring of ``pycbc.io.InferenceFile``.
