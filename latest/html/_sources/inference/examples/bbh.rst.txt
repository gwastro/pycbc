.. _inference_example_bbh:

.. |GaussianNoise| replace:: :py:class:`GaussianNoise <pycbc.inference.models.gaussian_noise.GaussianNoise>`
.. |MarginalizedPhaseGaussianNoise| replace:: :py:class:`MarginalizedPhaseGaussianNoise <pycbc.inference.models.marginalized_gaussian_noise.MarginalizedPhaseGaussianNoise>`

---------------------
Simulated BBH example
---------------------

This example recovers the parameters of a simulated binary black-hole (BBH)
that has similar parameters has GW150914.

^^^^^^^^^^^^^^^^^^^^^^^
1. Create the injection
^^^^^^^^^^^^^^^^^^^^^^^

First, we need to create an ``injection.hdf`` file that specifies the
parameters of the simulated signal. To do that we will use
``pycbc_create_injection``. Like ``pycbc_inference``,
``pycbc_create_injections`` uses a configuration file to set the parameters of
the injections it will create. To create a binary-black hole with parameters
similar to GW150914, use the following configuration file:

.. literalinclude:: ../../../examples/inference/bbh-injection/injection.ini
   :language: ini

:download:`Download <../../../examples/inference/bbh-injection/injection.ini>`

Note the similarity to the configuration file for ``pycbc_inference``: you must
have a ``[variable_params]`` section. If we wanted to randomize one or more
of the parameters, we would list them there, then add ``[prior]`` sections to
specify what distribution to draw the parameters from. In this case, however,
we want to fix the parameters, so we just put all of the necessary parameters
in the ``[static_params]`` section.

To create the injection file, run:

.. literalinclude:: ../../../examples/inference/bbh-injection/make_injection.sh
   :language: bash

:download:`Download <../../../examples/inference/bbh-injection/make_injection.sh>`

This will create the ``injection.hdf`` file, which we will give to
``pycbc_inference``. For more information on generating injection files, run
``pycbc_create_injections --help``.


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
2. Setup the configuration files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now we need to set up the configuration for ``pycbc_inference``. Since we
will be analyzing data, we will need to provide several additional options in a
``[data]`` section. To keep the configuration files easy to read, we will split
the data, sampler, and prior settings into their own configuration files.

Here are the model and prior settings we will use:

.. literalinclude:: ../../../examples/inference/priors/gw150914_like.ini
   :language: ini

:download:`Download <../../../examples/inference/priors/gw150914_like.ini>`

In the ``[model]`` section we have set the model to be ``gaussian_noise``.  As
described above, this is the standard model to use for CBC signals. It assumes
that the noise is wide-sense stationary Gaussian noise. Notice that we provide
a low frequency argument which is the lower bound for the likelihood integral.
(See the |GaussianNoise| docs for details.)

The |GaussianNoise| model will need to generate model waveforms in order to
evaluate the likelihood. This means that we need to provide it with a waveform
approximant to use. Which model to use is set by the ``approximant`` argument
in the ``[static_params]`` section. Here, we are using ``IMRPhenomPv2``. This
is a frequency-domain, precessing model that uses the dominant,
:math:`\ell=|m|=2` mode. For this reason, we are varying all three components
of each object's spin, along with the masses, location, orientation, and phase
of the signal.  We need to provide a lower frequency cutoff (``f_lower``),
which is the starting frequency of the waveform. This must be :math:`\leq` the
smallest low frequency cutoff set in the model section.

.. note::
   In this example we have to sample over a reference phase for the waveform
   (``coa_phase``). This is because we are using the |GaussianNoise| model. For
   dominant-mode only waveforms, it is possible to analytically marginalize
   over the phase using the |MarginalizedPhaseGaussianNoise| model. This can
   speed up the convergence of the sampler by a factor of 3 or faster. To use
   the marginalized phase model, change the model name to
   ``marginalized_phase``, and remove ``coa_phase`` from the
   ``variable_params`` and the prior. However, the marginalized phase model
   should not be used with fully precessing models or models that include
   higher modes. You can use it with ``IMRPhenomPv2`` due to some
   simplifications that the approximant makes.

We also need a prior for the coaslesence time ``tc``. We have done this by
setting a reference time in the ``static_params`` section, and varying a
+/-0.1s window around it with the ``delta_tc`` parameter. Notice that the
trigger time is not set to a value; instead, we reference the ``trigger-time``
option that is set in the ``[data]`` section. This way, we only need to set the
trigger time in one place; we can reuse this prior file for different BBH
events by simply providing a different data configuration file.



Here is the data configuration file we will use:

.. literalinclude:: ../../../examples/inference/bbh-injection/data.ini
   :language: ini

:download:`Download <../../../examples/inference/bbh-injection/data.ini>`

In this case, we are generating fake Gaussian noise (via the ``fake-strain``
option) that is colored by the `Advanced LIGO updated design sensitivity curve
<https://dcc.ligo.org/LIGO-T1800044/public>`_. (Note that this is ~3 times more
sensitive than what the LIGO detectors were when GW150914 was detected.) The
duration of data that will be analyzed is set by the
``analysis-(start|end)-time`` arguments. These values are measured with respect
to the ``trigger-time``. The analyzed data should be long enough such that it
encompasses the longest waveform admitted by our prior, plus our timing
uncertainty (which is determined by the prior on ``delta_tc``). Waveform
duration is approximately determined by the total mass of a system. The lowest
total mass (``= mass1 + mass2``) admitted by our prior is 20 solar masses. This
corresponds to a duration of ~6 seconds, so we start the analysis time 6
seconds before the trigger time. (See the :py:mod:`pycbc.waveform` module for
utilities to estimate waveform duration.) Since the trigger time is
approximately where we expect the merger to happen, we only need a small amount
of time afterward to account for the timing uncertainty and ringdown. Here, we
choose 2 seconds, which is a good safety margin.

We also have to provide arguments for estimating a PSD. Although we know the
exact shape of the PSD in this case, we will still estimate it from the
generated data, as this most closely resembles what you do with a real event.
To do this, we have set ``psd-estimation``  to ``median-mean`` and we have set
``psd-segment-length``, ``psd-segment-stride``, and ``psd-(start|end)-time``
(which are with respect to the trigger time). This means that a Welch-like
method will be used to estimate the PSD. Specifically, we will use 512s of data
centered on the trigger time to estimate the PSD. This time will be divided up
into 8s-long segments (the segment length) each overlapping by 4s (the segment
stride). The data in each segment will be transformed to the frequency domain.
Two median values will be determined in each frequency bin from across all
even/odd segments, then averaged to obtain the PSD.

The beginning and end of the analysis segment will be corrupted by the
convolution of the inverse PSD with the data. To limit the amount of time that
is corrupted, we set ``psd-inverse-length`` to ``8``. This limits the
corruption to at most the first and last four seconds of the data segment.  To
account for the corruption, ``psd-inverse-length/2`` seconds are
subtracted/added by the code from/to the analysis start/end times specified by
the user before the data are transformed to the frequency domain.
Consequently, our data will have a frequency resolution of :math:`\Delta f =
1/16\,` Hz. The 4s at the beginning/end of the segment are effectively ignored
since the waveform is contained entirely in the -6/+2s we set with the analysis
start/end time.

Finally, we will use the ``emcee_pt`` sampler with the following settings:

.. literalinclude:: ../../../examples/inference/samplers/emcee_pt-gw150914_like.ini
   :language: ini

:download:`Download <../../../examples/inference/samplers/emcee_pt-gw150914_like.ini>`

Here, we will use 200 walkers and 20 temperatures. We will checkpoint (i.e.,
dump results to file) every 2000 iterations. Since we have provided an
``effective-nsamples`` argument and a ``[sampler-burn_in]`` section,
``pycbc_inference`` will run until it has acquired 1000 independent samples
after burn-in, which is determined by a combination of the :py:meth:`nacl
<pycbc.inference.burn_in.MCMCBurnInTests.nacl>` and
:py:meth:`max_posterior
<pycbc.inference.burn_in.MCMCBurnInTests.max_posterior>` tests;
i.e., the sampler will be considered converged when both of these tests are
satisfied.

The number of independent samples is checked at each checkpoint: after dumping
the results, the burn-in test is applied and an autocorrelation length is
calculated. The number of independent samples is then
``nwalkers x (the number of iterations since burn in)/ACL``. If this number
exceeds ``effective-nsamples``, ``pycbc_inference`` will finalize the results
and exit.

^^^^^^
3. Run
^^^^^^

To perform the analysis, run:

.. literalinclude:: ../../../examples/inference/bbh-injection/run.sh
   :language: bash

:download:`Download <../../../examples/inference/bbh-injection/run.sh>`


Since we are generating waveforms and analyzing a 15 dimensional parameter
space, this run will be much more computationally expensive than the :ref:`analytic
example<inference_example_analytic>`. We recommend running this on a cluster or a computer with a
large number of cores. In the example, we have set the parallelization to use
10 cores. With these settings, it should checkpoint approximately every hour or
two. The run should complete in a few hours. If you would like to acquire more
samples, increase ``effective-nsamples``.

.. note::
   We have found that the ``emcee_pt`` sampler struggles to accumulate more
   than ~2000 independent samples with 200 walkers and 20 temps. The basic
   issue is that the ACT starts to grow at the same rate as new iterations,
   so that the number of independent samples remains the same. Increasing
   the number of walkers and decreasing the number of temperatures can help,
   but this sometimes leads to the sampler not fully converging before passing
   the burn in tests. If you want more than 2000 samples, we currently
   recommend doing multiple independent runs with different seed values, then
   combining posterior samples after they have finished.

