# Marginalized total mass example

This example analyzes a simulated binary black hole signal with the
`marginalized_mtotal` model, which brute force marginalizes the likelihood
over a grid of total masses.

Rather than generating a waveform at every point of that grid, the model
generates a single reference waveform at `fiducial_mtotal` and rescales it,
using the fact that at fixed mass ratio, spins and distance a waveform depends
on the total mass only through `t / M` with an amplitude proportional to `M`:

    h(f; M) = kappa**2 * h(kappa * f; M_ref),   kappa = M / M_ref

The rescaling is applied to the generated waveform, so the base model has to
accept a pre-generated waveform; at present only `marginalized_time` does.
The base model also applies the detector response, so the sky location, time,
polarization and distance must be marginalized over or held static.

Run it with:

    bash run.sh

`pycbc_inference_model_stats --reconstruct-parameters` draws a total mass from
the grid in proportion to its likelihood, so `mtotal` appears in the
reconstructed posterior along with the other marginalized parameters.

## Checking the model

`validation.ini` runs the same analysis with the base `marginalized_time`
model sampling the total mass directly instead of marginalizing over it. The
two runs should give the same evidence and the same posteriors:

    bash validate.sh
