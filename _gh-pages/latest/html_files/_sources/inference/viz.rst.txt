=============================================================
Plotting the posteriors (``pycbc_inference_plot_posterior``)
=============================================================

--------
Overview
--------

There is an executable that can plot the posteriors called
``pycbc_inference_plot_posterior``. You can use ``--plot-scatter``
to plot a each sample as a point or ``--plot-density`` to plot a density map.

By default the plotting executables will plot all the parameters in the input
file. In order to specify a different set of variables to plot, use the
``--parameters`` option. Examples for how to use this option are shown below.

By default the plotting executables will plot samples beginning at the end of
the burn in. If the burn-in was skipped, then it starts from the first sample.
It will then use a sample every autocorrelation length along the chain.
Examples on how to plot a specific iteration or change how the thinning is
performed are shown in the examples below.

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
                     mass1 mass2 spin1_a spin1_azimuthal spin1_polar \
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
                     mass1 mass2 spin1_a spin1_azimuthal spin1_polar \
                     spin2_a spin2_azimuthal spin2_polar \
                     "inclination*180/pi:$\iota$ (deg)" distance \
                     "coa_phase*180/pi:$\phi_0$ (deg)" tc
           
.. _inference_make_movie: 
         
===============================================
Making a movie (``pycbc_inference_plot_movie``)
===============================================       

``pycbc_inference_plot_movie`` is an executable similar to ``_plot_posterior`` that allows you to combine plots in to a small movie. Most options for ``_plot_movie`` are the same as ``_plot_posterior`` with a few differences. Again, the plotting executables will plot all the parameters in the input file unless the option ``--parameters`` is used to specify a set of parameters that you want to see. An example plotting every 20-th iteration into a directory called "movies"::

    INPUT_FILE=inference.hdf
    START_SAMPLE=1
    END_SAMPLE=12000
    FRAME_STEP=20
    OUTPUT_PREFIX=frame
    NPROCESSES=10
    MOVIE_FILE=~/src/pycbc/movies/movie.mp4
    DPI=100
    
    pycbc_inference_plot_movie \  
        --input-file ${INPUT_FILE} \
        --start-sample ${START_SAMPLE} \
        --end-sample ${END_SAMPLE} \
        --frame-step ${FRAME_STEP} \
        --output-prefix ${OUTPUT_PREFIX} \
        --nprocesses ${NPROCESSES} \
        --movie-file ${MOVIE_FILE} \
        --cleanup \
        --plot-scatter \
        --plot-marginal \
        --z-arg snr \
        --dpi ${DPI} \
        --parameters mass1 mass2 spin1_a spin1_azimuthal spin1_polar \
               	     spin2_a spin2_azimuthal spin2_polar \
            
This will create a 24-second movie for a selection of parameters. The option ``--cleanup`` deletes the individual frame files prefixed as specified by the variable ``OUTPUT_PREFIX``. This is optional. 
For a list of options use ``pycbc_inference_plot_movie --help``.

