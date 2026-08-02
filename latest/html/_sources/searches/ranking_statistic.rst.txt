####################################################################################
The ranking statistic used in PyCBC searches
####################################################################################

How triggers are ranked is defined by the ranking-statistic, sngl-ranking, statistic-features and statistic-keywords options.
- ``sngl-ranking`` = The ranking used for single-detector triggers, this is generally a re-weighting of the SNR.
- ``ranking-statistic`` = How the triggers from a set of detectors are ranked in order to calculate significance. This will take the form of an snr-like combination (``quadsum``, ``phasetd``, ``exp_fit_csnr``), or a log-rates-like statistic, ``exp_fit``. See Ranking Statistic table below for the options.
- ``statistic-features`` = If using ranking-statistic ``exp_fit``, then these are the features to add or subtract from the ranking statistic. These are described in the Statistic Features table below.
- ``statistic-keywords`` = Some statistics require keywords to modify the behavior of the statistic in certain situations. Keywords affecting the sngl-ranking calculation are also given here, starting with ``sngl_ranking_``. These are described in the Statistic Keywords table below.
- ``statistic-files`` = Files to be used in the statistic calculation, of particular note here are the files needed for DQ and KDE reranking.

.. list-table:: Ranking Statistic
   :widths: 25 75
   :header-rows: 1

   * - Statistic name
     - Description
   * - ``quadsum``
     - The quadrature sum of triggers in each detector in the triggered network. ``sngl_ranking_only`` can also be given and is exactly equivalent.
   * - ``phasetd``
     - The same as ``quadsum``, but reweighted by the coincident parameters.
   * - ``exp_fit_csnr``
     - This is a reworking of the exponential fit designed to resemble network SNR. Uses a monotonic function of the negative log noise rate density which approximates combined sngl-ranking for coincs with similar newsnr in each ifo
   * - ``exp_fit``
     - The ratio of signal-to-noise rates in the triggered network of detectors. The trigger density at a given sngl-ranking is approximated for each template, and this is combined for the triggered network.

.. list-table:: Statistic Features
   :widths: 25 75
   :header-rows: 1

   * - Feature name
     - Description
   * - ``phasetd``
     - Use a histogram of expected phase and time differences, and amplitude ratio, for signals to determine a factor to be added for the signal rate.
   * - ``sensitive_volume``
     - Signal rate is expected to be proportional to the cube of the sensitive distance. This feature adds a factor of :math:`log(\sigma^3)` minus a benchmark value, to make this zero in many cases.
   * - ``normalize_fit_rate``
     - Normalise the exponential fits to use a rate rather than an absolute count of triggers. This means that statistics should be comparable over differently-sized analyses.
   * - ``dq``
     - Apply a reweighting factor according to the rates of triggers during data-quality flags vs the rate outside this. Must supply a reranking file using ``statistic-files`` for each detector, with stat attribute '{detector}-dq_stat_info'
   * - ``kde``
     - Use a file to re-rank according to the signal and density rates calculated using a KDE approach. Must supply two reranking files using ``statistic-files`` with stat attributes 'signal-kde_file' and 'template-kde_file' respectively.
   * - ``chirp_mass``
     - Apply a factor of :math:`log((M_c / 20) ^{11 / 3})` to the statistic. This makes the signal rate uniform over chirp mass, as this factor cancels out the power of -11 / 3 caused by differences in the density of template placement.

.. list-table:: Statistic Keywords
   :widths: 25 75
   :header-rows: 1

   * - Keyword
     - Description
   * - ``benchmark_lograte``
     - This is a numerical factor to be subtracted from the log rate ratio in order to alter the dynamic range. Default -14.6.
   * - ``minimum_statistic_cutoff``
     - Cutoff for the statistic in order to avoid underflowing and very small statistic values. Default -30.
   * - ``alpha_below_thresh``
     - The fit coefficient (alpha) below the fit threshold (defined in the fit_by_template jobs) will be replaced by a standard value. This is as below this threshold, Gaussian noise can dominate over the glitch response that dominates above the threshold, and rates will be underestimated, boosting quiet things in noisy templates. For Gaussian noise, this will be approximately 6 (the default). To use whatever the fit value is, supply alpha_below_thresh:None.
   * - ``reference_ifos``
     - If using the ``sensitive_volume`` feature, these are the detectors used to determine the benchmark value by which the sensitive volume is compared. We use the median sensitive volume in the network of detectors supplied. Default H1,L1.
   * - ``max_chirp_mass``
     - If using the ``chirp_mass`` feature, this chirp mass defines a maximum weighting which can be applied to the statistic.
   * - ``sngl_ranking_*``
     - This is used to provide the keyword arguments to functions in `the events.ranking module <https://pycbc.org/pycbc/latest/html/_modules/pycbc/events/ranking.html>`_. For example, to use a different psdvar threshold in the newsnr_sgveto_psdvar_threshold function, we would use ``sngl_ranking_psd_var_val_threshold:10``.