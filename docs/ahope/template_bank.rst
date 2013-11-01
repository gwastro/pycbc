###############################
Template bank section
###############################

=============
Introduction
=============

The template bank section of ahope is responsible for gathering/generating the banks of waveforms that will be used to matched-filter against the data.

It should be able to run in a number of different modes

- With a pre-generated template bank (or banks) or generating banks on the fly
- Using the **same** template bank to analyse the same times in different detectors, or using a **different** template bank.
- Banks (if ~1 of them and it can be generated quickly) could be generated during ahope run time, or as part of the ahope workflow

The template bank section must be as independent of the matched-filtering section as possible (even though the two are obviously dependent on each other). It should be possible/easy for:

- The template bank analysis chunks to be longer/shorter than the inspiral chunks
- Allow for options sent to one job to be sent to the other **only** if desired. (No hardcoding of options that get sent to **both** matched-filter and template bank stage).

To do this I propose that the return from the template bank section of ahope is a list of template bank objects corresponding to each actual bank that will be (or are already) generated. Each template bank "object" would contain 4 pieces of information:

- The location of the template bank (which may not yet exist)
- The ifo that the bank is valid for
- The time span that the bank is valid for
- The dax job that will generate the bank (if appropriate)

This will be the **only** thing that will be passed from the template bank section to the future sections (including the matched filter section).

More details of the code can be found below:

:mod:`pycbc.ahope.tmpltbank_utils` Module
------------------------------------------

.. automodule:: pycbc.ahope.tmpltbank_utils
    :noindex:
    :members:
    :undoc-members:
    :show-inheritance:
