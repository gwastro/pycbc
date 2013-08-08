###############################
Matched-filter section
###############################

=============
Introduction
=============

The matched-filter section of ahope is responsible for matched-filtering the
data against the template bank(s) from the template bank section and generating
a list of "triggers" for each interferometer. These triggers should be a list
of any event where the signal to noise ratio and any signal consistency test
are such that that point should be sent forward to check for coincidence in
other ifos.

Any single-ifo signal consistency tests (ie. chi-squared tests etc.) should
be computed in this section and stored within the lists of triggers. Ahope
does not make any specifications on the output format, but obviously code in
the next stages of the workflow must know how to process that input.

The matched-filtering section should be as independent of the other stages of
the workflow as possible. This means that we don't require the data read in by
matched-filter jobs to match that read in by template bank jobs (however this
may be desireable in some cases, so **should** be possible where sensible).
Options should also not be hardcoded (so there are no cases where an option
that gets sent to a template bank job also gets sent to a matched-filter job
without any way to stop that). However, it is possible to duplicate options
where this is desireable (see :ref:`~pycbc/ahope/initialization_inifile.rst`).

The return from the matched-filter section of ahope is a list of AhopeOutFile
objects corresponding to each actual file (one for each job) that will be 
generated within the workflow. Each object would contain 4 pieces of
information (see :class:`~pycbc/ahope/AhopeOutFileList` and 
:class:`~pycbc/ahope/AhopeOutFile` for more):

- The location of the output file (which may not yet exist)
- The ifo that the file is valid for
- The time span that the file is valid for (generally this will be the time
  span over which triggers were generated).
- The dax job that will generate the bank (if appropriate)

This will be the **only** thing that will be passed from the matched-filter
section to the future sections.

More details of the code can be found below:

:mod:`pycbc.ahope.matchedfltr_utils` Module
--------------------------------------------

.. automodule:: pycbc.ahope.matchedfltr_utils
    :noindex:
    :members:
    :undoc-members:
    :show-inheritance:
