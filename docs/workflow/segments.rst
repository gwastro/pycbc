.. _workflowsegmentsmod:

#######################################
The workflow segment generation module
#######################################

=============
Introduction
=============

This page is designed to give you an introduction to the capabilities of the
pycbc workflow segment generation module and how to use this as part of a pycbc
workflow.

This module is designed to be able to support multiple ways of obtaining these
segments (different codes/interfaces whatever), although new code will always
be needed to support some code/interface that is not currently supported.

This module will generate science segments and any appropriate veto segments
and combine these together to identify a set of segments to be used in the
analysis. The various files will also be returned for later use in the analysis
(ie. for vetoing triggers with data-quality vetoes). The module also supports
generating cumulative and multiple-detector veto files for easy use with
ligolw_thinca and pipedown. If other workflows require similar combined files
these can be added on request.

=======
Usage
=======

Using this module requires a number of things

* A configuration file (or files) containing the information needed to tell this module how to generate the segments (described below).
* An initialized instance of the pycbc Workflow class, containing the ConfigParser.

The module is then called according to

.. autofunction:: pycbc.workflow.setup_segment_generation
       :noindex:

-------------------------
Configuration file setup
-------------------------

Here we describe the options given in the configuration file used in the
workflow that will be needed in this section

$$$$$$$$$$$$$$$$$$$$$$$$$$$
[workflow-segments] section
$$$$$$$$$$$$$$$$$$$$$$$$$$$

The configuration file must have a [workflow-segments] section, which is used to
tell the workflow how to construct the segments. The first option to choose
and provide is

segments-method = VALUE

The choices here and their description are as described below

* AT_RUNTIME - Use the setup_segment_gen_mixed to generate segments and generate all segment files at runtime
* CAT2_PLUS_DAG - Use the setup_segment_gen_mixed to generate segments, generate all veto files up to CATEGORY_1 at runtime, and add jobs to produce the remaining files to the workflow.
* CAT3_PLUS_DAG - Use the setup_segment_gen_mixed to generate segments, generate all veto files up to CATEGORY_2 at runtime, and add jobs to produce the remaining files to the workflow.
* CAT4_PLUS_DAG - Use the setup_segment_gen_mixed to generate segments, generate all veto files up to CATEGORY_3 at runtime, and add jobs to produce the remaining files to the workflow.

Each of these options will describe which subfunction to use. These are described here

.. autofunction:: pycbc.workflow.setup_segment_gen_mixed
          :noindex:

When using the setup_segment_gen_mixed function the following additional options apply

* segments-X1-science-name = NAME - REQUIRED. Where X1 is replaced by the ifo name for each ifo. The NAME should be the full channel name corresponding to analysable times for e.g. H1:DMT-SCIENCE:4
* segments-database-url = URL - REQUIRED. The URL to the segment databse that will be used to obtain this information
* segments-veto-definer-url = PATH - REQUIRED. The location to the veto-definer file that is used to identify which channels are CAT_1, which are CAT_2 etc.
* segments-veto-categories = COMMA-SEPARATED LIST OF INTS - OPTIONAL. Generate veto files for veto categories given by the ints in the list. These ranged from 1 through 4 or 5 for S5/S6 veto definers. Standard results have used categories 2,3,4.
* segments-minimum-segment-length = INT - OPTIONAL. If given, any segments of analysable data shorter than INT will not be included in the list of analysable times returned by this module.
* segments-generate-coincident-segments - OPTIONAL. Option takes no value. If given the module will generate cumulative, multiple detector coincidence files for easy use in ligolw_thinca and pipedown.
* segments-generate-segment-files - OPTIONAL (DEFAULT='always'). This option can be used if the user wants to re-use segment files generated previously. It is not recommended to use this option unless necessary. Options are

  * generate_segment_files='always' : DEFAULT: All files will be generated even if they already exist.
  * generate_segment_files='if_not_present': Files will be generated if they do not already exist. Pre-existing files will be read in and used.
  * generate_segment_files='error_on_duplicate': Files will be generated if they do not already exist. Pre-existing files will raise a failure.
  * generate_segment_files='never': Pre-existing files will be read in and used. If no file exists the code will fail.

$$$$$$$$$$$$$$$
[executables]
$$$$$$$$$$$$$$$

The following executable paths must be provided in the [executables] section
when running this module::

* segment_query = /home/ahnitz/local/lalsuite/bin/ligolw_segment_query
* segments_from_cats = /home/ahnitz/local/lalsuite/bin/ligolw_segments_from_cats
* llwadd = /home/ahnitz/local/lalsuite/bin/ligolw_add
* ligolw_combine_segments = /home/ahnitz/local/lalsuite/bin/ligolw_combine_segments

segment_query is used to obtain the science segments. segments_from_cats is used to obtain the files containing the CAT_1,2,3,4,5 segments. ligolw_combine_segments produces cumulative veto-files. llwadd is used to add the cumulative veto-files from different ifos together when producing cumulative, multiple-detector veto lists.

$$$$$$$$$$$$$$$$$$$
Other sections
$$$$$$$$$$$$$$$$$$$

For other sub-modules in the pycbc workflow module we would see sections like [segment_query], [segments_from_cats] etc. which would provide the options provided to those jobs. In this case the codes require rather specific input so for now these are hardcoded in this module and any segment like [segment_query] would either be ignored or
could break the code.

If there is a reason to do so we could add these sections in.

==========================================
:mod:`pycbc.workflow.segment` Module
==========================================

This is complete documentation of this module's code

.. automodule:: pycbc.workflow.segment 
    :noindex:
    :members:
    :undoc-members:
    :show-inheritance:

