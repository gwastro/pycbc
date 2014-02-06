#####################################
The ahope segment generation module
#####################################

=============
Introduction
=============

This page is designed to give you an introduction to the capabilities of the
ahope segment generation module and how to use this as part of an ahope
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

Using this module requires a number of things::

* A configuration file (or files) containing the information needed to tell this module how to generate the segments (described below).
* An initialized instance of the ahope workflow class, containing the ConfigParser.

The module is then called according to

.. autofunction:: pycbc.ahope.setup_segment_generation
       :noindex:

-------------------------
Configuration file setup
-------------------------

Here we describe the options given in the configuration file used in the ahope
workflow that will be needed in this section

$$$$$$$$$$$$$$$$$$$$$$$$$
[ahope-segments] section
$$$$$$$$$$$$$$$$$$$$$$$$$

The configuration file must have an [ahope-segments] section, which is used to
tell the workflow how to construct the segments. The first option to choose
and provide is

segments-method = VALUE

The choices here and their description are as described below::

* AT_RUNTIME - Use the setup_segment_gen_mixed to generate segments and generate all segment files at runtime
* CAT2_PLUS_DAG - Use the setup_segment_gen_mixed to generate segments, generate all veto files up to CATEGORY_1 at runtime, and add jobs to produce the remaining files to the workflow.
* CAT3_PLUS_DAG - Use the setup_segment_gen_mixed to generate segments, generate
 all veto files up to CATEGORY_2 at runtime, and add jobs to produce the remaini
ng files to the workflow.
* CAT4_PLUS_DAG - Use the setup_segment_gen_mixed to generate segments, generate
 all veto files up to CATEGORY_3 at runtime, and add jobs to produce the remaini
ng files to the workflow.

Each of these options will describe which subfunction to use. These are described here

.. autofunction:: pycbc.ahope.setup_segment_gen_mixed
          :noindex:

When using the setup_segment_gen_mixed function the following additional options apply::

* segments-X1-science-name = NAME - REQUIRED. Where X1 is replaced by the ifo name for each ifo. The NAME should be the full channel name corresponding to analysable times for e.g. H1:DMT-SCIENCE:4
* segments-database-url = URL - REQUIRED. The URL to the segment databse that will be used to obtain this information
* segments-veto-definer-url = PATH - REQUIRED. The location to the veto-definer file that is used to identify which channels are CAT_1, which are CAT_2 etc.
* segments-maximum-veto-category = INT - REQUIRED. Generate veto files for all veto categories up to and including INT. 4 or 5 is the standard value for ihope runs.
* segments-minimum-segment-length = INT - OPTIONAL. If given, any segments of analysable data shorter than INT will not be included in the list of analysable times returned by this module.
* segments-generate-coincident-segments - OPTIONAL. Option takes no value. If given the module will generate cumulative, multiple detector coincidence files for easy use in ligolw_thinca and pipedown.

$$$$$$$$$$$$$$$
[executables]
$$$$$$$$$$$$$$$

The following executable paths must be provided in the [executables] section
when running this module::

* segment_query = /home/ahnitz/local/lalsuite/bin/ligolw_segment_query
* segments_from_cats = /home/ahnitz/local/lalsuite/bin/ligolw_segments_from_cats
* llwadd = /home/ahnitz/local/lalsuite/bin/ligolw_add
* ligolw_segments_compat = /home/ahnitz/local/lalsuite/bin/ligolw_segments_compat

segment_query is used to obtain the science segments. segments_from_cats is used to obtain the files containing the CAT_1,2,3,4,5 segments. llwadd is used to add veto-files together when producing cumulative, multiple-detector veto lists. ligolw_segments_compat is used to ensure that the output of these codes are compatible with other output xml files. We want to create one front end to do what llwadd and ligolw_segments_compat does in one code.

$$$$$$$$$$$$$$$$$$$
Other sections
$$$$$$$$$$$$$$$$$$$

For other modules in ahope we would see sections like [segment_query], [segments_from_cats] etc. which would provide the options provided to those jobs. In this case the codes require rather specific input so for now these are hardcoded in this module and any segment like [segment_query] would either be ignored or
could break the code.

If there is a reason to do so we could add these sections in.

