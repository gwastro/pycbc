.. _ahopedatafindmod:
###########################################
The ahope datafind and validation module
###########################################

=============
Introduction
=============

This page is designed to give you an introduction to the capabilities of the
ahope datafind and validation module and how to use this as part of an ahope
workflow.

This module is designed to be able to support multiple ways of obtaining
data (different codes/interfaces whatever). Currently we only support datafind
through the glue.datafind module (which is equivalent to using gw_data_find).

This module will run the necessary queries to the datafind server to obtain
locations for frame files at the specfied times for each interferometer.

Optionally, it can also run a set of tests to verify this output and act
accordingly. This includes

* A check that the all times in the input segment lists are covered with frames,
and methods for dealing with cases where this is not true.
* A check that all returned frame files actually exist and are accessible on the
cluster.
* A check that segment_summary flags are defined for all frames that have been
returend.

=======
Usage
=======

Using this module requires a number of things

* A configuration file (or files) containing the information needed to tell this module how to generate the segments (described below).
* An initialized instance of the ahope workflow class, containing the ConfigParser.
* An ifo-keyed dictionary of glue.segments.segmentlist instances containing the times that should be analysed for each ifo. See :ref:`ahopesegmentsmod` for documentation of the segments module, which in most cases should be used to obtain this input.

The module is then called according to

.. autofunction:: pycbc.ahope.setup_segment_generation
       :noindex:

-------------------------
Configuration file setup
-------------------------

Here we describe the options given in the configuration file used in the ahope
workflow that will be needed in this section

$$$$$$$$$$$$$$$$$$$$$$$$$
[ahope-datafind] section
$$$$$$$$$$$$$$$$$$$$$$$$$

The configuration file must have an [ahope-datafind] section, which is used to
tell the workflow how to generate the datafind calls. The first option to choose
and provide is

datafind-method = VALUE

The choices here and their description are as described below

* AT_RUNTIME_SINGLE_FRAMES - Find frame files at runtime using setup_datafind_runtime_frames_single_call_perifo, use one datafind query for every interferometer. This may be the quickest approach, but some of the frames returned will not be suitable for analysis. The output list will contain a single entry for every frame file returned by the datafind queries, which will enable pegasus to more easily track input files for jobs that need to read from frame files.
* AT_RUNTIME_MULTIPLE_FRAMES - Find frame files at runtime using setup_datafind_runtime_frames_multi_calls_perifo, use one datafind query for every science segment. This may be a little slower than the above, but only frames that overlap analysable data stretches will be returned. The output list will contain a single entry for every frame file returned by the datafind queries, which will enable pegasus to more easily track input files for jobs that need to read from frame files.
* AT_RUNTIME_SINGLE_CACHES - Find frame files at runtime using setup_datafind_runtime_cache_single_call_perifo, use one datafind query for every interferometer. This may be the quickest approach, but some of the frames returned will not be suitable for analysis. The output list will contain a single entry for every call made to the datafind server, which will correspond to a .lcf "frame cache" file in the output directory.
* AT_RUNTIME_MULTIPLE_CACHES - Find frame files at runtime using setup_datafind_runtime_cache_multi_calls_perifo, use one datafind query for every science segment. This may be a little slower than the above, but only frames that overlap analysable data stretches will be returned. The output list will contain a single entry for every call made to the datafind server, which will correspond to a .lcf "frame cache" file in the output directory.

Each of these options will describe which subfunction to use. These are described here

.. autofunction:: pycbc.ahope.setup_datafind_runtime_cache_multi_calls_perifo
          :noindex:

.. autofunction:: pycbc.ahope.setup_datafind_runtime_cache_single_call_perifo
          :noindex:

.. autofunction:: pycbc.ahope.setup_datafind_runtime_frames_multi_calls_perifo
          :noindex:

.. autofunction:: pycbc.ahope.setup_datafind_runtime_frames_single_call_perifo
          :noindex:


When using any of these sub-modules the following other configuration options apply in the [ahope-datafind] section

* datafind-X1-frame-type = NAME - REQUIRED. Where X1 is replaced by the ifo name for each ifo. The NAME should be the full frame type, which is used when querying the database.
* datafind-ligo-datafind-server = URL - OPTIONAL. If provided use this server when querying for frames. If not provided, which is recommended for most applications, then the LIGO_DATAFIND_SERVER environment variable will be used to determine this.
* datafind-check-segment-gaps = STRING - OPTIONAL (default = "no_test"). If this option takes any value other than 'no_test' ahope will check that the local datafind server has returned frames covering all of the listed science times. Its behaviour is then as follows
  * 'no_test': Do not perform this test. Any discrepancies will cause later failures.
  * 'warn': Perform the test, print warnings covering any discrepancies but do nothing about them. Discrepancies will cause failures later in the workflow.
  * 'update_times': Perform the test, print warnings covering any discrepancies and update the input science times to remove times that are not present on the host cluster.
  * 'raise_error': Perform the test. If any discrepancies occur, raise a ValueError.
* datafind-check-frames-exist = STRING - OPTIONAL (default = "no_test"). If this options takes any value other than 'no_test' ahope will check that the frames returned by the local datafind server are accessible from the machine that is running ahope. Its behaviour is then as follows
  * 'no_test': Do not perform this test. Any discrepancies will cause later failures.
  * 'warn': Perform the test, print warnings covering any discrepancies but do nothing about them. Discrepancies will cause failures later in the workflow.
  * 'update_times': Perform the test, print warnings covering any discrepancies and update the input science times to remove times that are not present on the host cluster.
  * 'raise_error': Perform the test. If any discrepancies occur, raise a ValueError.
* datafind-check-segment-summary = STRING - OPTIONAL (default = "no_test"). If this option takes any value other than 'no_test' ahope will check that all frames returned by datafind are covered by the segment_summary table (for the science flag). Its behaviour is then as follows:
  * 'no_test': Do not perform this test. 
  * 'warn': Perform the test, print warnings covering any discrepancies but do nothing about them.
  * 'raise_error': Perform the test. If any discrepancies occur, raise a ValueError.


$$$$$$$$$$$$$$$
[executables]
$$$$$$$$$$$$$$$

Currently no executables are needed in the datafind section. Ahope will use the
glue.datafind module to run the calls to the datafind server.

$$$$$$$$$$$$$$$$$$$
Other sections
$$$$$$$$$$$$$$$$$$$

_______________
[datafind]
_______________

The other section that can be used in the datafind module is the [datafind]
section. This section contains option,value pairs that will be send as key-word
arguments when calling the datafind server. Valid options here

* match=STRING - If given return only those frames matching the given regular expression.
* urltype=TYPE - If given restrict the returned frames to the given scheme (e.g. "file").

the on_gaps keyword argument is not supported as sanity checking is handled by ahope. This is always set to 'ignore' (this can be overwritten, we don't recommend this).

==========================================
:mod:`pycbc.ahope.datafind_utils` Module
==========================================

This is complete documentation of this module's code

.. automodule:: pycbc.ahope.datafind_utils
    :noindex:
    :members:
    :undoc-members:
    :show-inheritance:

