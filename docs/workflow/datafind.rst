.. _workflowdatafindmod:

###########################################
The workflow datafind and validation module
###########################################

=============
Introduction
=============

This page is designed to give you an introduction to the capabilities of the
workflow datafind and validation module and how to use this as part of a pycbc
workflow.

This module is designed to be able to support multiple ways of obtaining
data (different codes/interfaces whatever). Currently we only support datafind
through the glue.datafind module (which is equivalent to using gw_data_find).

This module will run the necessary queries to the datafind server to obtain
locations for frame files at the specfied times for each interferometer.

Optionally, it can also run a set of tests to verify this output and act
accordingly. This includes

* A check that the all times in the input segment lists are covered with frames, and methods for dealing with cases where this is not true.
* A check that all returned frame files actually exist and are accessible on the cluster.
* A check that segment_summary flags are defined for all frames that have been returned.

=======
Usage
=======

Using this module requires a number of things

* A configuration file (or files) containing the information needed to tell this module how to generate the segments (described below).
* An initialized instance of the pycbc workflow class, containing the ConfigParser.
* An ifo-keyed dictionary of glue.segments.segmentlist instances containing the times that should be analysed for each ifo. See :ref:`workflowsegmentsmod` for documentation of the segments module, which in most cases should be used to obtain this input.

The module is then called according to

.. autofunction:: pycbc.workflow.setup_segment_generation
       :noindex:

-------------------------
Configuration file setup
-------------------------

Here we describe the options given in the configuration file used in the
workflow that will be needed in this section

$$$$$$$$$$$$$$$$$$$$$$$$$$$
[workflow-datafind] section
$$$$$$$$$$$$$$$$$$$$$$$$$$$

The configuration file must have an [workflow-datafind] section, which is used to
tell the workflow how to generate the datafind calls. The first option to choose
and provide is

datafind-method = VALUE

The choices here and their description are as described below

* AT_RUNTIME_SINGLE_FRAMES - Find frame files at runtime using setup_datafind_runtime_frames_single_call_perifo, use one datafind query for every interferometer. This may be the quickest approach, but some of the frames returned will not be suitable for analysis. The output list will contain a single entry for every frame file returned by the datafind queries, which will enable pegasus to more easily track input files for jobs that need to read from frame files.
* AT_RUNTIME_MULTIPLE_FRAMES - Find frame files at runtime using setup_datafind_runtime_frames_multi_calls_perifo, use one datafind query for every science segment. This may be a little slower than the above, but only frames that overlap analysable data stretches will be returned. The output list will contain a single entry for every frame file returned by the datafind queries, which will enable pegasus to more easily track input files for jobs that need to read from frame files.
* AT_RUNTIME_SINGLE_CACHES - Find frame files at runtime using setup_datafind_runtime_cache_single_call_perifo, use one datafind query for every interferometer. This may be the quickest approach, but some of the frames returned will not be suitable for analysis. The output list will contain a single entry for every call made to the datafind server, which will correspond to a .lcf "frame cache" file in the output directory.
* AT_RUNTIME_MULTIPLE_CACHES - Find frame files at runtime using setup_datafind_runtime_cache_multi_calls_perifo, use one datafind query for every science segment. This may be a little slower than the above, but only frames that overlap analysable data stretches will be returned. The output list will contain a single entry for every call made to the datafind server, which will correspond to a .lcf "frame cache" file in the output directory.
* FROM_PREGENERATED_LCF_FILES - Supply a set of pregenerated .lcf files containing a list of frame files to use for analysis. This option is intended to be used in cases where a datafind server is not available. Be warned that data does move around on clusters so on standard LDG clusters the AT_RUNTIME options are recommended.

Each of these options will describe which subfunction to use. These are described here

.. autofunction:: pycbc.workflow.setup_datafind_runtime_cache_multi_calls_perifo
          :noindex:

.. autofunction:: pycbc.workflow.setup_datafind_runtime_cache_single_call_perifo
          :noindex:

.. autofunction:: pycbc.workflow.setup_datafind_runtime_frames_multi_calls_perifo
          :noindex:

.. autofunction:: pycbc.workflow.setup_datafind_runtime_frames_single_call_perifo
          :noindex:

.. autofunction:: pycbc.workflow.setup_datafind_from_pregenerated_lcf_files
          :noindex:


When using any of the AT_RUNTIME sub-modules the following other configuration options apply in the [workflow-datafind] section

* datafind-X1-frame-type = NAME - REQUIRED. Where X1 is replaced by the ifo name for each ifo. The NAME should be the full frame type, which is used when querying the database.
* datafind-ligo-datafind-server = URL - OPTIONAL. If provided use this server when querying for frames. If not provided, which is recommended for most applications, then the LIGO_DATAFIND_SERVER environment variable will be used to determine this.
* datafind-backup-datafind-server = URL - OPTIONAL. This option is only available when using AT_RUNTIME_SINGLE_FRAMES or AT_RUNTIME_MULTIPLE_FRAMES. If given it will query a second datafind server (ie. a remote server) using gsiftp urltypes. This will then allow frames to be associated with both a file:// and gsiftp:// url, in the case that your local site is missing a frame file, or the file is not accessible, pegasus will copy the file from gsiftp://. **NOTE** This will not catch the case that the frame file is available at the **start** of a workflow but goes missing later. Pegasus can copy **all** frame files around at the start of the workflow, but you may not want this (remove symlink option from the basic_pegasus.conf if you want this).

When using the PREGENERATED sub-module the following configuartion options apply in the [workflow-datafind] section:

* datafind-pregenerated-cache-file-x1 = Path/to/file.lcf. This should be specified independently for each ifo and points to the pregenerated files.

The following configuration options apply in the [workflow-datafind] section for all sub-modules and can be used as sanity checks:

* datafind-check-segment-gaps = STRING - OPTIONAL (default = "no_test"). If this option takes any value other than 'no_test' the workflow module will check that the local datafind server has returned frames covering all of the listed science times. Its behaviour is then as follows:

  * 'no_test': Do not perform this test. Any discrepancies will cause later failures.
  * 'warn': Perform the test, print warnings covering any discrepancies but do nothing about them. Discrepancies will cause failures later in the workflow.
  * 'update_times': Perform the test, print warnings covering any discrepancies and update the input science times to remove times that are not present on the host cluster.
  * 'raise_error': Perform the test. If any discrepancies occur, raise a ValueError.

* datafind-check-frames-exist = STRING - OPTIONAL (default = "no_test"). If this options takes any value other than 'no_test' the workflow module will check that the frames returned by the local datafind server are accessible from the machine that is running the workflow generation. Its behaviour is then as follows:

  * 'no_test': Do not perform this test. Any discrepancies will cause later failures.
  * 'warn': Perform the test, print warnings covering any discrepancies but do nothing about them. Discrepancies will cause failures later in the workflow.
  * 'update_times': Perform the test, print warnings covering any discrepancies and update the input science times to remove times that are not present on the host cluster.
  * 'raise_error': Perform the test. If any discrepancies occur, raise a ValueError.

* datafind-check-segment-summary = STRING - OPTIONAL (default = "no_test"). If this option takes any value other than 'no_test' the workflow module will check that all frames returned by datafind are covered by the segment_summary table (for the science flag). Its behaviour is then as follows:

  * 'no_test': Do not perform this test. 
  * 'warn': Perform the test, print warnings covering any discrepancies but do nothing about them.
  * 'raise_error': Perform the test. If any discrepancies occur, raise a ValueError.


$$$$$$$$$$$$$$$
[executables]
$$$$$$$$$$$$$$$

Currently no executables are needed in the datafind section. Workflow will use
the glue.datafind module to run the calls to the datafind server.

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

the on_gaps keyword argument is not supported as sanity checking is handled by the workflow module. This is always set to 'ignore' (this can be overwritten, we don't recommend this).

==========================================
:mod:`pycbc.workflow.datafind` Module
==========================================

This is complete documentation of this module's code

.. automodule:: pycbc.workflow.datafind
    :noindex:
    :members:
    :undoc-members:
    :show-inheritance:

