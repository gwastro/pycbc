############################################################
HDF files within the PyCBC workflow
############################################################

.. note::

    Format specifications are provided here to aid in development. The canonical
    definition, as always, lives within the code itself. 

=========================
single inspiral triggers
=========================

*****************
Executables
*****************

 * pycbc_inspiral

*****************
Specification
*****************

All keys in the inspiral output are prefixed with the IFO name, e.g. H1, L1. Currently,
only a single ifo is present in each file, but at a future date, multiple may
be allowed.

The following table consists of columns of trigger data. Each column is of the same length
and an index into one column corresponds to the same trigger in each of the other columns.

.. csv-table:: Column vectors of trigger data
   :header: "path", "description"

   "IFO/snr", "The mangitude of the complex SNR"
   "IFO/coa_phase", "The phase of the complex SNR"
   "IFO/end_time", "The gps time of the trigger"
   "IFO/chisq", "Value of the bruce power chisq"
   "IFO/chisq_dof", "Not DOF. The number of bins in the chisq. DOF = 2 * (num_bins -1)"
   "IFO/bank_chisq", "Value of the bank chisq"
   "IFO/bank_chisq_dof", "Number of templates used to construct the bank chisq"
   "IFO/cont_chisq", "Value of the autochisq"
   "IFO/cont_chisq_dof", "Number of dof for the auto chisq"
   "IFO/template_duration", "Duration of the template approximant used for this trigger"
   "IFO/sigmasq", "The weighted power of the template, placed at 1Mpc, used for this trigger"
   
.. csv-table:: Additional Data
   :header: "path", "description"
   
   "IFO/search/start_time", "Array of gps times which denote the start of a valid period of triggers"
   "IFO/search/end_time", "Array of gps times which denote the corresponding end of a vlid period of triggers"


The following are columns that exists in the file, but should not be used by any user.
Its definition or existence is subject to change without notice.

.. csv-table:: reserved columns
   :header: "path"
   
   "IFO/template_hash"
   
==================================
combined single inspiral triggers
==================================

*****************
Executables
*****************

 * pycbc_coinc_mergetrigs

*****************
Specification
*****************

All keys in the inspiral output are prefixed with the IFO name, e.g. H1, L1. Currently,
only a single ifo is present in each file, but at a future date, multiple may
be allowed.

The following table consists of columns of trigger data. Each column is of the same length
and an index into one column corresponds to the same trigger in each of the other columns.


.. csv-table:: Column vectors of trigger data
   :header: "path", "description"

   "IFO/snr", "The mangitude of the complex SNR"
   "IFO/coa_phase", "The phase of the complex SNR"
   "IFO/end_time", "The gps time of the trigger"
   "IFO/chisq", "Value of the bruce power chisq"
   "IFO/chisq_dof", "Not DOF. The number of bins in the chisq. DOF = 2 * (num_bins -1)"
   "IFO/bank_chisq", "Value of the bank chisq"
   "IFO/bank_chisq_dof", "Number of templates used to construct the bank chisq"
   "IFO/cont_chisq", "Value of the autochisq"
   "IFO/cont_chisq_dof", "Number of dof for the auto chisq"
   "IFO/template_duration", "Duration of the template approximant used for this trigger"
   "IFO/sigmasq", "The weighted power of the template, placed at 1Mpc, used for this trigger"
   "IFO/template_id", "The unique template id value. This is the index into the hdf template file format"

The key feature that the combined trigger format adds is the convenience of precalculated region
references to access only the data produced by a given template. These values are stored in region 
reference arrays. The length of each array is the same as the number of templates, and an index
into the array matches the template_id number. Each array directly maps to a single column.

.. csv-table:: region reference arrays
   :header: "path"

    "IFO/bank_chisq_dof_template"
    "IFO/bank_chisq_template"
    "IFO/chisq_dof_template"
    "IFO/chisq_template"
    "IFO/coa_phase_template"
    "IFO/cont_chisq_dof_template"
    "IFO/cont_chisq_template"
    "IFO/end_time_template"
    "IFO/sigmasq_template"
    "IFO/snr_template"
    "IFO/template_boundaries"
    "IFO/template_duration_template"

.. csv-table:: Additional Data
   :header: "path", "description"
   
   
   "IFO/search/start_time", "Array of gps times which denote the start of a valid period of triggers"
   "IFO/search/end_time", "Array of gps times which denote the corresponding end of a vlid period of triggers"

*********************
Example uses
*********************

Accessing triggers by template

.. code-block:: python

        import h5py
        f = h5py.File('H1-testdata.hdf')
        snr_regs = f['H1/snr_template']
        snr_template_0 = f['H1/snr'][snr_regs[0]]



