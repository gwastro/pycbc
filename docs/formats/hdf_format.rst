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







