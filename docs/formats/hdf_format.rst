############################################################
HDF files within the PyCBC workflow
############################################################

.. note::

    Format specifications are provided here to aid in development. The canonical
    definition, as always, lives within the code itself. 

=========================
pycbc_inspiral
=========================

All keys in the inspiral output are prefixed with the IFO name, e.g. H1, L1. Currently,
only a single ifo is present in each file, but at a future date, multiple may
be allowed.

The following table consists of columns of trigger data. Each column is of the same length
and an index into one column corresponds to the same trigger in each of the other columns.

.. csv-table:: Column vectors of trigger data
   :header: "path", "description"

   "IFO/snr", ""
   "IFO/coa_phase", ""
   "IFO/end_time", ""
   "IFO/chisq", ""
   "IFO/chisq_dof", ""
   "IFO/bank_chisq", ""
   "IFO/bank_chisq_dof", ""
   "IFO/cont_chisq", ""
   "IFO/cont_chisq_dof", ""
   "IFO/template_duration", ""
   "IFO/sigmasq", ""
   
.. csv-table:: Additional Data
   :header: "path", "description"
   
   "IFO/search/start_time", ""
   "IFO/search/end_time", ""


The following are columns that exists in the file, but should not be used by any user.
Its definition or existence is subject to change without notice.

.. csv-table:: reserved columns
   :header: "path"
   
   "IFO/template_hash"





