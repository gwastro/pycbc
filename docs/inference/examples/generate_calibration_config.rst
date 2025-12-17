.. _generate_calibration_config:

----------------------------------------------------------------------------
Example: Generating calibration configuration file for GW150914 and GW170817
----------------------------------------------------------------------------

When analyzing real gravitational-wave data, it is often desirable to
marginalize over detector calibration uncertainties. PyCBC supports this through
a calibration configuration file that specifies frequency-dependent amplitude
and phase uncertainties for each interferometer.

This section describes how to download the official LIGO–Virgo calibration
uncertainty files and generate a calibration configuration file for GW150914 and
GW170817.

Downloading calibration uncertainty files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, we need to download the data from the `Gravitational Wave Open Science
Center <https://www.gwosc.org>`_. Run:

Calibration uncertainty files are provided by the LIGO–Virgo KAGRA Collaboration
and are available from the `LIGO DCC<https://dcc.ligo.org/T2100313/public>`.

First, create a directory to store the calibration files:

.. code-block:: bash

   OUT_DIR=calibration_uncertainty_files
   mkdir -p ${OUT_DIR}
   cd ${OUT_DIR}

Download the LIGO calibration uncertainty files. You need to do it once and then
use them for generation of any calibration configuration file in O1, O2, or O3
observation run of LIGO and Virgo detector network:

.. code-block:: bash

   wget https://dcc.ligo.org/public/0177/T2100313/003/LIGO_O1_cal_uncertainty.tgz
   wget https://dcc.ligo.org/public/0177/T2100313/003/LIGO_O2_cal_uncertainty.tgz
   wget https://dcc.ligo.org/public/0177/T2100313/003/LIGO_O3_cal_uncertainty.tgz

Download the Virgo calibration uncertainty files:

.. code-block:: bash

   wget https://dcc.ligo.org/public/0177/T2100313/003/Virgo_O2_cal_uncertainty.tgz
   wget https://dcc.ligo.org/public/0177/T2100313/003/Virgo_O3_cal_uncertainty.tgz

Extract all downloaded archives:

.. code-block:: bash

   tar -xzvf LIGO_O1_cal_uncertainty.tgz
   tar -xzvf LIGO_O2_cal_uncertainty.tgz
   tar -xzvf LIGO_O3_cal_uncertainty.tgz

   tar -xzvf Virgo_O2_cal_uncertainty.tgz
   tar -xzvf Virgo_O3_cal_uncertainty.tgz

After extraction, the directory will contain subdirectories for each observing
run and interferometer, with frequency-dependent calibration uncertainty
envelopes.

Generating a calibration configuration file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once the calibration uncertainty files are available locally, use
``pycbc_inference_create_calibration_config`` to generate a calibration
configuration file for GW150914.

Set the path to the extracted calibration files:

.. code-block:: bash

   CALIB_ENV_FILE_PATH=./calibration_uncertainty_files

Generate the calibration configuration file for the Hanford (H1) and Livingston
(L1) detectors:

.. code-block:: bash

   pycbc_inference_create_calibration_config \
       --calibration-files-path ${CALIB_ENV_FILE_PATH} \
       --ifos H1 L1 \
       --minimum-frequency H1:20 L1:20 \
       --maximum-frequency H1:1000 L1:1000 \
       --gps-time 1126259462.43 \
       --correction-type H1:data L1:data \
       --tag GW150914_095045

This command produces a calibration configuration file that can be supplied to
``pycbc_inference`` to enable calibration marginalization.

.. note::
   The calibration uncertainty files and frequency ranges should match the
   observing run, detectors, and frequency range used in the analysis.

Example: GW170817 (H1–L1–V1)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

For GW170817, calibration uncertainties are included for the Hanford (H1),
Livingston (L1), and Virgo (V1) detectors. The calibration configuration file
can be generated as follows:

.. code-block:: bash

   pycbc_inference_create_calibration_config \
       --calibration-files-path ${CALIB_ENV_FILE_PATH} \
       --ifos H1 L1 V1 \
       --minimum-frequency H1:21.6 L1:21.6 V1:21.6 \
       --maximum-frequency H1:2048 L1:2048 V1:2048 \
       --gps-time 1187008882.42 \
       --correction-type H1:data L1:data V1:template \
       --tag GW170817_124104

.. note::

   Unlike the GW150914 example, the Virgo (V1) detector for GW170817 uses
   ``template`` calibration corrections, meaning that calibration uncertainties
   are applied to the waveform model or template rather than directly to the
   data stream. The Hanford (H1) and Livingston (L1) detectors continue to use
   ``data`` calibration type.


