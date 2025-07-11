
# Provide a path to downloaded calibration files for O1-O2-O3 
# These files can be downloaded from
# https://dcc.ligo.org/T2100313/public

CALIB_ENV_FILE_PATH=./calibration_uncertainty_files

# Generate calibration configuration file for GW150914
pycbc_inference_create_calibration_config --calibration-files-path ${CALIB_ENV_FILE_PATH} --ifos H1 L1 --minimum-frequency H1:20 L1:20 --maximum-frequency H1:1000 L1:1000  --gps-time 1126259462.43 --correction-type H1:data L1:data
