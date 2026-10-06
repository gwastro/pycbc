
# Provide a path to downloaded calibration files for O4 
# These files can be downloaded from
# https://dcc.ligo.org/LIGO-T2500288/public

CALIB_ENV_FILE_PATH=./calibration_uncertainty_files

# Generate calibration configuration file for GW231123
pycbc_inference_create_calibration_config --calibration-files-path ${CALIB_ENV_FILE_PATH} --ifos H1 L1 --minimum-frequency H1:20 L1:20 --maximum-frequency H1:448 L1:448  --gps-time 1384782888.7 --correction-type H1:data L1:data
