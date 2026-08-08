
# Provide the path to download calibration files

OUT_DIR=calibration_uncertainty_files
mkdir -p $OUT_DIR
cd $OUT_DIR

# Download LIGO calibration uncertainties
wget https://dcc.ligo.org/public/0202/T2500288/001/LIGO_O4a_cal_uncertainty.tgz .

# Download Virgo files

# Extract the files
tar -xzvf LIGO_O4a_cal_uncertainty.tgz
 
