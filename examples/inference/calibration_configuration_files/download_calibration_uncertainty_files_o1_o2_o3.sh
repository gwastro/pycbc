
# Provide the path to download calibration files

OUT_DIR=$1

if [[ -z $OUT_DIR ]]
  then
    echo "ERROR: output-directory is mandatory argument. \n";
    exit 1;
  fi

cd $OUT_DIR

# Download LIGO calibration uncertainties
wget https://dcc.ligo.org/public/0177/T2100313/003/LIGO_O1_cal_uncertainty.tgz .
wget https://dcc.ligo.org/public/0177/T2100313/003/LIGO_O2_cal_uncertainty.tgz .
wget https://dcc.ligo.org/public/0177/T2100313/003/LIGO_O3_cal_uncertainty.tgz .

# Download Virgo files
wget https://dcc.ligo.org/public/0177/T2100313/003/Virgo_O2_cal_uncertainty.tgz .
wget https://dcc.ligo.org/public/0177/T2100313/003/Virgo_O3_cal_uncertainty.tgz

#
tar -xzvf LIGO_O1_cal_uncertainty.tgz 
tar -xzvf LIGO_O2_cal_uncertainty.tgz
tar -xzvf LIGO_O3_cal_uncertainty.tgz

tar -xzvf Virgo_O2_cal_uncertainty.tgz
tar -xzvf Virgo_O3_cal_uncertainty.tgz
