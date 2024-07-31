set -e
export OMP_NUM_THREADS=1
sh ../../examples/inference/lisa_smbhb_ldc/get.sh
sh ../../examples/inference/lisa_smbhb_ldc/run.sh
python ../../examples/inference/lisa_smbhb_ldc/advanced_plot.py
