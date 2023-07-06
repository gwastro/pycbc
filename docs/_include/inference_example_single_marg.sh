set -e
export OMP_NUM_THREADS=1
sh ../../examples/inference/single/get.sh
sh ../../examples/inference/single/run_marg.sh
sh ../../examples/inference/single/run_instant.sh
