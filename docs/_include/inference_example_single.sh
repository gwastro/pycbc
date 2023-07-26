set -e
export OMP_NUM_THREADS=1
sh ../../examples/inference/single/get.sh
sh ../../examples/inference/single/run.sh
sh ../../examples/inference/single/plot.sh
