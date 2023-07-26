set -e
export OMP_NUM_THREADS=1
sh ../../examples/inference/margtime/get.sh
sh ../../examples/inference/margtime/run.sh
