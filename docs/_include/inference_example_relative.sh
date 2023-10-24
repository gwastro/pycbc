set -e
export OMP_NUM_THREADS=1
sh ../../examples/inference/relative/get.sh
sh ../../examples/inference/relative/run.sh
sh ../../examples/inference/relative/plot.sh
