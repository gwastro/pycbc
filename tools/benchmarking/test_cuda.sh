#!/bin/bash

source /opt/intel/bin/compilervars.sh intel64
source profiling_utils.sh

# To run a single-threaded MKL instance on every core:
# schemes=$(get_schemes 0 1)
#
# To run a single-threaded MKL instance on every core except one,
# and an CUDA instance on that remaining one:
# schemes=$(get_schemes 1 1)
#
# Run a single cuda instance (on CPU 1)
# schemes='cuda|1'

# Run two cuda instances on the same GPU and different CPUs
schemes='cuda:1|1 cuda:2|2'

for data in clean loud grumbly
do
	run_tests $data "${schemes}" `hostname`_${data}_cuda_1.dat
done

