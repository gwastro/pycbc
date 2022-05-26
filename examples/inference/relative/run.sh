OMP_NUM_THREADS=1 python -m cProfile -o log `which pycbc_inference` \
--config-file `dirname "$0"`/relative.ini \
--nprocesses=1 \
--output-file relative.hdf \
--seed 0 \
--force \
--verbose
