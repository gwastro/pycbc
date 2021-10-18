#!/bin/sh
for f in cpnest_stub.ini emcee_stub.ini emcee_pt_stub.ini dynesty_stub.ini ultranest_stub.ini epsie_stub.ini; do
        echo $f
	pycbc_inference \
        --config-files `dirname $0`/simp.ini `dirname $0`/$f \
        --output-file $f.hdf \
        --nprocesses 2 \
        --seed 10 \
        --force
done

pycbc_inference_plot_posterior --input-file \
emcee_stub.ini.hdf:emcee \
emcee_pt_stub.ini.hdf:emcee_pt \
dynesty_stub.ini.hdf:dynesty \
ultranest_stub.ini.hdf:ultranest \
epsie_stub.ini.hdf:espie \
cpnest_stub.ini.hdf:cpnest \
--output-file sample.png
