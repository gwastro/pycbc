
N_INJ=10

pycbc_create_injections \
       --ninjections $N_INJ \
       --output-file injection_hlv_10_bbh.hdf\
       --config-files injection_hlv_10_bbh.ini\
       --force \
       --seed 233
