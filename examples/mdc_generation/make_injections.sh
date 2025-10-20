
N_INJ=10

pycbc_create_injections \
       --ninjections $N_INJ \
       --output-file injection_10inj.hdf\
       --config-files injections.ini\
       --force \
       --seed 233
