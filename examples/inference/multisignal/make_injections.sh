pycbc_create_injections --verbose --force \
    --ninjections 1 \
    --config-file event1_inj.ini \
    --output-file event1_inj.hdf

pycbc_create_injections --verbose --force \
    --ninjections 1 \
    --config-file event2_inj.ini \
    --output-file event2_inj.hdf

pycbc_merge_inj_hdf --injection-files event1_inj.hdf event2_inj.hdf --output-file combined_inj.hdf
