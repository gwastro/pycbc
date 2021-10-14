pycbc_submit_dax --dax gw.dax --no-grid --no-create-proxy \
--local-dir ./ \
--no-query-db \
--execution-sites condorpool_shared,condorpool_symlink,condorpool_copy \
--staging-sites condorpool_symlink=local,condorpool_copy=local
