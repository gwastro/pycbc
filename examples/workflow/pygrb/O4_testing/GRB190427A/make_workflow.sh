pycbc_make_offline_grb_workflow \
	--config-files \
	executables.ini \
	analysis_o3a.ini \
	injections_o3a.ini \
	postprocessing_o3a.ini \
	offline_o3a.ini \
	GRB190427A.ini \
	--config-overrides \
	workflow:output-directory:${PWD}/output \
	--workflow-name test_wf \
	--output output 
