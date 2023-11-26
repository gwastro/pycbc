# In this example, we run the workflow generation script and immediately
# submit it to the compute cluster
export PATH=$PATH:$PWD
python simple.py \
--workflow-name test \
--config-files simple.ini \
--submit-now
